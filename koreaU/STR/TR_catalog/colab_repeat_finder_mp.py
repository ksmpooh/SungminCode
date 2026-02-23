'''
# 전 컨티그 병렬 처리 (코어 자동 감지)
python3 perfect_repeat_finder_mp.py \
  --min-span 9 --min-repeats 3 --min-motif-size 2 --max-motif-size 6 \
  --show-progress-bar \
  -o KPPD019.1 \
  -t 16 \
  KPPD019.1.fa
# → KPPD019.1.bed

# 특정 구간만 단일 잡으로
python3 perfect_repeat_finder_mp.py \
  -i chr1:100000-200000 \
  -o chr1_100k_200k \
  KPPD019.1.fa
# → chr1_100k_200k.bed
'''
###
#!/usr/bin/env python3
import argparse
import os
import re
import tempfile
import shutil
from types import SimpleNamespace
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed

import pyfastx
import tqdm

from utils.perfect_repeat_tracker import PerfectRepeatTracker
from utils.plot_utils import plot_results


def detect_repeats(input_sequence, filter_settings, verbose=False, show_progress_bar=False, debug=False):
    """Detect repeats in a given input sequence and return list of (start_0based, end, motif)."""
    # basic validations
    if not getattr(filter_settings, "min_motif_size") or filter_settings.min_motif_size < 1:
        raise ValueError(f"min_motif_size is set to {filter_settings.min_motif_size}. It must be at least 1.")
    if not getattr(filter_settings, "max_motif_size") or filter_settings.max_motif_size < filter_settings.min_motif_size:
        raise ValueError(f"max_motif_size is set to {filter_settings.max_motif_size}. It must be at least min_motif_size.")
    if not getattr(filter_settings, "min_repeats") or filter_settings.min_repeats < 1:
        raise ValueError(f"min_repeats is set to {filter_settings.min_repeats}. It must be at least 1.")
    if not getattr(filter_settings, "min_span") or filter_settings.min_span < 1:
        raise ValueError(f"min_span is set to {filter_settings.min_span}. It must be at least 1.")

    seq_original = input_sequence.upper()

    # interval defaults
    interval_start_0based = getattr(filter_settings, "interval_start_0based", 0)
    interval_end = getattr(filter_settings, "interval_end", len(seq_original))
    if interval_end is None:
        interval_end = len(seq_original)

    # clamp
    interval_start_0based = max(0, interval_start_0based)
    interval_end = min(len(seq_original), interval_end)
    if interval_start_0based >= interval_end:
        return []

    # trim Ns at ends inside the interval
    left = interval_start_0based
    right = interval_end
    while left < right and seq_original[left] == "N":
        left += 1
    while right > left and seq_original[right - 1] == "N":
        right -= 1

    if right <= left:
        return []

    subseq = seq_original[left:right]

    output_intervals = {}
    repeat_trackers = {}
    for motif_size in range(filter_settings.min_motif_size, filter_settings.max_motif_size + 1):
        repeat_tracker = PerfectRepeatTracker(
            motif_size=motif_size,
            min_repeats=filter_settings.min_repeats,
            min_span=filter_settings.min_span,
            input_sequence=subseq,
            output_intervals=output_intervals,
        )
        repeat_trackers[motif_size] = repeat_tracker

    end_position = len(subseq)
    position_iter = range(len(subseq))
    if show_progress_bar:
        position_iter = tqdm.tqdm(position_iter, unit=" bp", unit_scale=True, total=end_position)

    for _pos in position_iter:
        any_in_middle_of_repeat = False
        for repeat_tracker in repeat_trackers.values():
            repeat_tracker.advance()
            any_in_middle_of_repeat = repeat_tracker.is_in_middle_of_repeat() or any_in_middle_of_repeat

    for repeat_tracker in repeat_trackers.values():
        assert not repeat_tracker.advance(), f"{repeat_tracker.motif_size}bp motif RepeatTracker did not reach end of the sequence"
        repeat_tracker.done()

    translated = []
    for (start0, end), motif in sorted(output_intervals.items()):
        translated.append((start0 + left, end + left, motif))
    return translated


def _process_one_contig(task):
    """Worker: reopen FASTA, run detect_repeats on one contig (optionally interval), write to a temp bed file."""
    fasta_path = task["fasta_path"]
    chrom = task["chrom"]
    min_motif_size = task["min_motif_size"]
    max_motif_size = task["max_motif_size"]
    min_repeats = task["min_repeats"]
    min_span = task["min_span"]
    show_progress = task["show_progress"]
    interval = task.get("interval", None)
    tmp_dir = task["tmp_dir"]

    # reopen FASTA inside worker
    fa = pyfastx.Fasta(fasta_path)
    seq = fa[chrom].seq
    seq_len = len(seq)

    if interval:
        start0, end = interval
        start0 = max(0, start0)
        end = min(seq_len, end)
        per_filters = SimpleNamespace(
            min_motif_size=min_motif_size,
            max_motif_size=max_motif_size,
            min_repeats=min_repeats,
            min_span=min_span,
            interval_start_0based=start0,
            interval_end=end,
        )
        announced_len = end - start0
    else:
        per_filters = SimpleNamespace(
            min_motif_size=min_motif_size,
            max_motif_size=max_motif_size,
            min_repeats=min_repeats,
            min_span=min_span,
            interval_start_0based=0,
            interval_end=seq_len,
        )
        announced_len = seq_len

    # run detection
    output_intervals = detect_repeats(
        seq,
        per_filters,
        verbose=False,
        show_progress_bar=False,  # per-process progress off; we show global bar in main
        debug=False,
    )

    # write to temp file
    tmp_path = os.path.join(tmp_dir, f"{re.sub('[^A-Za-z0-9_.-]', '_', chrom)}.bed.tmp")
    with open(tmp_path, "wt") as out:
        for s, e, motif in output_intervals:
            out.write(f"{chrom}\t{s}\t{e}\t{motif}\n")

    return {"chrom": chrom, "count": len(output_intervals), "announced_len": announced_len, "tmp_path": tmp_path}


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_argument_group("Repeat Filters")
    group.add_argument("-min", "--min-motif-size", default=1, type=int, help="Minimum motif size in base pairs.")
    group.add_argument("-max", "--max-motif-size", default=50, type=int, help="Maximum motif size in base pairs.")
    group.add_argument("--min-repeats", default=3, type=int, help="The minimum number of repeats to look for.")
    group.add_argument("--min-span", default=9, type=int, help="The repeats should span at least this many consecutive bases in the input sequence.")
    parser.add_argument("-i", "--interval", help="Only consider sequence from this interval (chrom:start_0based-end).")
    parser.add_argument("-p", "--plot", help="Write out a plot with this filename (only for direct sequence mode).")
    parser.add_argument("-o", "--output-prefix", help="Output filename prefix. For FASTA input, a BED will be generated.")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("--debug", action="store_true", help="Print debugging output.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show progress bar.")
    parser.add_argument("-t", "--threads", type=int, default=cpu_count(), help="Number of parallel processes for FASTA mode.")
    parser.add_argument("input_sequence", help="The nucleotide sequence, or a FASTA file path")

    args = parser.parse_args()

    # validations
    if args.min_motif_size < 1:
        parser.error(f"--min-motif-size is set to {args.min_motif_size}. It must be at least 1.")
    if args.max_motif_size < args.min_motif_size:
        parser.error(f"--max-motif-size is set to {args.max_motif_size}. It must be at least --min-motif-size.")
    if args.min_repeats < 1:
        parser.error(f"--min-repeats is set to {args.min_repeats}. It must be at least 1.")
    if args.min_span < 1:
        parser.error(f"--min-span is set to {args.min_span}. It must be at least 1.")

    interval_sequence = None

    # FASTA path input → parallel by contig
    if os.path.isfile(args.input_sequence):
        fasta_path = args.input_sequence
        if not args.output_prefix:
            args.output_prefix = re.sub(r"\.(fa|fasta)(\.gz)?$", "", os.path.basename(fasta_path), flags=re.IGNORECASE)
        output_bed_path = f"{os.path.basename(args.output_prefix)}.bed"

        fa = pyfastx.Fasta(fasta_path)
        contigs_in_order = list(fa.keys())

        # interval parsing (restrict to a single contig slice)
        interval_spec = None
        if args.interval:
            parts = re.split("[:-]", args.interval)
            if len(parts) != 3:
                parser.error("Invalid --interval format. Must be chrom:start_0based-end")
            chrom_i, s_i, e_i = parts
            try:
                s_i = int(s_i); e_i = int(e_i)
            except Exception:
                parser.error("Invalid --interval numbers. Must be integers: chrom:start_0based-end")
            if chrom_i not in fa:
                parser.error(f"Chromosome {chrom_i} not found in the input FASTA")
            # clamp
            clen = len(fa[chrom_i])
            s_i = max(0, s_i); e_i = min(clen, e_i)
            if s_i >= e_i:
                parser.error("Invalid --interval range: start must be < end within the contig length")
            contigs = [chrom_i]
            interval_spec = (s_i, e_i)
        else:
            contigs = contigs_in_order

        # prepare temp dir
        tmp_dir = tempfile.mkdtemp(prefix="prf_tmp_")

        # build tasks
        tasks = []
        for chrom in contigs:
            tasks.append({
                "fasta_path": fasta_path,
                "chrom": chrom,
                "min_motif_size": args.min_motif_size,
                "max_motif_size": args.max_motif_size,
                "min_repeats": args.min_repeats,
                "min_span": args.min_span,
                "show_progress": False,
                "interval": interval_spec if args.interval else None,
                "tmp_dir": tmp_dir,
            })

        total_jobs = len(tasks)
        print(f"Processing {total_jobs} contig(s) with {min(args.threads, total_jobs)} worker(s) ...")

        # run in parallel
        results = {}
        with ProcessPoolExecutor(max_workers=max(1, min(args.threads, total_jobs))) as ex:
            futures = {ex.submit(_process_one_contig, t): t["chrom"] for t in tasks}
            for fut in tqdm.tqdm(as_completed(futures), total=total_jobs, unit="contig", disable=(not args.show_progress_bar)):
                chrom = futures[fut]
                try:
                    res = fut.result()
                    results[chrom] = res
                    print(f"[{chrom}] repeats: {res['count']:,d} (len {res['announced_len']:,d} bp)")
                except Exception as e:
                    print(f"[{chrom}] ERROR: {e}")
                    raise

        # write final BED in original contig order (or single interval contig)
        with open(output_bed_path, "wt") as out:
            if args.interval:
                chrom = contigs[0]
                if chrom in results:
                    with open(results[chrom]["tmp_path"], "rt") as fin:
                        shutil.copyfileobj(fin, out)
            else:
                for chrom in contigs_in_order:
                    if chrom in results:
                        with open(results[chrom]["tmp_path"], "rt") as fin:
                            shutil.copyfileobj(fin, out)

        # cleanup
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print(f"Wrote results to {output_bed_path}")

    # direct sequence string (ACGTN) → single process
    elif set(args.input_sequence.upper()) <= set("ACGTN"):
        if args.interval:
            parser.error("The --interval option is only supported for FASTA files.")
        interval_sequence = args.input_sequence
        if not args.output_prefix:
            args.output_prefix = "repeats"
        output_tsv_path = f"{args.output_prefix}.tsv"

        per_filters = SimpleNamespace(
            min_motif_size=args.min_motif_size,
            max_motif_size=args.max_motif_size,
            min_repeats=args.min_repeats,
            min_span=args.min_span,
            interval_start_0based=0,
            interval_end=len(interval_sequence),
        )

        output_intervals = detect_repeats(interval_sequence, per_filters,
                                          verbose=args.verbose,
                                          show_progress_bar=args.show_progress_bar,
                                          debug=args.debug)
        print(f"Found {len(output_intervals):,d} repeats")
        with open(output_tsv_path, "wt") as tsv:
            tsv.write("\t".join(["start_0based", "end", "motif"]) + "\n")
            for s, e, m in output_intervals:
                tsv.write(f"{s}\t{e}\t{m}\n")
        print(f"Wrote results to {output_tsv_path}")

    else:
        parser.error(f"Invalid input: {args.input_sequence}. This should be a FASTA file path or a string of nucleotides.")

    # plotting (only for short direct sequences)
    if args.plot and interval_sequence:
        if len(interval_sequence) > 5_000:
            print(f"Warning: The input sequence is too long ({len(interval_sequence):,d} bp). Skipping plot...")
        else:
            plot_results(interval_sequence, output_intervals, args.max_motif_size, args.plot)


if __name__ == "__main__":
    main()