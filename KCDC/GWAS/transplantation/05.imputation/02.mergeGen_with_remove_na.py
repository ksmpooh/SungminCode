######
#tool -no_amf_align -buffer 1000 -int 1 12345 -h reference.hap.gz -l reference.legend.gz -m genetic_map -g phasing.haps.gz -o_gz -o chr1_1_12345

import glob,os
#05.imputation/OUTPUTs/02.merge

Dir = "/ADATA/smkim/JG/05.imputation/"
#Dir = "/ADATA/smkim/Gastric/Imputation/"
inDir = Dir + "INPUTs/"
outDir = Dir + "OUTPUTs/"

shDir = Dir + "SCRIPTs/02.mergeGen/"

phasingDir = "/ADATA/smkim/JG/04.phasing/OUTPUTs/03.5Ksplit/"
refDir = inDir + "KGP3KRG/"
mapDir = inDir + "map/"

tool ="/ADATA/smkim/JG/TOOLs/qctool"

def mk_sh(samples,chr,front,tail,ifront,itail):
#       os.system("mkdir "+shDir+chr)
        with open(shDir+"JG.imputed.%s.%s_%s.sampleMerge.sh"%(chr,ifront,itail),'w') as shwrite:
                shwrite.write("%s "%(tool))
                for sample in samples:
                        sample_input = phasingDir+"JG.phasing.%s.%s.sample"%(chr,sample)
                        gen_input = outDir + "01.imputation/JG.imputation.%s.%s.%s_%s.gen.gz"%(chr,sample,ifront,itail)
                        shwrite.write(" -g %s -s %s"%(gen_input,sample_input))
                output = outDir + "02.mergeGen/JG.imputation.mergeGen.%s.%s_%s.gen.gz"%(chr,ifront,itail)
                outsample = outDir + "02.mergeGen/JG.imputation.mergeGen.%s.%s_%s.sample"%(chr,ifront,itail)
                shwrite.write(" -og %s -os %s\n"%(output,outsample))
		gen = output
		out = output.replace(".mergeGen.",".mergeGen.processing.")
		shwrite.write("zcat %s |cut -d' ' -f 2- |gzip  -c > %s"%(gen,out))
#Gastric.phasing.chr9.5001_10000.sampleID.sample

def main():
        print("main : ...")
        samples = glob.glob("/ADATA/smkim/JG/04.phasing/INPUTs/5KsplitSample/*.sampleID")
        sample_index = [s.replace("/ADATA/smkim/JG/04.phasing/INPUTs/5KsplitSample/","") for s in samples]
#       sample_index = samples.replace("/DATA/smkim/Gastric/Phasing/INPUTs/5KsplitSample/","")

        ref_file = open(inDir + "imputation.POS.auto.txt","r")
        ref_list = [x.replace('\n','').split('\t') for x in ref_file]

        for ref_pos,imputation_pos in ref_list[1:]:
                chr,front,tail = ref_pos.split('_')
                ichr,ifront,itail = imputation_pos.split('_')
                mk_sh(sample_index,chr,front,tail,ifront,itail)

main()


