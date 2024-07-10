
'''
'ref.txt'
STR	MOTIFS	chrom	start	end
chr1_31555_31570	AAAAT	chr1	31555	31570
chr1_44835_44867	AAAT	chr1	44835	44867
chr1_50481_50513	GT	chr1	50481	50513
chr1_86229_86241	AATG	chr1	86229	86241


'rmsk.txt.gz'
585	463	13	6	17	chr1	10000	10468	-248945954	+	(TAACCC)n	Simple_repeat	Simple_repeat	1471	0	1
585	3612	114	215	13	chr1	10468	11447	-248944975	-	TAR1	Satellite	telo	-399	1712	483	2
585	484	251	132	0	chr1	11504	11675	-248944747	-	L1MC5a	LINE	L1	-2382	395	199	3

위 두개 파일을 비교하여 output 파일을 만들예정 (STR, MOTIF, geonName,genoStart,genoEed,repName,'repClass', 'repFamily')
아래 두가지 조건을 만족하는 것 뽑기.
1. ref.txt의 STR과 resk.txt.gz의 새로운 아이디를 만들어서 (geonName_genoStart_genoEed) 매치가 되는것
2. 위에서 안된것 중에서, chrom이 같고, start, end가 +-10 정도 차이나는것도 추출

'''
import gzip

# 파일 경로 설정
ref_file_path = '/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/STR.type.pass.txt'
rmsk_file_path = '/Users/ksmpooh/Desktop/KU/@research/STR/db/rmsk.txt.gz'
output_file_path = '/Users/ksmpooh/Desktop/KU/@research/STR/02.compare/STR_type/strtype.check.withrmsk.output.txt'

# ref.txt 파일 읽기
ref_data = []
with open(ref_file_path, 'r') as ref_file:
    next(ref_file)  # 헤더 스킵
    for line in ref_file:
        fields = line.strip().split('\t')
        ref_entry = {
            'STR': fields[0],
            'MOTIFS': fields[1],
            'chrom': fields[2],
            'start': int(fields[3]),
            'end': int(fields[4])
        }
        ref_data.append(ref_entry)

# output 파일 열기
with open(output_file_path, 'w') as output_file:
    # 헤더 쓰기
    output_file.write('STR\tMOTIFS\tgenoName\tgenoStart\tgenoEnd\trepName\trepClass\trepFamily\n')

    # rmsk.txt.gz 파일을 한 줄씩 읽어가며 처리
    with gzip.open(rmsk_file_path, 'rt') as rmsk_file:
        for line in rmsk_file:
            fields = line.strip().split('\t')
            rmsk_entry = {
                'genoName': fields[5],
                'genoStart': int(fields[6]),
                'genoEnd': int(fields[7]),
                'repName': fields[10],
                'repClass': fields[11],
                'repFamily': fields[12]
            }

            # 1. STR과 새로운 아이디가 일치하는 경우
            rmsk_id = f"{rmsk_entry['genoName']}_{rmsk_entry['genoStart']}_{rmsk_entry['genoEnd']}"
            matched = False
            for ref in ref_data:
                ref_id = f"{ref['chrom']}_{ref['start']}_{ref['end']}"
                if ref_id == rmsk_id:
                    output_file.write(f"{ref['STR']}\t{ref['MOTIFS']}\t{rmsk_entry['genoName']}\t{rmsk_entry['genoStart']}\t{rmsk_entry['genoEnd']}\t{rmsk_entry['repName']}\t{rmsk_entry['repClass']}\t{rmsk_entry['repFamily']}\n")
                    matched = True
                    break
            
            if not matched:
                # 2. chrom이 같고, start, end가 +-10 정도 차이나는 경우
                for ref in ref_data:
                    if ref['chrom'] == rmsk_entry['genoName'] and \
                       abs(ref['start'] - rmsk_entry['genoStart']) <= 10 and \
                       abs(ref['end'] - rmsk_entry['genoEnd']) <= 10:
                        output_file.write(f"{ref['STR']}\t{ref['MOTIFS']}\t{rmsk_entry['genoName']}\t{rmsk_entry['genoStart']}\t{rmsk_entry['genoEnd']}\t{rmsk_entry['repName']}\t{rmsk_entry['repClass']}\t{rmsk_entry['repFamily']}\n")
                        break