#20211008 by sm
#python3 01.split_digit_usingALLgen.py [intput.gen] [output Dir]
# input.gen : after HLA imputation using HanREF

import os,glob,sys



def split_digit(datain_path,dataoutDir):
    datain = open(datain_path,'r')
    two_digit = open(dataoutDir + "2digit.gen","w")
    four_digit = open(dataoutDir+ "4digit.gen","w")
    while 1:
        tmp = datain.readline()
        if not tmp:
            break
        tmp2 = tmp.split()
        chrom,ID,A,P = tmp2[0].split(":")
        HLA,gene,typing = ID.split("_")
        if len(typing) < 3:
            two_digit.write(tmp)
        else:
            four_digit.write(tmp)
    two_digit.close()
    four_digit.close()

def main():
    datain = sys.argv[1]
    dataout = sys.argv[2]
    if len(datain) == 0 or len(dataout) == 0:
        print("python split_digit_usingALL.py [input.gen] [outputDir]")
        return
    print("Input data : " + datain)
    print("Output Dir = "+dataout)
    split_digit(datain, dataout)
    print("Done, do next step")

main()


