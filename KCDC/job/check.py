#chunk check

chunks = ["101000001_106000000","126000001_131000000","171000001_171053406","20000001_25000000","25000001_30000000","35000001_40000000","45000001_50000000","55000001_58780059","66000001_71000000","71000001_76000000","91000001_96000000"]

import os,glob
def main():
    os.system("mkdir re")
    for chunk in chunks:
        df = glob.glob("*chr6*%s*sh"%chunk)
        for i in df:
            os.system("cp %s re/"%i)

main()

chunks = ["101000001_106000000","126000001_131000000","171000001_171053406","20000001_25000000","25000001_30000000","35000001_40000000","45000001_50000000","55000001_58780059","66000001_71000000","71000001_76000000","91000001_96000000"]
import os,glob
def main():
    df = glob.glob("*chr6*.sh")
    for chunk in chunks:
        ref = glob.glob("*chr6*%s%.sh"%chunk)
        

    #os.system("mkdir re")
    for chunk in chunks:
        df = glob.glob("*chr6*%s*"%chunk)
        for i in df:
            os.system("cp %s re/"%i)

main()