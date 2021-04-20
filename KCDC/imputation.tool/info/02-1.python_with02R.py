import os,glob

def main():
    dfs = glob.glob("*txt")
    outDir = "intersect/"
    os.system("mkdir "+outDir)
    for df in dfs:
        os.system("Rscript --vanilla intersect.R intersect.ref %s %s"%(df,outDir+df))


main()