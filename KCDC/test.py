
out = open("merge.list.txt","w")
for i in range(2,22+1):
    tmp = "JG.KR.rep.imputation.chr%s_MAF0.01_INFO0.8.filter"%(str(i))
    out.write("%s.bed\t%s.bim\t%s.fam\n"%(tmp,tmp,tmp))

out.close()

    
#JG.onlyKR.NODAT.imputation.filter.ALLchr.nosex

#plink --bfile JG.KR.rep.imputation.chr1_MAF0.01_INFO0.8.filter --merge-list merge.list.txt --allow-no-sex --make-bed --out JG.onlyKR.rep.imputation.filter.ALLchr

#awk '{print $1}' JG.onlyKR.rep.imputation.filter.ALLchr.bim |sort |uniq -c |wc -l