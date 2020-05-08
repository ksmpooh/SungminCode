#for i in $(seq 1 22);do
#	cat *chr$i.* | grep -v ^# | grep -v alternate_ids |  awk '{print $4"\t"$5"\t"$6}'| sort|uniq -d > chr$i.duplicated.txt
#	done



epacts='/backup/smkim/JG/08.asso/OUTPUTs/01.b.firth'
merge='/backup/smkim/JG/08.asso/OUTPUTs/01_1.merge'



for i in $(seq 1 22);do
	echo "chr[$i]...."
	cp $merge/header.txt $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt
	zcat $epacts/*chr$i.*.gz | grep -v ^# >> $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt
	gzip -c $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt > $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt.gz
	done




epacts='/epacts/Result/'
merge='/epacts/Merge/'


for i in $(seq 1 22);do
        echo "chr[$i]...."
        cp $merge/header.txt $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt
        zcat $epacts/*chr$i.*.gz | grep -v ^# >> $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt
        gzip -c $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt > $merge/chr$i.liver.organfailure.epacts.b.firth.result.txt.gz
        done



