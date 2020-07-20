open REF, "<HAN.MHC.reference.panel.bgl";
open OUT, ">HAN.MHC.reference.panel.fixed.bgl.phased";

$headoff = <REF>;
chomp($headoff);
my @head = split /\s+/, $headoff;
$phead = $headoff;
$phead =~ s/I id/P pedigree/;
my $zeros = "";
my $genders = "";
for(my $i=2; $i<=$#head; $i++){
	$zeros .= " 0";
	$genders .= " 2";
}

print OUT $phead."\n";
print OUT $headoff."\n";
print OUT "fID father".$zeros."\n";
print OUT "mID mother".$zeros."\n";
print OUT "C gender".$genders."\n";

while (my $line = <REF>){
	if($line !~ /\+/ && $line !~ /\-/){
		print OUT $line;
	}
	else {
		chomp($line);
		my @T = split /\s+/, $line;
		my %allele = ();

		for(my $i=2; $i<=$#T; $i++){
			$allele{$T[$i]}++;
		}

		my @AL = keys(%allele);
		my $new_allele = "";
		if($AL[0] =~ /[\+\-]/){
			$new_allele = $AL[0];
			$new_allele =~ s/[\+\-]/$AL[1]/;
			$AL[0] =~ s/([\+\-])/\\$1/;
			$line =~ s/$AL[0]/$new_allele/g;
		}
		elsif($AL[1] =~ /[\+\-]/){
			$new_allele = $AL[1];
			$new_allele =~ s/[\+\-]/$AL[0]/;
			$AL[1] =~ s/([\+\-])/\\$1/;
			$line =~ s/$AL[1]/$new_allele/g;
		}
		print OUT $line."\n";
	}
}
