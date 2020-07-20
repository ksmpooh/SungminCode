open INPUT, "<HAN.MHC.reference.panel.markers";
open OUT, ">HAN.MHC.reference.panel.fixed.markers";

while (my $line = <INPUT>){
	if($line !~ /\+/ && $line !~ /\-/){
		print OUT $line;
	}
	else {
		chomp($line);
		my @T = split /\s+/, $line;
		if($T[2] =~ /\+/ || $T[2] =~ /\-/){
			$T[2] =~ s/[\+\-]/$T[3]/;
		}
		elsif($T[3] =~ /\+/ || $T[3] =~ /\-/){
			$T[3] =~ s/[\+\-]/$T[2]/;
		}
		print OUT "$T[0] $T[1] $T[2] $T[3]\n";
	}
}
