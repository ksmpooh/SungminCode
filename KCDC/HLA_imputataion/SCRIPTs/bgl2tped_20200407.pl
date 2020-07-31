open REF, "<HAN.MHC.reference.panel.fixed.bgl.phased";
open TPED, ">HAN.MHC.reference.panel.fixed.tped";
open TFAM, ">HAN.MHC.reference.panel.fixed.tfam";
open POS, "<HAN.MHC.reference.panel.fixed.markers";

%POS_hash = ();
while (my $line = <POS>){
	chomp($line);
	my @T = split /\s+/, $line;
	$POS_hash{$T[0]} = $T[1];
}

$rmheadoff = <REF>;
$headoff = <REF>;
$rmheadoff = <REF>;
$rmheadoff = <REF>;
$rmheadoff = <REF>;
chomp($headoff);

my @ID = split /\s+/, $headoff;

for(my $i=2; $i<=$#ID; $i = $i + 2){
	print TFAM "$ID[$i]\t$ID[$i]\t0\t0\t0\t0\n";
}

while (my $line = <REF>){
	chomp($line);
	my @T = split /\s+/, $line;
	my $ID = $T[1];

	$line =~ s/$T[0] $T[1]/6 $ID 0 $POS_hash{$T[1]}/;
	print TPED $line."\n";
}
