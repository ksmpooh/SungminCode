#!usr/bin/perl

$geno_file = pop(@ARGV);
$ref_file = pop(@ARGV);

if(!-e $ref_file || !-e $geno_file){
	die "Please check your file name $geno_file, $ref_file\n";
}
if($ref_file !~ /bim/ || $geno_file !~ /bim/){
	die "Please check the extension of your files. Input files should be .bim format\n";
}

system("cp $geno_file $geno_file.org");

open INPUT, "<$geno_file.org";
open OUT, ">$geno_file";
open REF, "<$ref_file";

%REF_hash = ();
while (my $line = <REF>){
	chomp($line);
	my @T = split /\t/, $line;
	my $id = $T[0]."\t".$T[3];
	$REF_hash{$id} = $T[1];
}

while (my $line = <INPUT>){
	chomp($line);
	my @T = split /\t/, $line;
	my $id = $T[0]."\t".$T[3];
	if(exists $REF_hash{$id}){
		$line =~ s/$T[1]/$REF_hash{$id}/;
		print OUT $line."\n";
	}
	else {
		print OUT $line."\n";
	}
}
