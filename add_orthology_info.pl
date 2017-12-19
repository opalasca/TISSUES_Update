
#################################################################################################################################
##               Generate file containing scored gene-tissue pairs and 1:1 orthologs across the other three organisms, if they exist. Needed in computing the correlation between datasets across different organisms.                   #
##                       Comprehensive comparison of large-scale tissue expression datasets                                      #
##       Alberto Santos, Kalliopi Tsafou, Christian Stolte, Sune Pletscher-Frankild, Se??n I. O???Donoghue and Lars Juhl Jensen  #
##################################################################################################################################

#!/usr/bin/ perl

use strict;
use POSIX;

our $org1_id = "9606";
our $org2_id = "10090";
our $org3_id = "10116";
our $org4_id = "9823";

my $data1_human = "gnf";
my $data2_human = "exon";
my $data3_human = "hpa_rna";
my $data4_human = "rna";
my $data5_human = "human_text_mining";

my $data1_mouse = "mouse_gnf";
my $data2_mouse = "mouse_gnfv3";
my $data3_mouse = "mouse_rnaseq_encode";
my $data4_mouse = "mouse_rnaseq_mit";
my $data5_mouse = "mouse_text_mining";

my $data1_rat = "rat_array";
my $data2_rat = "rat_rnaseq_mit";
my $data3_rat = "rat_rnaseq_bodymap";
my $data4_rat = "rat_text_mining";

my $data1_pig = "pig_array";
my $data2_pig = "pig_rnaseq_aarhus";
my $data3_pig = "pig_rnaseq_wur";
my $data4_pig = "pig_text_mining";

my %human_data = ($data1_human=> 1,$data2_human=> 1,$data3_human=> 1, $data4_human=> 1,$data5_human=> 1);
my %mouse_data = ($data1_mouse => 1,$data2_mouse => 1,$data3_mouse => 1,$data3_mouse_bis => 1, $data4_mouse => 1,$data5_mouse => 1);
my %rat_data = ($data1_rat => 1,$data2_rat => 1,$data3_rat => 1,$data4_rat => 1);
my %pig_data = ($data1_pig => 1,$data2_pig => 1,$data3_pig => 1,$data4_pig => 1);

our %org1_org2 = ();
our %org2_org1 = ();
our %org1_org3 = ();
our %org3_org1 = ();
our %org1_org4 = ();
our %org4_org1 = ();
our %org2_org3 = ();
our %org3_org2 = ();
our %org2_org4 = ();
our %org4_org2 = ();
our %org3_org4 = ();
our %org4_org3 = ();

my @organisms = ($org1_id, $org2_id, $org3_id, $org4_id);	       


open IN, "< data/orthology/maNOG.members.tsv";
open INR, "< data/orthology/roNOG.members.tsv";

while (<IN>) {
        s/\r?\n//;
        my @fields = split("\t");
	my @items = split(/,/,$fields[5]);
	my %found;
	foreach my $i (@items){
		my ($org, $id) = split(/\./, $i);
		if ($org ~~ @organisms){
			$found{$org}{$id}=1;
		}
	}
	my @org1_orth = keys %{$found{$org1_id}};
	my @org2_orth = keys %{$found{$org2_id}};
	my @org3_orth = keys %{$found{$org3_id}};
	my @org4_orth = keys %{$found{$org4_id}};
		
	if (scalar(@org1_orth)==1 and scalar(@org2_orth)==1){ $org1_org2{$org1_orth[0]}=$org2_orth[0];$org2_org1{$org2_orth[0]}=$org1_orth[0];}
	if (scalar(@org1_orth)==1 and scalar(@org3_orth)==1){ $org1_org3{$org1_orth[0]}=$org3_orth[0];$org3_org1{$org3_orth[0]}=$org1_orth[0];}
	if (scalar(@org1_orth)==1 and scalar(@org4_orth)==1){ $org1_org4{$org1_orth[0]}=$org4_orth[0];$org4_org1{$org4_orth[0]}=$org1_orth[0];}
	if (scalar(@org2_orth)==1 and scalar(@org4_orth)==1){ $org2_org4{$org2_orth[0]}=$org4_orth[0];$org4_org2{$org4_orth[0]}=$org2_orth[0];}
	if (scalar(@org3_orth)==1 and scalar(@org4_orth)==1){ $org3_org4{$org3_orth[0]}=$org4_orth[0];$org4_org3{$org4_orth[0]}=$org3_orth[0];}
}
close IN;

while (<INR>) {
        s/\r?\n//;
        my @fields = split("\t");
	my @items = split(/,/,$fields[5]);
	my %found;
	foreach my $i (@items){
		my ($org, $id) = split(/\./, $i);
		if ($org ~~ @organisms){
			$found{$org}{$id}=1;
		}
	}
	my @org2_orth = keys %{$found{$org2_id}};
	my @org3_orth = keys %{$found{$org3_id}};	
	if (scalar(@org2_orth)==1 and scalar(@org3_orth)==1){ $org2_org3{$org2_orth[0]}=$org3_orth[0];$org3_org2{$org3_orth[0]}=$org2_orth[0];}
		
}
close INR;



open IN, "< data/pairs_major_tissues.tsv";
open OUT, "> data/pairs_major_tissues_orthologs.tsv";
open OUTF, "> data/pairs_major_tissues_orthologs_filtered.tsv";


while (<IN>) {
	s/\r?\n//;
	my ($dataset,$score,$stars,$tissue,$ens) = split("\t");
	if (defined $human_data{$dataset}){
		print OUT "$dataset\t$score\t$stars\t$tissue\t$ens\t$ens\t",&hval($ens,\%org1_org2),"\t",&hval($ens,\%org1_org3),"\t",&hval($ens,\%org1_org4),"\n";
		if ($stars>0){
			print OUTF "$dataset\t$score\t$stars\t$tissue\t$ens\t$ens\t",&hval($ens,\%org1_org2),"\t",&hval($ens,\%org1_org3),"\t",&hval($ens,\%org1_org4),"\n";
		}
	}
	if (defined $mouse_data{$dataset}){
		print OUT "$dataset\t$score\t$stars\t$tissue\t$ens\t",&hval($ens,\%org2_org1),"\t$ens\t",&hval($ens,\%org2_org3),"\t",&hval($ens,\%org2_org4),"\n";
		if ($stars>0){
			print OUTF "$dataset\t$score\t$stars\t$tissue\t$ens\t",&hval($ens,\%org2_org1),"\t$ens\t",&hval($ens,\%org2_org3),"\t",&hval($ens,\%org2_org4),"\n";
		}
	}
	if (defined $rat_data{$dataset}){
		print OUT "$dataset\t$score\t$stars\t$tissue\t$ens\t",&hval($ens,\%org3_org1),"\t",&hval($ens,\%org3_org2),"\t$ens\t",&hval($ens,\%org3_org4),"\n";
		if ($stars>0){
			print OUTF "$dataset\t$score\t$stars\t$tissue\t$ens\t",&hval($ens,\%org3_org1),"\t",&hval($ens,\%org3_org2),"\t$ens\t",&hval($ens,\%org3_org4),"\n";
		}
	}
	if (defined $pig_data{$dataset}){
		print OUT "$dataset\t$score\t$stars\t$tissue\t$ens\t",&hval($ens,\%org4_org1),"\t",&hval($ens,\%org4_org2),"\t",&hval($ens,\%org4_org3),"\t$ens\n";
		if ($stars>0){
			print OUTF "$dataset\t$score\t$stars\t$tissue\t$ens\t",&hval($ens,\%org4_org1),"\t",&hval($ens,\%org4_org2),"\t",&hval($ens,\%org4_org3),"\t$ens\n";
		}
	}
}	

close IN;
close OUT;
close OUTF;


sub hval(){
	my $ens = $_[0];
	my $orth_hash = $_[1];
	my $value;
	if (defined ${$orth_hash}{$ens}){
		$value = ${$orth_hash}{$ens};	
	}
	else {
		$value = "na";
	}
	return $value;
}

close STDERR;
close STDOUT;
POSIX::_exit(0);
