#################################################################################################################################
#		This code generates the files necessary to reproduce the figures shown in the manuscript			#
#				TISSUES: An integrative web resource on mammalian tissue expression                             #
#			   Oana Palasca, Alberto Santos, Christian Stolte, Jan Gorodkin, Lars Juhl Jensen                       #
#################################################################################################################################
	
#!/usr/bin/ perl
use warnings;
use strict;
use Data::Dumper;
#Data directory
my $data_dir = "data/";

#Input files: these files need to be downloaded from FigShare (https://figshare.com/s/84433e3f30d7b120fc24)
# and stored in the data/datasets/ directory

my $data1_human = "gnf";
my $data2_human = "exon";
my $data3_human = "hpa_rna";
my $data4_human = "rna";
my $data5_human = "tm";
my $data6_human = "human_text_mining";

my $data1_mouse = "mouse_gnf";
my $data2_mouse = "mouse_gnfv3";
my $data3_mouse = "mouse_rnaseq_encode";
my $data3_mouse_bis = "mouse_rnaseq_processed_rnaseq_encode_tpm";
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

my %options = ($data1_human => $data_dir."datasets/".$data1_human.".tsv",
		$data2_human => $data_dir."datasets/".$data2_human.".tsv",
		$data3_human => $data_dir."datasets/".$data3_human.".tsv",
		$data4_human => $data_dir."datasets/rna_seq.tsv",
		$data5_human => $data_dir."datasets/text_mining.tsv",
		$data6_human => $data_dir."datasets/".$data6_human.".tsv",
		$data1_mouse => $data_dir."datasets/".$data1_mouse.".tsv",
		$data2_mouse => $data_dir."datasets/".$data2_mouse.".tsv",
		$data3_mouse => $data_dir."datasets/".$data3_mouse.".tsv",
		$data3_mouse_bis => $data_dir."datasets/".$data3_mouse_bis.".tsv",
		$data4_mouse => $data_dir."datasets/".$data4_mouse.".tsv",
		$data5_mouse => $data_dir."datasets/".$data5_mouse.".tsv",
		$data1_rat => $data_dir."datasets/".$data1_rat.".tsv",
		$data2_rat => $data_dir."datasets/".$data2_rat.".tsv",
		$data3_rat => $data_dir."datasets/".$data3_rat.".tsv",
		$data4_rat => $data_dir."datasets/".$data4_rat.".tsv",
		$data1_pig => $data_dir."datasets/".$data1_pig.".tsv",
		$data2_pig => $data_dir."datasets/".$data2_pig.".tsv",
		$data3_pig => $data_dir."datasets/".$data3_pig.".tsv",
		$data4_pig => $data_dir."datasets/".$data4_pig.".tsv"
	);

my %human = ($data1_human=> 1,$data2_human=> 1,$data3_human=> 1, $data4_human=> 1,$data5_human=> 1,$data6_human=> 1);
my %mouse = ($data1_mouse => 1,$data2_mouse => 1,$data3_mouse => 1,$data3_mouse_bis => 1,$data4_mouse => 1,$data5_mouse => 1);
my %rat = ($data1_rat => 1,$data2_rat => 1,$data3_rat => 1,$data4_rat => 1);
my %pig = ($data1_pig => 1,$data2_pig => 1,$data3_pig => 1,$data4_pig => 1 );

#dictionary files: These files are used for backtracking tissues to their respective parent tissues
# to identify which genes are expressed in the tissues of interest (more details in the README file)
my $labels_file = $data_dir."dictionary/labels.tsv";
my $bto_entities = $data_dir."dictionary/tissues_entities.tsv";
my $bto_groups = $data_dir."dictionary/tissues_groups.tsv";

#common tissues: Lists of common tissues considered between transcriptomic datasets and the gold standard
my %common_tissues_1 =("heart"=>1, "liver"=>1, "intestine"=>1, "nervous system"=>1, "kidney"=>1);
my %common_tissues_2 =("heart"=>1, "liver"=>1, "muscle"=>1, "nervous system"=>1, "kidney"=>1);
my %common_tissues_3 =("heart"=>1, "muscle"=>1, "intestine"=>1, "nervous system"=>1, "kidney"=>1);
my %common_tissues_4 =("liver"=>1, "muscle"=>1, "lung"=>1, "nervous system"=>1, "spleen"=>1);

my %common_tissues_hash = ($data1_human => \%common_tissues_1,
		$data2_human => \%common_tissues_1,
		$data3_human => \%common_tissues_1,
		$data4_human => \%common_tissues_1,
		$data5_human => \%common_tissues_1,
		$data6_human => \%common_tissues_1,
		$data1_mouse => \%common_tissues_1,
		$data2_mouse => \%common_tissues_1,
		$data3_mouse => \%common_tissues_1,
		$data3_mouse_bis => \%common_tissues_1,
		$data4_mouse => \%common_tissues_1,
		$data5_mouse => \%common_tissues_1,
		$data1_rat => \%common_tissues_3,
		$data2_rat => \%common_tissues_1,
		$data3_rat => \%common_tissues_2,
		$data4_rat => \%common_tissues_1,
		$data1_pig => \%common_tissues_1,
		$data2_pig => \%common_tissues_2,
		$data3_pig => \%common_tissues_4,
		$data4_pig => \%common_tissues_1
	);

#Get the dataset name we are going to analyze   
my $dataset_file = "";

#Gold standards
my %goldstandards = ("human" => {"uniprot" => $data_dir."datasets/human_knowledge.tsv",
				"uniprot_orth" => $data_dir."datasets/uniprot_human_reduced_orthology_10116.tsv"},
		"mouse" => {"uniprot" => $data_dir."datasets/mouse_knowledge.tsv",
				"uniprot_orth" => $data_dir."datasets/uniprot_orthology_10090.tsv"},
		"rat" => {"uniprot" => $data_dir."datasets/rat_knowledge.tsv",
				"uniprot_orth" => $data_dir."datasets/uniprot_orthology_10116.tsv"},
		"pig" => {"uniprot" => $data_dir."datasets/pig_knowledge.tsv",
				"uniprot_orth" => $data_dir."datasets/uniprot_orthology_9823.tsv"}
	);


#sliding window
my $window_size = 100;
#output files
my $major_tissues_file = $data_dir."datasets_major_tissues.tsv"; #Specifies which datasets have which major tissues within the ones specified in the labels.tsv file
my $pairs_major_tissues_file = "pairs_major_tissues.tsv"; #Gene-tissue associations from all datasets for the 21 major tissues
##############

#  Execution #
##############
my ($labels) = &parse_labels_file();
my ($entities) = &parse_bto_entities_file();
my ($groups) = &parse_bto_groups_file();
my ($child_labels_hash) = &get_child_label($labels,$entities,$groups);

#dictionary with the number of proteins per tissue
my %major_tissues = ();

my %goldstandard_data = ();
my %gold_btos = ();
my %goldstandard_labels = ();
#dictionary with the protein-tissue pairs 
my %data_major_tissues = ();

&parse_goldstandard_file(\%goldstandards,\%goldstandard_data,\%gold_btos);

#print Dumper($child_labels_hash);

foreach my $org (keys %goldstandards){
	#print "$org  ";
	foreach my $gold (keys %{$goldstandards{$org}}){
		#print "$gold  ";
		#print Dumper($goldstandard_data{$org}{$gold});
		my ($goldstd_labels) = &convert_btos_labels(\%{$goldstandard_data{$org}{$gold}}, $child_labels_hash, \%major_tissues, $gold);
		$goldstandard_labels{$org}{$gold}=$goldstd_labels;
		&combine_data_major_tissues($goldstd_labels, \%data_major_tissues, $gold);
	}
}	

#dictionary with all the data available
my %dataset_scored_pairs = ();

foreach my $dataset_name (keys %options){
	$dataset_file = $options{$dataset_name};	
	my ($dataset) = &parse_studied_dataset_file($dataset_name);
	my $dataset_labels = undef;
	my $major_tissues = undef;
	($dataset_labels) = &convert_btos_labels($dataset,$child_labels_hash, \%major_tissues, $dataset_name);
	if (defined $common_tissues_hash{$dataset_name}){
		&format_score_pairs($dataset_labels, $common_tissues_hash{$dataset_name}, \%dataset_scored_pairs, $dataset_name);
	}
	&combine_data_major_tissues($dataset_labels, \%data_major_tissues, $dataset_name);

}

#calculate fold enrichment for the different datasets
foreach my $dataset_name (keys %dataset_scored_pairs){
	if (defined $human{$dataset_name}){	
		foreach my $goldstandard (keys $goldstandards{"human"}){
			my ($goldstandard_filtered) = &filter_common_tissues($goldstandard_labels{"human"}{$goldstandard}, $common_tissues_hash{$dataset_name});
			my $output_file = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis.tsv";
			my $output_file_full = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_full.tsv";	
			my $output_file_roc = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_roc.tsv";	
			&calculate_fold_enrichment(\%{$dataset_scored_pairs{$dataset_name}}, $goldstandard_filtered, $output_file, $window_size, $output_file_full,$output_file_roc);
		}
	}
	if (defined $mouse{$dataset_name}){	
		foreach my $goldstandard (keys $goldstandards{"mouse"}){
			my ($goldstandard_filtered) = &filter_common_tissues($goldstandard_labels{"mouse"}{$goldstandard}, $common_tissues_hash{$dataset_name});
			my $output_file = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis.tsv";
			my $output_file_full = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_full.tsv";		
			my $output_file_roc = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_roc.tsv";
			&calculate_fold_enrichment(\%{$dataset_scored_pairs{$dataset_name}}, $goldstandard_filtered, $output_file, $window_size, $output_file_full,$output_file_roc);
		}
	}
	if (defined $rat{$dataset_name}){	
		foreach my $goldstandard (keys $goldstandards{"rat"}){
			my ($goldstandard_filtered) = &filter_common_tissues($goldstandard_labels{"rat"}{$goldstandard}, $common_tissues_hash{$dataset_name});
			my $output_file = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis.tsv";
			my $output_file_full = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_full.tsv";	
			my $output_file_roc = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_roc.tsv";	
			&calculate_fold_enrichment(\%{$dataset_scored_pairs{$dataset_name}}, $goldstandard_filtered, $output_file, $window_size, $output_file_full,$output_file_roc);
		}
	}
	if (defined $pig{$dataset_name}){	
		foreach my $goldstandard (keys $goldstandards{"pig"}){
			my ($goldstandard_filtered) = &filter_common_tissues($goldstandard_labels{"pig"}{$goldstandard}, $common_tissues_hash{$dataset_name});
			my $output_file = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis.tsv";
			my $output_file_full = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_full.tsv";		
			my $output_file_roc = $data_dir."$dataset_name"."_$goldstandard"."_fold_enrichment_analysis_roc.tsv";
			&calculate_fold_enrichment(\%{$dataset_scored_pairs{$dataset_name}}, $goldstandard_filtered, $output_file, $window_size, $output_file_full,$output_file_roc);
		}
	}
}


#Print out the file with: Dataset Tissue NumberProteins
&print_2dim(\%major_tissues, $major_tissues_file);
#Print out the file with gene-tissue pairs
&print_3dim(\%data_major_tissues, $pairs_major_tissues_file, $data_dir);

##############
#  Functions #
##############
sub print_2dim(){
	my ($data) = $_[0];
	my ($output_file) = $_[1];
	
	open(OUT,">$output_file") or die "Unable to open the Output file $output_file\n";
	foreach my $dataset (keys %{$data}){
		foreach my $tissue (keys %{${$data}{$dataset}}){
			print OUT $dataset."\t".$tissue."\t".${$data}{$dataset}{$tissue}."\n";
		}
	}
	close(OUT);
}

sub print_3dim(){
    my ($data) = $_[0];
    my ($output_file) = $_[1];
    my ($data_dir) = $_[2];
    my $full_path = $data_dir.$output_file;
    open(OUT,">$full_path") or die "Unable to open the Output file $full_path\n";
    foreach my $dataset (keys %{$data}){
        foreach my $id (keys %{${$data}{$dataset}}){
            print OUT $dataset."\t".${$data}{$dataset}{$id}."\t".$id."\n";
        }
    }
    
    close(OUT);
}


sub parse_labels_file(){
    print "- Parsing the Labels file $labels_file with the tissues to be studied\n";
    
    open(LABELS,"$labels_file") or die "Unable to open $labels_file file\n";
    
    my %labels = ();
    my $tissues_type_code = "-25";
    
    while(<LABELS>){
        s/\r?\n//;
        my ($type,$name,$code) = split("\t",$_);
        if($type eq $tissues_type_code){
            $labels{$code} = $name;
        }
    }
    print "- The labels file has been correctly parsed\n";
    return \%labels;
}

sub parse_bto_entities_file(){
    print "- Parsing the BTO entities file $bto_entities\n";
    
    open(ENTITIES,"$bto_entities") or die "Unable to open $bto_entities file\n";
    
    my %entities = ();
    my $tissues_type_code = "-25";
    
    while(<ENTITIES>){
        s/\r?\n//;
        my ($id,$type,$bto) = split("\t",$_);
        if($type eq $tissues_type_code){
            $entities{$id} = $bto;
        }
    }
    print "- The BTO entities file has been correctly parsed\n";
    return \%entities;
}

sub parse_bto_groups_file(){
    print "- Parsing the BTO groups file $bto_groups\n";
    
    open(GROUPS,"$bto_groups") or die "Unable to open $bto_groups file\n";
    
    my %groups = ();
    
    while(<GROUPS>){
        s/\r?\n//;
        my ($id,$parent_id) = split("\t",$_);
        $groups{$id}{$parent_id} = 1;
    }
    print "- The BTO groups file has been correctly parsed\n";
    return \%groups;
}

sub get_child_label(){
    my ($labels) = $_[0];
    my($entities) = $_[1];
    my ($groups) = $_[2];
    
    my %child_labels = ();
    foreach my $id (keys %{$entities}){
        my $bto = ${$entities}{$id};
        foreach my $parent_id (keys %{${$groups}{$id}}){
            my $parent_bto = ${$entities}{$parent_id};
            if(exists ${$labels}{$parent_bto}){
                my $label = ${$labels}{$parent_bto};
                $child_labels{$bto}{$label} = 1;
            }
        }
        if(exists ${$labels}{$bto}){
            $child_labels{$bto}{${$labels}{$bto}} =1;
        }
    }
	return \%child_labels;
}

sub parse_studied_dataset_file(){
    	my $dataset_name = $_[0];
    	print "- Parsing the Dataset file $dataset_file\n";
    	open (DATASET,"$dataset_file") or die "Unable to open the Dataset file $dataset_file\n";
    	my %dataset = ();
    	while(<DATASET>){
		s/\r?\n//;
		my (undef,$ensp,undef,$bto,$txt_aux,$score_scoretype,$stars,undef) = split("\t",$_);
		my $score = $stars;
		if($dataset_name =~ /tm|mining/){
	    		$score = $txt_aux;
	    		$stars = $score_scoretype;
		}
		else{ 
	    		($score,undef) = split('\s',$score_scoretype);
		}
		$dataset{$ensp}{$bto} = "$score\t$stars";         
    	}
    	
    	print "- The Dataset file has been parsed\n";
    
	return \%dataset;
}

sub parse_goldstandard_file(){
	my ($goldstandards) = $_[0];
	my ($goldstandard_data) = $_[1];
	my ($gold_btos) = $_[2];
   
    	foreach my $org (keys %{$goldstandards}){
		foreach my $gold (keys %{${$goldstandards}{$org}}){
			my $file = ${$goldstandards}{$org}{$gold};
			print "- Parsing the goldstandard files $file\n";
			open (IN, "$file") or die "Unable to open the UniProt file $file\n";
			while(<IN>){
				s/\r?\n//;
				my (undef,$ensp,undef,$bto,undef) = split("\t",$_);
        			${$goldstandard_data}{$org}{$gold}{$ensp}{$bto} = 1;
     				${$gold_btos}{$org}{$gold}{$bto} = 1;
			}
    			print "- The gold standard file $file has been parsed\n";
			close IN;
		}	
	}		
}


sub get_label_values(){
    	my ($dataset) = $_[0];
    	my ($labels) = $_[1];
	my ($major_tissues) = $_[2];
	my ($dataset_name) = $_[3];
	#print "- Getting the labeled btos of the UniGene data\n";
    	my %dataset_labels = ();
    	foreach my $ensp (keys %{$dataset}){
		foreach my $bto (keys %{${$dataset}{$ensp}}){
	    		if(exists ${$labels}{$bto}){
				my $label = ${$labels}{$bto};
				my ($score,$star) = split("\t",${$dataset}{$ensp}{$bto});
				if(exists $dataset_labels{$ensp}{$label}){
		    			my ($p_score,$p_stars) = split("\t",$dataset_labels{$ensp}{$label});
		    			$dataset_labels{$ensp}{$label} = ($dataset_labels{$ensp}{$label}, ${$dataset}{$ensp}{$bto})[$p_score<$score];
					
				}
				else{
					$dataset_labels{$ensp}{$label} = ${$dataset}{$ensp}{$bto};
				}
				${$major_tissues}{$dataset_name}{$label}++;
	    		}
		}
    	}
    	print "- The parent btos have been associated to each ENSP\n";
    
	return \%dataset_labels;
}

sub convert_btos_labels(){
    my ($data) = $_[0];
    my ($btos_labels) = $_[1];
    my ($major_tissues) = $_[2];
    my ($dataset_name) = $_[3];
    
    my %ensp_labels = ();
    my %labels = ();
    foreach my $ensp (keys %{$data}){
        foreach my $bto (keys %{${$data}{$ensp}}){
            foreach my $label (keys %{${$btos_labels}{$bto}}){
                my ($score,$star) = split("\t",${$data}{$ensp}{$bto});
                if(exists $ensp_labels{$ensp}{$label}){
                    my ($p_score,$p_stars) = split("\t",$ensp_labels{$ensp}{$label});
                    $ensp_labels{$ensp}{$label} = ($ensp_labels{$ensp}{$label},${$data}{$ensp}{$bto})[$p_score<$score];
                }
                else{
                    $ensp_labels{$ensp}{$label} = ${$data}{$ensp}{$bto};
                }
                $labels{$label}= 1;
		${$major_tissues}{$dataset_name}{$label}++;
            }
        }
    }
    
    return \%ensp_labels;
}

sub format_score_pairs(){
    my ($dataset_labels) = $_[0];
	my ($common_tissues) = $_[1];
	my ($score_pairs) = $_[2];
	my ($dataset_name) = $_[3];

    	foreach my $ensp (keys %{$dataset_labels}){
		foreach my $label (keys %{${$dataset_labels}{$ensp}}){
			if (exists ${$common_tissues}{$label}){
				my ($score,$stars) = split("\t",${$dataset_labels}{$ensp}{$label});
				${$score_pairs}{$dataset_name}{$score}{$ensp}{$label} = $stars;
			}
		}
    	}
}

sub filter_common_tissues(){
	my ($data) = $_[0];
	my ($common_tissues) = $_[1];
	
	my %filtered = ();

	foreach my $ensp (keys %{$data}){
		foreach my $label (keys %{${$data}{$ensp}}){
			if(exists ${$common_tissues}{$label}){
				$filtered{$ensp}{$label} = 1;
			}
		}
	}

	return \%filtered;
}

sub combine_data_major_tissues(){
	my ($dataset_labels) = $_[0];
	my ($data_major_tissues) = $_[1];
	my ($dataset) = $_[2];

	foreach my $ensp (keys %{$dataset_labels}){
		foreach my $label (keys %{${$dataset_labels}{$ensp}}){
			my ($score, $stars) = split("\t",${$dataset_labels}{$ensp}{$label});
			if ($stars){
                #${$data_major_tissues}{"all"}{$dataset}{$label."\t".$ensp} = ${$dataset_labels}{$ensp}{$label};
                ${$data_major_tissues}{$dataset}{$label."\t".$ensp} = ${$dataset_labels}{$ensp}{$label};
            }
		}
	}
}


#Calculate the fold enrichment for each dataset using a moving average window over the score
#This allowed making the datasets comparable
sub calculate_fold_enrichment(){
    my ($dataset) = $_[0];
    my ($gold) = $_[1];
    my ($output_file) = $_[2];
    my ($window_size) = $_[3];
    my ($output_file_full) = $_[4];
    my ($output_file_roc) = $_[5];
    	
	print "- Calculating precision at the level of score for the dataset, $output_file\n";
    	my @pos_neg_array = ();
    	my @scores = ();
    	my @stars = ();
    	my $dataset_total_pairs = 0;
    	my %gold_valid_prots = ();
	foreach my $score (sort {$a <=> $b} keys %{$dataset}){
		foreach my $ensp (keys %{${$dataset}{$score}}){
	    		if(exists ${$gold}{$ensp}){ 
				$gold_valid_prots{$ensp} = 1;
		 		foreach my $label (keys %{${$dataset}{$score}{$ensp}}){
					$dataset_total_pairs++;
					my $star = ${$dataset}{$score}{$ensp}{$label};
					if(exists ${$gold}{$ensp}{$label}){
						push @pos_neg_array,1;
					}
					else{
						push @pos_neg_array,0;
					}
					push @scores,$score;
					push @stars,$star;
				}
	    		}
		}
    }
    print "- The dataset has $dataset_total_pairs pairs in total\n";
    	
	
    #calculate gold_standard pairs and proteins
	my $gold_num_pairs = 0;
 	my %gold_tissues = ();
    	foreach my $ensp (keys %gold_valid_prots){
		foreach my $label (keys %{${$gold}{$ensp}}){
			$gold_tissues{$label} = 1;
	    		$gold_num_pairs++;	
		}
    	}
    	my $gold_num_prots = scalar keys %gold_valid_prots;
    	my $gold_num_tissues = scalar keys %gold_tissues;
    	my $gold_prob = $gold_num_pairs/($gold_num_tissues*$gold_num_prots);
     
	
    	print "- Gold standard number pairs: $gold_num_pairs\n";
    	print "- Gold standard number proteins: $gold_num_prots\n";
    	print "- Gold standard number tissues: $gold_num_tissues\n";
    	print "- Gold standard probability: $gold_prob\n";
    	
    	use List::Util qw(sum);
    	
  open(OUT,">$output_file_roc") or die "Unable to open the Output file $output_file_roc\n";
	
	my $len = scalar @pos_neg_array-1;
  	my $k = $len;
	my $pos = $len-$k;
	my $sum = $pos_neg_array[$k];
	print OUT "$pos\t$sum\t$scores[$k]\n";
	$k = $k - 1;
	$pos = $pos + 1;

	while($k>=0){
		$sum = $sum + $pos_neg_array[$k];
		print OUT "$pos\t$sum\t$scores[$k]\n";
		$k = $k-1;
		$pos = $pos + 1;
	}
	close(OUT);

    
	open(OUT,">$output_file_full") or die "Unable to open the Output file $output_file_full\n";
	
	my $i = 0;
	my $j = $window_size -1;
	my $last = 0;

	$sum = sum(@pos_neg_array[$i..$j]);
	my $fold_enrichment = ($sum/$window_size)/$gold_prob;
	my $mean_stars = sum(@stars[$i..$j])/$window_size;
	my $mean_score = sum(@scores[$i..$j])/$window_size;

	print OUT "$mean_stars\t$fold_enrichment\t$mean_score\n";
	$i = $i + 1;
	if ($j + 1>(scalar @pos_neg_array-1)){
		$j = scalar @pos_neg_array-1;
		$last = 1;
	}
	else{
		$j = $j + 1;
	}

	while($j<(scalar @pos_neg_array-1)){
		$sum = sum(@pos_neg_array[$i..$j]);
		$fold_enrichment = ($sum/($j-$i+1))/$gold_prob;
		$mean_stars = sum(@stars[$i..$j])/($j-$i+1);
		$mean_score = sum(@scores[$i..$j])/($j-$i+1);
		$i=$i+1;
		if ($j + 1 >(scalar @pos_neg_array-1) and !$last){
			$j = scalar @pos_neg_array-1;
			$last = 1;
		}
		else{
			$j = $j + 1;
		}
		print OUT "$mean_stars\t$fold_enrichment\t$mean_score\n";
	}
	close(OUT);



	open(OUT,">$output_file") or die "Unable to open the Output file $output_file\n";
    	
    	$i = 0;
    	$j = $window_size -1;
    	$last = 0;
    	
    	$sum = sum(@pos_neg_array[$i..$j]);
    	$fold_enrichment = ($sum/$window_size)/$gold_prob;
    	$mean_stars = sum(@stars[$i..$j])/$window_size;
    	$mean_score = sum(@scores[$i..$j])/$window_size;
    	
    	print OUT "$mean_stars\t$fold_enrichment\t$mean_score\n";
    	
    	$i = $j + 1;
    	if ($j + $window_size>(scalar @pos_neg_array-1)){
		$j = scalar @pos_neg_array-1;
		$last = 1;
    	}
    	else{
		$j = $j + $window_size;
    	}    
    	while($j<(scalar @pos_neg_array-1)){
		$sum = sum(@pos_neg_array[$i..$j]);
		$fold_enrichment = ($sum/($j-$i+1))/$gold_prob;
		$mean_stars = sum(@stars[$i..$j])/($j-$i+1);
		$mean_score = sum(@scores[$i..$j])/($j-$i+1);
		$i=$j+1;
		if ($j + $window_size>(scalar @pos_neg_array-1) and !$last){
	    		$j = scalar @pos_neg_array-1;
	    		$last = 1;
		}
		else{
	    		$j = $j + $window_size;
		}   
		print OUT "$mean_stars\t$fold_enrichment\t$mean_score\n";
    	}
    	close(OUT);
	
    	print "- The fold_enrichment values have been calculated and stored in the $output_file file\n";
}


