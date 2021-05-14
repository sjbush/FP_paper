use strict;
use warnings;

# REQUIREMENTS
my $root    	     = '/home/ndm.local/steveb';
my $path    	     = "$root/test_suite";
my $progs   	     = "$root/programs";
my $in_dir  	     = "$path/FDA_ARGOS_assembly_vs_ref_genome_VCFs"; # from 1.download_fda_argos_reads_and_assemblies.pl
my $ref_dir 	     = "$path/RefSeq_reference_genomes"; # from 1.download_fda_argos_reads_and_assemblies.pl
my $blastn_path 	 = "$progs/ncbi-blast-2.10.0+/bin/blastn";
my $picard_path		 = "$progs/picard.jar";
my $samtools_path    = "$progs/samtools-1.7/bin/samtools";
my $bedtools_path	 = "$progs/bedtools2/bin/bedtools";
my $maskfasta_path   = "$progs/bedtools2/bin/maskFastaFromBed";
my $fatal   	     = 0;
if (!(-d($in_dir)))  		  { $fatal++; print "ERROR: cannot find $in_dir\n"; 		  }
if (!(-d($ref_dir))) 		  { $fatal++; print "ERROR: cannot find $ref_dir\n"; 		  }
if (!(-e($blastn_path))) 	  { $fatal++; print "ERROR: cannot find $blastn_path\n"; 	  }
if (!(-e($picard_path))) 	  { $fatal++; print "ERROR: cannot find $picard_path\n"; 	  }
if (!(-e($samtools_path)))    { $fatal++; print "ERROR: cannot find $samtools_path\n";    }
if (!(-e($bedtools_path)))    { $fatal++; print "ERROR: cannot find $bedtools_path\n";    }
if (!(-e($maskfasta_path)))   { $fatal++; print "ERROR: cannot find $maskfasta_path\n";   }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;

# OUTPUT
my $out_dir = "$path/masked_RefSeq_reference_genomes";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_file = "$path/summary_of_variants_called_in_FDA_ARGOS_samples_when_aligned_to_RefSeq_reference_genomes.tsv";
my $sh_file  = "$path/mask_non_consensus_positions.sh";
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";
open(OUT,'>',$out_file) or die $!;
print OUT "FDA-ARGOS sample ID\tReference genome (RefSeq assembly accession | assembly name)\tTotal no. of SNPs called, across all instances of nucmer and paftools\tTotal no. of non-SNPs (indels & complex variants) called, across all instances of nucmers and paftools\tNo. of SNPs identically called by all instances of nucmer and paftools\t% of SNPs identically called by all instances of nucmer and paftools\tNo. of non-SNPs identically called by all instances of nucmer and paftools\t% of non-SNPs identically called by all instances of nucmer and paftools\n";

opendir(DIR,$in_dir) or die $!;
my @sample_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
foreach my $sample_id (@sorted_sample_ids)
	{ next if (($sample_id eq '.') or ($sample_id eq '..'));
	  next if (!(-d("$in_dir/$sample_id")));
	  print "$sample_id...\n";
	  opendir(DIR,"$in_dir/$sample_id") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  my %consensus_vars = (); my %non_consensus_vars = ();
	  my $seen_file1 = 0; my $seen_file2 = 0;
	  foreach my $file (@files)
		{ next if (($file eq '.') or ($file eq '..'));
		  if (($file =~ /^$sample_id\.varcall\_relative\_to\.(.+?)\.consensus\_variants.vcf$/) or ($file =~ /^$sample_id\.varcall\_relative\_to\.(.+?)\.non\_consensus\_variants.vcf$/)) # \_incl\_short\_read\_confirmation
			{ my $ref_genome = $1;
			  open(IN,"$in_dir/$sample_id/$file") or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  next if ($line =~ /^\#/);
				  my @line = split(/\t/,$line);
				  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4];
				  if ($file =~ /^$sample_id\.varcall\_relative\_to\.(.+?)\.consensus\_variants.vcf$/) # \_incl\_short\_read\_confirmation
					{ $consensus_vars{$ref_genome}{"$chr:$pos:$ref/$alt"}++; }
				  if ($file =~ /^$sample_id\.varcall\_relative\_to\.(.+?)\.non\_consensus\_variants.vcf$/) # \_incl\_short\_read\_confirmation
					{ $non_consensus_vars{$ref_genome}{"$chr:$pos:$ref/$alt"}++; }
				}
			  close(IN) or die $!;
			  if ($file =~ /^$sample_id\.varcall\_relative\_to\.(.+?)\.consensus\_variants.vcf$/) 	   { $seen_file1++; } # \_incl\_short\_read\_confirmation
			  if ($file =~ /^$sample_id\.varcall\_relative\_to\.(.+?)\.non\_consensus\_variants.vcf$/) { $seen_file2++; } # \_incl\_short\_read\_confirmation
			}
		}
	  if (($seen_file1 == 0) or ($seen_file2 == 0))
		{ print "unable to process: cannot find consensus and non-consensus VCFs in $in_dir/$sample_id...\n"; }
	  next if (($seen_file1 == 0) or ($seen_file2 == 0));
	  if ($seen_file1 == 0) { print "ERROR: cannot find VCF of consensus variants for $sample_id\n"; 	 exit 1; }
	  if ($seen_file2 == 0) { print "ERROR: cannot find VCF of non-consensus variants for $sample_id\n"; exit 1; }
	  
	  # CREATE AN 'AMBIGUOUS POSITION' BED, FOR SUBSEQUENT MASKING...
	  # NOTE THAT WE CAN ONLY DO THIS IF THERE ARE NON-CONSENSUS VARIANTS TO BEGIN WITH!
	  if (scalar keys %non_consensus_vars == 0)
		{ print "note: there are no non-consensus vars in $sample_id, so no 'ambiguous' BED will be created\n"; }
	  while((my $ref_genome,my $irrel)=each(%non_consensus_vars))
		{ my %positions_to_mask = (); my %all_pos = ();
		  while((my $var,my $irrel)=each(%{$non_consensus_vars{$ref_genome}}))
			{ if ($var =~ /^(.*?)\:(.*?)\:(.*?)\/(.*?)$/)
				{ my $chr = $1; my $pos = $2; my $ref = $3; my $alt = $4;
				  my $len_ref = length($ref);
				  for(my $x=$pos;$x<=($pos+($len_ref-1));$x++)
					{ $positions_to_mask{$chr}{$x}++;
					  $all_pos{"$chr:$x"}++;
					}
				}
			}			
		  my $ambiguous_bed 	   = "$out_dir/ambiguous_positions_for_$sample_id.bed"; # i.e., these positions are ambiguous should you align the assembled $sample_id (the Illumina reads) to the reference genome
		  my $confident_bed 	   = "$out_dir/confident_positions_for_$sample_id.bed";
		  my $ref_genome_fa 	   = "$ref_dir/$ref_genome.fa";
		  my $ref_genome_size_file = "$out_dir/ref_genome_for_$sample_id.chr_lengths";
		  my $masked_ref_genome_fa = "$out_dir/ref_genome_for_$sample_id.fa";
		  if (!(-e($ref_genome_fa))) { print "ERROR: unable to find $ref_genome_fa\n"; exit 1; }
		  next if ( (-e($ambiguous_bed)) and (-e($confident_bed)) and (-e($masked_ref_genome_fa)) and (-e($ref_genome_size_file)) ); # CHECKPOINT: we've created this masked genome before
		  
		  # report chromosome sizes; this is necessary to run BEDtools complement
		  open(SIZE,'>',$ref_genome_size_file) or die $!;
		  my $chr = ''; my %seqs = ();
		  open(IN,$ref_genome_fa) or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  if ($line =~ /^\>(.*?) .*?$/)
				{ $chr = $1; }
			  next if ($line =~ /^\>/);
			  $seqs{$chr} .= $line unless ($chr eq '');
			}
		  close(IN) or die $!;
		  my @chrs = ();
		  while((my $chr,my $irrel)=each(%seqs))
			{ push(@chrs,$chr); }
		  my @sorted_chrs = sort {$a cmp $b} @chrs;
		  foreach my $chr (@sorted_chrs)
			{ my $length = length($seqs{$chr});
			  print SIZE "$chr\t$length\n";
			}
		  close(SIZE) or die $!;
		  
		  if (!(-e($ambiguous_bed)))
			{ open(BED,'>',$ambiguous_bed) or die $!;
			  my @chrs = ();
			  while((my $chr,my $irrel)=each(%positions_to_mask))
				{ push(@chrs,$chr); }
			  my @sorted_chrs = sort {$a cmp $b} @chrs;
			  my %printed_to_bed = ();
			  foreach my $chr (@sorted_chrs)
				{ my @pos = ();
				  while((my $pos,my $irrel)=each(%{$positions_to_mask{$chr}}))
					{ push(@pos,$pos); }
				  my @sorted_pos = sort {$a <=> $b} @pos;
				  foreach my $pos (@sorted_pos)
					{ my $pos_minus1 = $pos-1;
					  print BED "$chr\t$pos_minus1\t$pos\n" unless (exists($printed_to_bed{$pos}));
					  $printed_to_bed{$pos}++;
					}
				}
			  close(BED) or die $!;
			  # collapse the contents of $ambiguous_bed: combine overlapping or "book-ended" features in an interval file into a single feature (https://bedtools.readthedocs.io/en/latest/content/tools/merge.html)
			  my $ambiguous_bed_sorted = "$out_dir/$sample_id".".vs.$ref_genome.ambiguous_positions.sorted_bed";
			  print SH "sort -k1,1 -k2,2n $ambiguous_bed > $ambiguous_bed_sorted\n";
			  print SH "$bedtools_path merge -i $ambiguous_bed_sorted > $ambiguous_bed\n";
			  print SH "rm $ambiguous_bed_sorted\n";
			}
		  
		  if ( (-e($ambiguous_bed)) and (!(-e($confident_bed))) )
			{ # we will now create a complement to $ambiguous_bed: a file that represents "all intervals in a genome that are not covered by at least one interval in the input" (https://bedtools.readthedocs.io/en/latest/content/tools/complement.html)
			  # we need this file because when evaluating pipelines we will use it as input to hap.py. It will oblige hap.py to discard all calls made outside of these intervals
			  print SH "$bedtools_path complement -i $ambiguous_bed -g $ref_genome_size_file > $confident_bed\n";
			}
		  
		  # (SOFT) MASK REFERENCE GENOME TO HIDE NON-CONSENSUS VARIANTS
		  if ( (!(-e($masked_ref_genome_fa))) and (!(-e("$masked_ref_genome_fa.fai"))) )
			{ print SH "$maskfasta_path -fi $ref_genome_fa -fo $masked_ref_genome_fa -bed $ambiguous_bed -fullHeader -soft\n";
			  print SH "$samtools_path faidx $masked_ref_genome_fa\n";
			}
		}
	  
	  # ... ALTHOUGH IF THERE AREN'T ANY AMBIGUOUS POSITIONS, JUST MAKE A COPY OF THE ORIGINAL GENOME
	  while((my $ref_genome,my $irrel)=each(%consensus_vars))
		{ next if (exists($non_consensus_vars{$ref_genome})); # i.e. we would have seen this above
		  my $ref_genome_fa = "$ref_dir/$ref_genome.fa";
		  if (!(-e($ref_genome_fa))) { print "ERROR: unable to find $ref_genome_fa\n"; exit 1; }
		  my $masked_ref_genome_fa = "$out_dir/ref_genome_for_$sample_id.fa";
		  print SH "cp $ref_genome_fa $masked_ref_genome_fa\n"    unless (-e($masked_ref_genome_fa));
		  print SH "$samtools_path faidx $masked_ref_genome_fa\n" unless (-e("$masked_ref_genome_fa.fai"));
		  print SH "java -jar $picard_path CreateSequenceDictionary REFERENCE=$masked_ref_genome_fa OUTPUT=$out_dir/ref_genome_for_$sample_id.dict\n" unless (-e("$out_dir/ref_genome_for_$sample_id.dict"));
		}
	  
	  # OUTPUT SUMMARY STATISTICS
	  while((my $ref_genome,my $irrel)=each(%consensus_vars))
		{ my $ref_genome_acc = ''; my $assembly_name = '';
		  if ($ref_genome =~ /^(.+)\_(ASM.*?)$/)
			{ $ref_genome_acc = $1; $assembly_name = $2; }
		  # how many SNPs & non-SNPs (indels and complex variants) are there in total?
		  my $total_no_of_snps = 0; my $total_no_of_non_snps = 0;
		  my $no_of_consensus_snps = 0; my $no_of_consensus_non_snps = 0;
		  while((my $var,my $irrel)=each(%{$consensus_vars{$ref_genome}}))
			{ if ($var =~ /^(.*?)\:(.*?)\:(.*?)\/(.*?)$/)
				{ my $chr = $1; my $pos = $2; my $ref = $3; my $alt = $4;
				  my $len_ref = length($ref); my $len_alt = length($alt);
				  if (($len_ref == 1) && ($len_alt == 1))
					{ $total_no_of_snps++;
					  $no_of_consensus_snps++;
					}
				  else
					{ $total_no_of_non_snps++;
					  $no_of_consensus_non_snps++;
					}
				}
			}
		  while((my $var,my $irrel)=each(%{$non_consensus_vars{$ref_genome}}))
			{ if ($var =~ /^(.*?)\:(.*?)\:(.*?)\/(.*?)$/)
				{ my $chr = $1; my $pos = $2; my $ref = $3; my $alt = $4;
				  my $len_ref = length($ref); my $len_alt = length($alt);
				  if (($len_ref == 1) && ($len_alt == 1))
					{ $total_no_of_snps++; }
				  else
					{ $total_no_of_non_snps++; }
				}
			}
		  my $pc_of_consensus_snps 	   = sprintf("%.2f",(($no_of_consensus_snps/$total_no_of_snps)*100));
		  my $pc_of_consensus_non_snps = sprintf("%.2f",(($no_of_consensus_non_snps/$total_no_of_non_snps)*100));
		  print OUT "$sample_id\t$ref_genome_acc | $assembly_name\t$total_no_of_snps\t$total_no_of_non_snps\t$no_of_consensus_snps\t$pc_of_consensus_snps\t$no_of_consensus_non_snps\t$pc_of_consensus_non_snps\n";
		}
	}
close(SH) or die $!; close(OUT) or die $!;
exit 1;