use strict;
use warnings;

# REQUIREMENTS
my $root     = '/home/ndm.local/steveb';
my $path     = "$root/test_suite";
my $in_file1 = "$path/FDA_ARGOS_samples_used_for_testing.tsv"; # from 1.download_fda_argos_reads_and_assemblies.pl
my $in_file2 = "$path/summary_of_variants_called_in_FDA_ARGOS_samples_when_aligned_to_RefSeq_reference_genomes.tsv"; # from 2.create_masked_ref_genomes.pl
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file = "$path/summary_of_FDA_ARGOS_samples_used.tsv";
open(OUT,'>',$out_file) or die $!;
print OUT "FDA-ARGOS sample name\tSpecies\tTaxonomy ID\tSRA run accession(s) for associated Illumina reads\tRefSeq assembly accession for associated Illumina/PacBio hybrid assembly\tNo. of RefSeq reference genomes for this species (not including the FDA-ARGOS assembly)\tRefSeq assembly accessions for reference genomes\tMash distance between FDA-ARGOS assembly and reference genome\tLength of FDA-ARGOS assembly (bp)\t";
print OUT "Number of low-quality bases masked as N (those where, after re-mapping Illumina reads to the FDA-ARGOS assembly, there was no coverage or where the most common nucleotide was represented by < 99% of the aligned bases)\t% of low-quality bases masked as N\t";
print OUT "Number of discordant bases masked as N (those bases called as variants in the FDA-ARGOS assembly when using Snippy with the Illumina reads as input; these positions represent discordant calls between the Illumina and PacBio sequence)\t% of discordant bases masked as N\t";
print OUT "Total number of bases masked as N\tTotal % of bases masked as N\t";
print OUT "Total no. of SNPs called after aligning the FDA-ARGOS assembly to the reference genome, across all instances of nucmer and paftools\tTotal no. of non-SNPs (indels & complex variants) called after aligning the FDA-ARGOS assembly to the reference genome, across all instances of nucmers and paftools\tNo. of SNPs identically called by all instances of nucmer and paftools\t% of SNPs identically called by all instances of nucmer and paftools\tNo. of non-SNPs identically called by all instances of nucmer and paftools	% of non-SNPs identically called by all instances of nucmer and paftools\n";

# HOW MANY VARIANTS HAVE BEEN CALLED PER SAMPLE?
my %variant_line = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $fda_sample_id = $line[0];
	  $variant_line{$fda_sample_id} = "$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]";
	}
close(IN) or die $!;

# CREATE A SUMMARY TABLE DETAILING THE AMOUNT OF MASKING PERFORMED ON EACH SAMPLE, AND HOW MANY VARIANTS ARE AVAILBLE TO CALL.
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $fda_sample_id = $line[0]; my $species = $line[1]; my $taxon_id = $line[2]; my $sra_run_id = $line[3]; my $fda_assembly_acc = $line[4]; my $no_of_ref_genomes = $line[5]; my $ref_genome_accs = $line[6];
	  
	  # what is the mash distance between the FDA-ARGOS assembly and the reference genome(s)?
	  my $mash_dist = '';
	  my @ref_genome_accs = split(/\, /,$ref_genome_accs);
	  foreach my $ref_genome_acc (@ref_genome_accs)
		{ my $mash_dist_file = "$path/FDA_ARGOS_assemblies/$fda_sample_id/distance_to_$ref_genome_acc.tsv";
		  if (!(-e($mash_dist_file))) { print "WARNING: cannot find $mash_dist_file\n"; }
		  next if (!(-e($mash_dist_file)));
		  open(TMP,$mash_dist_file) or die $!;
		  while(<TMP>)
			{ my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  my $dist = $line[2];
			  $mash_dist .= "$dist, ";
			}
		  close(TMP) or die $!;
		}
	  $mash_dist =~ s/\, $//;
	  if ($mash_dist eq '') { $mash_dist = 'unknown'; }
	  
	  # how many bases in the FDA-ARGOS assembly have been masked as low-quality, both due to low coverage and due to miscalls between the Illumina & PacBio?
	  my $proportion_masked_file_lowcov  = "$path/FDA_ARGOS_assemblies/$fda_sample_id/proportion_of_bases_masked_in_$fda_assembly_acc.lowcov.tsv";
	  my $proportion_masked_file_miscall = "$path/FDA_ARGOS_assemblies/$fda_sample_id/proportion_of_bases_masked_in_$fda_assembly_acc.miscall.tsv";
	  if (!(-e($proportion_masked_file_lowcov)))  { print "WARNING: cannot find $proportion_masked_file_lowcov\n";  }
	  if (!(-e($proportion_masked_file_miscall))) { print "WARNING: cannot find $proportion_masked_file_miscall\n"; }
	  my $assembly_length_lowcov = ''; my $no_of_n_lowcov = ''; my $pc_of_n_lowcov = '';
	  if (-e($proportion_masked_file_lowcov))
		{ open(TMP,$proportion_masked_file_lowcov) or die $!;
		  while(<TMP>)
			{ next if ($. == 1);
			  my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  $assembly_length_lowcov = $line[0]; $no_of_n_lowcov = $line[1]; $pc_of_n_lowcov = $line[2];
			}
		  close(TMP) or die $!;
		}
	  my $assembly_length_miscall = ''; my $no_of_n_miscall = ''; my $pc_of_n_miscall = '';
	  if (-e($proportion_masked_file_miscall))
		{ open(TMP,$proportion_masked_file_miscall) or die $!;
		  while(<TMP>)
			{ next if ($. == 1);
			  my $line = $_; chomp($line);
			  my @line = split(/\t/,$line);
			  $assembly_length_miscall = $line[0]; $no_of_n_miscall = $line[1]; $pc_of_n_miscall = $line[2];
			}
		  close(TMP) or die $!;
		}
	  if ($assembly_length_lowcov != $assembly_length_miscall)
		{ print "WARNING: assembly length discrepancy between the 'low cov' and 'miscall' files for $fda_sample_id: $assembly_length_lowcov != $assembly_length_miscall\n"; }
	  my $total_no_of_n = $no_of_n_lowcov+$no_of_n_miscall;
	  my $total_pc_of_n = sprintf("%.4f",(($total_no_of_n/$assembly_length_lowcov)*100));
	  
	  next if (!(exists($variant_line{$fda_sample_id})));
	  
	  print OUT "$fda_sample_id\t$species\t$taxon_id\t$sra_run_id\t$fda_assembly_acc\t$no_of_ref_genomes\t$ref_genome_accs\t$mash_dist\t$assembly_length_lowcov\t$no_of_n_lowcov\t$pc_of_n_lowcov\t$no_of_n_miscall\t$pc_of_n_miscall\t$total_no_of_n\t$total_pc_of_n\t$variant_line{$fda_sample_id}\n";
	}
close(IN) or die $!;

close(OUT) or die $!;
exit 1;