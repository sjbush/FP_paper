use strict;
use warnings;

### REQUIREMENTS ###
# ESSENTIAL INPUT FILES: SRA RUN ACCESSIONS FOR THE FDA-ARGOS DATA, PLUS THE NCBI REFSEQ ASSEMBLY INDEX
my $root    = '/home/ndm.local/steveb';
my $path    = "$root/test_suite";
my $progs   = "$root/programs";
my $sra_ids = "$path/PRJNA231221.tsv"; # see https://www.nature.com/articles/s41467-019-11306-6#data-availability
my $ass_sum = "$path/assembly_summary.13.02.20.txt"; # ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt; updates daily
my $fatal 	= 0;
if (!(-d($path)))    { $fatal++; print "error: cannot find $path\n";    }
if (!(-d($progs)))   { $fatal++; print "error: cannot find $progs\n";   }
if (!(-e($sra_ids))) { $fatal++; print "error: cannot find $sra_ids\n"; }
if (!(-e($ass_sum))) { $fatal++; print "error: cannot find $ass_sum\n"; }
# PATHS TO ALIGNERS, CALLERS, AND OTHER SUBSIDARY TOOLS
my $vt_dir    	  		= "$progs/vt";
my $bwa_dir 	 	    = "$progs/bwa-0.7.17";
my $mash_path			= "$progs/mash-Linux64-v2.1/mash";
my $seqtk_dir    	  	= "$progs/seqtk";
my $bgzip_path     		= "$progs/bin/bgzip";
my $tabix_path     		= "$progs/bin/tabix";
my $wgsim_path	  		= "$progs/samtools-1.7/bin/wgsim";
my $nucmer_dir 			= "$progs/mummer-4.0.0beta2";
my $prokka_dir			= "$progs/prokka/bin";
my $seqtk_path    	    = "$progs/seqtk/seqtk";
my $prepy_path			= "$progs/hap.py-build/bin/pre.py";
my $picard_path   		= "$progs/picard.jar";
my $nucmer_path  		= "$progs/mummer-4.0.0beta2/nucmer";
my $snippy_path			= "$progs/snippy/bin/snippy";
my $spandx_path			= "$progs/SPANDx/SPANDx.sh";
my $bcftools_dir  		= "$progs/bcftools";
my $dnadiff_path 		= "$progs/mummer-4.0.0beta2/dnadiff";
my $samtools_dir		= "$progs/samtools-1.7/bin";
my $samtools_path 		= "$progs/samtools-1.7/bin/samtools";
my $bcftools_path  		= "$progs/bcftools/bcftools";
my $minimap2_path	    = "$progs/minimap2/minimap2";
my $paftools_path		= "$progs/minimap2/misc/paftools.js";
my $showsnps_path 		= "$progs/mummer-4.0.0beta2/show-snps";
my $maskfasta_path		= "$progs/bedtools2/bin/maskFastaFromBed";
my $MUMmerSNPs2VCF_path = "$progs/MUMmerSNPs2VCF.py"; # from https://github.com/liangjiaoxue/PythonNGSTools
if (!(-d($vt_dir)))    			 { $fatal++; print "error: cannot find $vt_dir\n";    			}
if (!(-d($bwa_dir)))    		 { $fatal++; print "error: cannot find $bwa_dir\n";    			}
if (!(-e($mash_path)))  		 { $fatal++; print "error: cannot find $mash_path\n";  			}
if (!(-d($seqtk_dir)))    		 { $fatal++; print "error: cannot find $seqtk_dir\n";    		}
if (!(-d($mauve_dir)))  		 { $fatal++; print "error: cannot find $mauve_dir\n";  			}
if (!(-e($mauve_jar)))  		 { $fatal++; print "error: cannot find $mauve_jar\n";  			}
if (!(-e($mauve_path)))  		 { $fatal++; print "error: cannot find $mauve_path\n";  		}
if (!(-e($bgzip_path)))    		 { $fatal++; print "error: cannot find $bgzip_path\n";   		}
if (!(-e($tabix_path)))    		 { $fatal++; print "error: cannot find $tabix_path\n";    		}
if (!(-e($wgsim_path)))    		 { $fatal++; print "error: cannot find $wgsim_path\n";    		}
if (!(-d($nucmer_dir)))  		 { $fatal++; print "error: cannot find $nucmer_dir\n";  		}
if (!(-d($prokka_dir))) 		 { $fatal++; print "error: cannot find $prokka_dir\n";  	    }
if (!(-e($seqtk_path)))  		 { $fatal++; print "error: cannot find $seqtk_path\n";  		}
if (!(-e($prepy_path)))  		 { $fatal++; print "error: cannot find $prepy_path\n";  		}
if (!(-e($picard_path))) 		 { $fatal++; print "error: cannot find $picard_path\n";  	    }
if (!(-e($nucmer_path)))  		 { $fatal++; print "error: cannot find $nucmer_path\n";  		}
if (!(-e($snippy_path))) 		 { $fatal++; print "error: cannot find $snippy_path\n";  	  	}
if (!(-e($spandx_path))) 		 { $fatal++; print "error: cannot find $spandx_path\n";  	  	}
if (!(-e($bcftools_dir)))  		 { $fatal++; print "error: cannot find $bcftools_dir\n";  		}
if (!(-e($dnadiff_path))) 		 { $fatal++; print "error: cannot find $dnadiff_path\n"; 		}
if (!(-e($samtools_dir)))  		 { $fatal++; print "error: cannot find $samtools_dir\n";  		}
if (!(-e($samtools_path))) 		 { $fatal++; print "error: cannot find $samtools_path\n"; 		}
if (!(-e($bcftools_path))) 		 { $fatal++; print "error: cannot find $bcftools_path\n";  		}
if (!(-e($minimap2_path))) 		 { $fatal++; print "error: cannot find $minimap2_path\n"; 		}
if (!(-e($paftools_path))) 		 { $fatal++; print "error: cannot find $paftools_path\n"; 		}
if (!(-e($showsnps_path))) 		 { $fatal++; print "error: cannot find $showsnps_path\n"; 		}
if (!(-e($maskfasta_path))) 	 { $fatal++; print "error: cannot find $maskfasta_path\n"; 		}
if (!(-e($MUMmerSNPs2VCF_path))) { $fatal++; print "error: cannot find $MUMmerSNPs2VCF_path\n"; }
exit if ($fatal > 0);

### PARAMETERS ###
my $num_procs = 12;
my @shortlist = ('FDAARGOS_338','FDAARGOS_536','FDAARGOS_598','FDAARGOS_687','FDAARGOS_700','FDAARGOS_751');
my %shortlist = map {$_ => 1} @shortlist;

### OUTPUT ###
my $out_dir1 = "$path/FDA_ARGOS_illumina_reads";
my $out_dir2 = "$path/FDA_ARGOS_assemblies";
my $out_dir3 = "$path/RefSeq_reference_genomes";
my $out_dir4 = "$path/FDA_ARGOS_negative_controls";
my $out_dir5 = "$path/FDA_ARGOS_assembly_vs_ref_genome_VCFs";
if (!(-d($out_dir1))) { mkdir $out_dir1 or die $!; }
if (!(-d($out_dir2))) { mkdir $out_dir2 or die $!; }
if (!(-d($out_dir3))) { mkdir $out_dir3 or die $!; }
if (!(-d($out_dir4))) { mkdir $out_dir4 or die $!; }
if (!(-d($out_dir5))) { mkdir $out_dir5 or die $!; }
my $sh_file  = "$path/download_FDA_ARGOS_data.sh";
my $err_file = "$path/download_FDA_ARGOS_data.err";
my $out_file = "$path/FDA_ARGOS_samples_used_for_testing.tsv";
open(SH,'>',$sh_file) or die $!; open(ERR,'>',$err_file) or die $!; open(OUT,'>',$out_file) or die $!;
print OUT "FDA-ARGOS sample name\tSpecies\tTaxonomy ID\tSRA run accession(s) for associated Illumina reads\tRefSeq assembly accession for associated Illumina/PacBio hybrid assembly\tNo. of RefSeq reference genomes for this species (not including the FDA-ARGOS assembly)\tRefSeq assembly accessions for reference genomes\n";
print SH "#!/bin/bash\n";
print SH "export PATH=:$progs/bin:\$PATH\n"; # required for Snippy (specifically, to obtain GNU Parallel)
print SH "export PATH=:$bwa_dir:\$PATH\n"; # required for Snippy
print SH "export PATH=:$bcftools_dir:\$PATH\n"; # required for Snippy
print SH "PATH=\$PATH:$samtools_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$prokka_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$seqtk_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$vt_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$progs\n"; # required for k8 (necessary to use minimap2's paftools)

# STORE ENA RUN ACCESSION IDs FOR EACH FDA-ARGOS (BACTERIAL) SAMPLE
my %data = (); my %orgs = (); my %biosample_ids = (); my %fda_argos_biosamples = ();
open(IN,$sra_ids) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $assay_type = $line[0]; my $biosample_id = $line[2]; my $library_layout = $line[10]; my $library_source = $line[12]; my $species = $line[17]; my $platform = $line[18]; my $run_id = $line[20]; my $sample_name = $line[22];
	  next if ($assay_type ne 'WGS');
	  next if ($library_source ne 'GENOMIC');
	  if ( (($platform eq 'ILLUMINA') and ($library_layout eq 'PAIRED')) or (($platform eq 'PACBIO_SMRT') and ($library_layout eq 'SINGLE')) ) # CHECKPOINT: only bacterial samples in the FDA-ARGOS dataset have paired Illumina & PacBio reads
		{ $data{$sample_name}{$platform}{$run_id}++;
		  $orgs{$sample_name} = $species;
		  $biosample_ids{$sample_name} = $biosample_id;
		  $fda_argos_biosamples{$biosample_id}++;
		}
	}
close(IN) or die $!;

# IDENTIFY THE SET OF COMPLETE NCBI BACTERIAL GENOMES. WE WILL USE THIS DATA TO FIND THE ILLUMINA/PACBIO HYBRID ASSEMBLY ASSOCIATED WITH EACH $sample_name, ABOVE, *AND* TO FIND THE REPRESENTATIVE/REFERENCE GENOMES OF THAT SPECIES.
my %fda_argos_assemblies = (); my %other_assemblies = (); my %species_taxid_for_assembly_accessions = ();
open(IN,$ass_sum) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line =~ /^\#/);
	  my @line = split(/\t/,$line);
	  my $assembly_accession = $line[0]; my $biosample_id = $line[2]; my $refseq_category = $line[4]; my $taxid = $line[5]; my $species_taxid = $line[6]; my $organism_name = $line[7]; my $intraspecific_name = $line[8]; my $version = $line[10]; my $assembly_level = $line[11]; my $genome_rep = $line[13]; my $asm_name = $line[15]; my $ftp_path = $line[19];
	  $organism_name =~ s/Mycobacterium tuberculosis variant bovis/Mycobacterium bovis/; # to avoid downstream retention of M. bovis as a reference genome of "M. tuberculosis", equivalent to H37Rv
	  my $genus = ''; my $species = '';
	  if ($organism_name =~ /^(\w+) (.+?) .+$/)
		{ $genus   = $1;
		  if ($2 !~ /sp\./)
			{ $species = "$1 $2"; }
		  else
			{ $species = $organism_name; }
		}
	  elsif ($organism_name =~ /^(\w+) (.+)$/)
		{ $genus   = $1;
		  $species = "$1 $2";
		}
	  next if (($species eq '') or ($genus eq ''));
	  next if ($species !~ /^[A-Z].*?$/);
	  next if ($species =~ /^.*? sp\.$/);
	  next if (($taxid !~ /\d+/) or ($species_taxid !~ /\d+/));
	  next if ($version ne 'latest');
	  next if ($assembly_level ne 'Complete Genome');
	  next if ($genome_rep ne 'Full'); # CHECKPOINTS: we require unambiguous classification at the levels of both genus and species, and are interested only in the latest versions of complete genomes
	  
	  # STORE THE COMPLETE GENOME ASSEMBLIES AVAILABLE FOR THE FDA-ARGOS SAMPLES...
	  if (exists($fda_argos_biosamples{$biosample_id}))
		{ if ($intraspecific_name =~ /^strain\=(.*?)$/)
			{ my $fda_argos_name = $1;
			  $fda_argos_assemblies{$species}{$biosample_id}{$ftp_path}{acc} 			= $assembly_accession;
			  $fda_argos_assemblies{$species}{$biosample_id}{$ftp_path}{species_taxid}  = $species_taxid;
			  $fda_argos_assemblies{$species}{$biosample_id}{$ftp_path}{fda_argos_name} = $fda_argos_name;
			}
		}
	  else
		{ # ... AND STORE ALL OTHER COMPLETE GENOME ASSEMBLIES AVAILABLE FOR THIS SPECIES, PROVIDED THEY ARE CONSIDERED A 'REFERENCE GENOME' (DEFINED AT ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt AS: "a manually selected high quality genome assembly that NCBI and the community have identified as being important as a standard against which other data are compared")
		  if ($refseq_category eq 'reference genome')
			{ next if (($species =~ /^Salmonella enterica/) and ($assembly_accession ne 'GCF_000195995.1'));
			  next if (($species =~ /^Bacillus anthracis/) and ($assembly_accession ne 'GCF_000007845.1'));
			  next if (($species =~ /^Escherichia coli/) and ($assembly_accession ne 'GCF_000005845.2'));
			  $other_assemblies{$species}{$assembly_accession} = $ftp_path;
			  $species_taxid_for_assembly_accessions{$assembly_accession} = $species_taxid;
			}
		}
	}
close(IN) or die $!;

# DOWNLOAD THE SET OF FDA-ARGOS ILLUMINA READS FOR THOSE BACTERIA THAT HAVE A CORRESPONDING ASSEMBLY (UNIQUELY NAMED AS BEING OF FDA-ARGOS ORIGIN AND CONSIDERED 'COMPLETE'), PLUS THE FDA-ARGOS ASSEMBLY ITSELF *AND* ONE OR MORE NCBI REFERENCE GENOMES FOR THE SAME SPECIES.
# WE WILL THEN RUN A NEGATIVE CONTROL: ALIGNING THE ILLUMINA READS AGAINST THE FDA-ARGOS ASSEMBLY AND CALLING VARIANTS. AS THE FORMER IS IN PRINCIPLE SOURCED FROM THE LATTER, WE EXPECT NO VARIANTS TO BE FOUND.
my @sample_names = ();
while((my $sample_name,my $irrel)=each(%data))
	{ push(@sample_names,$sample_name); }
my @sorted_sample_names = sort {$a cmp $b} @sample_names;
my %downloaded_already = ();
foreach my $sample_name (@sorted_sample_names)
	{ next if (!(exists($data{$sample_name}{'PACBIO_SMRT'}))) or (!(exists($data{$sample_name}{'ILLUMINA'}))); # CHECKPOINT: we require that each sample be associated with both Illumina and PacBio reads (although we're only going to obtain the former)
	  next if ( (!(exists($biosample_ids{$sample_name}))) or (!(exists($orgs{$sample_name}))) );
	  next if (!(exists($shortlist{$sample_name}))); # CHECKPOINT: restrict to a shortlist
	  my $biosample_id = $biosample_ids{$sample_name};
	  my $species = $orgs{$sample_name};
	  next if (!(exists($fda_argos_assemblies{$species}{$biosample_id})));
	  next if ((scalar keys %{$fda_argos_assemblies{$species}{$biosample_id}}) != 1); # CHECKPOINT: we require that each FDA-ARGOS sample be associated with one, and only one, complete genome assembly (lower levels of assembly are available for many FDA-ARGOS samples, however)
	  my $failure = 0;
	  my $ftp_path_FDA_ASSEMBLY = ''; my $assembly_accession_FDA_ASSEMBLY = ''; my $species_taxid_FDA_ASSEMBLY = '';
	  while((my $ftp_path,my $irrel)=each(%{$fda_argos_assemblies{$species}{$biosample_id}}))
		{ $ftp_path_FDA_ASSEMBLY = $ftp_path;
		  my $assembly_accession = $fda_argos_assemblies{$species}{$biosample_id}{$ftp_path}{acc};
		  my $species_taxid      = $fda_argos_assemblies{$species}{$biosample_id}{$ftp_path}{species_taxid};
		  my $fda_argos_name     = $fda_argos_assemblies{$species}{$biosample_id}{$ftp_path}{fda_argos_name};
		  $species_taxid_FDA_ASSEMBLY      = $species_taxid;
		  $assembly_accession_FDA_ASSEMBLY = $assembly_accession;
		  if ($sample_name ne $fda_argos_name) # CHECKPOINT: is the FDA-ARGOS sample name the same for the Illumina reads as for the assembly these reads were supposed to contribute to? If not, abort.
			{ $failure++; }
		}
	  next if ($failure > 0);
	  next if (!(exists($other_assemblies{$species}))); # CHECKPOINT: there are no other reference genomes for this species
	  my $no_of_reference_genomes = scalar keys %{$other_assemblies{$species}};
	  my @ref_accs = ();
	  while((my $ref_acc,my $irrel)=each(%{$other_assemblies{$species}}))
		{ my $species_taxid = $species_taxid_for_assembly_accessions{$ref_acc};
		  next if ($species_taxid != $species_taxid_FDA_ASSEMBLY); # CHECKPOINT: is the taxonomy ID of the FDA-ARGOS assembly the same as that assigned to the reference genome (the species name will be)?
		  push(@ref_accs,$ref_acc);
		}
	  my @sorted_ref_accs = sort {$a cmp $b} @ref_accs;
	  my $ref_accs = join(", ",@sorted_ref_accs);
	  
	  # (1) DOWNLOAD ILLUMINA READS
	  if (!(-d("$out_dir1/$sample_name"))) { print SH "mkdir $out_dir1/$sample_name\n"; }
	  print SH "cd $out_dir1/$sample_name\n";
	  my $fq_1_list = ''; my $fq_2_list = '';
	  while((my $platform,my $irrel)=each(%{$data{$sample_name}}))
		{ next if ($platform ne 'ILLUMINA'); # CHECKPOINT: restrict download to Illumina reads only; we don't need the PacBio if we've already got the assemblies
		  my @run_ids = ();
		  while((my $run_id,my $irrel)=each(%{$data{$sample_name}{$platform}}))
			{ my $first_3; my $first_6; my $digits;
			  if ($run_id =~ /^(.{3}).*?$/) { $first_3 = $1; }
			  if ($run_id =~ /^(.{6}).*?$/) { $first_6 = $1; }
			  if ($run_id =~ /^.*?(\d+)$/)  { $digits  = $1; }
			  my $number_of_digits = length($digits);
			  my $ena1 = ''; my $ena2 = ''; my $ena_single = '';
			  if ($number_of_digits == 6)
				{ $ena1 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
				  $ena2 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
				}
			  elsif ($number_of_digits == 7)
				{ if ($digits =~ /^.+?(\d{1})$/)
					{ my $last_digit = $1;
					  $ena1 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
					  $ena2 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
					}
				}
			  elsif ($number_of_digits == 8)
				{ if ($digits =~ /^.+?(\d{2})$/)
					{ my $last_two_digits = $1;
					  $ena1 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
					  $ena2 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
					}
				}
			  elsif ($number_of_digits == 9)
				{ if ($digits =~ /^.+?(\d{3})$/)
					{ my $last_three_digits = $1;
					  $ena1 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
					  $ena2 = "ftp://ftp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
					}
				}
			  else
				{ print "ERROR: unable to determine download URL for $run_id\n"; exit 1; }
			  my $fq_1 = "$out_dir1/$sample_name/$run_id"."_1.fastq.gz";
			  if (!(-e($fq_1)))
				{ print SH "if curl --head --fail --silent \"$ena1\" >/dev/null; then wget $ena1; else echo \"ERROR: unable to find file at URL $ena1\" >> $err_file && exit 1; fi\n"; # CHECKPOINT: fail if cannot find the file at this URL
				  print SH "if [ ! -e $fq_1 ]; then echo \"ERROR: unable to find $fq_1 after supposedly downloading\" >> $err_file && exit 1; fi\n"; # CHECKPOINT: fail if we have cannot find the file we've supposedly downloaded
				}
			  my $fq_2 = "$out_dir1/$sample_name/$run_id"."_2.fastq.gz";
			  if (!(-e($fq_2)))
				{ print SH "if curl --head --fail --silent \"$ena2\" >/dev/null; then wget $ena2; else echo \"ERROR: unable to find file at URL $ena2\" >> $err_file && exit 1; fi\n";
				  print SH "if [ ! -e $fq_2 ]; then echo \"ERROR: unable to find $fq_2 after supposedly downloading\" >> $err_file && exit 1; fi\n";
				}
			  push(@run_ids,$run_id);
			  $fq_1_list .= "$fq_1 ";
			  $fq_2_list .= "$fq_2 ";
			}
		  my @sorted_run_ids = sort {$a cmp $b} @run_ids;
		  my $run_ids = join(", ",@sorted_run_ids);
		  print OUT "$sample_name\t$species\t$species_taxid_FDA_ASSEMBLY\t$run_ids\t$assembly_accession_FDA_ASSEMBLY\t$no_of_reference_genomes\t$ref_accs\n";
		}
	  $fq_1_list =~ s/\s$//; $fq_2_list =~ s/\s$//;
	
	  # (2) DOWNLOAD FDA-ARGOS ASSEMBLY
	  if (!(-d("$out_dir2/$sample_name"))) { print SH "mkdir $out_dir2/$sample_name\n"; }
	  print SH "cd $out_dir2/$sample_name\n";
	  my $last_bit = '';
	  if ($ftp_path_FDA_ASSEMBLY =~ /^.+\/(.*?)$/) { $last_bit = $1; }
	  $ftp_path_FDA_ASSEMBLY .= "/$last_bit"."_genomic.fna.gz";
	  my $fda_assembly 		  = "$out_dir2/$sample_name/$last_bit.fa";
	  my $fda_assembly_idx    = "$out_dir2/$sample_name/$last_bit";
	  if ( (!(-e($fda_assembly))) and (!(exists($downloaded_already{$fda_assembly}))) )
		{ print SH "if curl --head --fail --silent \"$ftp_path_FDA_ASSEMBLY\" >/dev/null; then wget $ftp_path_FDA_ASSEMBLY -O $fda_assembly.gz; else echo \"ERROR: unable to find file at URL $ftp_path_FDA_ASSEMBLY\" >> $err_file && exit 1; fi\n";
		  print SH "if [ ! -e $fda_assembly.gz ]; then echo \"ERROR: unable to find $fda_assembly.gz after supposedly downloading\" >> $err_file && exit 1; fi\n";
		  print SH "gunzip $fda_assembly.gz\n";
		  print SH "if [ ! -e $fda_assembly ]; then echo \"ERROR: unable to find $fda_assembly after supposedly unzipping\" >> $err_file && exit 1; fi\n";
		  $downloaded_already{$fda_assembly}++;
		}
	  
	  # (3) DOWNLOAD REFERENCE GENOME ASSEMBLIES
	  while((my $ref_acc,my $ftp_path)=each(%{$other_assemblies{$species}}))
		{ my $last_bit = '';
		  if ($ftp_path =~ /^.+\/(.*?)$/) { $last_bit = $1; }
		  $ftp_path .= "/$last_bit"."_genomic.fna.gz";
		  my $ref_assembly = "$out_dir3/$last_bit.fa";
		  if ( (!(-e($ref_assembly))) and (!(exists($downloaded_already{$ref_assembly}))) )
			{ print SH "if curl --head --fail --silent \"$ftp_path\" >/dev/null; then wget $ftp_path -O $ref_assembly.gz; else echo \"ERROR: unable to find file at URL $ftp_path\" >> $err_file && exit 1; fi\n";
			  print SH "if [ ! -e $ref_assembly.gz ]; then echo \"ERROR: unable to find $ref_assembly.gz after supposedly downloading\" >> $err_file && exit 1; fi\n";
			  print SH "gunzip $ref_assembly.gz\n";
			  print SH "if [ ! -e $ref_assembly ]; then echo \"ERROR: unable to find $ref_assembly after supposedly unzipping\" >> $err_file && exit 1; fi\n";
			  $downloaded_already{$ref_assembly}++;
			}
		}
		
      # (4) DETERMINE THE MASH DISTANCE BETWEEN THE FDA-ARGOS ASSEMBLY AND THE REFERENCE GENOME.
	  while((my $ref_acc,my $ftp_path)=each(%{$other_assemblies{$species}}))
		{ my $last_bit = '';
		  if ($ftp_path =~ /^.+\/(.*?)$/) { $last_bit = $1; }
		  $ftp_path .= "/$last_bit"."_genomic.fna.gz";
		  my $ref_assembly = "$out_dir3/$last_bit.fa";
		  if (!(-e("$out_dir2/$sample_name/distance_to_$ref_acc.tsv")))
			{ print SH "cd $out_dir2/$sample_name\n";
			  print SH "$mash_path sketch -p $num_procs -o $sample_name $fda_assembly\n";
			  print SH "$mash_path sketch -p $num_procs -o $ref_acc $ref_assembly\n";
			  print SH "$mash_path dist -p $num_procs $out_dir2/$sample_name/$sample_name.msh $out_dir2/$sample_name/$ref_acc.msh > $out_dir2/$sample_name/distance_to_$ref_acc.tsv\n";
			  print SH "rm $out_dir2/$sample_name/$sample_name.msh $out_dir2/$sample_name/$ref_acc.msh\n";
			}
		}

	  # (5) RE-MAP THE ILLUMINA READS BACK TO THE FDA-ARGOS ASSEMBLY IN ORDER TO IDENTIFY REGIONS WHERE THERE IS POOR MAPPING.
	  if (!(-e("$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bam")))
		{ print SH "cd $out_dir2/$sample_name\n";
		  print SH "$minimap2_path -d $fda_assembly_idx.mmi $fda_assembly\n" unless (-e("$fda_assembly_idx.mmi")); # first index the FDA-ARGOS assembly
		  print SH "$minimap2_path -ax sr $fda_assembly -R '\@RG\\tID:group\\tSM:sample\\tPL:Illumina\\tLIB:lib\\tPU:unit' -t $num_procs <(zcat $fq_1_list) <(zcat $fq_2_list) | $samtools_path view -Shb - > $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.unsorted.bam\n";
		  print SH "java -jar $picard_path CleanSam INPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.unsorted.bam OUTPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.cleaned.bam TMP_DIR=$out_dir2/$sample_name\n";
		  print SH "rm $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.unsorted.bam\n";
		  print SH "mv $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.cleaned.bam $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.unsorted.bam\n";
		  print SH "java -jar $picard_path SortSam INPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.unsorted.bam OUTPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.sorted.bam SORT_ORDER=coordinate TMP_DIR=$out_dir2/$sample_name\n";
		  print SH "java -jar $picard_path MarkDuplicates INPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.sorted.bam OUTPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bam METRICS_FILE=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.metrics ASSUME_SORTED=true\n";
		  print SH "java -jar $picard_path BuildBamIndex INPUT=$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bam\n";
		  print SH "rm $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.unsorted.bam $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.sorted.bam $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.metrics\n";
		}
		
	  # (6) IDENTIFY LOW-QUALITY REGIONS FROM THE BAM - THOSE WHERE THERE IS EITHER NO COVERAGE, OR WHERE THE MOST COMMON NUCLEOTIDE REPRESENTS <$min_percent_depth% OF THE TOTAL DEPTH AT THAT POSITION. THIS INFORMATION IS OUTPUT IN A BED FILE, LATER USED TO CREATE A REPEAT-MASKED VERSION OF THE FDA-ARGOS ASSEMBLY.
	  # Code from https://raw.githubusercontent.com/martinghunt/bioinf-scripts/master/perl/bam_to_low_qual_mask.pl
	  if ( (!(-e("$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bed"))) and (-e("$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bam")) )
		{ open my $f_in, "samtools mpileup -aa $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bam |" or die $!;
		  open my $f_out, "> $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bed" or die $!;
		  
		  my $min_percent_depth = 99;
		  my $current_contig = "";
		  my $current_start = -1;
		  my $current_end = -1;
		  my $i = 0;
		  print STDERR "Gathering depths from mpileup.\n";
		  while (<$f_in>)
			{ my @a = split;
			  print STDERR "$a[0] $a[1]\n" if ($a[1] == 1 or $a[1] % 100000 == 0);
			  my $is_good = pileup_array_line_is_good($min_percent_depth, \@a);
			  if ($is_good and $current_contig eq "")
				{ next; }
			  elsif ($is_good and $current_contig ne "")
				{ print $f_out "$current_contig\t" . ($current_start - 1) . "\t$current_end\n";
				  $current_contig = "";
				  $current_start = -1;
				  $current_end = -1;
				}
			  elsif ($current_contig eq $a[0] and $current_end + 1 == $a[1])
				{ $current_end++;
				}
			  else
				{ unless ($current_contig eq "")
					{ print $f_out "$current_contig\t" . ($current_start - 1) . "\t$current_end\n"; }
				  $current_contig = $a[0];
				  $current_start = $a[1];
				  $current_end = $a[1];
				}
			}
		  if ($current_contig ne "")
			{ print $f_out "$current_contig\t" . ($current_start - 1) . "\t$current_end\n";
			}
		  print STDERR "finished\n";
		  close $f_in or die $!;
		  close $f_out or die $!;
		}
	  
	  # (7) MASK LOW-QUALITY REGIONS WITHIN THE FDA-ARGOS ASSEMBLY - THOSE WITH POOR COVERAGE.
	  my $masked_fda_assembly = "$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.masked_lowcov.fa";
	  if ( (!(-e($masked_fda_assembly))) and (-e("$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bed")) )
		{ print SH "$maskfasta_path -fi $fda_assembly -fo $masked_fda_assembly -bed $out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.bed -fullHeader\n";
		}
		
	  # (8) HOW MANY BASES HAVE BEEN MASKED AS LOW COVERAGE?
	  if ( (-e($masked_fda_assembly)) and (!(-e("$out_dir2/$sample_name/proportion_of_bases_masked_in_$assembly_accession_FDA_ASSEMBLY.lowcov.tsv"))) )
		{ my $masked_seq = '';
		  open my $f_out, "> $out_dir2/$sample_name/proportion_of_bases_masked_in_$assembly_accession_FDA_ASSEMBLY.lowcov.tsv" or die $!;
		  print $f_out "total bases\tnumber N\t% N\n";
		  open(IN,$masked_fda_assembly) or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  next if ($line =~ /^\>/);
			  $masked_seq .= uc($line);
			}
		  close(IN) or die $!;
		  my $len  = length($masked_seq);
		  my $no_n = $masked_seq =~ s/N/N/g;
		  my $pc_n = sprintf("%.4f",(($no_n/$len)*100));
		  print $f_out "$len\t$no_n\t$pc_n\n";
		  close $f_out or die $!;
		}
	  
	  # (8) RUN A NEGATIVE CONTROL: ALIGN THE ILLUMINA READS AGAINST THE MASKED FDA-ARGOS ASSEMBLY AND CALL VARIANTS. AS THE FORMER IS SOURCED FROM THE LATTER, WE EXPECT NO VARIANTS TO BE FOUND. ANY VARIANTS ARE THEREFORE EITHER INFORMATIC SNAFUS OR SEQUENCING ERRORS. THESE POSITIONS ALSO NEED TO BE MASKED.
	  if ( (-e($masked_fda_assembly)) and (!(-e("$out_dir4/$sample_name.bed"))) and (!(-e("$out_dir4/$sample_name.vcf"))) )
		{ print SH "$snippy_path --cpus $num_procs --outdir $out_dir4/$sample_name --prefix $sample_name --cleanup --ref $masked_fda_assembly --R1 $fq_1_list --R2 $fq_2_list\n";
		  print SH "mv $out_dir4/$sample_name/$sample_name.vcf $out_dir4/$sample_name.vcf\n";
		  print SH "mv $out_dir4/$sample_name/$sample_name.bed $out_dir4/$sample_name.bed\n";
		  print SH "rm -r $out_dir4/$sample_name\n";
		}

	  # (9) MASK SEQUENCING ERRORS IDENTIFIED AFTER CALLING VARIANTS AGAINST THE (LOW COVERAGE) MASKED FDA-ARGOS ASSEMBLY USING THE ILLUMINA READS.
	  my $masked_fda_assembly2 = "$out_dir2/$sample_name/$assembly_accession_FDA_ASSEMBLY.masked_lowcov_and_miscall.fa";
	  if ( (!(-e($masked_fda_assembly2))) and (-e("$out_dir4/$sample_name.bed")) )
		{ print SH "$maskfasta_path -fi $masked_fda_assembly -fo $masked_fda_assembly2 -bed $out_dir4/$sample_name.bed -fullHeader\n";
		}
	  
	  # (10) HOW MANY BASES HAVE BEEN MASKED DUE TO DISCORDANT ILLUMINA CALLS (THAT IS, THE CALL AT THAT POSITION IS DISCORDANT RELATIVE TO THE ONE MADE BY PACBIO)?
	  if ( (-e($masked_fda_assembly)) and (-e($masked_fda_assembly2)) and (-e("$out_dir2/$sample_name/proportion_of_bases_masked_in_$assembly_accession_FDA_ASSEMBLY.lowcov.tsv")) and (!(-e("$out_dir2/$sample_name/proportion_of_bases_masked_in_$assembly_accession_FDA_ASSEMBLY.miscall.tsv"))) )
		{ open my $f_out, "> $out_dir2/$sample_name/proportion_of_bases_masked_in_$assembly_accession_FDA_ASSEMBLY.miscall.tsv" or die $!;
		  print $f_out "total bases\tnumber N\t% N\n";
		  my $masked_seq_lowcov = '';
		  open(IN,$masked_fda_assembly) or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  next if ($line =~ /^\>/);
			  $masked_seq_lowcov .= uc($line);
			}
		  close(IN) or die $!;
		  my $masked_seq_miscall = '';
		  open(IN,$masked_fda_assembly2) or die $!;
		  while(<IN>)
			{ my $line = $_; chomp($line);
			  next if ($line =~ /^\>/);
			  $masked_seq_miscall .= uc($line);
			}
		  close(IN) or die $!;
		  my $len   = length($masked_seq_miscall);
		  my $no_n1 = $masked_seq_lowcov  =~ s/N/N/g;
		  my $no_n2 = $masked_seq_miscall =~ s/N/N/g;
		  my $no_n  = $no_n2-$no_n1;
		  my $pc_n  = sprintf("%.4f",(($no_n/$len)*100));
		  print $f_out "$len\t$no_n\t$pc_n\n";
		  close $f_out or die $!;
		}
	  
	  # (11) ALIGN THE (DOUBLY-MASKED) FDA-ARGOS ASSEMBLY TO THE REFERENCE GENOME, SO CALLING VARIANTS RELATIVE TO IT. WE WILL DO THIS TWICE, USING NUCMER AND PAFTOOLS (WITH A RANGE OF PARAMETERS FOR EACH), STANDARDISING THE REPRESENTATION OF EACH VCF USING PRE.PY. WE WILL THEN TAKE THE SET OF CONSENSUS CALLS TO BE THE 'TRUTH SET'. THESE TRUTH POSITIONS WILL BE WHAT WE WILL LOOK FOR WHEN ALIGNING THE ORIGINAL ILLUMINA READS TO THE REFERENCE GENOME.
	  # IMPORTANT: before using pre.py we must manually change "convert_gvcf_to_vcf=args.convert_gvcf" to "convert_gvcf_to_vcf=False" to avoid an error. See https://github.com/Illumina/hap.py/issues/106.
	  if (!(-d("$out_dir5/$sample_name"))) { print SH "mkdir $out_dir5/$sample_name\n"; }
	  print SH "cd $out_dir5/$sample_name\n";
	  while((my $ref_acc,my $ftp_path)=each(%{$other_assemblies{$species}}))
		{ my $last_bit = '';
		  if ($ftp_path =~ /^.+\/(.*?)$/) 	 { $last_bit = $1; }
		  $ftp_path .= "/$last_bit"."_genomic.fna.gz";
		  my $ref_assembly = "$out_dir3/$last_bit.fa";
		  my $ref_assembly_root_name = "$out_dir3/$last_bit";
		  
		  if (!(-e("$ref_assembly.fai"))) { print SH "$samtools_path faidx $ref_assembly\n"; }
		  
		  my @vcfs = ();

		  # ALIGN USING NUCMER, VARYING THE -c (MIN. CLUSTER LENGTH), -g (MAX. GAP SIZE) AND -b (BREAK LENGTH) PARAMETERS
		  for(my $c=25;$c<=200;$c+=10) # default: 65
			{ for(my $g=90;$g<=900;$g+=90) # default: 90
				{ for(my $b=200;$b<=400;$b+=100) # default: 200
					{ my $aligner = "nucmer.c$c.g$g.b$b";
						{ if ( (-e($masked_fda_assembly2)) and (!(-e("$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.regularised.vcf.gz"))) )
							{ print SH "$nucmer_path -c $c -g $g -p $sample_name $ref_assembly $masked_fda_assembly2\n";
							  print SH "$dnadiff_path -p $sample_name -d $out_dir5/$sample_name/$sample_name.delta\n";
							  print SH "$showsnps_path -Clr -x 1 -T $out_dir5/$sample_name/$sample_name.delta > $out_dir5/$sample_name/$sample_name.snps.filter\n";
							  print SH "rm $out_dir5/$sample_name/$sample_name.1coords $out_dir5/$sample_name/$sample_name.delta $out_dir5/$sample_name/$sample_name.1delta $out_dir5/$sample_name/$sample_name.mcoords $out_dir5/$sample_name/$sample_name.mdelta $out_dir5/$sample_name/$sample_name.qdiff $out_dir5/$sample_name/$sample_name.rdiff\n";
							  print SH "rm $out_dir5/$sample_name/$sample_name.snps $out_dir5/$sample_name/$sample_name.report\n";
							  print SH "rm $out_dir5/$sample_name/$sample_name.unref $out_dir5/$sample_name/$sample_name.unqry\n"; # NOTE: these files do not always get created by $dnadiff_path, so this step may report a warning ("unable to delete non-existent file...")
							  print SH "python $MUMmerSNPs2VCF_path $out_dir5/$sample_name/$sample_name.snps.filter $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf\n";
							  print SH "rm $out_dir5/$sample_name/$sample_name.snps.filter\n";
							  print SH "$bgzip_path $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf\n";
							  print SH "$tabix_path $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz\n";
							  print SH "$prepy_path --leftshift --decompose -r $ref_assembly $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.regularised.vcf.gz\n";
							  print SH "rm $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz.tbi\n";
							}
						  if (-e($masked_fda_assembly2))
							{ push(@vcfs,"$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.regularised.vcf.gz"); }
						}
					}
				}
			}

		  # ALIGN USING PAFTOOLS, VARYING THE -l (MIN. ALIGNMENT LENGTH TO COMPUTE COVERAGE & CALL VARIANTS) PARAMETER
		  for(my $l=25;$l<=200;$l+=25)
			{ my $aligner = "paftools_l$l";
			  if ( (-e($masked_fda_assembly2)) and (!(-e("$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.regularised.vcf.gz"))) )
				{ print SH "$minimap2_path -t $num_procs -c -cx asm5 --cs $ref_assembly $masked_fda_assembly2 | sort -k6,6 -k8,8n | $paftools_path call -f $ref_assembly -s $sample_name -q 20 -l $l -L $l - > $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf\n";
				  print SH "$bgzip_path $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf\n";
				  print SH "$tabix_path $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz\n";
				  print SH "$prepy_path --leftshift --decompose -r $ref_assembly $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.regularised.vcf.gz\n";
				  print SH "rm $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.vcf.gz.tbi\n";
				}
			  if (-e($masked_fda_assembly2))
				{ push(@vcfs,"$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.$aligner.regularised.vcf.gz"); }
			}
			
		  # INTERSECT NUCMER AND PAFTOOLS VCFs TO CREATE AN INTERSECT (CONSENSUS) VCF OF 'TRUTH' SNPs (THOSE IDENTICALLY CALLED BY ALL INSTANCES OF NUCMER AND PAFTOOLS) AND A VCF OF 'AMBIGUOUS' POSITIONS (THOSE CALLED ONLY BY SOME COMBINATION OF ALIGNERS AND/OR PARAMETERS, BUT NOT ALL, AND WHICH WILL LATER BE MASKED).
		  my $total_no_of_vcfs = @vcfs; my $all_but_one_vcf = $total_no_of_vcfs-1;
		  my $vcfs_line = join(" ",@vcfs);
		  if (-e($masked_fda_assembly2))
			{ if (!(-e("$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.consensus_variants.vcf")))
				{ print SH "$bcftools_path isec --threads $num_procs -n=$total_no_of_vcfs -w1 -o $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.consensus_variants.vcf $vcfs_line\n";
				}
			  if (!(-e("$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.non_consensus_variants.vcf")))
				{ print SH "$bcftools_path isec --threads $num_procs -n-$all_but_one_vcf -w1 -o $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.non_consensus_variants.vcf $vcfs_line\n";
				}
			}
			
		  # CREATE A UNION VCF FROM THE SET OF NUMCER AND PAFTOOLS OUTPUT
		  if (-e($masked_fda_assembly2))
			{ if (!(-e("$out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.union_variants.vcf")))
				{ print SH "$bcftools_path isec --threads $num_procs -n+1 -w1 -o $out_dir5/$sample_name/$sample_name.varcall_relative_to.$last_bit.union_variants.vcf $vcfs_line\n";
				}
			}
		}
	}
close(SH) or die $!; close(ERR) or die $!; close(OUT) or die $!;
exit 1;

sub pileup_array_line_is_good
	{ my $min_depth_percent = shift;
	  my $array = shift;
	  return 0 if ($array->[3] == 0);
	  my %nuc_counts = (A => 0, G => 0, C => 0, T => 0);
	  my $total = 0;
	  for my $n (split(//, $array->[4]))
		{ if (exists $nuc_counts{uc $n})
			{ $nuc_counts{uc $n}++;
			  $total++;
			}
		}
	  return 0 if $total < 1;
	  return 1 if 1 == scalar keys %nuc_counts;
	  my @counts = reverse sort {$a <=> $b} values %nuc_counts;
	  return (100 * $counts[0] / $total) >= $min_depth_percent;
	}