use strict;
use warnings;

# REQUIREMENTS
my $root    = '/home/ndm.local/steveb';
my $path    = "$root/test_suite";
my $vcf_dir = "$path/pipeline_VCFs"; # from 4.run_pipelines.pl
my $fatal   = 0;
if (!(-d($vcf_dir))) { $fatal++; print "error: cannot find $vcf_dir\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my @aligners  = (qw/bowtie2 bwa-mem bwa+stampy bwa-sw gassst gem hisat2 minimap2 mosaik ngm smalt snap stampy yara breseq snippy spandx speedseq/);
my @callers   = (qw/deepvariant freebayes gatk lofreq mpileup octopus pilon platypus snver snvsniffer solsnp strelka varscan breseq snippy spandx speedseq/);
my %aligners  = map {$_ => 1} @aligners;
my %callers   = map {$_ => 1} @callers;
my @shortlist = ('FDAARGOS_338','FDAARGOS_536','FDAARGOS_598','FDAARGOS_687','FDAARGOS_700','FDAARGOS_751');
my %shortlist = map {$_ => 1} @shortlist;

# OUTPUT
my $out_file1 = "$path/summary_of_pipeline_performance.tsv";
my $out_file2 = "$path/characteristics_of_FP_SNPs_using_real_data.tsv";
my $out_file3 = "$path/characteristics_of_TP_SNPs_using_real_data.tsv";
my $out_file4 = "$path/characteristics_of_FP_indels_using_real_data.tsv";
my $out_file5 = "$path/characteristics_of_TP_indels_using_real_data.tsv";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!; open(OUT3,'>',$out_file3) or die $!; open(OUT4,'>',$out_file4) or die $!; open(OUT5,'>',$out_file5) or die $!;
print OUT1 "Sample ID\tPipeline\tAligner\tCaller\tNo. of SNPs in truth set\tPrecision (SNPs)\tRecall (SNPs)\tF-score (SNPs)\tNo. of true positives (SNPs)\tNo. of false positives (SNPs)\tNo. of false negatives (SNPs)\tNo. of indels in truth set\tPrecision (indels)\tRecall (indels)\tF-score (indels)\tNo. of true positives (indels)\tNo. of false positives (indels)\tNo. of false negatives (indels)\n";
print OUT2 "Sample ID\tAligner\tCaller\tAligner/caller\tChr\tPos\tReference base\tVariant base\tVariant call quality\tRead depth (no. of reads mapping to this position)\tNo. of reads mapping to forward strand\tNo. of reads mapping to reverse strand\tNo. of reads supporting reference allele\tNo. of reads supporting variant allele\tProportion of reads supporting the variant allele\tNo. of reads supporting the variant allele and mapped to its left\tNo. of reads supporting the variant allele and mapped to its right\tAverage quality of variant-supporting reads\tDistance to nearest SNP (bp)\tDistance to nearest indel (bp)\n";
print OUT3 "Sample ID\tAligner\tCaller\tAligner/caller\tChr\tPos\tReference base\tVariant base\tVariant call quality\tRead depth (no. of reads mapping to this position)\tNo. of reads mapping to forward strand\tNo. of reads mapping to reverse strand\tNo. of reads supporting reference allele\tNo. of reads supporting variant allele\tProportion of reads supporting the variant allele\tNo. of reads supporting the variant allele and mapped to its left\tNo. of reads supporting the variant allele and mapped to its right\tAverage quality of variant-supporting reads\tDistance to nearest SNP (bp)\tDistance to nearest indel (bp)\n";
print OUT4 "Sample ID\tAligner\tCaller\tAligner/caller\tChr\tPos\tReference base\tVariant base\tVariant call quality\tRead depth (no. of reads mapping to this position)\tNo. of reads mapping to forward strand\tNo. of reads mapping to reverse strand\tNo. of reads supporting reference allele\tNo. of reads supporting variant allele\tProportion of reads supporting the variant allele\tNo. of reads supporting the variant allele and mapped to its left\tNo. of reads supporting the variant allele and mapped to its right\tAverage quality of variant-supporting reads\tDistance to nearest SNP (bp)\tDistance to nearest indel (bp)\n";
print OUT5 "Sample ID\tAligner\tCaller\tAligner/caller\tChr\tPos\tReference base\tVariant base\tVariant call quality\tRead depth (no. of reads mapping to this position)\tNo. of reads mapping to forward strand\tNo. of reads mapping to reverse strand\tNo. of reads supporting reference allele\tNo. of reads supporting variant allele\tProportion of reads supporting the variant allele\tNo. of reads supporting the variant allele and mapped to its left\tNo. of reads supporting the variant allele and mapped to its right\tAverage quality of variant-supporting reads\tDistance to nearest SNP (bp)\tDistance to nearest indel (bp)\n";

# DETERMINE SUMMARY STATISTICS FOR EACH FDA-ARGOS SAMPLE, AND PIPELINE
opendir(DIR,$vcf_dir) or die $!;
my @sample_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
my $sample_ids_seen = 0; my $sample_ids_total = @sample_ids;
foreach my $sample_id (@sorted_sample_ids)
	{ $sample_ids_seen++;
	  print "$sample_ids_seen of $sample_ids_total\n";
	  next if (($sample_id eq '.') or ($sample_id eq '..'));
	  next if (!(-d("$vcf_dir/$sample_id")));
	  next if (!(exists($shortlist{$sample_id})));
	  opendir(DIR,"$vcf_dir/$sample_id");
	  my @contents = readdir(DIR);
	  closedir(DIR) or die $!;
	  my @sorted_contents = sort {$a cmp $b} @contents;
	  my $pipelines_seen = 0; my $pipelines_total = @sorted_contents;
	  foreach my $content (@sorted_contents)
		{ next if (($content eq '.') or ($content eq '..'));
		  $pipelines_seen++;
		  my $aligner = ''; my $caller = '';
		  if ($content =~ /^happy\-(.*?)\.(.*?)$/)
			{ $aligner = $1; $caller = $2; }
		  next if ( (!(exists($aligners{$aligner}))) or (!(exists($callers{$caller}))) );
		  
		  if (-d("$vcf_dir/$sample_id/$content"))
			{ if ($content =~ /^happy\-(.*?)\.(.*?)$/)
				{ my $aligner = $1; my $caller = $2;
				  my $summary = "$vcf_dir/$sample_id/$content/$sample_id.summary.csv";
				  my $hap_vcf = "$vcf_dir/$sample_id/$content/$sample_id.vcf.gz";
				  my $var_vcf = "$vcf_dir/$sample_id/$sample_id.$aligner.$caller.vcf.gz";
				  next if ( (!(-e($summary))) or (!(-e($hap_vcf))) or (!(-e($var_vcf))) );
				  
				  # OBTAIN SUMMARY STATISTICS FOR EACH PIPELINE'S PERFORMANCE
				  my $truth_set_total_IND = 0; my $precision_IND = 0; my $recall_IND = 0; my $f_score_IND = 0; my $no_tp_IND = 0; my $no_fp_IND = 0; my $no_fn_IND = 0;
				  my $truth_set_total_SNP = 0; my $precision_SNP = 0; my $recall_SNP = 0; my $f_score_SNP = 0; my $no_tp_SNP = 0; my $no_fp_SNP = 0; my $no_fn_SNP = 0;
				  open(IN,$summary) or die $!;
				  while(<IN>)
					{ my $line = $_; chomp($line);
					  my @line = split(/\,/,$line);
					  my $truth_set_total = $line[2]; my $precision = $line[11]; my $recall = $line[10]; my $f_score = $line[13]; my $no_tp = $line[3]; my $no_fp = $line[6]; my $no_fn = $line[4];
					  if ($. == 3) # INDEL PASS
						{ $truth_set_total_IND = $truth_set_total; $precision_IND = $precision; $recall_IND = $recall; $f_score_IND = $f_score; $no_tp_IND = $no_tp; $no_fp_IND = $no_fp; $no_fn_IND = $no_fn; }
					  elsif ($. == 5) # SNP PASS
						{ $truth_set_total_SNP = $truth_set_total; $precision_SNP = $precision; $recall_SNP = $recall; $f_score_SNP = $f_score; $no_tp_SNP = $no_tp; $no_fp_SNP = $no_fp; $no_fn_SNP = $no_fn; }
					}
				  close(IN) or die $!;
				  $precision_IND = sprintf("%.4f",$precision_IND) unless ($precision_IND !~ /\d+/); $recall_IND = sprintf("%.4f",$recall_IND) unless ($recall_IND !~ /\d+/); $f_score_IND = sprintf("%.4f",$f_score_IND) unless ($f_score_IND !~ /\d+/);
				  $precision_SNP = sprintf("%.4f",$precision_SNP) unless ($precision_SNP !~ /\d+/); $recall_SNP = sprintf("%.4f",$recall_SNP) unless ($recall_SNP !~ /\d+/); $f_score_SNP = sprintf("%.4f",$f_score_SNP) unless ($f_score_SNP !~ /\d+/);
				  if ($f_score_SNP !~ /\d+/) { $f_score_SNP = 0; }
				  if ($f_score_IND !~ /\d+/) { $f_score_IND = 0; }
				  print OUT1 "$sample_id\t$aligner/$caller\t$aligner\t$caller\t$truth_set_total_SNP\t$precision_SNP\t$recall_SNP\t$f_score_SNP\t$no_tp_SNP\t$no_fp_SNP\t$no_fn_SNP\t$truth_set_total_IND\t$precision_IND\t$recall_IND\t$f_score_IND\t$no_tp_IND\t$no_fp_IND\t$no_fn_IND\n";
				  
				  # FOR EACH HIGH-PERFORMING PIPELINE, IDENTIFY, FROM THE HAPPY VCF, THE LOCATION OF TRUE AND FALSE POSITIVE CALLS.
				  my %fps_IND = (); my %tps_IND = ();
				  my %fps_SNP = (); my %tps_SNP = ();
				  open(IN,"gunzip -c $hap_vcf | ") or die $!;
				  while(<IN>)
					{ my $line = $_; chomp($line);
					  next if ($line =~ /^#/);
					  my @line = split(/\t/,$line);
					  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4]; my $filter = $line[6]; my $target = $line[9]; my $query = $line[10];
					  next if ($filter ne '.'); # i.e. this is a pipeline-specific filter code, carried through into this hap.py VCF; it was never a PASS value in the original VCF to begin with
					  if (($target =~ /^.*?\:.*?\:.*?$/) and ($query =~ /^.*?\:FP\:.*?$/)) # middle position of $target is usually NOCALL; this might not get all of them, because there can be multiple calls per line?
						{ if (((length($ref)==1)) and (length($alt)==1))
							{ $fps_SNP{"$chr|$pos|$ref|$alt"}++; } # unless ($f_score_SNP < 0.95); }
						  else
							{ $fps_IND{"$chr|$pos|$ref|$alt"}++; } # unless ($f_score_IND < 0.80); }
						  print "$sample_id\t$aligner/$caller\t$aligner\t$caller\t$chr $pos $ref/$alt\n";
						}
					  if (($target =~ /^.*?\:TP\:.*?$/) and ($query =~ /^.*?\:TP\:.*?$/))
						{ if (((length($ref)==1)) and (length($alt)==1))
							{ $tps_SNP{"$chr|$pos|$ref|$alt"}++; } # unless ($f_score_SNP < 0.95); }
						  else
							{ $tps_IND{"$chr|$pos|$ref|$alt"}++; } # unless ($f_score_IND < 0.80); }
						}
					}
				  close(IN) or die $!;
				  
				  # DETERMINE THE LOCATIONS OF EACH SNP AND INDEL, SO THAT WE CAN DETERMINE THE DISTANCE TO THE NEAREST
				  my %snps = (); my %indels = ();
				  open(IN,"gunzip -c $var_vcf | ") or die $!;
				  while(<IN>)
					{ my $line = $_; chomp($line);
					  next if ($line =~ /^#/);
					  my @line = split(/\t/,$line);
					  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4]; my $filter = $line[6]; my $info = $line[7];
					  next if ($alt eq '.');
					  my $len_ref = length($ref); my $len_alt = length($alt);
					  if ($caller eq 'breseq')
						{ if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'deepvariant')
						{ next if ($filter ne 'PASS');
						  if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'freebayes')
						{ if ($info =~ /^.*?\;TYPE\=(.*?)\;.*?$/)
							{ my $type = $1;
							  if (($type eq 'del') or ($type eq 'ins'))
								{ push(@{$indels{$chr}},$pos); }
							  elsif (($type eq 'snp') or ($type eq 'mnp') or ($type eq 'complex')) # 'complex' vars are multiple SNPs in a row where length($ref) == length($alt)
								{ push(@{$snps{$chr}},$pos); }
							}
						}
					  elsif ($caller eq 'gatk')
						{ if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'lofreq') # only reports SNPs
						{ next if ($filter ne 'PASS');
						  push(@{$snps{$chr}},$pos);
						}
					  elsif ($caller eq 'mpileup')
						{ if ($info =~ /INDEL/)
							{ push(@{$indels{$chr}},$pos); }
						  else
							{ push(@{$snps{$chr}},$pos); }
						}
					  elsif ($caller eq 'octopus')
						{ next if ($filter ne 'PASS');
						  if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'pilon')
						{ next if ($filter ne 'PASS');
						  if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'platypus')
						{ next if ($filter ne 'PASS');
						  if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'snippy')
						{ if ($info =~ /^AB\=\d+\;AO\=\d+\;DP\=\d+\;QA\=\d+\;QR\=\d+\;RO\=\d+\;TYPE\=(.*?)$/)
							{ my $type = $1;
							  if (($type eq 'del') or ($type eq 'ins'))
								{ push(@{$indels{$chr}},$pos); }
							  elsif (($type eq 'snp') or ($type eq 'mnp') or ($type eq 'complex')) # 'complex' vars are multiple SNPs in a row where length($ref) == length($alt)
								{ push(@{$snps{$chr}},$pos); }
							}
						}
					  elsif ($caller eq 'snver') # only reports SNPs
						{ push(@{$snps{$chr}},$pos);
						}
					  elsif ($caller eq 'snvsniffer')
						{ if ($info =~ /^DP\=\d+\;VT\=(.*?)$/)
							{ my $type = $1;
							  if (($type =~ 'DEL') or ($type =~ 'INS'))
								{ push(@{$indels{$chr}},$pos); }
							  elsif ($type =~ 'SNP')
								{ push(@{$snps{$chr}},$pos); }
							}
						}
					  elsif ($caller eq 'solsnp') # only reports SNPs
						{ push(@{$snps{$chr}},$pos);
						}
					  elsif ($caller eq 'spandx')
						{ next if ($filter ne 'PASS');
						  if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'speedseq')
						{ if ($info =~ /^.*?\;TYPE\=(.*?)\;.*?$/)
							{ my $type = $1;
							  if (($type eq 'del') or ($type eq 'ins'))
								{ push(@{$indels{$chr}},$pos); }
							  elsif (($type eq 'snp') or ($type eq 'mnp') or ($type eq 'complex')) # 'complex' vars are multiple SNPs in a row where length($ref) == length($alt)
								{ push(@{$snps{$chr}},$pos); }
							}
						}
					  elsif ($caller eq 'strelka')
						{ next if ($filter ne 'PASS');
						  if (($len_ref == 1) and ($len_alt == 1))
							{ push(@{$snps{$chr}},$pos); }
						  else
							{ push(@{$indels{$chr}},$pos); }
						}
					  elsif ($caller eq 'varscan') # only reports SNPs
						{ next if ($filter ne 'PASS');
						  push(@{$snps{$chr}},$pos);
						}
					}
				  close(IN) or die $!;
				  
				  # SORT THE SNP AND INDEL POSITION ARRAYS
				  my %sorted_snps = ();
				  while((my $chr,my $irrel)=each(%snps))
					{ my @snp_pos = @{$snps{$chr}};
					  my @sorted_snp_pos = sort {$a <=> $b} @snp_pos;
					  $sorted_snps{$chr} = \@sorted_snp_pos;
					}
				  my %sorted_indels = ();
				  while((my $chr,my $irrel)=each(%indels))
					{ my @indel_pos = @{$indels{$chr}};
					  my @sorted_indel_pos = sort {$a <=> $b} @indel_pos;
					  $sorted_indels{$chr} = \@sorted_indel_pos;
					}
				  
				  # DETERMINE THE INFO LINES FOR THE TRUE AND FALSE POSITIVE CALLS - SO THAT WE CAN START LOOKING AT THEIR COMMON CHARACTERISTICS.
				  open(IN,"gunzip -c $var_vcf | ") or die $!;
				  while(<IN>)
					{ my $line = $_; chomp($line);
					  next if ($line =~ /^#/);
					  my @line = split(/\t/,$line);
					  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4]; my $qual = $line[5]; my $filter = $line[6]; my $info = $line[7]; my $format = $line[9];
					  my $len_ref = length($ref); my $len_alt = length($alt);
					  if ($qual =~ /^\.$/) { $qual = 'NA'; }
					  next if (($ref =~ /\,/) or ($alt =~ /\,/));
					  for(my $x=0;$x<=1;$x++)
						{ my $skipped = 0;
						  my %hash = ();
						  for(my $y=0;$y<=1;$y++)
							{ if ($x == 0) # SNPs
								{ if 	($y == 0) { %hash = %fps_SNP; }
								  elsif ($y == 1) { %hash = %tps_SNP; }
								}
							  elsif ($x == 1) # indels
								{ if 	($y == 0) { %hash = %fps_IND; }
								  elsif ($y == 1) { %hash = %tps_IND; }
								}
							  print "reading $sample_id ($sample_ids_seen of $sample_ids_total) $aligner/$caller ($pipelines_seen of $pipelines_total) $chr|$pos|$ref|$alt\n";
							  next if (!(exists($hash{"$chr|$pos|$ref|$alt"})));
							  if ($x == 0) # SNPs
								{ my $allele_balance = 'NA'; my $dp = 'NA'; my $sum_ref_qual = 'NA'; my $sum_alt_qual = 'NA';
								  my $reads_mapped_to_ref = 'NA'; my $reads_mapped_to_alt = 'NA';
								  my $reads_mapped_to_fwd = 'NA'; my $reads_mapped_to_rev = 'NA';
								  my $reads_mapped_left = 'NA'; my $reads_mapped_right = 'NA';
								  if ($caller eq 'breseq')
									{ if ($info =~ /^.*?\;DP\=(.*?)$/) { $dp 			 = $1; }
									  if ($info =~ /^AF\=(.*?)\;.*?$/) { $allele_balance = $1; }
									  if ($info =~ /^.*?\;AD\=(.*?)\;.*?$/)
										{ my $ad = $1;
										  $reads_mapped_to_ref = $dp-$ad;
										  $reads_mapped_to_alt = $ad;
										}
									}
								  elsif ($caller eq 'deepvariant')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref != $len_alt) and ($len_ref != 1)) { $skipped++; }
									  next if (($len_ref != $len_alt) and ($len_ref != 1)); # i.e. not a SNP
									  if ($format =~ /^.*?\:.*?\:(\d+)\:(.*?)\:(.*?)\:.*?$/)
										{ $dp = $1; my $ad = $2; $allele_balance = $3;
										  if ($ad =~ /^(\d+)\,(\d+)$/)
											{ $reads_mapped_to_ref = $1;
											  $reads_mapped_to_alt = $2;
											}
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $format\n"; exit 1; }
									}
								  elsif ($caller eq 'freebayes')
									{ my $type = ''; if ($info =~ /^.*?TYPE\=(.*?)\;.*?$/) { $type = $1; }
									  if ($type ne 'snp') { $skipped++; }
									  next if ($type ne 'snp');
									  my $saf = ''; my $srf = ''; my $sar = ''; my $srr = '';
									  if ($info =~ /^.*?\;DP\=(\d+)\;.*?$/)   { $dp  = $1; }
									  if ($info =~ /^AB\=(.*?)\;.*?$/)   	  { $allele_balance = $1; $allele_balance = 1-$allele_balance; }
									  if ($info =~ /^.*?\;AO\=(\d+)\;.*?$/)   { $reads_mapped_to_alt = $1; }
									  if ($info =~ /^.*?\;QA\=(\d+)\;.*?$/)   { $sum_alt_qual = $1; }
									  if ($info =~ /^.*?\;QR\=(\d+)\;.*?$/)   { $sum_ref_qual = $1; }
									  if ($info =~ /^.*?\;RO\=(\d+)\;.*?$/)   { $reads_mapped_to_ref = $1; }
									  if ($info =~ /^.*?\;RPR\=(\d+)\;.*?$/)  { $reads_mapped_right  = $1; }
									  if ($info =~ /^.*?\;RPL\=(\d+)\;.*?$/)  { $reads_mapped_left   = $1; }
									  if ($info =~ /^.*?\;SAF\=(\d+)\;.*?$/)  { $saf = $1; }
									  if ($info =~ /^.*?\;SRF\=(\d+)\;.*?$/)  { $srf = $1; }
									  if ($info =~ /^.*?\;SAR\=(\d+)\;.*?$/)  { $sar = $1; }
									  if ($info =~ /^.*?\;SRR\=(\d+)\;.*?$/)  { $srr = $1; }							  
									  $reads_mapped_to_fwd = $saf+$srf;
									  $reads_mapped_to_rev = $sar+$srr;
									}
								  elsif ($caller eq 'gatk')
									{ if (($len_ref != $len_alt) and ($len_ref != 1)) { $skipped++; }
									  next if (($len_ref != $len_alt) and ($len_ref != 1)); # i.e. not a SNP
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($format =~ /^.*?\:(.*?)\:\d+\:.*?\:.*?$/)
										{ my $ad = $1;
										  if ($ad =~ /^(\d+)\,(\d+)$/)
											{ $reads_mapped_to_ref = $1;
											  $reads_mapped_to_alt = $2;
											  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt) unless (($reads_mapped_to_ref+$reads_mapped_to_alt)==0);
											}
										}
									}
								  elsif ($caller eq 'lofreq')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if ($info =~ /^DP\=(\d+)\;.*?$/) 		{ $dp 			  = $1; }
									  if ($info =~ /^.*?\;AF\=(.*?)\;.*?$/) { $allele_balance = $1; }
									  my $dp4 = '';
									  if ($info =~ /^.*?\;DP4\=(.*?)$/) { $dp4 = $1; }
									  if ($dp4 =~ /^(\d+)\,(\d+)\,(\d+)\,(\d+)$/)
										{ my $fwd_ref = $1; my $rev_ref = $2; my $fwd_alt = $3; my $rev_alt = $4;
										  $reads_mapped_to_alt = $fwd_alt+$rev_alt;
										  $reads_mapped_to_ref = $fwd_ref+$rev_ref;
										  $reads_mapped_to_fwd = $fwd_ref+$fwd_alt;
										  $reads_mapped_to_rev = $rev_ref+$rev_alt;
										}
									}
								  elsif ($caller eq 'mpileup')
									{ if ($info =~ /INDEL/) { $skipped++; }
									  next if ($info =~ /INDEL/);
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/) { $dp = $1; }
									  my $dp4 = '';
									  if ($info =~ /^.*?DP4\=(.*?)\;.*?$/) { $dp4 = $1; }
									  if ($dp4 =~ /^(\d+)\,(\d+)\,(\d+)\,(\d+)$/)
										{ my $fwd_ref = $1; my $rev_ref = $2; my $fwd_alt = $3; my $rev_alt = $4;
										  $reads_mapped_to_alt = $fwd_alt+$rev_alt;
										  $reads_mapped_to_ref = $fwd_ref+$rev_ref;
										  $reads_mapped_to_fwd = $fwd_ref+$fwd_alt;
										  $reads_mapped_to_rev = $rev_ref+$rev_alt;
										  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt);
										}
									}
								  elsif ($caller eq 'octopus')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref != $len_alt) and ($len_ref != 1)) { $skipped++; }
									  next if (($len_ref != $len_alt) and ($len_ref != 1)); # i.e. not a SNP
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/) { $dp = $1; }
									}
								  elsif ($caller eq 'pilon')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref != $len_alt) and ($len_ref != 1)) { $skipped++; }
									  next if (($len_ref != $len_alt) and ($len_ref != 1)); # i.e. not a SNP
									  if ($info =~ /^DP\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($info =~ /^.*?\;AF\=(.*?)$/) { $allele_balance = $1; }
									  my $num_a = 0; my $num_c = 0; my $num_g = 0; my $num_t = 0;
									  if ($info =~ /^.*?BC\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;.*?$/) { $num_a = $1; $num_c = $2; $num_g = $3; $num_t = $4; }
									  if ($ref eq 'A') { $reads_mapped_to_ref = $num_a; }
									  if ($ref eq 'C') { $reads_mapped_to_ref = $num_c; }
									  if ($ref eq 'G') { $reads_mapped_to_ref = $num_g; }
									  if ($ref eq 'T') { $reads_mapped_to_ref = $num_t; }
									  if ($alt eq 'A') { $reads_mapped_to_alt = $num_a; }
									  if ($alt eq 'C') { $reads_mapped_to_alt = $num_c; }
									  if ($alt eq 'G') { $reads_mapped_to_alt = $num_g; }
									  if ($alt eq 'T') { $reads_mapped_to_alt = $num_t; }
									}
								  elsif ($caller eq 'platypus')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref != $len_alt) and ($len_ref != 1)) { $skipped++; }
									  next if (($len_ref != $len_alt) and ($len_ref != 1)); # i.e. not a SNP
									  if ($info =~ /^.*?TC\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($info =~ /^.*?NF\=(\d+)\;.*?$/) { $reads_mapped_to_fwd = $1; }
									  if ($info =~ /^.*?NR\=(\d+)\;.*?$/) { $reads_mapped_to_rev = $1; }
									  if ($info =~ /^.*?TR\=(\d+)\;.*?$/) { $reads_mapped_to_alt = $1; }
									  $reads_mapped_to_ref = $dp-$reads_mapped_to_alt;
									  $allele_balance = $reads_mapped_to_alt/$dp;
									}
								  elsif ($caller eq 'snippy')
									{ if ($info =~ /^AB\=(\d+)\;AO\=(\d+)\;DP\=(\d+)\;QA\=(\d+)\;QR\=(\d+)\;RO\=(\d+)\;TYPE\=(.*?)$/)
										{ my $type = $7;
										  if ($type ne 'snp') { $skipped++; }
										  next if ($type ne 'snp');
										  $allele_balance = $1; $reads_mapped_to_alt = $2; $dp = $3; $sum_alt_qual = $4; $sum_ref_qual = $5; $reads_mapped_to_ref = $6; 
										  $allele_balance = 1-$allele_balance;								  
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $info\n"; exit 1; }
									}
								  elsif ($caller eq 'snver')
									{ if ($format =~ /^.*?\:.*?\:(\d+)\:(\d+)\:(\d+)\:(\d+)$/)
										{ my $ac1 = $1; my $ac2 = $2; my $rc1 = $3; my $rc2 = $4;
										  if (($rc1+$rc2)>0) { $allele_balance = ($ac1+$ac2)/($ac1+$ac2+$rc1+$rc2); } else { $allele_balance = 1; }
										  $dp = $ac1+$ac2+$rc1+$rc2;
										  $reads_mapped_to_alt = $ac1+$ac2;
										  $reads_mapped_to_ref = $rc1+$rc2;
										  $reads_mapped_to_fwd = $ac1+$rc1;
										  $reads_mapped_to_rev = $ac2+$rc2;
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $format\n"; exit 1; }
									}
								  elsif ($caller eq 'snvsniffer')
									{ if ($info =~ /^DP\=(\d+)\;VT\=(.+?)\;.*?$/) # this will only identify SNPs
										{ my $type = $2;
										  if ($type ne 'SNP') { $skipped++; }
										  next if ($type ne 'SNP');
										  $dp = $1;
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $info\n"; exit 1; }
									}
								  elsif ($caller eq 'solsnp')
									{ if ($info =~ /^CL\=.*?\;DP\=(.*?)\;AR\=(.*?)\;PL\=.*?$/)
										{ $dp = $1; $allele_balance = $2; }
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $info\n"; exit 1; }
									}
								  elsif ($caller eq 'spandx')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if ($info =~ /^.*?\;DP\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($format =~ /^.*?\:(.*?)\:.+?$/)
										{ my $ad = $1;
										  if ($ad =~ /^(\d+)\,(\d+)$/)
											{ $reads_mapped_to_ref = $1;
											  $reads_mapped_to_alt = $2;
											  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt);
											}
										}
									}
								  elsif ($caller eq 'speedseq')
									{ my $type = ''; if ($info =~ /^.*?TYPE\=(.*?)\;.*?$/) { $type = $1; }
									  if ($type ne 'snp') { $skipped++; }
									  next if ($type ne 'snp');
									  my $saf = ''; my $srf = ''; my $sar = ''; my $srr = '';
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/)    { $dp  = $1; }
									  if ($info =~ /^AB\=(.*?)\;.*?$/)   	 { $allele_balance = $1; $allele_balance = 1-$allele_balance; }
									  if ($info =~ /^.*?\;AO\=(.*?)\;.*?$/)  { $reads_mapped_to_alt = $1; }
									  if ($info =~ /^.*?\;RO\=(.*?)\;.*?$/)  { $reads_mapped_to_ref = $1; }
									  if ($info =~ /^.*?\;RPR\=(\d+)\;.*?$/) { $reads_mapped_right  = $1; }
									  if ($info =~ /^.*?\;RPL\=(\d+)\;.*?$/) { $reads_mapped_left   = $1; }
									  if ($info =~ /^.*?\;SAF\=(\d+)\;.*?$/) { $saf = $1; }
									  if ($info =~ /^.*?\;SRF\=(\d+)\;.*?$/) { $srf = $1; }
									  if ($info =~ /^.*?\;SAR\=(\d+)\;.*?$/) { $sar = $1; }
									  if ($info =~ /^.*?\;SRR\=(\d+)\;.*?$/) { $srr = $1; }							  
									  $reads_mapped_to_fwd = $saf+$srf;
									  $reads_mapped_to_rev = $sar+$srr;
									}
								  elsif ($caller eq 'strelka')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref != $len_alt) and ($len_ref != 1)) { $skipped++; }
									  next if (($len_ref != $len_alt) and ($len_ref != 1) and ($alt ne '.')); # i.e. not a SNP
									  if ($format =~ /^.*?\:.*?\:.*?\:(\d+)\:.*?\:.*?\:(.*?)\:(.*?)\:.*?\:.*?\:.*?$/)
										{ $dp = $1; my $adf = $2; my $adr = $3;
										  my $fwd_ref; my $rev_ref; my $fwd_alt; my $rev_alt;
										  if ($adf =~ /^(\d+)\,(\d+)$/) { $fwd_ref = $1; $fwd_alt = $2; }
										  if ($adr =~ /^(\d+)\,(\d+)$/) { $rev_ref = $1; $rev_alt = $2; }
										  $reads_mapped_to_alt = $fwd_alt+$rev_alt;
										  $reads_mapped_to_ref = $fwd_ref+$rev_ref;
										  $reads_mapped_to_fwd = $fwd_ref+$fwd_alt;
										  $reads_mapped_to_rev = $rev_ref+$rev_alt;
										  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt);
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $format\n"; exit 1; }
									}
								  elsif ($caller eq 'varscan')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if ($format =~ /^.*?\:.*?\:.*?\:(\d+)\:.*?\:.*?\:(.*?)\:.*?\:.*?\:.*?\:(\d+)\:(\d+)\:(\d+)\:(\d+)$/)
										{ $dp = $1; $allele_balance = $2; my $rdf = $3; my $rdr = $4; my $adf = $5; my $adr = $6;
										  $allele_balance =~ s/\%//; $allele_balance = $allele_balance/100;
										  $reads_mapped_to_alt = $adf+$adr;
										  $reads_mapped_to_ref = $rdf+$rdr;
										  $reads_mapped_to_fwd = $rdf+$adf;
										  $reads_mapped_to_rev = $rdr+$adr;
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $format\n"; exit 1; }
									}
								  next if ($skipped > 0);
							      my @sorted_snp_pos = @{$sorted_snps{$chr}};
								  my @dist_to_snp = ();
								  for(my $x=0;$x<@sorted_snp_pos;$x++)
									{ my $this_pos = $sorted_snp_pos[$x];
									  next if ($this_pos != $pos);
									  if (($x != 0) and (defined($sorted_snp_pos[$x-1])))
										{ my $prev_pos = $sorted_snp_pos[$x-1];
										  my $dist_to_prev = $this_pos-$prev_pos;
										  push(@dist_to_snp,$dist_to_prev);
										}
									  if (defined($sorted_snp_pos[$x+1]))
										{ my $next_pos = $sorted_snp_pos[$x+1];
										  my $dist_to_next = $next_pos-$this_pos;
										  push(@dist_to_snp,$dist_to_next);
										}
									}
								  my $distance_to_nearest_indel = 'NA';
								  if (exists($indels{$chr}))
									{ my $pushed = 0;
									  my @sorted_indel_pos_plus_snp = ();
									  my @sorted_indel_pos = @{$sorted_indels{$chr}};
									  for(my $x=0;$x<@sorted_indel_pos;$x++) # note that $pos (which is a SNP) is not present in this array, and so we need to insert it
										{ my $this_pos = $sorted_indel_pos[$x];
										  if (($x == 0) and ($pos < $this_pos) and ($pushed == 0))
											{ push(@sorted_indel_pos_plus_snp,$pos);
											  $pushed++;
											}
										  if (defined($sorted_indel_pos[$x+1]))
											{ my $next_pos = $sorted_indel_pos[$x+1];
											  if (($pos > $this_pos) and ($pos <= $next_pos) and ($pushed == 0))
												{ push(@sorted_indel_pos_plus_snp,$this_pos,$pos);
												  $pushed++;
												}
											  else
												{ push(@sorted_indel_pos_plus_snp,$this_pos); }
											}
										  else
											{ push(@sorted_indel_pos_plus_snp,$this_pos); }
										}
									  if ($pushed == 0)
										{ push(@sorted_indel_pos_plus_snp,$pos); }									  
									  my @dist_to_indel = ();
									  for(my $x=0;$x<@sorted_indel_pos_plus_snp;$x++)
										{ my $this_pos = $sorted_indel_pos_plus_snp[$x];
										  next if ($this_pos != $pos);
										  if (($x != 0) and (defined($sorted_indel_pos_plus_snp[$x-1])))
											{ my $prev_pos = $sorted_indel_pos_plus_snp[$x-1];
											  my $dist_to_prev = $this_pos-$prev_pos;
											  push(@dist_to_indel,$dist_to_prev);
											}
										  if (defined($sorted_indel_pos_plus_snp[$x+1]))
											{ my $next_pos = $sorted_indel_pos_plus_snp[$x+1];
											  my $dist_to_next = $next_pos-$this_pos;
											  push(@dist_to_indel,$dist_to_next);
											}
										}
									  my @sorted_dist_to_indel   = sort {$a <=> $b} @dist_to_indel;
									  $distance_to_nearest_indel = $sorted_dist_to_indel[0] unless (!(defined($sorted_dist_to_indel[0])));
									}
								  my @sorted_dist_to_snp  = sort {$a <=> $b} @dist_to_snp;
								  my $dist_to_nearest_snp = $sorted_dist_to_snp[0];
								  my $avg_qual_of_variant_containing_reads = 'NA';
								  if (($qual !~ /NA/) and ($reads_mapped_to_alt !~ /NA/))
									{ if ($reads_mapped_to_alt > 0)
										{ $avg_qual_of_variant_containing_reads = sprintf("%.3f",($qual/$reads_mapped_to_alt)); }
									}
								  if ($allele_balance !~ /NA/)
									{ $allele_balance = sprintf("%.3f",$allele_balance); }
							      if ($y == 0) # FPs
									{ print OUT2 "$sample_id\t$aligner\t$caller\t$aligner/$caller\t$chr\t$pos\t$ref\t$alt\t$qual\t$dp\t$reads_mapped_to_fwd\t$reads_mapped_to_rev\t$reads_mapped_to_ref\t$reads_mapped_to_alt\t$allele_balance\t$reads_mapped_left\t$reads_mapped_right\t$avg_qual_of_variant_containing_reads\t$dist_to_nearest_snp\t$distance_to_nearest_indel\n"; }
								  elsif ($y == 1) # TPs
									{ print OUT3 "$sample_id\t$aligner\t$caller\t$aligner/$caller\t$chr\t$pos\t$ref\t$alt\t$qual\t$dp\t$reads_mapped_to_fwd\t$reads_mapped_to_rev\t$reads_mapped_to_ref\t$reads_mapped_to_alt\t$allele_balance\t$reads_mapped_left\t$reads_mapped_right\t$avg_qual_of_variant_containing_reads\t$dist_to_nearest_snp\t$distance_to_nearest_indel\n"; }
								}
							  elsif ($x == 1) # indels
								{ my $allele_balance = 'NA'; my $dp = 'NA'; my $sum_ref_qual = 'NA'; my $sum_alt_qual = 'NA';
								  my $reads_mapped_to_ref = 'NA'; my $reads_mapped_to_alt = 'NA';
								  my $reads_mapped_to_fwd = 'NA'; my $reads_mapped_to_rev = 'NA';
								  my $reads_mapped_left = 'NA'; my $reads_mapped_right = 'NA';
								  if ($caller eq 'breseq')
									{ if ($info =~ /^.*?\;DP\=(.*?)$/) { $dp 			 = $1; }
									  if ($info =~ /^AF\=(.*?)\;.*?$/) { $allele_balance = $1; }
									  if ($info =~ /^.*?\;AD\=(.*?)\;.*?$/)
										{ my $ad = $1;
										  $reads_mapped_to_ref = $dp-$ad;
										  $reads_mapped_to_alt = $ad;
										}
									}
								  elsif ($caller eq 'deepvariant')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref == $len_alt) and ($len_ref == 1)) { $skipped++; }
									  next if (($len_ref == $len_alt) and ($len_ref == 1)); # i.e. not an indel
									  if ($format =~ /^.*?\:.*?\:(\d+)\:(.*?)\:(.*?)\:.*?$/)
										{ $dp = $1; my $ad = $2; $allele_balance = $3;
										  if ($ad =~ /^(\d+)\,(\d+)$/)
											{ $reads_mapped_to_ref = $1;
											  $reads_mapped_to_alt = $2;
											}
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $format\n"; exit 1; }
									}
								  elsif ($caller eq 'freebayes')
									{ my $type = ''; if ($info =~ /^.*?TYPE\=(.*?)\;.*?$/) { $type = $1; }
									  if (($type ne 'ins') and ($type ne 'del')) { $skipped++; }
									  next if (($type ne 'ins') and ($type ne 'del'));
									  my $saf = ''; my $srf = ''; my $sar = ''; my $srr = '';
									  if ($info =~ /^.*?\;DP\=(\d+)\;.*?$/)   { $dp  = $1; }
									  if ($info =~ /^AB\=(.*?)\;.*?$/)   	  { $allele_balance = $1; $allele_balance = 1-$allele_balance; }
									  if ($info =~ /^.*?\;AO\=(\d+)\;.*?$/)   { $reads_mapped_to_alt = $1; }
									  if ($info =~ /^.*?\;QA\=(\d+)\;.*?$/)   { $sum_alt_qual = $1; } # inexplicably, is 0 in all cases
									  if ($info =~ /^.*?\;QR\=(\d+)\;.*?$/)   { $sum_ref_qual = $1; } # inexplicably, is 0 in all cases
									  if ($info =~ /^.*?\;RO\=(\d+)\;.*?$/)   { $reads_mapped_to_ref = $1; }
									  if ($info =~ /^.*?\;RPR\=(\d+)\;.*?$/)  { $reads_mapped_right  = $1; }
									  if ($info =~ /^.*?\;RPL\=(\d+)\;.*?$/)  { $reads_mapped_left   = $1; }
									  if ($info =~ /^.*?\;SAF\=(\d+)\;.*?$/)  { $saf = $1; }
									  if ($info =~ /^.*?\;SRF\=(\d+)\;.*?$/)  { $srf = $1; }
									  if ($info =~ /^.*?\;SAR\=(\d+)\;.*?$/)  { $sar = $1; }
									  if ($info =~ /^.*?\;SRR\=(\d+)\;.*?$/)  { $srr = $1; }							  
									  $reads_mapped_to_fwd = $saf+$srf;
									  $reads_mapped_to_rev = $sar+$srr;
									}
								  elsif ($caller eq 'gatk')
									{ if (($len_ref == $len_alt) and ($len_ref == 1)) { $skipped++; }
									  next if (($len_ref == $len_alt) and ($len_ref == 1)); # i.e. not an indel
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($format =~ /^.*?\:(.*?)\:\d+\:.*?\:.*?$/)
										{ my $ad = $1;
										  if ($ad =~ /^(\d+)\,(\d+)$/)
											{ $reads_mapped_to_ref = $1;
											  $reads_mapped_to_alt = $2;
											  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt) unless (($reads_mapped_to_ref+$reads_mapped_to_alt)==0);
											}
										}
									}
								  elsif ($caller eq 'mpileup')
									{ if ($info !~ /INDEL/) { $skipped++; }
									  next if ($info !~ /INDEL/);
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/) { $dp = $1; }
									  my $dp4 = '';
									  if ($info =~ /^.*?DP4\=(.*?)\;.*?$/) { $dp4 = $1; }
									  if ($dp4 =~ /^(\d+)\,(\d+)\,(\d+)\,(\d+)$/)
										{ my $fwd_ref = $1; my $rev_ref = $2; my $fwd_alt = $3; my $rev_alt = $4;
										  $reads_mapped_to_alt = $fwd_alt+$rev_alt;
										  $reads_mapped_to_ref = $fwd_ref+$rev_ref;
										  $reads_mapped_to_fwd = $fwd_ref+$fwd_alt;
										  $reads_mapped_to_rev = $rev_ref+$rev_alt;
										  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt);
										}
									}
								  elsif ($caller eq 'octopus')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref == $len_alt) and ($len_ref == 1)) { $skipped++; }
									  next if (($len_ref == $len_alt) and ($len_ref == 1)); # i.e. not an indel
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/) { $dp = $1; }
									}
								  elsif ($caller eq 'pilon')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref == $len_alt) and ($len_ref == 1)) { $skipped++; }
									  next if (($len_ref == $len_alt) and ($len_ref == 1)); # i.e. not an indel
									  if ($info =~ /^DP\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($info =~ /^.*?\;AF\=(.*?)$/) { $allele_balance = $1; }
									  my $num_a = 0; my $num_c = 0; my $num_g = 0; my $num_t = 0;
									  if ($info =~ /^.*?BC\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;.*?$/) { $num_a = $1; $num_c = $2; $num_g = $3; $num_t = $4; }
									  if ($ref eq 'A') { $reads_mapped_to_ref = $num_a; }
									  if ($ref eq 'C') { $reads_mapped_to_ref = $num_c; }
									  if ($ref eq 'G') { $reads_mapped_to_ref = $num_g; }
									  if ($ref eq 'T') { $reads_mapped_to_ref = $num_t; }
									  if ($alt eq 'A') { $reads_mapped_to_alt = $num_a; }
									  if ($alt eq 'C') { $reads_mapped_to_alt = $num_c; }
									  if ($alt eq 'G') { $reads_mapped_to_alt = $num_g; }
									  if ($alt eq 'T') { $reads_mapped_to_alt = $num_t; }
									}
								  elsif ($caller eq 'platypus')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref == $len_alt) and ($len_ref == 1)) { $skipped++; }
									  next if (($len_ref == $len_alt) and ($len_ref == 1)); # i.e. not an indel
									  if ($info =~ /^.*?TC\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($info =~ /^.*?NF\=(\d+)\;.*?$/) { $reads_mapped_to_fwd = $1; }
									  if ($info =~ /^.*?NR\=(\d+)\;.*?$/) { $reads_mapped_to_rev = $1; }
									  if ($info =~ /^.*?TR\=(\d+)\;.*?$/) { $reads_mapped_to_alt = $1; }
									  $reads_mapped_to_ref = $dp-$reads_mapped_to_alt;
									  $allele_balance = $reads_mapped_to_alt/$dp;
									}
								  elsif ($caller eq 'snippy')
									{ if ($info =~ /^AB\=(\d+)\;AO\=(\d+)\;DP\=(\d+)\;QA\=(\d+)\;QR\=(\d+)\;RO\=(\d+)\;TYPE\=(.*?)$/)
										{ my $type = $7;
										  if (($type ne 'ins') and ($type ne 'del')) { $skipped++; }
										  next if (($type ne 'ins') and ($type ne 'del'));
										  $allele_balance = $1; $reads_mapped_to_alt = $2; $dp = $3; $sum_alt_qual = $4; $sum_ref_qual = $5; $reads_mapped_to_ref = $6; 
										  $allele_balance = 1-$allele_balance;								  
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $info\n"; exit 1; }
									}
								  elsif ($caller eq 'snvsniffer')
									{ if (($info =~ /^DP\=(\d+)\;VT\=INS$/) or ($info =~ /^DP\=(\d+)\;VT\=DEL$/)) # this will only identify indels
										{ $dp = $1;
										}
									}
								  elsif ($caller eq 'spandx')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if ($info =~ /^.*?\;DP\=(\d+)\;.*?$/) { $dp = $1; }
									  if ($format =~ /^.*?\:(.*?)\:.+?$/)
										{ my $ad = $1;
										  if ($ad =~ /^(\d+)\,(\d+)$/)
											{ $reads_mapped_to_ref = $1;
											  $reads_mapped_to_alt = $2;
											  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt);
											}
										}
									}
								  elsif ($caller eq 'speedseq')
									{ my $type = ''; if ($info =~ /^.*?TYPE\=(.*?)\;.*?$/) { $type = $1; }
									  if (($type ne 'ins') and ($type ne 'del')) { $skipped++; }
									  next if (($type ne 'ins') and ($type ne 'del'));
									  my $saf = ''; my $srf = ''; my $sar = ''; my $srr = '';
									  if ($info =~ /^.*?DP\=(\d+)\;.*?$/)    { $dp  = $1; }
									  if ($info =~ /^AB\=(.*?)\;.*?$/)   	 { $allele_balance = $1; $allele_balance = 1-$allele_balance; }
									  if ($info =~ /^.*?\;AO\=(.*?)\;.*?$/)  { $reads_mapped_to_alt = $1; }
									  if ($info =~ /^.*?\;RO\=(.*?)\;.*?$/)  { $reads_mapped_to_ref = $1; }
									  if ($info =~ /^.*?\;RPR\=(\d+)\;.*?$/) { $reads_mapped_right  = $1; }
									  if ($info =~ /^.*?\;RPL\=(\d+)\;.*?$/) { $reads_mapped_left   = $1; }
									  if ($info =~ /^.*?\;SAF\=(\d+)\;.*?$/) { $saf = $1; }
									  if ($info =~ /^.*?\;SRF\=(\d+)\;.*?$/) { $srf = $1; }
									  if ($info =~ /^.*?\;SAR\=(\d+)\;.*?$/) { $sar = $1; }
									  if ($info =~ /^.*?\;SRR\=(\d+)\;.*?$/) { $srr = $1; }							  
									  $reads_mapped_to_fwd = $saf+$srf;
									  $reads_mapped_to_rev = $sar+$srr;
									}
								  elsif ($caller eq 'strelka')
									{ if ($filter ne 'PASS') { $skipped++; }
									  next if ($filter ne 'PASS');
									  if (($len_ref == $len_alt) and ($len_ref == 1)) { $skipped++; }
									  next if (($len_ref == $len_alt) and ($len_ref == 1)); # i.e. not an indel
									  if ($format =~ /^.*?\:.*?\:.*?\:.*?\:.*?\:(.*?)\:(.*?)\:.*?\:.*?$/)
										{ my $adf = $1; my $adr = $2;
										  my $fwd_ref; my $rev_ref; my $fwd_alt; my $rev_alt;
										  if ($adf =~ /^(\d+)\,(\d+)$/) { $fwd_ref = $1; $fwd_alt = $2; }
										  if ($adr =~ /^(\d+)\,(\d+)$/) { $rev_ref = $1; $rev_alt = $2; }
										  if ( (!(defined($fwd_ref))) or (!(defined($fwd_alt))) or (!(defined($rev_ref))) or (!(defined($rev_alt))) )
											{ $skipped++; }
										  next if ( (!(defined($fwd_ref))) or (!(defined($fwd_alt))) or (!(defined($rev_ref))) or (!(defined($rev_alt))) );
										  $reads_mapped_to_alt = $fwd_alt+$rev_alt;
										  $reads_mapped_to_ref = $fwd_ref+$rev_ref;
										  $reads_mapped_to_fwd = $fwd_ref+$fwd_alt;
										  $reads_mapped_to_rev = $rev_ref+$rev_alt;
										  $allele_balance = $reads_mapped_to_alt/($reads_mapped_to_ref+$reads_mapped_to_alt);
										}
									  else
										{ print "error with $sample_id, $chr, $pos, $caller - unable to parse $format\n"; exit 1; }
									}
								  next if ($skipped > 0);
								  my $dist_to_nearest_snp = 'NA';
								  if (exists($sorted_snps{$chr}))
									{ my $pushed = 0;
									  my @sorted_snp_pos_plus_indels = ();
									  my @sorted_snp_pos = @{$sorted_snps{$chr}};
									  for(my $x=0;$x<@sorted_snp_pos;$x++) # note that $pos (which is a SNP) is not present in this array, and so we need to insert it
										{ my $this_pos = $sorted_snp_pos[$x];
										  if (($x == 0) and ($pos < $this_pos) and ($pushed == 0))
											{ push(@sorted_snp_pos_plus_indels,$pos);
											  $pushed++;
											}
										  if (defined($sorted_snp_pos[$x+1]))
											{ my $next_pos = $sorted_snp_pos[$x+1];
											  if (($pos > $this_pos) and ($pos <= $next_pos) and ($pushed == 0))
												{ push(@sorted_snp_pos_plus_indels,$this_pos,$pos);
												  $pushed++;
												}
											  else
												{ push(@sorted_snp_pos_plus_indels,$this_pos); }
											}
										  else
											{ push(@sorted_snp_pos_plus_indels,$this_pos); }
										}
									  if ($pushed == 0)
										{ push(@sorted_snp_pos_plus_indels,$pos); }									  
									  my @dist_to_snp = ();
									  for(my $x=0;$x<@sorted_snp_pos_plus_indels;$x++)
										{ my $this_pos = $sorted_snp_pos_plus_indels[$x];
										  next if ($this_pos != $pos);
										  if (($x != 0) and (defined($sorted_snp_pos_plus_indels[$x-1])))
											{ my $prev_pos = $sorted_snp_pos_plus_indels[$x-1];
											  my $dist_to_prev = $this_pos-$prev_pos;
											  push(@dist_to_snp,$dist_to_prev);
											}
										  if (defined($sorted_snp_pos_plus_indels[$x+1]))
											{ my $next_pos = $sorted_snp_pos_plus_indels[$x+1];
											  my $dist_to_next = $next_pos-$this_pos;
											  push(@dist_to_snp,$dist_to_next);
											}
										}
									  my @sorted_dist_to_snp  = sort {$a <=> $b} @dist_to_snp;
									  $dist_to_nearest_snp = $sorted_dist_to_snp[0] unless (!(defined($sorted_dist_to_snp[0])));
									}
								  my $distance_to_nearest_indel = 'NA';
								  if (exists($indels{$chr}))
									{ my @sorted_indel_pos = @{$sorted_indels{$chr}};
									  my @dist_to_indel = ();
									  for(my $x=0;$x<@sorted_indel_pos;$x++)
										{ my $this_pos = $sorted_indel_pos[$x];
										  next if ($this_pos != $pos);
										  if (($x != 0) and (defined($sorted_indel_pos[$x-1])))
											{ my $prev_pos = $sorted_indel_pos[$x-1];
											  my $dist_to_prev = $this_pos-$prev_pos;
											  push(@dist_to_indel,$dist_to_prev);
											}
										  if (defined($sorted_indel_pos[$x+1]))
											{ my $next_pos = $sorted_indel_pos[$x+1];
											  my $dist_to_next = $next_pos-$this_pos;
											  push(@dist_to_indel,$dist_to_next);
											}
										}
									  my @sorted_dist_to_indel   = sort {$a <=> $b} @dist_to_indel;
									  $distance_to_nearest_indel = $sorted_dist_to_indel[0] unless (!(defined($sorted_dist_to_indel[0])));
									}
								  my $avg_qual_of_variant_containing_reads = 'NA';
								  if (($qual !~ /NA/) and ($reads_mapped_to_alt !~ /NA/) and ($reads_mapped_to_alt !~ /\,/))
									{ if ($reads_mapped_to_alt > 0)
										{ $avg_qual_of_variant_containing_reads = sprintf("%.3f",($qual/$reads_mapped_to_alt)); }
									}
								  if ($allele_balance !~ /NA/)
									{ $allele_balance = sprintf("%.3f",$allele_balance); }
								  if ($y == 0) # FPs
									{ print OUT4 "$sample_id\t$aligner\t$caller\t$aligner/$caller\t$chr\t$pos\t$ref\t$alt\t$qual\t$dp\t$reads_mapped_to_fwd\t$reads_mapped_to_rev\t$reads_mapped_to_ref\t$reads_mapped_to_alt\t$allele_balance\t$reads_mapped_left\t$reads_mapped_right\t$avg_qual_of_variant_containing_reads\t$dist_to_nearest_snp\t$distance_to_nearest_indel\n"; }
								  elsif ($y == 1) # TPs
									{ print OUT5 "$sample_id\t$aligner\t$caller\t$aligner/$caller\t$chr\t$pos\t$ref\t$alt\t$qual\t$dp\t$reads_mapped_to_fwd\t$reads_mapped_to_rev\t$reads_mapped_to_ref\t$reads_mapped_to_alt\t$allele_balance\t$reads_mapped_left\t$reads_mapped_right\t$avg_qual_of_variant_containing_reads\t$dist_to_nearest_snp\t$distance_to_nearest_indel\n"; }
								}
							}
						}
					}
				  close(IN) or die $!;
				}
			}
		}
	}
close(OUT1) or die $!; close(OUT2) or die $!; close(OUT3) or die $!; close(OUT4) or die $!; close(OUT5) or die $!;
exit 1;