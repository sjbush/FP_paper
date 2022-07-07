use strict;
use warnings;

# REQUIREMENTS
my $in_file = $ARGV[0]; my $out_file = $ARGV[1];
if (!(-e($in_file))) { die "ERROR: cannot find $in_file\n"; }

# OUTPUT
open(OUT,'>',$out_file) or die $!;

# IDENTIFY WHETHER THIS PIPELINE POPULATES THE FILTER FIELD; IF IT DOESN'T, EACH ROW WILL BE MARKED ONLY BY "."
my $filter_field_in_use = 0;
open(IN,$in_file) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line =~ /^\#/);
	  my @line = split(/\t/,$line);
	  my $filter = $line[6];
	  if ($filter ne '.')
		{ $filter_field_in_use++; }
	}
close(IN) or die $!;

# FILTER VCF TO REMOVE REFERENCE CALLS, NON-PASSING CALLS (IF RELEVANT), AND NULL ENTRIES (THOSE THAT REVISE A MASKED, I.E. "N", POSITION)
open(IN,$in_file) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\#/) { print OUT "$line\n"; }
	  next if ($line =~ /^\#/);
	  my @line = split(/\t/,$line);
	  my $ref = $line[3]; my $alt = $line[4]; my $filter = $line[6]; my $vals = $line[$#line];
	  next if (($filter_field_in_use > 0) and ($filter ne 'PASS'));
	  next if (($ref eq '.') or ($alt eq '.'));
	  next if (($ref eq 'N') or ($alt eq 'N'));
	  next if (($vals =~ /^0\/0\:.*?$/) || ($vals =~ /^0\|0\:.*?$/));
	  print OUT "$line\n";
	}
close(IN) or die $!;

close(OUT) or die $!;
exit 1;