#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::MoreUtils qw(firstidx);

=head1 program description
    
    filter VCF in genotype fields of records,eg: GQ,DP... 
    
=head1 command line options
    
    --input [required]  the raw vcf file
	--dp	[required] the depth cutoff
	--dp_up	[required] the up level of depth 2*sequence coverage
	--gq	[requried] the GQ cutoff
    --out   [required]   output file
    --help          illustration of this script
    
=head1 Usage
    
    perl filter_genotype.pl --help
    
=cut

my ($input,$dp,$dp_up,$gq,$out,$help);

GetOptions(
    "input=s"=>\$input,
	"dp=s"=>\$dp,
	"dp_up=s"=>\$dp_up,
	"gq=s"=>\$gq,
    "out=s"=>\$out,
    "help"=>\$help,
);
die `pod2text $0` if($help);
die `pod2text $0` if(!defined $input);
die `pod2text $0` if(!defined $dp);
die `pod2text $0` if(!defined $dp_up);
die `pod2text $0` if(!defined $gq);
die `pod2text $0` if(!defined $out);

open IN,"<$input";
open OUT,">$out";
while(<IN>)
{
    if (/^#/)
    {
        print OUT $_;
    }
    else
    {
		chomp;
		my @all=split(/\t/);
		my $flag = 0;
		if($all[4]=~/,/)
		{
			next;
		}
		else{
			if($all[8] =~ /:GQ:/)
			{
				my @b=split(/:/,$all[8]);
				my $GQidx = firstidx { $_ eq 'GQ' } @b;
				my $DPidx = firstidx { $_ eq 'NR' } @b;
				for(my $i=9;$i<=$#all;$i++)
				{
					if ($all[$i] =~ /^\./)
					{
						$flag += 1;
					}
					else
					{
						my @tmp = split(/:/,$all[$i]);
						if($tmp[$GQidx] < $gq || $tmp[$DPidx] < $dp || $tmp[$DPidx] > $dp_up)
						{
							$flag += 1;
						}
					}
				}
			}
			else
			{
				print "no GQ\n";
			}
			if($flag == 0)
			{
				print OUT $_,"\n";
			}
			else
			{
				next;
			}
		}
    } 
}
close IN;
close OUT;