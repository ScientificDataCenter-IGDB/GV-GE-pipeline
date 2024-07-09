###################################################
# Author     : Fang Lu                            #
# E-mail     : l.fang@genetics.ac.cn              #
# Date       : Thu Mar 08 2024                    #
# Version    : V1.0                               #
# You are using the program scripted by Fang Lu   #
# Please let me know if you have any suggestions  #
# Thank you !                                     #
# Best wishes !                                   #
###################################################
#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Math::BigFloat;
use List::MoreUtils qw(firstidx);

=head1 program description

    Defined offspring denovo mutations with vaf cutoff;

=head1 command line options

    --input [required]  the raw vcf file
    --samplename [required] sample name
    --vaf   [required] vaf default 0.2
    --out   [required]   output file
    --help          illustration of this script

=head1 Usage

    perl defined_DNM_platypus.pl --help

=cut

my ($input,$samplename,$vaf,$out,$help);

GetOptions(
    "input=s"=>\$input,
    "samplename=s"=>\$samplename,
    "vaf=s"=>\$vaf,
    "out=s"=>\$out,
    "help"=>\$help,
);

die `pod2text $0` if($help);
die `pod2text $0` if(!defined $input);
die `pod2text $0` if(!defined $samplename);
die `pod2text $0` if(!defined $vaf);
die `pod2text $0` if(!defined $out);

open IN,"<$input";
open OUT,">$out";

my $offidx;
while(<IN>)
{
    chomp;
    if (/^#/)
    {
        if(/^#CHROM/)
        {
            my @tmp=split(/\t/);
            $offidx = firstidx { $_ eq  $samplename} @tmp;
            print OUT $_,"\n";
        }
        else{
            print OUT $_,"\n";
        }
    }
    else
    {
        my @all=split(/\t/);
        my @c=split(":",$all[8]);
        my $GTidx = firstidx { $_ eq "GT" } @c;
        my $NRidx = firstidx { $_ eq 'NR' } @c;
        my $NVidx = firstidx { $_ eq 'NV' } @c;

        ### revise genotype
        my (@off_GT,$off_ref,$off_var,$off_vaf,$off_geno,$off_ref_vaf);
        my (@fa_GT,$fa_geno,@fa,$fa_ref,$fa_var);
        my (@mo_GT,$mo_geno,@mo,$mo_ref,$mo_var);
        my $f = 0;
        for (my $i=9;$i<=11;$i++)
        {
            if ($i == $offidx)
            {
                my @off = split(":",$all[$offidx]);
                @off_GT = split("[\/|\|]",$off[$GTidx]);
                $off_ref = $off[$NRidx] - $off[$NVidx];
                $off_var = $off[$NVidx];
            }
            else{
                if($f == 0)
                {
                    @fa = split(":",$all[$i]);
                    @fa_GT = split("[\/|\|]",$fa[$GTidx]);
                    $fa_ref = $fa[$NRidx] - $fa[$NVidx];
                    $fa_var = $fa[$NVidx];
                    $f = 1;
                }
                else{
                    @mo = split(":",$all[$i]);
                    @mo_GT = split("[\/|\|]",$mo[$GTidx]);
                    $mo_ref = $mo[$NRidx] - $mo[$NVidx];
                    $mo_var = $mo[$NVidx];
                }
            }
        }
        my $fa_alt_reads = $fa_var;
        my $mo_alt_reads = $mo_var;
        my $fa_ref_reads = $fa_ref;
        my $mo_ref_reads = $mo_ref;
        my $de = join("_",@off_GT,@fa_GT,@mo_GT);
        if($de =~ /\./)
        {
            next;
        }
        
        $off_geno = get_genotype($off_GT[0],$off_GT[1]);
        $fa_geno = get_genotype($fa_GT[0],$fa_GT[1]);
        $mo_geno = get_genotype($mo_GT[0],$mo_GT[1]);

        $off_vaf = $off_var/($off_ref + $off_var);
        $off_ref_vaf = $off_ref/($off_ref + $off_var);
        #### denovo mutation
        if(($off_ref + $off_var == 0) or ($fa_ref + $fa_var == 0) or ($mo_ref + $mo_var == 0))
        {
            next;
        }
        else{
            if($fa_geno eq "0/0" and $mo_geno eq "0/0")
            {
                if(($off_geno eq "0/1") or ($off_geno eq "1/1"))
                {
                    if($fa_alt_reads >= 2 or $mo_alt_reads >= 2)
                    {
                        next;
                    }
                    else{
                        if($off_vaf > $vaf)
                        {
                            print OUT $_,"\n";
                        }
                    }   
                }
                else
                {
                    next;
                }
            }
            elsif(($fa_geno eq "0/1" and $mo_geno eq "0/1"))
            {
                next;
            }
            elsif(($fa_geno eq "0/0" and $mo_geno eq "1/1") or ($mo_geno eq "0/0" and $fa_geno eq "1/1"))
            {
                if($off_geno eq "0/0")
                {
                    if(($fa_ref_reads >= 2 or $mo_ref_reads >= 2))
                    {
                        next;
                    }
                    else{
                        if($off_ref_vaf > $vaf)
                        {
                            print OUT $_,"\n";
                        }
                        else{
                            next;
                        }
                    }
                }
                elsif($off_geno eq "1/1")
                {
                    if(($fa_alt_reads >= 2 or $mo_alt_reads >= 2))
                    {
                        next;
                    }
                    else{
                        if($off_vaf > $vaf)
                        {
                            print OUT $_,"\n";
                        }
                        else{
                            next;
                        }
                    }
                }
                else{
                    next;
                }
            }
            elsif(($fa_geno eq "0/1" and $mo_geno eq "1/1") or ($mo_geno eq "0/1" and $fa_geno eq "1/1"))
            {
                if($off_geno eq "0/0")
                {
                    if(($fa_ref_reads >= 2 or $mo_ref_reads >= 2))
                    {
                        next;
                    }
                    else{
                        if($off_ref_vaf > $vaf)
                        {
                            print OUT $_,"\n";
                        }
                        else{
                            next;
                        }
                    }
                }
                else{
                    next;
                }
            }
            elsif(($fa_geno eq "0/0" and $mo_geno eq "0/1") or ($mo_geno eq "0/0" and $fa_geno eq "0/1"))
            {
                if($off_geno eq "1/1")
                {
                    if(($fa_alt_reads >= 2 or $mo_alt_reads >= 2))
                    {
                        next;
                    }
                    else{
                        if($off_vaf > $vaf)
                        {
                            print OUT $_,"\n";
                        }
                        else{
                            next;
                        }
                    }
                }
                else{
                    next;
                }
            }
            elsif(($fa_geno eq "1/1" and $mo_geno eq "1/1"))
            {
                if(($off_geno eq "0/1") or ($off_geno eq "0/0"))
                {
                    if(($fa_ref_reads >= 2) or ($mo_ref_reads >= 2))
                    {
                        next;
                    }
                    else{
                        if ($off_ref_vaf > $vaf)
                        {
                            print OUT $_,"\n";
                        }
                        else
                        {
                            next;
                        }
                    } 
                }
                else
                {
                    next;
                }
            }
        }
    }
}
close IN;
close OUT;


sub get_genotype{
    my ($a,$b) = @_;
    my $geno;
    if( $a != $b )
    {
        $geno = "0/1";
    }
    else{
        $geno = $a."/".$b;
    }
    return $geno;
}