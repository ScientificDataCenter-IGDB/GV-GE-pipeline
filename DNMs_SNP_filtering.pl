###################################################
# Author     : Fang Lu                            #
# E-mail     : l.fang@genetics.ac.cn              #
# Date       : Thu Mar 08 2024                    #
# Version    : V2.0                               #
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

    filter VCF in off fa mo bam files

=head1 command line options

    --input [required]  the raw vcf file
    --off   [required] offspring sample name
    --fa    [required] father sample name
    --mo    [required] mother sample name
	--spe	[required] reference sequences, format: fasta
    --out   [required]   output file
    --help  illustration of this script

=head1 Usage

    perl home_check_INDEL_v2.pl --help

=cut

my ($input,$off,$fa,$mo,$spe,$out,$help);

GetOptions(
    "input=s"=>\$input,
    "off=s"=>\$off,
    "fa=s"=>\$fa,
    "mo=s"=>\$mo,
	"spe=s"=>\$spe,
    "out=s"=>\$out,
    "help"=>\$help,
);
die `pod2text $0` if($help);
die `pod2text $0` if(!defined $input);
die `pod2text $0` if(!defined $off);
die `pod2text $0` if(!defined $fa);
die `pod2text $0` if(!defined $mo);
die `pod2text $0` if(!defined $spe);
die `pod2text $0` if(!defined $out);

open IN,"<$input";
open OUT,">$out";

my ($offidx,$faidx,$moidx);
while(<IN>)
{
    chomp;
    if (/^#/)
    {
        print OUT $_,"\n";
        if(/^#CHROM/)
        {
            my @tmp=split(/\t/);
            $offidx = firstidx { $_ eq $off} @tmp;
            $faidx = firstidx { $_ eq  $fa} @tmp;
            $moidx = firstidx { $_ eq  $mo} @tmp;
        }
    }
    else{
        my $info = $_;
        my @all = split(/\t/);

        my ($off_geno, $fa_geno, $mo_geno);
        $off_geno = get_genotype($all[$offidx]);
        $fa_geno = get_genotype($all[$faidx]);
        $mo_geno = get_genotype($all[$moidx]);

        my $pos = $all[0].":".$all[1];
        my $ref = $all[3];
        my $alt = $all[4];
       
        system("samtools tview -p $pos /public-de6000/fanglu/Project/Zhangyongqing_Lab/02_GATK_GVCF/${off}/${off}_sorted_merge_dedup.bam $spe -d T > off_snp.log");
        system("samtools tview -p $pos /public-de6000/fanglu/Project/Zhangyongqing_Lab/02_GATK_GVCF/${fa}/${fa}_sorted_merge_dedup.bam $spe -d T > fa_snp.log");
        system("samtools tview -p $pos /public-de6000/fanglu/Project/Zhangyongqing_Lab/02_GATK_GVCF/${mo}/${mo}_sorted_merge_dedup.bam $spe -d T > mo_snp.log");

        my $off_var_num = get_var_reads_number("off_snp.log", $alt);
        my $fa_var_num = get_var_reads_number("fa_snp.log", $alt);
        my $mo_var_num = get_var_reads_number("mo_snp.log", $alt);

        if(($off_geno eq "0/1") or ($off_geno eq "1/1"))
        {
            if (($fa_var_num > 1) or ($mo_var_num > 1))
            {
                next;
            }
            else{
                print OUT $info,"\n";
            }
        }
        else
        {
            if (($fa_var_num > 1) and ($mo_var_num > 1) and ($off_var_num < 2))
            {
                print OUT $info,"\n";
            }
            else{
                next;
            }
        }
    }
}
close IN;
close OUT;

sub get_genotype{
    my $value = shift;
    my $geno;
    my @sam = split(":",$value);
    my @sam_GT = split("[\/|\|]",$sam[0]);
    if( $sam_GT[0] != $sam_GT[1] )
    {
        $geno = "0/1";
    }
    else{
        $geno = $sam_GT[0]."/".$sam_GT[1];
    }
    return $geno;
}

sub get_var_reads_number{
    my ($file, $alt) = @_;
    open R,$file;
    my $first_line = <R>;
    my $second_line = <R>;
    my $third_line = <R>;

    my $var_num=0;
    while(<R>)
    {
        chomp;
        my @b = split(//);
        my $var = uc($b[0]);
        if($var eq $alt)
        {
            $var_num+=1;
        }
    }
    close R;
    return $var_num;
}