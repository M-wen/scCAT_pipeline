#!/usr/bin/perl -w
use strict;
my $atac=shift;
my $hex=shift;

my %hash;
open IN,"$atac";
while(<IN>){
    chomp;
    my @tmp=split /\t/;
    $hash{$tmp[0]}=$tmp[-1];


}close IN;


open IN,"$hex";
while(<IN>){
    chomp;
    my @tmp=split /\t/;
    if($hash{$tmp[-1]}){
        print $tmp[0]."\t".$hash{$tmp[-1]}."\n";
    }

}close IN;

