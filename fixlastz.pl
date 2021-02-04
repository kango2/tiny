#!/usr/bin/perl
use strict;
use warnings;

my ($input, $output, $sizes) = @ARGV;

my %sizes = ();
open (S, "<$sizes") or die $!;
while (<S>){
	chomp $_;
	my @a = split ("\t", $_);
	$sizes{$a[0]} = $a[1];
}
close S;

open (O, ">$output") or die $!;
open (I, "gunzip -c $input |") or die $!;
while (<I>){
	chomp $_;
	if ($_ =~ /^\d+/) { 
		my @a = split (" ", $_);
		if ($a[7] eq "-") {
			my $base = (int($a[5]/1000000) * 1000000);
			if ($base < int($sizes{$a[4]}/1000000) * 1000000) {
				$a[5] = $sizes{$a[4]} - $base + $a[5] - $base - 1000000;
				$a[6] = $sizes{$a[4]} - $base + $a[6] - $base - 1000000;
			}
			elsif ($base == int($sizes{$a[4]}/1000000) * 1000000) {
        $a[5] = $a[5] - $base;
        $a[6] = $a[6] - $base;
			}
		}
		print O join(" ", @a)."\n";
	}
	else {
		print O "$_\n";
	}	
}
close I;
close O;
