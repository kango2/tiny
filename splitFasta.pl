#!/usr/bin/perl -l
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

##process arguments
my $config = configure(scalar @ARGV);
##read fasta file
my ($seqInfo, $totalLength) = readFasta($config->{'inputfasta'});

##intialize variables
my $lenperchunk = int($totalLength / $config->{'chunks'});
my $currentchunksize = 0;
my $chunkID = 0;
my $suffixLen = length($config->{'chunks'});

##process sequences into desired number of chunks, windowsize, overlap
open (CHUNK, ">$config->{'outputbase'}".sprintf("%0${suffixLen}d", $chunkID).".fa") or die $!;
foreach my $seq (@$seqInfo){
  for (my $i=0;$i<length($$seq{'seq'});$i+=$config->{'windowsize'} - $config->{'overlap'}){
    if ($currentchunksize > $lenperchunk) {
      $chunkID++;
      $currentchunksize = 0;
      close CHUNK;
      open (CHUNK, ">$config->{'outputbase'}".sprintf("%0${suffixLen}d", $chunkID).".fa") or die $!;
    }
    my $subseq = substr($$seq{'seq'}, $i, $config->{'windowsize'});
    if ($subseq =~ /[ACGTacgt]/) {
      print CHUNK ">$$seq{'header'}.s$i.e".(($i + $config->{'windowsize'} > length($$seq{'seq'})) ? length($$seq{'seq'}) : ($i + $config->{'windowsize'})).(($$seq{'description'} ne "NoDescription") ? " $$seq{'description'}" : "");
      for (my $i = 0; $i < length($subseq); $i += 80) {
        print CHUNK substr($subseq, $i, 80);
      }
    }
    $currentchunksize += (length($subseq) > $config->{'windowsize'} - $config->{'overlap'}) ? $config->{'windowsize'} - $config->{'overlap'} : length($subseq);
  }
}
close CHUNK;

sub readFasta {
  my $file = shift;
  my @seq = ();
  my ($header, $seq, $description);
  my $length = 0;
  if ($file =~ /\.gz$/){
    open (F, "gunzip -c $file |") or die $!;
  }
  else{
    open (F, "<$file") or die $!;
  }
  while (<F>){
    chomp $_;
    if ($_ =~ />(\S+)\s*(.*)/){
      if (defined $seq && length($seq)>0) {
        push (@seq, {'header'=>$header, 'description'=>$description, 'seq'=>$seq});
        $length += length($seq);
      }
      $header = $1;
      $description = $2;
      $description = "NoDescription" unless (length($description) > 0);
      $seq = "";
    }
    else {
      $seq .= $_;
    }
  }
  close F;
  if (defined $seq && length($seq)>0) {
    push (@seq, {'header'=>$header, 'description'=>$description, 'seq'=>$seq});
    $length += length($seq);
  }
  return (\@seq, $length);
}



sub configure {
  my $args = shift;
  my $config = {};
  $config->{'overlap'}    = 0;
  $config->{'chunks'}     = 1;
  $config->{'windowsize'} = 1000000;
  GetOptions(
  $config,
  'inputfasta|i=s',
  'outputbase|o=s',
  'windowsize|w=i',
  'overlap|s:i',
  'chunks|c:i',
  ) or usage(1);

  if ($config->{'outputbase'}) {
    runcommand("mkdir -p ". dirname($config->{'outputbase'}));
    $config->{'outputbase'} = abs_path(dirname($config->{'outputbase'}))."/".basename($config->{'outputbase'});
  }
  else {
    print "\nERROR: Provide output base for writing split files.";
    usage(1);
  }

  if ($config->{'inputfasta'} && -e $config->{'inputfasta'}) {
    $config->{'inputfasta'} = abs_path($config->{'inputfasta'});
  }
  else{
    print "\nERROR: Provide path to input fasta file.";
    usage(1);
  }

  unless ($config->{'windowsize'} >= 0) {
    print "\nERROR: Provide non-negative value for the windowsize.";
    usage(1);
  }

  unless ($config->{'overlap'} < $config->{'windowsize'} && $config->{'overlap'} >= 0) {
    print "\nERROR: Provide non-negative value less than windowsize of $config->{'windowsize'}.";
    usage(1);
  }
  unless ($config->{'chunks'} >= 1) {
    print "\nERROR: Provide value >=1 for the chunks parameters.";
    usage(1);
  }
  return $config;
}

sub runcommand {
  my ($command, $errormsg) = @_;
  unless (system($command) == 0){
    print "\nERROR: Follwing command failed\n-----------------";
    print $command;
    print "-----------------";
    print $errormsg if ($errormsg);
    print "-----------------";
    exit(1);
  }
}


sub usage {
  my $exit_code = shift;
  print <<USAGEMSG;

USAGE:

  splitFasta.pl -inputfasta input.fa.gz -outputbase outputdir/chunkname -windowsize {1..n} -overlap {0..n} -chunks {1..n}

  This utility can divide up each sequence into a fixed size sequence with a predefined overlap if desired.

Options:
  -inputfasta/-i      Input fasta file for aligned regions. It can be gzipped.
  -outputbase/-o      Path to output directory and basename for output files.
  -windowsize/-w      Maximum size of the sequence in the output. Default is 1000000.
  -overlap/-s         Number of overlapping basepairs between two consecutive windows (sliding window). Default is 0bp.
  -chunks/-c          Number of output files to divide the output into. If desired, output can be split across multiple files. Default is 1 output file.
USAGEMSG
  exit($exit_code);
}

