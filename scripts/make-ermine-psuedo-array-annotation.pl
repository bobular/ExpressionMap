#!/usr/bin/perl -w

use lib './';
use ExpressionMap::Map;
use Config::General;


# same config file is also loaded by ExpressionMap::Map
my $conf = new Config::General("config.txt");
my %config = $conf->getall;

my $map_name = shift || die "must give map name as first argument\n";

my $biomart = shift || 'biomart.txt';

my %goterms; # geneid->[goterms]

#
# get gene ids from map
#

my $map = ExpressionMap::Map->retrieve(name=>$map_name, bootstrap=>0);
die "can't get map: $map_name\n" unless (defined $map);

foreach my $mapnode ($map->mapnodes) {
  my @mappings = $map->mappings(mapnode=>$mapnode);
  foreach $mapping (@mappings) {
    my $id = $mapping->gene_id;
    $goterms{$id} = [];
  }
}

open(BIOMART, $biomart) || die "can't open $biomart\n";
while (<BIOMART>) {
  my ($geneid) = split;
  my @goterms = /(GO:\d+)/g;
  if (defined $geneid && defined $goterms{$geneid} && @goterms) {
    push @{$goterms{$geneid}}, @goterms;
  }
}
close(BIOMART);

print "Probe\tGene\tDesc\tGOterms\n";
foreach my $geneid (sort keys %goterms) {
  print "$geneid\t$geneid\t\t@{$goterms{$geneid}}\n";
}
