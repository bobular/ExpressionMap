#!/usr/bin/perl -w
#                           -*- mode: cperl -*-
use CGI qw(:standard);

use Config::General;
use lib './';
use ExpressionMap::Map;

my $conf = new Config::General("config.txt");
my %config = $conf->getall;

my $map_name = $config{map_name} || die;
my $map = ExpressionMap::Map->retrieve(name=>$map_name, bootstrap=>0);

my ($gene_id1, $gene_id2, $radius) = param('args');

print header();
print $map->bootstrap_pair($gene_id1, $gene_id2, $radius);
