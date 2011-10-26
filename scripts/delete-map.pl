#!/usr/bin/perl -w

use Getopt::Long;
use PDL;
use lib './';
use PDL::Kohonen;
use ExpressionMap::Map;


my $mapname = "test map";
my $bootstrap;

GetOptions("mapname=s"=>\$mapname,
	   "bootstrap=i"=>\$bootstrap,
	  );

die unless ($mapname);

warn "this takes a while...\n";

if (defined $bootstrap) {
  ExpressionMap::Map->search(name => $mapname,
			     bootstrap=>$bootstrap)->delete_all;
} else {
  ExpressionMap::Map->search(name => $mapname)->delete_all;
}
warn "attempted to delete $mapname\n";
