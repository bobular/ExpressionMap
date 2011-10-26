#!/usr/bin/perl -w
#                           -*- mode: cperl -*-

#
# usage:
#
# edit config.txt
#
# ./run-go-ora.pl "VectorBase AgamP3.4 1.0.3" ~/ermineJ.data/AgamP3.4-1.0.3.ermine.txt
#


use lib './';
use ExpressionMap::Map;
use Config::General;

# same config file is also loaded by ExpressionMap::Map
my $conf = new Config::General("config.txt");
my %config = $conf->getall;

my $map_name = shift || die "must give map name as first argument\n";
my $ermineJ_annotation = shift || die "must give ermineJ-style annotation as second argument\n";
die unless (-s $ermineJ_annotation);

my $nonrandom_test_max = $config{nonrandom_test_max} || 500;
my $nonrandom_test_reps = $config{nonrandom_test_reps} || 100;
my $mapheight = $config{map_height} || 300;
my $maxmapheight = 1200;
my $max_e_selected = $config{max_e_selected} || 3;
my $gene_char = $config{gene_char} || 'O';
my $data_type = $config{data_type} || 'binary';

my $ermineJ = $config{ermineJ} || 'ermineJ.sh';
my $ermineJ_goterms = $config{ermineJ_goterms} || '/home/maccallr/ermineJ.data/go_2009-03-02-termdb.rdf-xml';
my $ermineJ_max_as_fraction = $config{ermineJ_max_as_fraction} || 0.25;
my $ermineJ_min = $config{ermineJ_min} || 10;

#
# get the map object
#

my $map = ExpressionMap::Map->retrieve(name=>$map_name, bootstrap=>0);
die "can't get map: $map_name\n" unless (defined $map);

my @geneids;
my %mapnodegeneid;

foreach my $mapnode ($map->mapnodes) {
  my @mappings = $map->mappings(mapnode=>$mapnode);
  foreach $mapping (@mappings) {
    my $id = $mapping->gene_id;
    $mapnodegeneid{$mapnode}{$id} = 1;
    push @geneids, $id;
  }
}

my $ermineJ_max = int($ermineJ_max_as_fraction*@geneids);

foreach my $mapnode ($map->mapnodes) {
  my @results = $mapnode->go_ora_results();
  if (@results) {
    warn "skipping mapnode ".$mapnode->x.",".$mapnode->y."\n";
    next;
  }
  warn "running mapnode ".$mapnode->x.",".$mapnode->y."\n";
  my $tmpscores = "/tmp/tmpejscores$$";
  my $tmpresults = "/tmp/tmpejresults$$";
  open(TMP, ">$tmpscores") || die;
  foreach my $id (@geneids) {
    printf TMP "$id\t%d\n", $mapnodegeneid{$mapnode}{$id} ? 1 : 0;
  }
  close(TMP);
  open(ERMINEJ, "$ermineJ -b -a $ermineJ_annotation -c $ermineJ_goterms -s $tmpscores -t 0.5 -x $ermineJ_max -y $ermineJ_min -e 2 -n 0 -o $tmpresults && cat $tmpresults|") || die;
  while (<ERMINEJ>) {
    if (/^\!/) {
      chomp;
      my ($bang, $name, $acc, $numprobes, $numgenes, $score, $pvalue, $corrected_pvalue, $same_as) = split /\t/, $_;
      my @same_as = $same_as =~ /(GO:\d+)/g;
      my $serial_number = 1;
      foreach my $go_acc ($acc, @same_as) {
	my $go_ora_result = ExpressionMap::GoOraResult->insert({
							      mapnode=>$mapnode,
							      term=>$go_acc,
							      score=>$score,
							      numgenes=>$numgenes,
							      pvalue=>$pvalue,
							      corrected_pvalue=>$corrected_pvalue,
							      serial_number=>$serial_number++,
							      });
      }
    }
  }
  close(ERMINEJ);
  unlink $tmpscores, $tmpresults;
}
