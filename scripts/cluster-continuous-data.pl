#!/usr/bin/perl -w

use Getopt::Long;
use PDL;
use lib './';
use PDL::Kohonen;
use ExpressionMap::Map;



my $commandline = "$0 @ARGV";
my $mapdims = "25x20";
my $epochs = 20;
my $debug;
my $mapname = "test map continuous";
my $seed;
my $numbootstraps = 0;
my $alpha = 0.1;

my $printnumberofgenesandexit;

# custom experiment/condition selection
my $includeconditions; # regexp
my $excludeconditions; # regexp

# continuous options
my $median_shift;     # shift all gene vectors to median
my $min_range = 0;    # filter out all with top-bottom <= min_range
my $min_absolute = 0; # filter out all with max(abs(values)) <= min_absolute
my $pearson;
my $mingenes = 1;     # conditions with fewer values than this are removed
my $minconds = 1;     # genes with fewer values ...as above... (filter is applied after genes)


GetOptions("mapdims=s"=>\$mapdims,
	   "epochs=i"=>\$epochs,
	   "alpha=s"=>\$alpha,
	   "debug"=>\$debug,
	   "mapname=s"=>\$mapname,
	   "seed=i"=>\$seed,
           "bootstrap=i"=>\$numbootstraps,

	   "median_shift|medianshift"=>\$median_shift,
	   "min_range|minrange=s"=>\$min_range,
	   "min_absolute|minabsolute=s"=>\$min_absolute,
	   "pearson"=>\$pearson,


	   "printnumberofgenesandexit"=>\$printnumberofgenesandexit,

	   "mingenes=i"=>\$mingenes,
	   "minconds=i"=>\$minconds,

	   "includeconditions=s"=>\$includeconditions, # regexp, e.g. tissues|development (case insensitive)
	   "excludeconditions=s"=>\$excludeconditions, # regexp, e.g. tissues|development (case insensitive)

	  );


my @genes;
my %averages; # gene = perl array of expression values
my $ndims = 0;
my $headers = <>;
chomp($headers);
my @headers = split /\t/, $headers;
shift @headers;


while (<>) {
  next if (/^#/);
  chomp;
  my ($gene, @averages) = split /\t/, $_;
  my @good = grep { defined $_ && $_ ne 'NaN' && $_ ne ''; } @averages;
  if (@good) {
    if ($min_range > 0) {
      my ($min, $max) = minmax(pdl(@good));
      next unless ($max - $min > $min_range);
    }
    if ($min_absolute > 0) {
      my $maxabs = max(abs(pdl(@good)));
      next unless ($maxabs > $min_absolute);
    }

    $averages{$gene} = \@averages;
    my $n = scalar @averages;
    $ndims = $n if ($n>$ndims);
    push @genes, $gene;
  }
}


splice @genes, 100 if ($debug);

my %count;

# the following must be floating point so we can have NaN
my $data = zeroes $ndims, scalar @genes;
$data .= NaN;

for (my $i=0; $i<@genes; $i++) {
  my $gene = $genes[$i];

  for (my $j=0; $j<$ndims; $j++) {
    my $avg = $averages{$gene}[$j];
    if (defined $avg && $avg ne 'NaN' && $avg ne '') {
      $data->set($j, $i, $avg);
      $count{$headers[$j]}++;
    }
  }
  if ($median_shift) {
    my $this = $data->slice(":,$i");
    $this -= $this->where($this->isfinite)->median;
  }
}

my @ok_conditions;
my @bad_conditions;
for (my $j=0; $j<$ndims; $j++) {
  if ($count{$headers[$j]} >= $mingenes &&
      (!$includeconditions || $headers[$j] =~ /$includeconditions/i) &&
      (!$excludeconditions || $headers[$j] !~ /$excludeconditions/i)
      ) {
    push @ok_conditions, $j;
  } else {
    push @bad_conditions, $j;
    warn "removing $headers[$j]\n";
  }
}
if (@bad_conditions) {
  # remove the bad data (keep good)
  $data = $data->dice(\@ok_conditions);
  # remove the bad headings
  my @old_headers = splice @headers; # headers is now empty
  foreach my $j (@ok_conditions) {
    push @headers, $old_headers[$j];
  }
}

# now check for genes which now have too little data
my @ok_genes;
my @bad_genes;
for (my $i=0; $i<@genes; $i++) {
  if ($data->slice(":,$i")->isfinite->sum >= $minconds ) {
    push @ok_genes, $i;
  } else {
    push @bad_genes, $i;
    warn "removing $genes[$i]\n";
  }
}
if (@bad_genes) {
  # remove the bad data
  $data = $data->dice(X,\@ok_genes);
  # remove the bad genes;
  my @old_genes = splice @genes;
  foreach my $i (@ok_genes) {
    push @genes, $old_genes[$i];
  }
}



# theoretically one should repeat gene and condition filtering until convergence...
# but we won't...


if ($printnumberofgenesandexit) {
  print scalar(@genes), " genes passed filters\n";
  print scalar(@headers), " conditions passed filters\n";
  exit;
}


my $index = 0;
foreach my $header (@headers) {
  printf "# %-2d %-5d %s\n", $index++, $count{$header}, $header;
}



# for (my $i=0; $i<@genes; $i++) {
#   my $row = $data->slice(":,$i");
#   printf "%-10s %4d %s\n", $genes[$i], sum($row==1), join '', map { sprintf "%4s", $_ } $row->list;
# }


for (my $bootstrap=0; $bootstrap<=$numbootstraps; $bootstrap++) {

  my $som = new PDL::Kohonen();
  my @mapdims = split /\D/, $mapdims;
  die "sorry, map must be two-dimensional\n" unless (@mapdims == 2);
  srand($seed + $bootstrap) if (defined $seed);
  $som->init($data, @mapdims);

  $som->train($data, { epochs => $epochs, progress=>'on',
		       alpha => $alpha,
		       winnerfunc=>$pearson ? 'pearson' : 'euclidean'
		     });


  #
  # save map to database
  #

  my $savedmap = ExpressionMap::Map->insert({
					     name => $mapname,
					     xsize => $mapdims[0],
					     ysize => $mapdims[1],
					     # seed => $seed,
					     bootstrap => $bootstrap,
					     commandline => $commandline,
					    });

  my $i = 0;
  foreach $tsumm (@headers) {
    ExpressionMap::VectorDimension->insert({
					    map => $savedmap,
					    idx => $i++,
					    name => $tsumm,
					    count_pos => $count{$tsumm}  });
  }


  for (my $x=0; $x<$mapdims[0]; $x++) {
    for (my $y=0; $y<$mapdims[1]; $y++) {
      my $mapnode = ExpressionMap::MapNode->insert({
						    map=>$savedmap,
						    'x'=>$x, 'y'=>$y });
    }
  }
  my @vecdims = $savedmap->vecdims;

  # only save weights for "main" map
  if ($bootstrap == 0) {
    foreach my $mapnode ($savedmap->mapnodes) {
      foreach my $vecdim (@vecdims) {
	ExpressionMap::Weight->insert({
				     mapnode=>$mapnode,
				     vecdim=>$vecdim,
				     weight=>$som->at($vecdim->idx,
						      $mapnode->x,$mapnode->y)
				    });
      }
    }
  }

  print STDERR "wrote map '".$savedmap->name."' to database (bootstrap = $bootstrap)\n";


  #
  # apply map to data and save that
  #

  my ($locs, $errors) = $som->apply($data,
				    {progress=>'on',
				     winnerfunc=>$pearson ? 'pearson' : 'euclidean'
				    },
				   );


  # wcols "%10g", $data->mv(-1,0)->dog,  $locs->mv(-1,0)->dog, $errors;

  my @out;

  for (my $i=0; $i<@genes; $i++) {
    my $row = $data->slice(":,($i)");
    my $loc = $locs->slice(":,$i");
    my ($x, $y) = $loc->list;

    push @out, sprintf "%03d %03d %-10s %2d %.4f : %s\n", $x, $y, $genes[$i], sum($row==1), $errors->at($i), join('', map { sprintf "%7.4f", $_ } $row->list) if ($numbootstraps == 0);

    my $mapping = ExpressionMap::Mapping->insert({
						  map=>$savedmap,
						  mapnode=>$savedmap->mapnodes('x'=>$x, 'y'=>$y)->first,
						  gene_id=>$genes[$i],
						  error=>$errors->at($i)
						 });

    # only save bit vector info for main mapping
    if ($bootstrap == 0) {
      foreach my $vecdim (@vecdims) {
	my $value = $row->at($vecdim->idx);
	$value = undef if ($value eq 'nan' || $value eq 'inf');
	ExpressionMap::ExpressionValue->insert({
						mapping=>$mapping,
						vecdim=>$vecdim,
						value=>$value
					       });
      }
    }

  }

  print sort @out if ($numbootstraps == 0);
}
