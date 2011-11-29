#!/usr/bin/perl -w
#                           -*- mode: cperl -*-
use CGI qw(:standard);
use CGI::Ajax;

#use GO::AppHandle;
use PDL;

use lib './';
use ExpressionMap::Map;
use GeneOntology::Term;
use Statistics::Distributions;

use Config::General;

# same config file is also loaded by ExpressionMap::Map
# make sure this is protected by .htaccess
my $conf = new Config::General("config.txt");
my %config = $conf->getall;

my $map_name = $config{map_name} || die;
my $nonrandom_test_max = $config{nonrandom_test_max} || 500;
my $nonrandom_test_reps = $config{nonrandom_test_reps} || 100;
my $mapheight = $config{map_height} || 300;
my $maxmapheight = 1200;
my $max_e_selected = $config{max_e_selected} || 3;
my $gene_char = $config{gene_char} || 'O';

my %comparator_text = ( gt => 'greater than', lt => 'less than',
			gte => 'greater than or equal to', lte => 'less than or equal to');

#
# process CGI params
#

if (defined param('r')) {
  my ($id) = param('r') =~ /(\w+)/;
  if ($id) {
    print redirect("http://funcgen.vectorbase.org/ExpressionData/gene/$id");
    exit;
  }
}




my %q_genes;
if (defined param('g')) {
  map { $q_genes{$_} = 1 if ($_) } map { split /\W+/, $_ } param('g');
}

my @f_ids;
my $f_query;
if (defined param('f')) {
  my %f_ids;
  foreach my $f (param('f')) {
    map { $f_ids{$_} = 1 } $f =~ /(IPR\d+|GO:\d+)/g;
  }
  @f_ids = sort keys %f_ids;

  if (@f_ids == 0) {
    ($f_query) = param('f') =~ /(\b[\w ]{3,}\b)/;
  }
}

my @e_indices;
my @e_comps;
my @e_thresholds;

foreach my $i (1 ..3) {
  if (defined param("e$i")) {
    my ($index) = param("e$i") =~ /(-?\d+)/;
    if (defined $index && $index >= 0) {
      push @e_indices, $index;
      push @e_comps, param("e_comp$i") || 'gt';
      my $thresh = 0;
      if (defined param("e_thresh$i") && param("e_thresh$i") =~ /(-?[0-9.]+)/) {
	$thresh = $1;
      }
      push @e_thresholds, $thresh;
    }
  }
}

my $e_logic = 'AND';
if (defined param('e_logic') && param('e_logic') eq 'OR') {
  $e_logic = 'OR';
}
my $f_logic = 'OR';
if (defined param('f_logic') && param('f_logic') eq 'AND') {
  $f_logic = 'AND';
}

my ($x, $y);
if (defined param('x') && defined param('y')) {
  ($x) = param('x') =~ /(\d+)/;
  ($y) = param('y') =~ /(\d+)/;

  if (defined param('e')) {
    @e_indices = grep { $_ >= 0 } map { split /[^0-9-]/, $_ } param('e');
    if (@e_indices > $max_e_selected) {
      @e_indices = ();
    }
  }
} elsif (@e_indices == 0 && param('e')) {
# http://base.vectorbase.org/emtest/vb105p/?Highlight+checked+expression%2Ffunction+on+map...=Highlight+checked+expression%2Ffunction+on+map...&e=0+gt++0.796&e=60+gt++0.160&g=&tt=&h=400
  my $i = 0;
  foreach my $e (param('e')) {
    $i++;
    my ($index, $comp, $thresh) = split " ", $e;
    if (defined $thresh) {
      param("e$i", $index);
      param("e_comp$i", $comp);
      param("e_thresh$i", $thresh);
      push @e_indices, $index;
      push @e_comps, $comp;
      push @e_thresholds, $thresh;
    }
  }
}



if (defined param('h') && param('h') =~ /(\d+)/ && $1 <= $maxmapheight) {
  $mapheight = $1;
}

#
# get the map object
#

my $map = ExpressionMap::Map->retrieve(name=>$map_name, bootstrap=>0);
my ($xmax, $ymax) = ($map->xsize-1, $map->ysize-1);

#
# get the basic map stats
#

my %index2tsumm;
my %tsummcount;
my %index2vecdim;

foreach my $vecdim ($map->vecdims) {
  my $index = $vecdim->idx;
  my $count = $vecdim->count_pos;
  my $tsumm = $vecdim->name;
  $tsummcount{$index} = $count;
  $index2tsumm{$index} = $tsumm;
  $index2vecdim{$index} = $vecdim;
}

#
# read in biomart annotations
#

open MART, "biomart.txt";
my $headers = <MART>;
chomp($headers);
my @headers = split /\t/, $headers;
my %header_idx;
for (my $i=0; $i<@headers; $i++) {
  $header_idx{$headers[$i]} = $i;
}
my %biomart;
my %iprdesc;
my %godesc;

my %genedesc;
my %genesymbol;

my $hi_iprid = $header_idx{"Interpro ID"};
my $hi_iprdesc = $header_idx{"Interpro Description"};
my $hi_goid = $header_idx{"GO ID"};
my $hi_godesc = $header_idx{"GO description"};

#my $apph = GO::AppHandle->connect(-dbname=>'go_200807',
#				  -dbuser=>'maccallr',
#				  -dbauth=>'ca66age',
#				  -dbhost=>'bio-iisrv1.bio.ic.ac.uk',
#				  );

while (<MART>) {
  chomp;
  my @cols = split /\t/, $_;
  my $id = $cols[$header_idx{"Ensembl Gene ID"}];
  die "no 'Ensembl Gene ID' column in biomart data\n" unless (defined $id);

  my $goacc = $cols[$header_idx{'GO ID'}];

  if ($goacc) {
    $biomart{$id}{'GO ID'}{$goacc} = 1;
#    my $term = $apph->get_term({acc=>$goacc});
#    if (defined $term) {
#      my $parents = $apph->get_parent_terms($term);
#      foreach $parent (@$parents) {
#	$biomart{$id}{'GO ID'}{$parent->acc} = 1;
#      }
#    }
  }

  my $ipracc = $cols[$header_idx{'Interpro ID'}];
  $biomart{$id}{'Interpro ID'}{$ipracc} = 1 if ($ipracc);

  $genedesc{$id} = $cols[$header_idx{Description}] || '' unless ($genedesc{$id});
  $genesymbol{$id} = $cols[$header_idx{'Anopheles symbol'}] || '' unless ($genesymbol{$id});

  $iprdesc{$cols[$hi_iprid]} = $cols[$hi_iprdesc] if (defined $cols[$hi_iprid]);
  $godesc{$cols[$hi_goid]} = $cols[$hi_godesc] if (defined $cols[$hi_goid]);
}
close(MART);

my $tt = 'whole';

if ($f_query) {
  my ($perlregexp, $mysqlregexp);

  if (param('tt') && param('tt') eq 'partial') {
    $tt = 'partial word';
    $perlregexp = qr/$f_query/;
    $mysqlregexp = $f_query;
  } else {
    $perlregexp = qr/\b$f_query\b/i;
    $mysqlregexp = '[[:<:]]'.$f_query.'[[:>:]]';
  }

  foreach my $iprid (keys %iprdesc) {
    if ($iprdesc{$iprid} =~ /$perlregexp/) {
      push @f_ids, $iprid;
    }
  }
  push @f_ids, map { $_->acc } GeneOntology::Term->search_name_regexp($mysqlregexp);
}

my @orig_f_ids = @f_ids; # before child expansion

#
# expand to include child terms
#
if (@f_ids) {
  my %f_ids;
  foreach my $f (@f_ids) {
    $f_ids{$f} = 1;
    my $go = GeneOntology::Term->retrieve(acc=>$f);
    $godesc{$f} = $go->name if (defined $go);

    if ($f =~ /^GO:\d+$/) {
      foreach my $child (GeneOntology::Term->search_descendants($f)) {
	# if we've seen the child in the biomart file
	# then we can add it as a query term
	if ($godesc{$child->acc}) {
	  $f_ids{$child->acc} = 1;
	}
	## $godesc{$child->acc} = $child->name;
      }
    }
  }
  @f_ids = sort keys %f_ids;
}


#
# heading
#
print header(), start_html( -style=>{-src=>"stylesheet.css"},
			     -title => "VectorBase Expression Map (BETA)");

#
# some AJAX magic
#

my $pjx = new CGI::Ajax( 'bootstrap_pair' => 'bootstrap_pair.cgi' );
print $pjx->show_javascript();

my $pjx2 = new CGI::Ajax( 'bootstrap_neighbours' => 'bootstrap_neighbours.cgi' );
print $pjx2->show_javascript();

print div({-style=>'float: left;'}, h1("<a class=vblink href=\"http://www.vectorbase.org\">VectorBase</a> Gene Expression Map<sub style='color: blue; font-size: 66%;'>BETA</sub> for <i>Anopheles gambiae</i>"));


print div({-style=>'float: right; margin-right: 20px; '},
div( {-style=>"text-align: right"}, a({-href=>"http://www.vectorbase.org/Help/Expression_Map_BETA"}, "Help")),br,
"Name: ",b($map->name()),br,
"Dimensions: ",b("W=".$map->xsize()." H=".$map->ysize()),

	 );


print br({-style=>'clear: both;'});

if (@e_thresholds) {
  print div({-id=>'patient', -style=>"text-align: center; color: red;"}, "Please be patient while waiting for this page to load.  Thousands of expression values are being processed.");
}

my @seen_genes;

if (defined $x && defined $y) {
  #
  # print node details
  #

  #
  # first load SOM vector for this node
  #
  my @vector; # fill this in arbitrary order, hope there are no gaps!

  my $mapnode = $map->mapnodes('x'=>$x, 'y'=>$y)->first;
  foreach my $weight ($mapnode->weights) {
    $vector[$weight->vecdim->idx] = $weight->weight;
  }

  my @mappings = $mapnode->mappings(); # map->mappings(mapnode=>$mapnode);
  my $n = @mappings;

  #
  # some javascript
  #

  print <<"EOJS";
<script type="text/javascript">
function clearAll() {
  form = document.forms[0];
  nelem = form.elements.length;
  var i=0;
  for (i=0; i<nelem; i++) {
    if (form.elements[i].name == 'e' ||
        form.elements[i].name == 'f') {
      form.elements[i].checked = 0;
    }
  }
}

function limitNumberOfSelections(checkbox) {
  var count = 0;
  var form = document.forms[0];
  var nelem = form.elements.length;
  for (var i=0; i<nelem; i++) {
    if (form.elements[i].name == 'e' && form.elements[i].checked == true) count++;
  }
  if (count > $max_e_selected) {
    alert("Sorry, you are not allowed to make more than $max_e_selected selections.");
    checkbox.checked = false;
    return false;
  }
  return true;
}
</script>
EOJS


  #
  # now build the table
  #

  param('f', @orig_f_ids);

  my @e_checkbox_vals;
  foreach my $i (@e_indices) {
    my $val = sprintf "%6.3f", $vector[$i];
    my $comp = $val >= 0 ? 'gt' : 'lt';
    push @e_checkbox_vals, "$i $comp $val";
  }
  param('e', @e_checkbox_vals);

  my @gene_ids;
  my @rows;

  my @ths = ( th('&nbsp;'), th("Experiment and condition"), th("Weight") );
  foreach my $mapping (@mappings) {
    my $id = $mapping->gene_id();
    push @ths, th({-class=>'gene'}, gene_symbol_link($id));
    push @gene_ids, $id;
  }

  push @rows, Tr({-valign=>'top'}, @ths);

  my @sorted_indices = sort { $vector[$b] <=> $vector[$a] } 0 .. @vector-1;

  foreach my $i (@sorted_indices) {
    my $vecdim = $index2vecdim{$i};
    my @tds;
    my $style = 'font-weight: bold;';
    my $val = sprintf "%6.3f", $vector[$i];
    my $comp = $val >= 0 ? 'gt' : 'lt';
    push @tds, td(checkbox(-name=>'e', -value=>"$i $comp $val", -label=>'', -onchange=>'limitNumberOfSelections(this)'));
    push @tds, td(div({-style=>"width: 450px; min-height: 200%;"}, $index2tsumm{$i}));
    $val =~ s/ /&nbsp;/g;
    push @tds, td({-style=>$style, -class=>'c'}, $val);

    foreach my $mapping (@mappings) {
      my $id = $mapping->gene_id();
      my $val = $mapping->expression_values(vecdim=>$vecdim)->first->value;
      if (!defined $val) {
	$val = '?';
      } else {
	$val = sprintf "%5.2f", $val;
	$val =~ s/ /&nbsp;/g;
      }
      push @tds, td({-class=>'c'}, $val);
    }

    push @rows, Tr({-onmouseover=>"this.style.background = '#eee';", -onmouseout=>"this.style.background = 'white';"}, @tds)."\n";
  }

  #
  # now count GO and interpro occurrences per gene
  #

  my %f_count;  # interpro or GO ID => count
  foreach my $mapping (@mappings) {
    my $id = $mapping->gene_id;
    foreach my $f (keys %{$biomart{$id}{'Interpro ID'}}) {
      $f_count{$f}++ if ($f);
    }

    # now start with each assigned GO term and 
    my %go_accs;
    foreach my $f (keys %{$biomart{$id}{'GO ID'}}) {
      $go_accs{$f} = 1;
      foreach my $ancestor (GeneOntology::Term->search_ancestors($f)) {
	$go_accs{$ancestor->acc} = 1;
	$godesc{$ancestor->acc} = $ancestor->name;
      }
    }
    foreach my $go (keys %go_accs) {
      $biomart{$id}{"GO ID"}{$go} = 1;
      $f_count{$go}++ if ($go);
    }
  }

  my %f_pvalue; # interpro or GO ID => corrected pvalue from ermineJ
  foreach my $go_ora_result ($mapnode->go_ora_results()) {
    my $pval = sprintf "%5.1g", $go_ora_result->corrected_pvalue();
    $pval =~ s/\s//g;
    $f_pvalue{$go_ora_result->term} = $pval;
  }

  #
  # and add data to table
  #

  @ths = ( th('&nbsp;'), th("Interpro or GO ID",
 #span({-style=>'font-weight: normal;'}, ' ( checkbox combination: ',
 # radio_group(-name => 'f_logic', -values => [ 'OR', 'AND' ]), ')'
			   ), th("Count<br>(GO p-value)") );
  foreach my $mapping (@mappings) {
    my $id = $mapping->gene_id();
    push @ths, th({-class=>'gene'}, gene_symbol_link($id));
  }
  push @rows, Tr({-valign=>'top'}, @ths);

  foreach my $f (sort { $f_count{$b} <=> $f_count{$a} } keys %f_count) {
    my @tds;
    push @tds, td(checkbox(-name=>'f', -value=>$f, -label=>''));
    my $descrip = $godesc{$f} || $iprdesc{$f} || '?';
    push @tds, td(f_links($f)), td($f_count{$f}.(defined $f_pvalue{$f} && $f_count{$f} > 1 ? $f_pvalue{$f}<0.05 ? "<br><font color=red>($f_pvalue{$f})</font>" : "<br>($f_pvalue{$f})" : ''));

    foreach my $mapping (@mappings) {
      my $id = $mapping->gene_id;
      my $bit = '&nbsp;';
      if ($biomart{$id}{"GO ID"}{$f} ||
	  $biomart{$id}{"Interpro ID"}{$f}) {
	$bit = gene_link($id);
      }
      push @tds, td({-class=>'gene'}, $bit);
    }
    push @rows, Tr({-onmouseover=>"this.style.background = '#eee';", -onmouseout=>"this.style.background = 'white';"}, @tds);
  }

  print h2("Expression and function details for ".a({-href=>"#listofgenes"}, "$n genes")." in map node/cluster $x,$y");

  my $reset_button = span({-onclick=>'clearAll()', -onmouseover=>"this.style.background = '#eee';", -onmouseout=>"this.style.background = '#ddd';", -style=>'margin-left: 10px; background: #ddd; padding: 4px;'}, "clear all checkboxes");

  print start_form(-method=>'GET'),
    submit('Highlight checked expression/function on map...'),
      $reset_button, br, br,
      table({-class=>'nodeinfo'}, @rows),
	submit('Highlight checked expression/function on map...'),
	  $reset_button, hidden(-name=>'g'), hidden(-name=>'tt'), hidden(-name=>'h'),
	    end_form;

  print a({-name=>'listofgenes'}), h2("List of $n genes mapping to this node/cluster");
  print table({-class=>'geneidlist'}, map Tr(td($_).td($genesymbol{$_}).td($genedesc{$_})), map $_->gene_id, @mappings);

} else {

  #
  # print whole map
  #

  print <<"EOJS";
<script type="text/javascript">
function getElementsByCondition(condition,container) {
  container = container || document;
  var all = container.all || container.getElementsByTagName('*');
  var arr = [];
  for(var k=0;k<all.length;k++) {
    var elm = all[k];
    if(condition(elm,k));
    arr[arr.length] = elm;
  }
  return arr
}

function toggleNotSure() {
  var hidden = getElementsByCondition(
    function(el) {
      if(el.className=='n') {
        if (el.style.display=='none') {
          el.style.display='inline';
        } else {
          el.style.display='none';
        }
        return el;
      }
    }
  );
}

function limitNumberOfSelections(dropdown) {
  var count = 0;
  for (var i = 0; i < dropdown.options.length; i++) {
    if (dropdown.options[i].selected==true) count++;
  }
  if (count > $max_e_selected) {
    alert("Sorry, you are not allowed to make more than $max_e_selected selections.");
    for (var i = 0; i < dropdown.options.length; i++) {
      dropdown.options[i].selected = false;
    }
    return false;
  }
  return true;
}
</script>
EOJS

  #
  # write out table
  #


  my @trs;
  my ($have_hilited_gene, $have_hilited_function, $have_hilited_expression, $have_hilited_notsure) = (0,0,0,0);

  my @e_hit_coords;
  my @f_hit_coords;
  my @g_hit_coords;
  my @f_geneids;
  my @e_geneids;
  my @all_coords;
  my @all_e_avail_coords;
  my %seen_genes;

  my @e_vecdims = map $index2vecdim{$_}, @e_indices;
  my $tot_e = @e_indices;
  my $tot_f = @f_ids;

  for $y (0 .. $ymax) {
    my @tds;
    my @ths;
#    $y = sprintf "%03d", $y;
    for $x (0 .. $xmax) {
#      $x = sprintf "%03d", $x;

      my $mapnode = $map->mapnodes('x'=>$x, 'y'=>$y)->first;
      my @mappings = $mapnode->mappings(); # map->mappings(mapnode=>$mapnode);

      my @tdcontent;
      my %tdclasses;

      foreach $mapping (@mappings) {
	my $id = $mapping->gene_id;
	my @classes;
	push @all_coords, [ $x, $y ];

	# flag gene as 'selected'
	if ($q_genes{$id} || ($genesymbol{$id} && $q_genes{$genesymbol{$id}})) {
	  push @classes, 'g';
	  $tdclasses{g} = 1;
	  $have_hilited_gene++;
	  push @g_hit_coords, [$x, $y];
	  $seen_genes{$id} = 1;
	}
	# flag gene as having desired expression
	if (@e_indices) {
	  my ($e_hits, $e_nans) = (0,0);
	  for (my $i=0; $i<@e_vecdims; $i++) {
	    my $vecdim = $e_vecdims[$i];
	    my $e_comp = $e_comps[$i];
	    my $e_thresh = $e_thresholds[$i];
	    my $val = $mapping->expression_values(vecdim=>$vecdim)->first->value;
	    if (!defined $val) {
	      $e_nans++ }
	    else {
	      $e_hits++ if (($e_comp eq 'gt' && $val > $e_thresh) ||
			    ($e_comp eq 'lt' && $val < $e_thresh) ||
			    ($e_comp eq 'gte' && $val >= $e_thresh) ||
			    ($e_comp eq 'lte' && $val <= $e_thresh)
			   );
	    }
	  }
	  if ($e_hits == $tot_e) {
	    push @classes, 'e';
	    $have_hilited_expression++;
	    push @e_hit_coords, [ $x, $y ];
	    push @all_e_avail_coords, [ $x, $y ];
	    push @e_geneids, $id;
	  } elsif ($e_nans > 0) {
	    push @classes, 'n'; # not sure/nan
	    $have_hilited_notsure++;
	  } else {
	    push @all_e_avail_coords, [ $x, $y ];
	  }
	}

	# flag gene as having desired function

	if (@f_ids) {
	  my $f_hits = 0;
	  foreach $f (@f_ids) {
	    if ($biomart{$id}{"GO ID"}{$f} ||
			  $biomart{$id}{"Interpro ID"}{$f}) {
	      $f_hits++;
	    }
	  }
	  if ( ($f_logic eq 'AND' && $f_hits == $tot_f) ||
	       ($f_logic eq 'OR' && $f_hits > 0) ) {
	    push @classes, 'f';
	    $have_hilited_function++;
	    push @f_hit_coords, [ $x, $y ];
	    push @f_geneids, $id;
	  }
	}


#	push @tdcontent, a({-href=>"?r=$id", (@classes ? ('-class'=>"@classes") : ())}, "&#x25fc;");
	push @tdcontent, a({-href=>"?r=$id", (@classes ? ('-class'=>"@classes") : ())}, $gene_char);
      }

      push @tds, td({keys %tdclasses ? ('-class'=>join ' ', keys %tdclasses) : ()}, join ' ', @tdcontent);
      my $nodelabel = sprintf "%d,%d", $x, $y;
      my $g = join ',', keys %q_genes;
      my $e = join ',', @e_indices;
      push @ths, th(a({-href=>"?x=$x;y=$y;g=$g;f=__FIDS__;e=$e;h=$mapheight"}, $nodelabel));
    }
    push @trs, Tr(@ths), "\n";
    push @trs, Tr(@tds), "\n";
  }

  my $f = $f_query ? "$f_query;tt=$tt" : join ',', @orig_f_ids;
  grep { s/__FIDS__/$f/g } @trs;

  print div({-class=>"mapwindow", -style=>"height: ${mapheight}px;"}, table({-class=>'map'}, @trs));


  #
  # do some hit clustering stats
  #

  my $e_nonrandom_test = '';
  if ($have_hilited_expression > 1) {
    if ($have_hilited_expression <= $nonrandom_test_max) {
      $e_nonrandom_test = nonrandom_test(\@e_hit_coords, \@all_e_avail_coords, $nonrandom_test_reps);
    } else {
      $e_nonrandom_test = "Too many genes for non-random distribution score calculation (maximum $nonrandom_test_max)";
    }
  }

  my $f_nonrandom_test = '';
  if ($have_hilited_function > 1) {
    if ($have_hilited_function <= $nonrandom_test_max) {
      $f_nonrandom_test = nonrandom_test(\@f_hit_coords, \@all_coords, $nonrandom_test_reps);
    } else {
      $f_nonrandom_test = "Too many genes for non-random distribution score calculation (maximum $nonrandom_test_max)";
    }
  }

  my $g_nonrandom_test = '';
  if ($have_hilited_gene > 1) {
    if ($have_hilited_gene <= $nonrandom_test_max) {
      $g_nonrandom_test = nonrandom_test(\@g_hit_coords, \@all_coords, $nonrandom_test_reps);
    } else {
      $g_nonrandom_test = "Too many genes for non-random distribution score calculation (maximum $nonrandom_test_max)";
    }
  }


  #
  # key
  #

  print h2("Key to symbols:");

  my @tds;

  push @tds, td(table({-class=>'map inset'}, Tr(td({-class=>'inset'}, a({-href=>''}, $gene_char)))), "Gene", br, small(br. 'Each dot on the map is a gene.  Click for the gene\'s expression summary.'));

  push @tds, td(table({-class=>'map inset'}, Tr(th({-class=>'inset'}, a({-href=>''}, "13,8")))), "Map&nbsp;location", br, small(br, 'column,row - click on these headings for more information about the expression and function of genes in each map location'));

  @seen_genes = sort keys %seen_genes;
  if ($have_hilited_gene) {
    push @tds, td(table({-class=>'map inset'}, Tr(td({-class=>'g inset'}, a({-href=>'', -class=>'g'}, $gene_char)))), "Highlighted&nbsp;gene(s)", br, ($g_nonrandom_test ? br.$g_nonrandom_test.br : ''), small(br, join(" ", map gene_symbol_link($_), @seen_genes)));
  }

  if ($have_hilited_function) {
    my $yq = $f_query ? "matching your query '$f_query' " : '';
    push @tds, td(table({-class=>'map inset'}, Tr(td({-class=>'inset'}, a({-href=>'', -class=>'f'}, $gene_char)))),  "Highlighted Gene Ontology terms or InterPro domains $yq($have_hilited_function genes)<br>", ($f_nonrandom_test ? br.$f_nonrandom_test.br : ''), div({-class=>'keylist'}, join(" $f_logic ", map f_links($_), @orig_f_ids)));
  }

  if ($have_hilited_expression || $have_hilited_notsure) {

    my @not_available = (br, br, table({-class=>'map inset'}, Tr(td({-class=>'inset'}, a({-href=>'', -class=>'n x'}, $gene_char)))), "Expression data not available ($have_hilited_notsure genes)");
    push @not_available, br, small(checkbox(-name=>'click to hide/show these genes', -onclick=>'toggleNotSure()'), "(please wait a few seconds after clicking)");

    my @e_descrips;
    for (my $i=0; $i<@e_indices; $i++) {
      push @e_descrips, $index2tsumm{$e_indices[$i]}." <b>is $comparator_text{$e_comps[$i]}</b> $e_thresholds[$i]";
    }

    push @tds, td(table({-class=>'map inset'}, Tr(td({-class=>'inset'}, a({-href=>'', -class=>'e'}, $gene_char)))),  "Highlighted expression ($have_hilited_expression genes)<br>", ($e_nonrandom_test ? br.$e_nonrandom_test.br : ''), div({-class=>'keylist'}, join br.AND.br, @e_descrips),
		  $have_hilited_notsure ? @not_available : ()
		 );
  }

  print table({-class=>'key'}, Tr({-valign=>'top'}, @tds));

  #
  # print form below
  #

  param('f', $f_query ? $f_query : join(',', @orig_f_ids));
  param('e', @e_indices);

  print start_form(-method=>'GET');

  print h2("Customize map display:");

#  print submit("Submit");

  # make labels like "Blood meal time-series blah blah (1234 genes)"
  my %tsummlabels;
  grep { $tsummlabels{$_} = "$index2tsumm{$_} ($tsummcount{$_} genes)" } keys %index2tsumm;
  $tsummlabels{-1} = '-- no selection --';

  print "<table class=customizeform><tr>\n";
  print "<td class=submit style='width: 150px;'>\n";
  print submit("Redraw map...");
  print br,br,br,h3("Map display height: "), popup_menu(-name=>'h',
				       -values=>[ 300, 400, 500, 600, 700, 800 ],
				       -default=>$config{map_height});
  print "</td>\n<td style='width: 300px;'>\n";
  print h3("Highlight by gene:"), "Enter comma or space-separated gene IDs or symbols", br, textarea(-name=>'g',
				   -rows=>2,
				   -columns=>24,
				  );
  print br."[e.g.: ".join(', ', map a({-href=>"?h=$mapheight;g=$_"}, $_), qw(AGAP007000 AGAP000106), 'CTLMA2,CTL4')."]";

  print "</td><td>\n";
  print h3("Highlight by annotated function:"), "Enter comma or space-separated GO IDs and/or Interpro IDs, or <b>one</b> Interpro/GO text query", br, textarea(-name=>'f',
						     -rows=>2,
						     -columns=>24,
						    );
  print br."[e.g.: ".join(', ', map a({-href=>"?h=$mapheight;f=$_"}, $_), qw(IPR000488 GO:0005549 GO:0003676 helicase clathrin), 'histone modification')."]";
#  print br, radio_group(-name => 'f_logic',
#			-values => [ 'OR', 'AND'],
#		       );
  print br, "Treat text query as: ", radio_group(-name => 'tt',
			-values => [ 'whole', 'partial' ]), ' word';

  print "</td>\n</tr>\n";
  print "<tr><td colspan=3>\n";

  print h3("Highlight by expression value:"),

    join br.'AND'.br, map { popup_menu(-name=>"e$_",
				       -values => [ -1, sort { $index2tsumm{$a} cmp $index2tsumm{$b} } keys %index2tsumm ],
				       -labels => \%tsummlabels,
				       -size=>5,
				       -style=>'font-size: 9px;').br.popup_menu(-name=>"e_comp$_", -values=>[ qw(gt gte lt lte) ],
										-labels=> \%comparator_text,
										-style=>'font-size: 9px;').textfield(-name=>"e_thresh$_", -size=>5, -default=>0.5) } (1 .. 3);

  print br.br.submit("Redraw map...");
  print "</td></tr>\n";
  print "</table>\n";
  print hidden(-name=>'verbose');
  print end_form;

  if (param('verbose') && @f_geneids) {
    print p(h3("Genes highlighted by function"),
	    small(join ',<wbr>', @f_geneids));
  }

  if (param('verbose') && @e_geneids) {
    print p(h3("Genes highlighted by expression"),
	    small(join ',<wbr>', @e_geneids));


    print p(h3("Gene expression expression"),
	    join ' && ', map { "{$e_indices[$_]} $e_comps[$_] $e_thresholds[$_]" } ( 0 .. @e_indices-1 )
	   );
  }
}

print br, h2("Co-clustering bootstrap scores:");

print p('The mapping procedure involves a random initiation step, which means that the map shown above is just one of a very large number of possible mappings.  We have repeated the mapping many times and can tell you (below) how often a pair of genes end up within a certain (Euclidean) distance of each other on the map.  (A distance cutoff of zero means the same map node.)');

print p('Please use only gene IDs below  (e.g. AGAP001234).');

print p('Please allow up to 10 seconds for the results to appear.');

print <<"EOJS";
<script type="text/javascript">
  function clearbyid(id) {
    var obj = document.getElementById(id);
    if (obj != null) obj.innerHTML = '';
  }
</script>
EOJS

print "\n",  start_form, h3('Gene pair'),
  'Gene ID 1: ', textfield(-name=>'gene1', -size=>10, -default=>$seen_genes[0] || '', -onkeyup=>"clearbyid('bootstrap_result')"),
  ' Gene ID 2: ', textfield(-name=>'gene2', -size=>10, -default=>$seen_genes[1] || '', -onkeyup=>"clearbyid('bootstrap_result')"), ' Distance cutoff: ', textfield(-name=>'dist', -size=>3,  -onkeyup=>"clearbyid('bootstrap_result')", -default=>0), '&nbsp;&nbsp;',
  button(-name=>' = ', -onClick=>"bootstrap_pair(['gene1', 'gene2', 'dist'], ['bootstrap_result'])"), '&nbsp;'x2, span({-id=>'bootstrap_result', -style=>'font-weight: bold;' }, ''), end_form;

print "\n", start_form, h3('Gene neighbours'),
  'Gene ID &nbsp;: ', textfield(-name=>'ngene', -size=>10, -default=>$seen_genes[0] || '', -onkeyup=>"clearbyid('bootstrap_neighbours')"),
  ' Distance cutoff: ', textfield(-name=>'ndist', -size=>3,  -onkeyup=>"clearbyid('bootstrap_neighbours')", -default=>0), '&nbsp;&nbsp;',
  button(-name=>'calculate...', -onClick=>"bootstrap_neighbours(['ngene', 'ndist'], ['bootstrap_neighbours'])"), br, div({-id=>'bootstrap_neighbours'}, ''), end_form;

print br, br, a({href=>'.'}, "[Reset all filters / start again]"),
br, a({-href=>'map.cgi'}, "[Newer version]"),
br, a({-href=>'..'}, "[More maps (other versions and species)...]");

if (@e_thresholds) {
  print << 'EOJS';
<script type="text/javascript">
clearbyid('patient');
</script>
EOJS
}

print end_html;


sub f_links {
  my $id = shift;

  if ($iprdesc{$id}) {
    return a({-href=>"http://www.ebi.ac.uk/interpro/DisplayIproEntry?ac=$id"}, $id)." ($iprdesc{$id})";
  } elsif ($godesc{$id}) {
    return a({-href=>"http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=$id"}, $id)." ($godesc{$id})";
  } else {
    return $id;
  }
}

sub gene_link {
  my $id = shift;
  my $split = $id;
  $split =~ s/AGAP/AGAP<wbr>/;
  return a({-href=>"http://www.vectorbase.org/Search/Keyword/?quicksearch=true&term=$id"}, $split);
}

sub gene_symbol_link {
  my $id = shift;
  return gene_link($id).($genesymbol{$id} ? " ($genesymbol{$id})" : '');
}


sub nonrandom_test {
  my ($hits_aref, $allcoords_aref, $reps) = @_;
  my $numhits = @$hits_aref;
  my $numall = @$allcoords_aref;

  my $hits = pdl($hits_aref);
  my $realmeasure = measure($hits);

  my $oldseed = int rand 654321;
  srand($numhits);

  my $sumfracs = 0;

  for (1 .. $reps) {
    my %randcoords;
    for (1 .. $numhits) {
      my $index = int rand $numall;
      $randcoords{$index} = $allcoords_aref->[$index];
    }
    while (keys %randcoords < $numhits) {
      my $index = int rand $numall;
      $randcoords{$index} = $allcoords_aref->[$index];
    }

    my $rhits = pdl(values %randcoords);
    my $randmeasure = measure($rhits);
    $sumfracs++ if ($randmeasure < $realmeasure);
  }

  srand($oldseed);

  my $sigfigs = sprintf "%.0f", log($nonrandom_test_reps)/log(10);
  my $p = sprintf "%.${sigfigs}f", $sumfracs/$reps;
  $p = "&lt;".(1/$reps) if ($p == 0);
  return "Empirical non-random distribution p-value: <b>$p</b>";
}

#
# for each gene, calculate closest distance to other genes
# then return the average of all these
#
sub measure {
  my $hits = shift;
  my $n = $hits->getdim(1);
  my $cityblockdists = sumover(abs($hits->dummy(2,$n) - $hits->dummy(1,$n)));
  return $cityblockdists->pctover(1/($n-1))->avg();
}
