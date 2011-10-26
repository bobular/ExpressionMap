#!/usr/bin/perl -w
#                           -*- mode: cperl -*-
use CGI qw(:standard);

#use GO::AppHandle;
#use PDL;

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

#
# get the map object
#

my $map = ExpressionMap::Map->retrieve(name=>$map_name, bootstrap=>0);
my ($map_id, $map_height, $map_width) = ($map->id(), $map->ysize(), $map->xsize());

#
# heading
#
print header(), start_html( -style=>{-src=>"stylesheet.css"},
			    -script=>[
				      { -type => 'text/javascript',
					-src => 'https://ajax.googleapis.com/ajax/libs/prototype/1.7.0.0/prototype.js'
				      },
				      { -type => 'text/javascript',
					-src => 'jsonpath-0.8.0.js'
				      },
				      { -type => 'text/javascript',
					-src => 'raphael-min.js'
				      },
				      { -type => 'text/javascript',
					-src => 'scriptaculous.js?load=effects,controls,slider,dragdrop'
				      },
				      { -type => 'text/javascript',
					-src => 'config.js'
				      },
				      { -type => 'text/javascript',
					-src => 'expressionmap.js'
				      },
				     ],
			    -title => "VectorBase Expression Map");

print div({-style=>'float: left;'}, h1("<a class=vblink href=\"http://www.vectorbase.org\">VectorBase</a> Gene Expression Map for <i>Anopheles gambiae</i>"));


print div({-style=>'float: right; margin-right: 20px; '},
"Name: ",b($map->name()),br,
					span({-id=>'num-genes'}, '').br,
					span({-id=>'num-conditions'}, '')
				 # "Dimensions: ",b("W=$map_width H=$map_height"),
	 );

#
# build a hash of node sizes (number of genes) and record the maximum
#

my @nodesize; # [ x ] [ y ] -> num genes
my $maxsize = 0;
for (my $y=0; $y<$map_height; $y++) {
	for (my $x=0; $x<$map_width; $x++) {
		my $mapnode = $map->mapnodes('x'=>$x, 'y'=>$y)->first;
		my $size = $mapnode->mappings()->count();
		$nodesize[$x][$y] = $size;
		$maxsize = $size if ($size > $maxsize);
	}
}
#my $maxsizesqrt = sqrt($maxsize);

print "<br clear='both' />\n";

print "<div style='float: left;'>
<div id='please-wait'>Please wait a few seconds while the map is being initialised...
<br/>
<br/>
<span style='font-size: 70%'>If this message does not disappear, you may need to try a newer web browser.
<br/>
Internet Explorer 9 and most recent versions of Firefox, Safari and Chrome should work.
<br/>
If you continue to have problems, you can use the <a href='index.cgi'>beta version</a> instead.
</span></div>
<table class='newmap' id='map1'>\n";
for (my $y=0; $y<$map_height; $y++) {
	print "<tr id='row-$y'>\n";
	for (my $x=0; $x<$map_width; $x++) {
		my $numgenes = $nodesize[$x][$y];
#		my $radius = int(15*sqrt($numgenes)/$maxsizesqrt);
		print "<td id='node-$x-$y' numgenes=$numgenes></td>\n";
#		print "<td id='node-$x-$y' numgenes=$numgenes><svg title='node $x,$y ($numgenes genes, click for more)' xmlns='http://www.w3.org/2000/svg' version='1.1'><circle cx=15 cy=15 r=$radius fill='#ccc' /></svg></td>\n";
	}
	print "</tr>\n";
}
print "</table></div>\n";


# make some preparations for the control panel

foreach my $vecdim ($map->vecdims) {
  my $index = $vecdim->idx;
  my $count = $vecdim->count_pos;
  my $tsumm = $vecdim->name;
  $tsummcount{$index} = $count;
  $index2tsumm{$index} = $tsumm;
  $index2vecdim{$index} = $vecdim;
}

$index2tsumm{-1} = '-- select a condition or start typing in the box above --';

# print the control panel
print div({-style=>"float: left;", -class=>"control-panel"},
					start_form(-id=>'control-panel'),
					b("Highlight genes:").br."Enter gene ID(s) and/or GO ID(s) and/or Interpro ID(s) <b>or</b> <i>one</i> GO/Interpro text query (e.g. 'RNA processing', 'cuticle').".br,

					join(br, map {
						div({-class=>"search$_ search"}, b("Gene query $_:")).
								textarea(-name=>"search$_", -id=>"search$_", -class=>"search").
								button(-name=>'go', -class=>"search", -onclick=>"em_handle_search_submit($map_id, $_)")
					} (1..3)),
					# span({-id=>'mouse-over-help'}, 'Mouse-over yellow/orange/pink rectangles for gene info.'),
					br() x 3,
					b("Highlight areas by node-averaged expression:"),br,

					join(br, map {
						div({-class=>"efilter$_ efilter"}, b("Condition $_:")).
						textfield(-id=>"efiltertext$_", -class=>'text', -autocomplete=>'off').br.
						popup_menu(-class=>'ehighlight',
							   -name=>"e$_", -id=>"e$_",
							   -values => [ sort { $index2tsumm{$a} cmp $index2tsumm{$b} } keys %index2tsumm ],
							   -labels => \%index2tsumm,
							   # -size=>5,
							   -style=>'font-size: 9px;',
							   -onchange=>"em_handle_condition_change($map_id,this,$_)"
							  ).br.

						div({-id=>"track$_", -class=>'slider-track'},

								div({-id=>"handle$_", -class=>"efilter$_ slider-handle"}, '')
							 ).div({-class=>"efilterinfo"},
										 span({-id=>"efilterval$_"}, '0'),
										 span({-id=>"efilternodes$_"}, ''))

						 } (1 .. 3)),
					reset(-onclick=>'em_reset_all()'),
					end_form
)."\n";

print br({-clear=>"both"})."\n";

print h2("Tutorial (use case examples)");

print p("Note, these examples are based on the <a href=\"http://funcgen.vectorbase.org/ExpressionMap/Anopheles_gambiae/current/map.cgi\">Anopheles gambiae</a> expression map.");

print h3("Your favourite gene"),
  p("Enter the VectorBase gene ID (e.g. AGAP001234) or official gene symbol (e.g. CLIPA8) into the yellow 'Gene query 1' text box and hit 'go'.  The node (cluster) containing the gene will be highlighted with a yellow sausage.  You can then click on this node to see a popup summary of the genes, their expression and any enrichment of function in the cluster.  From here you could perform another query with an enriched GO term in the orange box, or you could explore the extent of genes with similar expression (if the cluster shows high expression in ovaries, enter ovary into the red box, select an experimental condition and move the slider to the right (and let go). ");

print h3("Your favourite expression pattern"),
  p("Let's say we're interested in genes which are highly expressed soon after blood feeding.   Type 'blood' into the red box and select the 'blood-fed 3h' condition.  Move the slider to +1 (2-fold upregulated).  You should see two regions of clusters (of genes) with high expression after blood feeding.  By clicking on the nodes you should be able to see that the smaller cluster predominantly contains cuticle proteins and also has an embryonic expression signature.  You could also highlight by high expression at 43h in embryos (enter '43h' in the green box...) and see how the two 'blood-fed 3h' regions differ in this second aspect of gene expression.  The map is telling us that these genes are probably involved in remodelling the exoskeleton to accommodate the blood meal.   Now, returning to the larger 'blood-fed 3h' region - we can set another expression highlighter for '>1' in 'Blood-fed adult female tissues (Marinotti et al., 2005) midgut [TGMA:0001036]' (this one is possibly easier to find via the drop-down menu).  We can now see that some of the 3h responding genes are in the midgut but not all.  We could go more fine-grained than this by setting a third highlighter for 'gastric caeca' or 'posterior midgut' etc.");

print h3("Your favourite protein domain or Gene Ontology term"),
  p("Just enter an Interpro domain ID, GO ID, or free text search (e.g. 'lectin', 'DNA repair') in one (or more) of the gene query boxes.   Clusters containing genes annotated with the domains and/or GO terms will be highlighted with a yellow, orange or pink sausage.  The widths of the sausages reflect the number of genes matching the query (mouse-over or click for more details).  You may find that genes with the function of interest are located in one or more regions of the map.  These regions can be explored via the node-popups and expression highlighting sliders as above.  You may enter multiple GO and domain IDs, but only one text query into each box.   With regard to GO terms, genes annotated with child terms of the query are highlighted.");

print h3("A gene list from your own experiment or analysis"),
  p("Simply paste them into one of the gene query boxes.  They will be highlighted on the map.  Up to three lists can be visualised.");

print h3("Genes from a particular pathway"),
  p("As above, you can paste in a list of genes annotated as belonging to a pathway.  One useful source is ".a({-target=>'_wikipathways', -href=>"http://www.wikipathways.org/index.php?title=Special%3ABrowsePathwaysPage&browse=Anopheles+gambiae&browseCat=All+Categories"}, 'WikiPathways').", from which you can conveniently copy-paste the entire contents of the 'External references' table into a gene query box.");

print h3("Visualise or analyse the data another way"),
  p("You can download the processed, normalised data used to generate this map from ".a({-href=>"http://www.vectorbase.org/GetData/Downloads/?&type=Expression&archive_status=current"}, 'the VectorBase download area').".  Look for the files named 'shifted-means'.  Be aware that missing values are represented by empty strings (tab-delimited).");

print <<"EOJS";
<script type="text/javascript">
function finish_page() {
  em_init_autocompleters($map_id);
  em_init_sliders.defer($map_id);
  em_draw_nodes.defer($map_id, $map_width, $map_height, $maxsize);
  em_activate_nodes.defer($map_id, $map_width, $map_height);
}
window.onload = finish_page.defer();
</script>
EOJS


print end_html;

