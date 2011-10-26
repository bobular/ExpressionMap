#!/usr/bin/perl -w
#                           -*- mode: cperl -*-
use CGI qw(:standard);
use JSON;
use lib './';
use ExpressionMap::Map;
use GeneOntology::Term;
# use feature 'switch'; # no decent perl on prod server!

#
# get the map object
#

print header('application/json');

my $map_id = param('id');
unless (defined $map_id && $map_id =~ /^\d+$/) {
	print to_json({ error => 'malformed id parameter' });
	exit;
}
my $map = ExpressionMap::Map->retrieve(id=>$map_id);
my ($map_height, $map_width) = ($map->ysize(), $map->xsize());

my $action = param('action');

#given ($action) {

	if ($action eq 'nodeinfo') {
		my $biomart = read_biomart_txt(); # WARNING: takes 1s!
		my ($x) = param('x') =~ /(\d+)/;
		my ($y) = param('y') =~ /(\d+)/;
		my ($fold_thresh) = param('ft') =~ /([0-9.]+)/;
		my ($num_conds) = param('nc') =~ /(\d+)/;
		if (defined $x && defined $y && $x >= 0 && $y >= 0 &&
				$x < $map_width && $y < $map_height) {
			my $mapnode = $map->mapnodes('x'=>$x, 'y'=>$y)->first;
			# warn $mapnode->weights;

			my @weights;
			if (defined $fold_thresh) {
				# conditions above a threshold
				my $log2_thresh = log($fold_thresh)/log(2);
				@weights = grep { abs($_->weight)>= $log2_thresh } $mapnode->weights;
			} elsif (defined $num_conds) {
				# top N conditions
				@weights = sort { abs($b->weight) <=> abs($a->weight) } $mapnode->weights;

				splice @weights, $num_conds;
			}

			print to_json({ # x => $x, y => $y, map_id => $map_id,
											weights => {
																	map { $_->vecdim->name, $_->weight } @weights
																 },
											genes => [
																sort map { symbolise($_->gene_id, $biomart->{gene2symbol}) } $mapnode->mappings
															 ],
										  go_ora_results => [
																				 map {
																					 my $type = $_->term->term_type;
																					 $type =~ s/_/ /;
																					 { term => $_->term->name,
																						 acc => $_->term->acc,
																						 type => $type,
																						 score => $_->score,
																						 pvalue => $_->corrected_pvalue,
																					 } } grep { $_->score > 1 && $_->corrected_pvalue < 0.05 } $mapnode->go_ora_results
																				],
										});

		} else {
			print to_json({ error => 'malformed x and y params' });
		}
	}	elsif ($action eq 'search') {
		my $biomart = read_biomart_txt();
		my $max_node_hits = 0;

		my @results;
		#
		# tokenise query (only GO:0012345 and IPR012123 and AGAP123123 and CLIP8A style allowed)
		#
		my @tokens = param('q') =~ /([A-Z0-9:]{3,})/g;
		# make lookup and uniqify
		my %tokens = map { $_ => 1 } @tokens;

		# maybe it was a text query (remember: the DNA of DNA synthase would have been tokenised above)
		my ($f_query) = param('q') =~ /\b([\w\s-]+)\b/;
		if ($f_query) {
			my $perlregexp = qr/\b$f_query\b/i;
			my $mysqlregexp = '[[:<:]]'.$f_query.'[[:>:]]';

			foreach my $iprid (keys %{$biomart->{interpro_descriptions}}) {
				if ($biomart->{interpro_descriptions}{$iprid} =~ /$perlregexp/) {
					$tokens{$iprid} = 1;
				}
			}
			map { $tokens{$_->acc}=1 } GeneOntology::Term->search_name_regexp($mysqlregexp);
		}

		#
		# expand to child terms of query GO terms
		#
		my @children;
		foreach my $token (keys %tokens) {
			if ($token =~ /^GO:\d+$/) {
				foreach my $child (GeneOntology::Term->search_descendants($token)) {
					push @children, $child->acc;
				}
			}
		}
		# add children to %tokens
		map { $tokens{$_} = 1 } @children;

		if (keys %tokens) {
			# turn the go/interpro/genesymbol tokens into gene ids.
			my @temp = keys %tokens;
			map { map { $tokens{$_}=1; } keys %{$biomart->{go}{$_}} } @temp;
			map { map { $tokens{$_}=1; } keys %{$biomart->{interpro}{$_}} } @temp;
			map { map { $tokens{$_}=1; } keys %{$biomart->{symbol2gene}{$_}} } @temp;

			my @mapnodes = $map->mapnodes;
			foreach my $mapnode (@mapnodes) {
				my $x = $mapnode->x;
				my $y = $mapnode->y;
				my %hits; # gene_id => 1
				foreach my $mapping ($mapnode->mappings) {
					my $gene_id = $mapping->gene_id;
					if ($tokens{$gene_id}) {
						$hits{symbolise($gene_id, $biomart->{gene2symbol})} = 1;
					}
				}
				my @hits = sort keys %hits;
				push @results, { x => $x, y => $y, hits => \@hits } if (@hits);
				$max_node_hits = @hits if (@hits > $max_node_hits);
			}
		}
		print to_json({ search_results => [ @results ],
										max_node_hits => $max_node_hits,
									});
	}	elsif ($action eq 'ehighlight') {
		my ($log2_thresh) = param('lt') =~ /(-?[0-9.]+)/;
		my ($vecdim_idx) = param('v') =~ /(\d+)/;
		if (defined $log2_thresh && defined $vecdim_idx) {
			my $vecdim = $map->vecdims(idx=>$vecdim_idx)->first;
			my $vecdim_id = $vecdim->id;
			my @mapnodes = $map->mapnodes;

			my @highlighted;
			if ($log2_thresh > 0) {  # UGLY CUT AND PASTE!!
				@highlighted = grep { $_->weights('vecdim'=>$vecdim_id)->first->weight > $log2_thresh } @mapnodes;
			} else {
				@highlighted = grep { $_->weights('vecdim'=>$vecdim_id)->first->weight < $log2_thresh } @mapnodes;
			}
			print to_json({ nodes => [ map { [ $_->x, $_->y ] } @highlighted ] });
		} else {
			print to_json({ error => 'malformed v and lt params' });
		}
	} elsif ($action eq 'conditions') {
    my $json = to_json({ conditions => { map { $_->idx, $_->name } $map->vecdims() } });
    $json =~ s/sailvary/salivary/g; ### bugfix alert ###
    print $json;
	}	else {
		print to_json({ error => 'unrecognised action' });
	}
#}; # given



# turn a gene_id into "gene_id (SYMBOL1/SYMBOL2/...)"
sub symbolise {
	my ($gene_id, $gene2symbol) = @_;
	my $retval = $gene_id;
	if ($gene2symbol->{$gene_id}) {
		$retval .= ' ('.join('/', sort keys %{$gene2symbol->{$gene_id}}).')';
	}
	return $retval;
}

sub read_biomart_txt {
	open MART, "biomart.txt";
	my $headers = <MART>;
	chomp($headers);
	my @headers = split /\t/, $headers;
	my %header_idx;
	for (my $i=0; $i<@headers; $i++) {
		$header_idx{$headers[$i]} = $i;
	}
	my %goterms; # go id => gene id => 1
	my %interpro; # interpro id => gene id => 1
	my %iprdesc;
	my %symbol2gene; # symbol => gene id => 1
	my %gene2symbol; # gene id => symbol => 1

	my $hi_iprid = $header_idx{"Interpro ID"};
	my $hi_iprdesc = $header_idx{"Interpro Description"};
	my $hi_goid = $header_idx{"GO ID"};
	my $hi_godesc = $header_idx{"GO description"};
	my $hi_symbol = $header_idx{'Anopheles symbol'};
	my $hi_geneid = $header_idx{"Ensembl Gene ID"};
	while (<MART>) {
		chomp;
		my @cols = split /\t/, $_;
		my $id = $cols[$hi_geneid];
		die "no 'Ensembl Gene ID' column in biomart data\n" unless (defined $id);

		my $goacc = $cols[$hi_goid];
		$goterms{$goacc}{$id} = 1 if ($goacc);

		my $ipracc = $cols[$hi_iprid];
		$interpro{$ipracc}{$id} = 1 if ($ipracc);

		my $symbol = $cols[$hi_symbol];
		$symbol2gene{$symbol}{$id} = 1 if ($symbol);
		$gene2symbol{$id}{$symbol} = 1 if ($symbol);

		my $iprdesc = $cols[$hi_iprdesc];
		if ($ipracc && $iprdesc) {
			$iprdesc =~ s/_/ /g; # bit of a hack to allow IPR short descrips to be searchable
			$iprdesc{$ipracc} = $iprdesc;
		}
	}
	close(MART);

	return { go => \%goterms,
					 interpro => \%interpro,
					 symbol2gene => \%symbol2gene,
					 gene2symbol => \%gene2symbol,
					 interpro_descriptions => \%iprdesc,
					 };
}
