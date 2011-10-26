package ExpressionMap::Map;

use CGI;
use base 'ExpressionMap::DBI';

__PACKAGE__->table('map');

__PACKAGE__->columns(All =>
    qw/id name xsize ysize bootstrap commandline/);

__PACKAGE__->has_many(vecdims => "ExpressionMap::VectorDimension",
                     { order_by => 'idx' });
__PACKAGE__->has_many(mapnodes => "ExpressionMap::MapNode");
__PACKAGE__->has_many(mappings => "ExpressionMap::Mapping");
__PACKAGE__->has_many(mappings_ids => [ "ExpressionMap::Mapping", 'gene_id' ]);

my $dbh = __PACKAGE__->db_Main();
$pair_support = $dbh->prepare(q{
select count(*) as support from mapping m left join map on (m.map=map.id) left join mapnode mn on (m.mapnode=mn.id), mapping n left join mapnode nn on (n.mapnode=nn.id) where m.map=n.map and (mn.x-nn.x)*(mn.x-nn.x)+(mn.y-nn.y)*(mn.y-nn.y) <= ? and bootstrap>0 and m.gene_id=? and n.gene_id=? and map.name=? group by m.gene_id, n.gene_id
});

$gene_neighbours = $dbh->prepare(q{
select m.gene_id, count(*) as support from mapping m left join map on (m.map=map.id) left join mapnode mn on (m.mapnode=mn.id), mapping n left join mapnode nn on (n.mapnode=nn.id) where m.map=n.map and (mn.x-nn.x)*(mn.x-nn.x)+(mn.y-nn.y)*(mn.y-nn.y) <= ? and bootstrap>0 and n.gene_id=? and map.name = ? group by m.gene_id, n.gene_id order by support desc;
});

# old without euclidean distance
# select count(*) as support from mapping m left join map on (m.map=map.id), mapping n where m.map=n.map and m.mapnode=n.mapnode and bootstrap>0 and m.gene_id=? and n.gene_id=? and map.name=? group by m.gene_id, n.gene_id


sub bootstrap_pair {
  my ($self, $gene1, $gene2, $radius) = @_;

  $radius = 0 unless (defined $radius);

  my $r2 = $radius*$radius;
  # check that the genes are on the map
  my $mapping1 = $self->mappings(map=>$self, gene_id=>$gene1)->first;
  my $mapping2 = $self->mappings(map=>$self, gene_id=>$gene2)->first;

  if (defined $mapping1 && defined $mapping2) {
    $pair_support->execute($r2, $gene1, $gene2, $self->name);
    my ($support) = $pair_support->fetchrow_array();
    $pair_support->execute($r2, $gene1, $gene1, $self->name);
    my ($total) = $pair_support->fetchrow_array();
    return sprintf("%.0f%%", (100*$support/$total || 0));
  }
  return "?";
}


sub bootstrap_neighbours {
  my ($self, $gene, $radius) = @_;

  $radius = 0 unless (defined $radius);

  my $r2 = $radius*$radius;
  # check that the genes are on the map
  my $mapping = $self->mappings(map=>$self, gene_id=>$gene)->first;

  if (defined $mapping) {
    $gene_neighbours->execute($r2, $gene, $self->name);

    my $array = $gene_neighbours->fetchall_arrayref();
    my @rows;
    push @rows, CGI::Tr(CGI::td("Neighbour gene"),
		       CGI::td("Bootstrap score"));
    my $max;
    foreach my $row (@$array) {
      my ($g, $s) = @$row;
      $max = $s if ($g eq $gene);
    }
    return 'nomax error!' unless ($max);
    foreach my $row (@$array) {
      my ($g, $s) = @$row;
      next if ($g eq $gene);
      my $pc = sprintf "%.0f%%", 100*$s/$max;
      push @rows, CGI::Tr(CGI::td({}, [$g, $pc]));
    }
    return CGI::table(@rows);
  }
  return "?";
}
1;
