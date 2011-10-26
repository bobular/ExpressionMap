package GeneOntology::Term;

use base 'GeneOntology::DBI';

__PACKAGE__->table('term');
__PACKAGE__->columns(All => qw/acc name term_type/);

__PACKAGE__->set_sql(descendants => qq{
SELECT DISTINCT descendant.id, descendant.name, descendant.term_type, descendant.acc
FROM term
INNER JOIN graph_path ON (term.id=graph_path.term1_id)
INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id)
WHERE term.acc=? AND distance > 0
});

__PACKAGE__->set_sql(ancestors => qq{
SELECT ancestor.id, ancestor.name, ancestor.term_type, ancestor.acc
FROM term child, graph_path, term ancestor
WHERE child.id=graph_path.term2_id
AND ancestor.id=graph_path.term1_id
AND distance > 0
AND ancestor.is_root = 0
AND child.acc=?
});

__PACKAGE__->set_sql(name_regexp => qq{SELECT __ESSENTIAL__ FROM __TABLE__ WHERE name REGEXP ?});

1;
