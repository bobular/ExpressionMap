package ExpressionMap::GoOraResult;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('go_ora_result');
__PACKAGE__->columns(Primary => qw/id/);
__PACKAGE__->columns(Essential => qw/term corrected_pvalue score/);
__PACKAGE__->columns(Others => qw/mapnode pvalue numgenes serial_number/);
__PACKAGE__->has_a(term => 'GeneOntology::Term');

#
# need to deal with"same as" redundant term
#

1;
