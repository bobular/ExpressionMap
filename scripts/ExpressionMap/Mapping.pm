package ExpressionMap::Mapping;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('mapping');
__PACKAGE__->columns(Primary => qw/id/);
__PACKAGE__->columns(Essential => qw/gene_id/);
__PACKAGE__->columns(Others => qw/map mapnode error/);
__PACKAGE__->has_many(bits => 'ExpressionMap::Bit');
__PACKAGE__->has_many(expression_values => 'ExpressionMap::ExpressionValue');

1;
