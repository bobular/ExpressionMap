package ExpressionMap::VectorDimension;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('vecdim');
__PACKAGE__->columns(All => qw/id map idx name count_pos/);

1;
