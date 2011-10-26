package ExpressionMap::Weight;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('weight');
__PACKAGE__->columns(All => qw/id mapnode vecdim weight/);

__PACKAGE__->has_a(vecdim => 'ExpressionMap::VectorDimension');


1;
