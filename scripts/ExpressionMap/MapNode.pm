package ExpressionMap::MapNode;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('mapnode');
__PACKAGE__->columns(All => qw/id map x y/);
__PACKAGE__->has_many(weights => 'ExpressionMap::Weight');
__PACKAGE__->has_many(mappings => 'ExpressionMap::Mapping');
__PACKAGE__->has_many(go_ora_results => 'ExpressionMap::GoOraResult');
1;
