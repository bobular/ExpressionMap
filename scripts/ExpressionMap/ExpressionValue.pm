package ExpressionMap::ExpressionValue;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('expression');
__PACKAGE__->columns(Primary => qw/id/);
__PACKAGE__->columns(Essential => qw/value/);
__PACKAGE__->columns(Others => qw/mapping vecdim/);

1;
