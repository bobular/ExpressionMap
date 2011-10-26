package ExpressionMap::Bit;

use base 'ExpressionMap::DBI';

__PACKAGE__->table('bit');
__PACKAGE__->columns(Primary => qw/id/);
__PACKAGE__->columns(Essential => qw/state/);
__PACKAGE__->columns(Others => qw/mapping vecdim/);

1;
