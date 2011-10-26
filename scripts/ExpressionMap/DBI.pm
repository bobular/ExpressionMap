package ExpressionMap::DBI;

use base 'Class::DBI';
use Config::General;

my $conf = new Config::General("config.txt"); # same config file as index.cgi
my %config = $conf->getall;

__PACKAGE__->connection("dbi:mysql:database=$config{dbname};host=$config{dbhost};mysql_read_default_file=/etc/my.cnf", $config{dbuser}, $config{dbpass});

1;
