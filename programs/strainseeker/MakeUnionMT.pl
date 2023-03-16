use strict;
use warnings;
use Getopt::Long;
my $glc = "./glistcompare"; # Glistcompare
my $glq = "./glistquery"; # Glistquery
if (!-e $glq) { die "GListquery not found!\n" }
if (!-e $glc) { die "GListcompare not found!\n" }

my $threads = 32; # Number of parallel processes
my $help;
GetOptions(
    't=i' => \$threads,
    'threads=i' => \$threads,
    'h' => \$help,
    'help' => \$help,

     ) or die printHelp()."\n";
if ($help || !@ARGV){
	die printHelp()."\n";
}
my $word;
my $i = 0;
foreach (@ARGV){
	if (system "$glq $_ -stat >/dev/null 2>&1"){
		next;
	} else {
		$word = (split ("\t", (qx/$glq $_ -stat/)[1]))[1];
		last;
	}
}
chomp $word;
my @nodes = @ARGV;
my %forCmd;
my $cmd;
my $rm = "rm ";
my $left = (scalar @nodes)-1;
my $unioncounter = 0;
my $last_union;

if (scalar @nodes == 1){
	die "Atleast 2 valid files must be given!\n";
} elsif(scalar @nodes == 0){
	die printHelp()."\n";
}

while (@nodes){
	if (scalar @nodes > 1){
		my $unionname = "MULTITHREADED_union_".$unioncounter;
		my $node1;
		my $node2;
		while ((!defined $node1) || (system "$glq $node1 -stat >/dev/null 2>&1")){
			$node1 = shift @nodes;
		}
		while ((!defined $node2) || (system "$glq $node2 -stat >/dev/null 2>&1")){
			$node2 = shift @nodes;
		}

		if(!defined $node1 || !defined $node2){
				die "Atleast 2 valid files must be given!\n";
		}
		@{$forCmd{$unionname}} = ($node1, $node2);
		$unioncounter++;
	}

	if ((scalar @nodes <= 1) || (scalar keys %forCmd >= $threads)){
		print "Making union of lists, remaining: $left\n";
		foreach (keys %forCmd){
			$cmd .= "$glc $forCmd{$_}[0] $forCmd{$_}[1] -u -o $_ & ";
			if ("$forCmd{$_}[0]" =~ /MULTITHREADED_union_/) {
				$rm .= "$forCmd{$_}[0] ";
			}
			if ("$forCmd{$_}[1]" =~ /MULTITHREADED_union_/) {
				$rm .= "$forCmd{$_}[1] ";
			}
			push (@nodes, $_."_$word\_union.list");
			$last_union = $_."_$word\_union.list";
		}

		if(system "$cmd wait"){
			die "Problem with main command: $cmd";
		}

		if(($rm ne "rm ") && (system "$rm")){
			die "Problem with main command: $rm";
		}

		$left-= scalar keys %forCmd;
		undef %forCmd;
		$cmd = "";
		$rm = "rm ";
		if (scalar @nodes == 1){
			last;
		}
	} elsif ((scalar @nodes == 1) && (scalar keys %forCmd == 1)){
		die "Atleast 2 valid files must be given!\n";
	}
}

if (system "mv $last_union union\_$word\_union.list"){
	die "Problem with command: mv $last_union\_$word\_union.list union\_$word\_union.list";
}
if (glob("MULTITHREADED_union_*") && (system "rm MULTITHREADED_union_*")){
	warn "Problem with removing: rm MULTITHREADED_union_*";
}

sub printHelp {
	print "Usage: $0 <file1.list> <file2.list> <fileN.list>\n";
	print "Options:
	-h, --help\t - Print this help
	-t, --threads\t - Number of cores used (default $threads)\n";
	return "";
}
print "Created file union\_$word\_union.list\n";