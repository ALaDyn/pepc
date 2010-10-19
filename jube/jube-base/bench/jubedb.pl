#!/usr/bin/perl -w

use strict;

my $patint="([\\+\\-\\d]+)";   # Pattern for Integer number
my $patfp ="([\\+\\-\\d.E]+)"; # Pattern for Floating Point number
my $patwrd="([\^\\s]+)";       # Pattern for Work (all noblank characters)
my $patnint="[\\+\\-\\d]+";    # Pattern for Integer number, no () 
my $patnfp ="[\\+\\-\\d.E]+";  # Pattern for Floating Point number, no () 
my $patnwrd="[\^\\s]+";        # Pattern for Work (all noblank characters), no () 
my $patbl ="\\s+";             # Pattern for blank space (variable length)

use Getopt::Long;
use File::stat qw(:FIELDS);
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use Time::HiRes qw ( time );
use Tk;
use Tk::Tree;

use jube_data_dbhash;
use jube_data_mysql;

my $pwd=`pwd`;
chomp($pwd);

my $opt_verbose=0;
my $opt_dump=0;
my $opt_gui=0;
my $opt_debug=undef;
my $opt_scan=undef;
my $opt_dbtype="mysql";

usage($0) if( ! GetOptions( 
			    'verbose:i'        => \$opt_verbose,
			    'dump'             => \$opt_dump,
			    'gui'              => \$opt_gui,
			    'scan=s'           => \$opt_scan,
			    'dbtype=s'         => \$opt_dbtype
			    ) );

 


my $startdate=localtime(time());

my ($benchxmlfile,$benchlogfile,$startspec,$enddate);

my $dbobj;


$dbobj=bench_data_dbhash->new() if($opt_dbtype eq "dbhash");
$dbobj=bench_data_mysql->new("localhost","nbench","nbench","")  if($opt_dbtype eq "mysql");

$dbobj->load_db();

if($opt_scan) {
    $dbobj->scan_for_new_benchmarks($opt_scan);
}

exit;


$dbobj->analyse();

my ($name);
foreach $name (keys %{$dbobj->{B_NAME}}) {
    print "WF Benchmark: $name\n";
}

my ($platform);
foreach $platform (keys %{$dbobj->{B_PLATFORM}}) {
    print "WF Benchmark: $platform\n";
}

if($opt_gui) {
    &gui();
}

$dbobj->save_db();

sub usage {
    die "Usage: $_[0] <options> <xml-file>
                -verbose           : verbose
                -dump              : dump XML-file structure
                -debug             : don't submit jobs
                -rmtmp             : remove temp directory directly

";
}

sub gui {
    my $top = new MainWindow( -title => "DirTree" );

    my $current_bench;
    my $tree = $top->Scrolled( qw/Tree
                               -width 35 -height 30
                               -selectmode browse -exportselection 1 
                               -scrollbars osoe/ );

    my $lab = $top->Label( -text => "Benchmark:" );
    my $ent = $top->Entry( -textvariable => \$current_bench );
 
    my $ok = $top->Button( qw/-text Ok -underline 0 -width 6/ );
    my $cancel = $top->Button( qw/-text Cancel -underline 0 -width 6/ );

    $tree->configure(   -browsecmd => sub { $current_bench = shift } );
    $tree->configure(   -command   => sub { do_it( $current_bench ) } );
    $ok->configure(     -command   => sub { do_it( $current_bench ) } );
    $cancel->configure( -command   => sub { exit } );
 
    $tree->pack(   qw/-expand yes -fill both -padx 10 -pady 10 -side top/ );
    $lab->pack(    qw/-anchor w/ );
    $ent->pack(    qw/-fill x/ );
    $ok->pack(     qw/-side left  -padx 10 -pady 10/ );
    $cancel->pack( qw/-side right -padx 10 -pady 10/ );

    MainLoop();
}









