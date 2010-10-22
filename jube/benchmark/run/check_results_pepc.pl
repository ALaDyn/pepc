#!/usr/local/bin/perl -w

use strict;
use Carp;

if((scalar @ARGV) != 3) {
    printf(STDERR "incorrect number of parameter (%d) of $0 (3 required)\n",scalar @ARGV);
    exit(-1);
}

my $xmloutfile = $ARGV[0];

my $vcheck=0;
my $vcomment="not tested";

my $limit=1e-6;

my $line1; 
my $line2;

my $L2norm = 0;

open(file1, "$ARGV[1]") or die "can not open file $ARGV[1]";
open(file2, "$ARGV[2]") or die "can not open file $ARGV[2]";

print "compare trajectories in $ARGV[1] and $ARGV[2]\n";

readline(file1);
readline(file2);

my $lines_compared = 0;
while( defined ($line1 = readline(file1)) && defined ($line2 = readline(file2)) )
{
    print "line1: $line1";
    print "line2: $line2";

    my @values1 = split(/\s+/, $line1);
    my @values2 = split(/\s+/, $line2);

    my $diffx = $values1[2] - $values2[2];
    my $diffy = $values1[3] - $values2[3];
    my $diffz = $values1[4] - $values2[4];

    print "linediff: $diffx $diffy $diffz\n";
    print "\n";

    $L2norm = $L2norm + $diffx*$diffx + $diffy*$diffy + $diffz*$diffz;
    
    $lines_compared += 1;
}

if( $lines_compared != 4)
{
    print "not enough lines compared";
    $L2norm = 10*$limit*$limit;
}

$L2norm = sqrt($L2norm);

print "L2norm difference $L2norm\n";

close(file1);
close(file2);

if( $L2norm < $limit )
{
    $vcomment = "verify ok";
    $vcheck   = 1;
}
else
{
    $vcomment = "verify ERROR";
    $vcheck   = 0;
}   


open(XMLOUT,"> $xmloutfile") || die "cannot open file $xmloutfile";
print XMLOUT "<verify>\n";
print XMLOUT " <parm name=\"vcheck\" value=\"$vcheck\" type=\"bool\" unit=\"\" />\n";
print XMLOUT " <parm name=\"vcomment\" value=\"$vcomment\" type=\"string\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"L2norm\" value=\"$L2norm\" type=\"float\" unit=\"\"/>\n";
print XMLOUT "</verify>\n";
print XMLOUT "\n";
close(XMLOUT);


exit(0);
