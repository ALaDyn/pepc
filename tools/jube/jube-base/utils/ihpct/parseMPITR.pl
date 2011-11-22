use strict;
use warnings;

my $toParseFileNumber = @ARGV;

my $fileCnt;
for($fileCnt=0; $fileCnt < $toParseFileNumber; $fileCnt++)
{
    open(FILE, "$ARGV[$fileCnt]") or die "can not open file $ARGV[$fileCnt]";
    
    my $line;
    while( defined ($line = readline(FILE)) )
    {
	print "IHPCT: MPITracer: $line";
    }

    close(FILE);
}
    
