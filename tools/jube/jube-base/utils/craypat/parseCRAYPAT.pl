use strict;
use warnings;

my $toParseFileNumber = @ARGV - 1;

my $mode = $ARGV[0];

my $proccnt=0;

my $fileCnt;
for($fileCnt=1; $fileCnt < ($toParseFileNumber+1); $fileCnt++)
{
    open(FILE, "$ARGV[$fileCnt]") or die "can not open file $ARGV[$fileCnt]";
    
    my $line;
    my $outstage = 0; #0: start; 1: change into output; 2: output; 3: change out; 4: end
    while( defined ($line = readline(FILE)) )
    {
	if( ($mode eq "HWC") || ($mode eq "HEAP"))
	{
	    if($outstage == 1) {$outstage = 2};
	    if($outstage == 3) {$outstage = 4};
	    if($line =~ /Table 1/ && $outstage == 0) {$outstage = 1;}
	    if($line =~ /Additional details/ && $outstage == 2) {$outstage = 3;}

	    if($outstage == 2) {print "JuBE: CRAYPAT: $mode: $line";}
	}
	elsif( ($mode eq "TIME") || ($mode eq "FLOPS"))
	{
	    if($outstage == 1) {$outstage = 2};
	    if($outstage == 3) {$outstage = 4};
	    if($line =~ /\|USER/ && $outstage == 0) {$outstage = 1;}
	    if($line =~ /\|\|---/ && ($outstage == 0 || $outstage == 2) ) {$outstage = 1;}
	    if($line =~ /\|\|===/ && $outstage == 2) {$outstage = 3;}

	    if($outstage == 2) {$proccnt++; print "JuBE: CRAYPAT: $mode: proc $proccnt: $line";}
	}

    }

    close(FILE);
}
    
