my $file = open(FILE,"$ARGV[2]");
my $line, $tmpline, $time, $proccnt;
my @splitted;

while (defined ($line = <FILE>)) {
    $tmpline = $line;
    chomp $tmpline;
    my $size = length $tmpline;
    if( $line =~ '  %   cumulative   self              self     total')
    {
#	print $line;
	$line = <FILE>;
#	print $line;
	$time = 1;
	$proccnt = 0;
	while($time > 0.1)
	{
	    $proccnt++;
	    $line = <FILE>;
#	    print "$line";
	    $tmpline = $line;
	    chomp $tmpline;
	    $size = length $tmpline;
	    @splitted = split(/\s+/, $line);
	    $time = $splitted[1];
#	    foreach $x (@splitted) {print("$x\n");}
	    print("JuBE: gprof: proc $proccnt: $splitted[$ARGV[0]] $splitted[$ARGV[1]]\n");
#	    print("\t\t$time\n");
	}
    }
}
