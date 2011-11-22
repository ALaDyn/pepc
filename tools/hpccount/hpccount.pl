#!/usr/bin/env perl

#*==========================================================*
#*                                                          *
#*  Wrapper Program for BGP_hpccount                        *
#*                                                          *
#*  This program sets environment variables and then calls  *
#*  the application                                         *
#*                                                          *
#*  Initial version: Oct 1, 2010                            *
#*  C. Pospiech ** ACTC  ** IBM Germany                     *
#*  © IBM 2010                                              *
#*                                                          *            
#*==========================================================*

sub print_helptext {
    my($i);
    my(@helptext) = (
    "\n", 
    " usage:\n",
    "   hpccount [-o <name>] [-u] [-n] [-x] [-g <group>] "
                ."[-a <plugin>] <program>\n",
    "   OR           \n",
    "   hpccount [-h]\n",
    " where: <program>                 program to be executed\n",
    "        -h, --help                displays this help message\n",
    "        -o <name>, --output       output file name\n",
    "        -u, --unique              make the file name <name> unique\n",
    "        -n, --no-stdout           no hpccount output to stdout\n",
    "        -x, --xtra-print          adds formulas for derived metrics\n",
    "        -g <group>, --group       BlueGene UPC group number(s)\n",
    "        -a <plugin>, --aggregate  aggregate counters using the"
                                           ." plugin <plugin>\n",
    "\n", 
    "        options can be specified in either of the following three ways.\n",
    "         -o <value>\n",
    "        --option <value>\n",
    "        --option=<value>\n",
    "\n"
		     );
    open(HELPFILE,"| more");
    for $i (@helptext) {
	print HELPFILE $i;
    }
    close(HELPFILE);
    exit;
}
#-------------------------------------------------------------------------
sub get_options() {
    use Getopt::Long;

    GetOptions ('output=s'    => \$out_fname,
		'unique'      => \$unique_fname,
                'no-stdout'   => \$no_stdout,
                'xtra-print'  => \$print_formula,
                'group=s'     => \$event_set,
                'aggregate=s' => \$aggregate,
                'export=s'    => \$export_str, 
		'help'        => \$help ) || print_helptext();

    # The option --export=opt_string is only used for MPI implementations
    # different from BlueGene MPI. It refers to the name of the option
    # to mpirun or mpiexec that is used to propagate environment variables 
    # to the MPI tasks.
    #
    # The default value is "env" because this is the value for BlueGene.

    print_helptext() if ( $help );
}
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#  main program
#-------------------------------------------------------------------------

# Step 1 - default values

$export_str = "env";

# Step 2 - get command line options

&get_options();

# Command line options take precedence over environment variables.
# So the latter are changed to meet the command line options.

$ENV{"HPM_OUTPUT_NAME"}      = $out_fname if ($out_fname);
$ENV{"HPM_UNIQUE_FILE_NAME"} = "yes"      if ($unique_fname);
$ENV{"HPM_STDOUT"}           = "no"       if ($no_stdout);
$ENV{"HPM_PRINT_FORMULA"}    = "yes"      if ($print_formula);
$ENV{"HPM_EVENT_SET"}        = $event_set if ($event_set);
$ENV{"HPM_AGGREGATE"}        = $aggregate if ($aggregate);

# Step 3 - find the right values for LD_LIBRARY_PATH and LD_PRELOAD

$pwd = `pwd`;
chomp($pwd);
$ld_preload = `which $0`;               # Start with the path to hpccount
chomp($ld_preload);
$ld_preload =~ s/[^\/]+$/lib/;          # replace "hpccount[.pl]" with "lib"
$ld_preload = $pwd."/".$ld_preload      # convert to absolute path name 
    if ( $ld_preload !~ m/^\// );       # if appropriate
while ( $ld_preload =~ m/\/\.+\// ) {
    $ld_preload =~ s/\/\.\//\//;        # clean up; change /./ -> /
    $ld_preload =~ s/\/[^\/]+\/\.\.\//\//;        # change /something/../ -> /
}

if ( $ENV{"LD_LIBRARY_PATH"} ) {
    $ld_library_path = $ld_preload.":".$ENV{"LD_LIBRARY_PATH"};
} else {
    $ld_library_path = $ld_preload;
}
$ld_library_path = $ENV{"IHPCT_BASE"}."/lib:$ld_library_path"
    if ( $ENV{"IHPCT_BASE"} );
$ld_preload .= "/libhpccnt.so";

# Step 4 - Stack all relevant environment variables into the command line

$command = shift @ARGV;

unshift @ARGV, "-$export_str LD_PRELOAD=$ld_preload";
unshift @ARGV, "-$export_str LD_LIBRARY_PATH=$ld_library_path";

for $key ( keys %ENV ) {
    if ( $key =~ m/HPM_/ ) {
	$val = $ENV{$key};
	unshift @ARGV, "-$export_str $key=$val";
    }
}

# fit the command back into the command line
unshift @ARGV, $command;

# Step 5 - call the command

exec(join(' ',@ARGV)) 
    or die "HPM ERROR - calling the command $command failed\n"
    ."OS error message: $!\n";

