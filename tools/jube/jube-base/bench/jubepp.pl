#!/usr/bin/perl -w
#
#####################################################################################
#                                                                                   #
#     JuBE Pretty Print: Tool for creating benchmark reports from JuBE output       #
#                                                                                   #
#####################################################################################
#
#
#
use FindBin;
use lib "$FindBin::RealBin/lib";
#use Find::File;
use File::Basename;
 
use strict;
use Carp;

use XML::Simple;
use XML::Writer;
use IO::File;

use Data::Dumper;
use Time::HiRes qw(time);
use Getopt::Long;

my $optLongOutput = undef;
my $optVersion = undef;
my $optAnalyse = undef;
my $optXMLCreate = undef;
my $optParamLong = undef;


usage($0) if(!GetOptions( 
			    'long'        => \$optLongOutput,
			    'Version|V'   => \$optVersion,
			    'analyse'     => \$optAnalyse,
			    'xmlcreate'   => \$optXMLCreate,
			    'param_long'  => \$optParamLong
			    ));

&usage($0) if (!$ARGV[0]);
&version() if ($optVersion);


# Patterns
my $compilerPattern = "\^\\s+(cc|cpp|cxx|f77|f90|f95|mpi_\\w+)";


# Open files and create appropriate name for output file
my $logfile = "xmllogs/" . shift;

my ($logfileName, $logfilePath, $logfileExtension) = fileparse($logfile);

my $xmlFileReport = undef;
my $XMLReport = undef;
my $xmlFileReportName = undef;


# XML file for writing the report in XML format; usable by XLST
if ($optXMLCreate) {
    $xmlFileReportName = $logfileName . ".xml";

# rename file if it already exists
    if (-f $xmlFileReportName) {
	my $cmd = "mv $xmlFileReportName $xmlFileReportName.$$";
	system($cmd);
    }
    $xmlFileReport = new IO::File("> $xmlFileReportName");
    $XMLReport = new XML::Writer(OUTPUT => $xmlFileReport,
				 NEWLINES => 1);
   }

my $report  = $logfileName . ".pp";
# rename file if it already exists
    if (-f $report) {
	my $cmd = "mv $report $report.$$";
	system($cmd);
    }

# detects pattern files
open(PATTERNFILES, '<', "pattern_files.dat") or die "Can't open pattern_files.dat: $!";
my $filesCounter = 0;
my @analyseFiles = undef;
while(<PATTERNFILES>){
    next if ($_ =~ /^\s/);
    chomp($_);
    $analyseFiles[$filesCounter] = $_;
    ++$filesCounter;
}

print "which files had been detected \n";
foreach (@analyseFiles) {
print "$_\n"
}

my @compileFlags = undef;
my @compileDirectories = undef;
my @compilerDefinitions = undef;
my @dirAndLibs = undef;
my @analysePatterns = undef;
my @analysePatternsIndex = undef;
my @params = undef;

&collector();

my $input = XML::Simple->new();
my $analyseXML = XML::Simple->new();

# read in log file
my $startTime = time;
my $benchrun = $input->XMLin($logfile, KeyAttr => {benchmark => "+platform"
			     },
			     ForceArray => 1,
			     ForceContent => 1,
			     SuppressEmpty => 1,
			     ContentKey => '-content');

my $compileCommandRef = $benchrun->{compile};
my $compileParamsRef = $compileCommandRef->[0]->{params}->[0];
my $analyseRef = $benchrun->{analyse};
my $stdoutRef = $benchrun->{stdoutfile};


my $timeDiff = time - $startTime;



&headline_print();
&benchmark_print();
&compile_print();
&params_print()if ($optParamLong);
&stdoutput_print() if ($optLongOutput); 
&analyse_print() if ($optAnalyse);
&footnote_print();
&result_print();



sub headline_print{

    open (REPORT, '>>', $report) or die $!;   
    print_separator();
    printf REPORT ("Report performed on " . localtime() . "\n");
    printf REPORT ("%s\n", &getversion());
    print_separator();
    printf REPORT ("%-30s:%s\n", "Parsed Logfile" ,$logfile);
    printf REPORT ("%-30s:%s\n", "Report written to", $report);
    print_separator();
    printf REPORT ("%-30s:%s\n", "Calculation performed on", $compileParamsRef->{platform});
    printf REPORT ("%-30s:   %s\n","Executable", $compileParamsRef->{execname});
    print_separator();
    printf REPORT ("%-30s:%s\n", "Number of CPUs", $compileParamsRef->{ncpus});
    printf REPORT ("%-30s:%s\n", "Number of nodes", $compileParamsRef->{nodes});
    printf REPORT ("%-30s:%s\n", "Threads per task", $compileParamsRef->{threadspertask});
    printf REPORT ("%-30s:%s\n", "Tasks per node", $compileParamsRef->{taskspernode});
    print_separator();
    close (REPORT);

    
    if ($optXMLCreate) {
	my $performingTime = localtime();
	my $jubeppVersion = &getversion();
	$XMLReport->pi('xml-stylesheet', 'href="jube_report.xsl" type="text/xsl"');
	$XMLReport->startTag("report");
	$XMLReport->startTag('specifications');
	$XMLReport->dataElement('inputfile', $logfile);
	$XMLReport->dataElement('outputfile', $report);
	$XMLReport->dataElement('jubepp_version', $jubeppVersion);
	$XMLReport->dataElement('date', $performingTime);
	$XMLReport->dataElement('platform', $compileParamsRef->{platform});
	$XMLReport->dataElement('executable', $compileParamsRef->{execname} );
	$XMLReport->dataElement( 'ncpus', $compileParamsRef->{ncpus});
	$XMLReport->dataElement('nnodes', $compileParamsRef->{nodes});
	$XMLReport->dataElement('threads_per_task', $compileParamsRef->{threadspertask});
	$XMLReport->dataElement('tasks_per_node', $compileParamsRef->{taskspernode});
	$XMLReport->endTag('specifications');
    }
}


sub collector {

    my $flagsCounter = 0;
    my $compilerCounter = 0;
    my $dirLibCounter = 0;
    my $paramsCounter = 0;

    open (LOGFILE, '<', $logfile) or die "Can't open $logfile. $!";
    while(<LOGFILE>)
    {
	if ($_ =~ /^\s+(\w+flags)/) {
	    $compileFlags[$flagsCounter] = $1;
	    ++$flagsCounter;
	}
	elsif ($_ =~/$compilerPattern/) {
	    $compilerDefinitions[$compilerCounter] = $1;
	    ++$compilerCounter;
	}
	elsif (($_ =~/^\s+(\w+_dir)/) ||($_ =~/^\s+(\w+_lib)/) || ($_ =~/^\s+(\w+_inc)/)  ) {
	    $dirAndLibs[$dirLibCounter] = $1;
	    ++$dirLibCounter;
	}
	elsif ($_ =~ /^\s+(\w+)/) {
	    $params[$paramsCounter] = $1;
	    ++$paramsCounter;
	}
    }

    my $arrayCounter = 0;  
    my $indexCounter = 0; 
    foreach (@analyseFiles) {
	open(ANALYSE, '<', $_) or die "Can't open $_: $!";
#	print "analysefile: $_\n";
	while(<ANALYSE>) {
	    if ($_ =~ /parm name\s*=\s*\"(\w+)\".+mode\s*=\s*\"line,\s*index.+\"/) {
		    $analysePatternsIndex[$indexCounter] = $1;
		    ++$indexCounter;
		}
		elsif ($_ =~ /parm name\s*=\s*\"(\w+)\".+mode\s*=\s*\"(line|span|derived)/){
		    $analysePatterns[$arrayCounter] = $1;
		    ++$arrayCounter;
		}
	    }
	close (ANALYSE);
    }
}

sub benchmark_print{

    my $varmvx = $benchrun->{platform};
    open (REPORT, '>>', $report) or die $!; 

    close (REPORT);

}

sub compile_print {
   
    open (REPORT, '>>', $report) or die $!; 
    printf REPORT ("%-30s: %s\n", "Compile command", $compileCommandRef->[0]->{command}->[0]->{content});

    my $content;
    foreach $content (@compileFlags) {
	printf REPORT ("%-30s: %s\n",$content, $compileParamsRef->{$content});
    }
    print_separator();

    foreach $content (@compilerDefinitions) {
	printf REPORT ("%-30s: %s\n",$content, $compileParamsRef->{$content});
    }

    print_separator();
    printf REPORT ("Directories, Libraries and Includes:\n");
    foreach $content (@dirAndLibs) {
	printf REPORT ("%-30s: %s\n",$content, $compileParamsRef->{$content});
    }
    close (REPORT);
    if ($optXMLCreate) {
	$XMLReport->startTag('compile');
	foreach (@compilerDefinitions) {
	    $XMLReport->dataElement('compiler' ,  $compileParamsRef->{$_},
				    alias => $_);
	}

	foreach (@compileFlags) {
	    $XMLReport->dataElement('flags' ,  $compileParamsRef->{$_},
				    alias => $_);
	}

	my (@libArray, @dirArray, @incArray) = LibDirIncSort(@dirAndLibs);
	
	foreach (@libArray) {
	    $XMLReport->dataElement('library' , $compileParamsRef->{$_},
				    alias => $_);
	}
	
	foreach (@dirArray) {
	    $XMLReport->dataElement('directory' , $compileParamsRef->{$_},
				    alias => $_);
	}
	
	foreach (@incArray) {
	    $XMLReport->dataElement('include' , $compileParamsRef->{$_},
				    alias => $_);
	}
	
	$XMLReport->endTag();
    }

}

sub params_print {
    open (REPORT, '>>', $report) or die $!; 
    print REPORT "\n\n---------------------------- PARAMS SECTION ---------------------------\n";
    foreach (@params) {
	printf REPORT ("%-30s: %s\n",$_, $compileParamsRef->{$_}) if ($compileParamsRef->{$_});
    }
    print REPORT "\n-------------------------- END OF PARAMS SECTION --------------------------\n";
    close (REPORT);

    if($optXMLCreate) {
	$XMLReport->startTag('params');
	foreach (@params) {
	    if ($compileParamsRef->{$_}) {
		$XMLReport->dataElement('param', $compileParamsRef->{$_},
					variable => $_);	
	    }		 
	}
	$XMLReport->endTag();
    }
}

sub LibDirIncSort {
    my @dirAndLibsArray = @_;
    my (@libArray, @dirArray, @incArray);
    foreach (@dirAndLibsArray) {
	if ($_ =~ /\w+_lib/) {
	    push @libArray, $_;
	} elsif ($_ =~ /\w+_dir/) {
	    push @dirArray, $_
	} elsif ($_ =~ /\w+_inc/) {
	    push @incArray, $_;
	}
    }
    (@libArray, @dirArray, @incArray);
}

sub analyse_print() {
    open (REPORT, '>>', $report) or die $!; 
    print REPORT "\n\n---------------------------- ANALYSE SECTION ---------------------------\n";
    my $analyseValue = undef;
    my $analyseUnit = undef;
    my $analyseMultipleValues = undef;
    my $analyseMultipleNames = undef;
    my $content;
    foreach $content (@analysePatterns) {
	$analyseValue = $analyseRef->[0]->{$content}->[0]->{value};
	$analyseUnit = $analyseRef->[0]->{$content}->[0]->{unit};
	printf REPORT ("%-30s: %s %s\n",$content, $analyseValue, $analyseUnit) if ($analyseValue);
    }
    foreach $content (@analysePatternsIndex) {
	my $counter = $analyseRef->[0]->{$content}->[0]->{count};
	if(!$counter) {next;}
	printf REPORT ("%-30s\n", $content);
	for (my $i = 0 ; $i < $counter ; ++$i){
	    $analyseMultipleValues = $analyseRef->[0]->{$content}->[0]->{values}->[$i]->{value};
	    $analyseMultipleNames = $analyseRef->[0]->{$content}->[0]->{values}->[$i]->{name};
	    if($analyseMultipleNames) {
		printf REPORT ("%-60s: %s\n",$analyseMultipleNames, $analyseMultipleValues);
	    }
	}
	print_separator();
    }
    print REPORT "\n-------------------------- END OF ANALYSE SECTION --------------------------\n";
    close (REPORT);

    if ($optXMLCreate) {
	$XMLReport->startTag('analyse');
	foreach $content (@analysePatterns) {
	 #   print "--->Analyse Patterns: $content\n";
	    $analyseValue = $analyseRef->[0]->{$content}->[0]->{value};
	    $analyseUnit = $analyseRef->[0]->{$content}->[0]->{unit};
	    if ($analyseValue) {
	#	print "\$analyseUnit: $analyseUnit\n";
	#	print "\$analyseValue: $analyseValue\n";
		$XMLReport->startTag('result', name => $content, unit => $analyseUnit);
#		$XMLReport->emptyTag('result', "name" => $content,
#				     "value" =>$analyseValue,
#				     "unit" => $analyseUnit);
		$XMLReport->dataElement('value', $analyseValue);
		$XMLReport->endTag();
	    }
      	}

	foreach $content (@analysePatternsIndex) {
	    my $counter = $analyseRef->[0]->{$content}->[0]->{count};
	    if (!$counter) {next;}
	    $XMLReport->startTag('result', name => $content) if ($counter > 0);
	    $XMLReport->startTag('values') if ($counter > 0);
	    for(my $i=0 ; $i<$counter ; ++$i) {
		$analyseMultipleValues = $analyseRef->[0]->{$content}->[0]->{values}->[$i]->{value};
		$analyseMultipleNames = $analyseRef->[0]->{$content}->[0]->{values}->[$i]->{name};
		if ($analyseMultipleNames) {
		    $XMLReport->startTag('data', name => $analyseMultipleNames);
		    $XMLReport->dataElement('value', $analyseMultipleValues);
		    $XMLReport->endTag();
		}
	    }
	    $XMLReport->endTag() if ($counter > 0);
	    $XMLReport->endTag() if ($counter > 0);
	}
	$XMLReport->endTag();
    }
}

sub stdoutput_print {
    open (REPORT, '>>', $report) or die $!;
    if($stdoutRef->[0]->{content}) {
	print REPORT "\n\n---------------------------- OUTPUT SECTION ---------------------------\n";
	print REPORT $stdoutRef->[0]->{content};
	print REPORT "\n-------------------------- END OF OUTPUT SECTION --------------------------\n";


#	if ($optXMLCreate) {
#	    $XMLReport->startTag('output');
#
#	    $XMLReport->endTag();
#	}
    }
    close (REPORT);
}

sub result_print {
    open (REPORT, '>>', $report) or die $!; 
#    printf REPORT (": %s\n", $analyseParamsRef);
    close (REPORT);
}

sub footnote_print {
    open (REPORT, '>>', $report) or die $!; 
#    print  REPORT Dumper($benchrun);
#    print  Dumper($analysePattern);
    print_separator();
    printf REPORT ("\n\tneeded time for parsing %s: %6.4f sec\n", $logfile, $timeDiff);
    print_separator();
    close (REPORT);

    if ($optXMLCreate) {
# close XML file
	$XMLReport->endTag();
	$XMLReport->end();
	$xmlFileReport->close();
    }
}

sub version {
    print &getversion(); print "\n";
    exit(0);
}


sub getversion {
    my $text = "jubepp version 0.1p00";
    return $text;
}

sub usage {
    die "Usage: $_[0] <options> <name.longlog>
                -analyse           : prints out the data defined and requested resp. in analyse.xml and analyse-pattern.xml
                -long              : includes stdout file in the report 
                -Version           : prints out the current version
                -xmlcreate         : prints out a XML version of the report
                -param_long        : prints out a full list of defined parameters
";
}


sub print_separator {
print REPORT ("------------------------------------------------------------------------------------------\n");
}
