# $Source: /cvsroot/esa4-t3/DEISA_BENCH/bench/bench_data_dbhash.pm,v $
# $Author: mahermanns $
# $Revision: 1.1.1.1 $
# $Date: 2007/08/07 06:49:12 $
#
#
package jube_data_dbhash;
use strict;
use File::stat qw(:FIELDS);
use XML::Simple;
use Slurp;                     # Slurpt ganze Dateien in Variablen
use Data::Dumper;
use File::Listing;
use Time::HiRes qw ( time );

my $patint="([\\+\\-\\d]+)";   # Pattern for Integer number
my $patfp ="([\\+\\-\\d.E]+)"; # Pattern for Floating Point number
my $patwrd="([\^\\s]+)";       # Pattern for Work (all noblank characters)
my $patnint="[\\+\\-\\d]+";    # Pattern for Integer number, no () 
my $patnfp ="[\\+\\-\\d.E]+";  # Pattern for Floating Point number, no () 
my $patnwrd="[\^\\s]+";        # Pattern for Work (all noblank characters), no () 
my $patbl ="\\s+";             # Pattern for blank space (variable length)

$Data::Dumper::Indent=0;

my($debug)=2;

sub new {
    my $self  = {};
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my($tstart,$tdiff);
    printf("\t\tBench_Data_Dbhash: new %s\n",ref($proto)) if($debug>=3);
    $self->{VERBOSE}=0;
    $self->{BENCHDBREF}=0;

    $tstart=time;
    $self->{XS}=XML::Simple->new();
    $tdiff=time-$tstart;
    printf("WF: init XML-Parser in  %6.4f sec\n",$tdiff);
    
    bless $self, $class;
    return $self;
}

sub load_db {
    my($self) = shift;
    my($tstart,$tdiff);
    my ($VAR1);
    if(-f "benchdb.dump" ) {
	$tstart=time;
	my $dump = slurp( "benchdb.dump" ) # Inhalt der Dumper-Datei einlesen
	    or die $!;
	eval( $dump );    # Macht eine Datenstruktur daraus
	$self->{BENCHDBREF} = $VAR1;
	$tdiff=time-$tstart;
	printf("WF: reading database in %6.4f sec\n",$tdiff);
    }
    return(1);
}

sub save_db {
    my($self) = shift;
    my($tstart,$tdiff);
    $tstart=time;
    open(DB,">benchdb.dump");
    print DB Dumper($self->{BENCHDBREF});
    $tdiff=time-$tstart;
    printf("WF: dumping database in %6.4f sec\n",$tdiff);
    close(DB);
}

sub load_db_from_xml {
    my($self) = shift;
    my($tstart,$tdiff);
    if(-f "benchdb.dat" ) {
	$tstart=time;
	$self->{BENCHDBREF} = XMLin("benchdb.dat");
	$tdiff=time-$tstart;
	printf("WF: parsing benchdb.dat in %6.4f sec\n",$tdiff);
    }
}

sub save_db_to_xml {
    my($self) = shift;
    my($tstart,$tdiff);
    my $xml = XMLout($self->{BENCHDBREF});
    open(DB,">benchdb.dat");
    print DB $xml;
    close(DB);
}

sub scan_for_new_benchmarks {
    my($self) = shift;
    my($dir)=@_;
    my($tstart,$tdiff);
    my($rc,$fn,$afn,$n,$scanfile,$st);
    $n=0;
    $rc=opendir(DIR,$dir);
    while($fn=readdir(DIR)) {
	$afn="$dir/$fn";
	next if($fn!~/\.log$/);
	next if($fn!~/benchlog.*i$patint\_/);
	stat($afn) or die "No $afn: $!";
	$scanfile=1;
	if(exists($self->{BENCHDBREF}->{$fn})) {
	    if(exists($self->{BENCHDBREF}->{$fn}{SIZE}) && exists($self->{BENCHDBREF}->{$fn}{MTIME})) {
		$scanfile=0;
		$scanfile=1 if($st_size != $self->{BENCHDBREF}->{$fn}{SIZE});
		$scanfile=1 if($st_mtime != $self->{BENCHDBREF}->{$fn}{MTIME});
	    } 
	}
	if($scanfile) {
	    $n++;
	    printf("WF[%4d]: %s/%s\n",$n,$dir,$fn);
	    $self->{BENCHDBREF}->{$fn}{SIZE}=$st_size;
	    $self->{BENCHDBREF}->{$fn}{MTIME}=$st_mtime;
	    $tstart=time;
	    $self->{BENCHDBREF}->{$fn}=$self->{XS}->XMLin("$dir/$fn", KeyAttr => { 'map' => "n",
									   'benchmark' => "+name"
									   },
						  ForceArray => 1);
	    $tdiff=time-$tstart;
#	    printf("WF: parsing $dir/$fn in %6.4f sec\n",$tdiff);
	}
    }
}


sub clean_statdata {
    my($self) = shift;
    my($tstart,$tdiff,$id);
    foreach $id (keys(%{$self->{B_NAME}}))       {	delete($self->{B_NAME}->{$id});          }
    foreach $id (keys(%{$self->{B_PLATFORM}}))       {	delete($self->{B_PLATFORM}->{$id});          }
}

sub analyse {
    my($self) = shift;
    my($tstart,$tdiff);
    my($key,$val);
    while (($key,$val)=each %{$self->{BENCHDBREF}}) {
	my $name=$val->{name};
	my $platform=$val->{platform};
	$self->{B_NAME}->{$name}++;
	$self->{B_PLATFORM}->{$platform}++;
    }
}


1;


