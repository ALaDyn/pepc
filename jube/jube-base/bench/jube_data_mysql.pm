# $Source: /cvsroot/esa4-t3/DEISA_BENCH/bench/bench_data_mysql.pm,v $
# $Author: mahermanns $
# $Revision: 1.1.1.1 $
# $Date: 2007/08/07 06:49:12 $
#
#
package jube_data_mysql;
use strict;
use Data::Dumper;
use File::stat qw(:FIELDS);
use XML::Simple;
use File::Listing;
use Time::HiRes qw ( time );
use DBI;
$Data::Dumper::Indent = 1;

my $patint="([\\+\\-\\d]+)";   # Pattern for Integer number
my $patfp ="([\\+\\-\\d.E]+)"; # Pattern for Floating Point number
my $patwrd="([\^\\s]+)";       # Pattern for Work (all noblank characters)
my $patnint="[\\+\\-\\d]+";    # Pattern for Integer number, no () 
my $patnfp ="[\\+\\-\\d.E]+";  # Pattern for Floating Point number, no () 
my $patnwrd="[\^\\s]+";        # Pattern for Work (all noblank characters), no () 
my $patbl ="\\s+";             # Pattern for blank space (variable length)


my($debug)=2;

sub new {
    my $self  = {};
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my($tstart,$tdiff);

    $self->{HOST} = shift;
    $self->{DATABASE} = shift;
    $self->{USER} = shift;
    $self->{PASSWD} = shift;
    $self->{CONNECT} = "DBI:mysql:hostname=$self->{HOST}:database=$self->{DATABASE}";
    printf("\t\tBench_Data_mysql: new %s\n",ref($proto)) if($debug>=3);
    $self->{VERBOSE}=0;
    $tstart=time;
    print "WF: >$self->{CONNECT}<>$self->{USER}<>$self->{PASSWD}<\n";
    $self->{DBH} = DBI->connect($self->{CONNECT},$self->{USER},$self->{PASSWD}) || die "Cannot connect to DB";
    $self->{XS}=XML::Simple->new();
    $tdiff=time-$tstart;
    printf("WF: connect to database in  %6.4f sec\n",$tdiff);
    
    bless $self, $class;
    return $self;
}

sub DESTROY {
    my($self) = shift;
    my($tstart,$tdiff);
    $tstart=time;
    $self->{DBH}->disconnect() || die "Cannot disconnect to DB";
    $tdiff=time-$tstart;
    printf("WF: disconnect to database in  %6.4f sec\n",$tdiff);
}

sub load_db {
    my($self) = shift;
   
    if(!$self->table_exists("STATUS")) {
	$self->create_table("STATUS");	
	$self->init_table("STATUS");	

    }
    if(!$self->table_exists("BENCHRUNS")) {
	$self->create_table("BENCHRUNS");
    }
    if(!$self->table_exists("BENCHMARKS")) {
	$self->create_table("BENCHMARKS");
    }

    return(1);
}

sub save_db {
    my($self) = shift;
}


sub scan_for_new_benchmarks {
    my($self) = shift;
    my($dir)=@_;
    my($tstart,$tdiff);
    my($rc,$fn,$afn,$n,$scanfile,$st,$xmltree,$id,$subid,$nr,$snr,$benchname);
    my($NR,$ID,$FN,$FNSIZE,$MTIME);
    $n=0;


    my $sth = $self->{DBH}->prepare( q{ 
	SELECT NR,SUBID,FN,FNSIZE,MTIME FROM BENCHMARKS WHERE FN=?; 
    }) or die "Can't prepare statement: $DBI::errstr";


    $rc=opendir(DIR,$dir);
    while($fn=readdir(DIR)) {
	$afn="$dir/$fn";
	next if($fn!~/\.log$/);
	next if($fn!~/benchlog.*i$patint\_/);
	stat($afn) or die "No $afn: $!";

	$scanfile=1;
	# check if file is already in DB
	my $rc = $sth->execute($fn)  or die "Can't execute statement: $DBI::errstr";
	if(($NR,$ID,$FN,$FNSIZE,$MTIME) = $sth->fetchrow()) {
	    $scanfile=0 if (($FNSIZE==$st_size) && ($MTIME==$st_mtime));
	}
	$sth->finish();

	if($scanfile) {
	    $n++;
#	    printf("WF[%4d]: %s/%s\n",$n,$dir,$fn);
	    $tstart=time;
	    $xmltree=$self->{XS}->XMLin("$dir/$fn", KeyAttr => { 'map' => "n",
								 'benchmark' => "+name"
								 },
					ForceArray => 1);
	    $tdiff=time-$tstart;
#	    print Dumper($xmltree);
	    ($benchname)=(keys(%{$xmltree->{'benchmark'}}));
	    $id=$xmltree->{'benchmark'}->{$benchname}->{'identifier'};
	    $subid=$xmltree->{'benchmark'}->{$benchname}->{'subid'};
	    $nr=$self->add_to_benchruns($id,$fn,$st_size,$st_mtime);
#	    $snr=$self->add_to_benchmarks($subid,$fn,$st_size,$st_mtime);
	    printf("WF[%4d]: parsing $dir/$fn in %6.4f sec\n",$n,$tdiff);
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


#####################################################################################################
# Utility functions
#####################################################################################################

sub table_exists {
    my($self) = shift;
    my($tablename)=@_;
    my($table,$tablespace,$found);
    my $sth = $self->{DBH}->prepare( q{ 
	show tables; 
    }) or die "Can't prepare statement: $DBI::errstr";
    
    my $rc = $sth->execute() or die "Can't execute statement: $DBI::errstr";
    
    $found=0;
    while(($table) = $sth->fetchrow()) {
	print "WF: table_exists: >$table<\n";
	$found=1 if($table eq $tablename);
    }
    $sth->finish();
    return($found);
}

sub create_table {
    my($self) = shift;
    my($tablename)=@_;
    my($createstmt);
    my($tstart,$tdiff);

    $createstmt=qq{
	CREATE table STATUS
	    (
	     NAME              VARCHAR(20)  NOT NULL,
	     VAL               INTEGER(12)  NOT NULL
	     )
	} if ($tablename eq "STATUS");
    $createstmt=qq{
	CREATE table BENCHRUNS
	    (
	     NR                  INTEGER(12)  NOT NULL,
	     ID                  VARCHAR(255)  NOT NULL,
	     COUNTBM             INTEGER(12)   NOT NULL
	     )
	} if ($tablename eq "BENCHRUNS");
    $createstmt=qq{
	CREATE table BENCHMARKS
	    (
	     NR                  INTEGER(12)  NOT NULL,
	     NR_IN_BENCHMARKS    INTEGER(12)  NOT NULL,
	     SUBID               VARCHAR(255)  NOT NULL,
	     FN                  VARCHAR(255)  NOT NULL,
	     FNSIZE              INTEGER(12)   NOT NULL,
	     MTIME               INTEGER(12)   NOT NULL
	     )
	} if ($tablename eq "BENCHMARKS");

    $tstart=time;
    my $sth = $self->{DBH}->prepare( $createstmt ) or die "Can't prepare statement: $DBI::errstr";
    my $rc = $sth->execute()                       or die "Can't execute statement: $DBI::errstr";
    $sth->finish();
    $tdiff=time-$tstart;
    printf("WF: table $tablename created ($rc) in %6.4f sec\n",$tdiff);

    return($rc);
}

sub init_table {
    my($self) = shift;
    my($tablename)=@_;
    my $createstmt=0;
    my($tstart,$tdiff,$rc);

    $createstmt=qq{
	INSERT INTO STATUS (NAME,VAL) VALUES('BENCHRUNS_NEXTID',0); 
	INSERT INTO STATUS (NAME,VAL) VALUES('BENCHMARKS_NEXTID',0); 
    } if ($tablename eq "STATUS");

    if($createstmt) {
	$tstart=time;
	my $sth = $self->{DBH}->prepare( $createstmt ) or die "Can't prepare statement: $DBI::errstr";
	$rc = $sth->execute()                       or die "Can't execute statement: $DBI::errstr";
	$sth->finish();
	$tdiff=time-$tstart;
	printf("WF: table $tablename created ($rc) in %6.4f sec\n",$tdiff);
    }
    return($rc);
}

sub set_next_id {
    my($self) = shift;
    my($table,$nextid)=@_;
    my($sth);
    my $entry=$table."_NEXTID";
    $sth = $self->{DBH}->prepare( q{ 
	UPDATE STATUS SET VAL=? WHERE NAME =?; 
    }) or die "Can't prepare statement: $DBI::errstr";
    my $rc = $sth->execute($nextid,$entry) or die "Can't execute statement: $DBI::errstr";
    $sth->finish();
    return($rc);
}

sub get_next_id {
    my($self) = shift;
    my($table) = shift;
    my $nextid=-1;
    my $rnextid;
    my $entry=$table."_NEXTID";
    my $sth = $self->{DBH}->prepare( q{ 
	SELECT VAL FROM STATUS WHERE NAME =?; 
    }) or die "Can't prepare statement: $DBI::errstr";
    my $rc = $sth->execute($entry) or die "Can't execute statement: $DBI::errstr";
    if(($rnextid) = $sth->fetchrow()) {
	$nextid=$rnextid;
    }
    $sth->finish();
    
#    printf("WF: get_next_id: %d\n",$nextid);
    
    return($nextid);
}

sub add_to_benchruns {
    my($self) = shift;
    my($id)=@_;
    my($nr,$rid,$countbm,$rc,$inserted,$rc,$sth,$nextid);

    # check if already inserted
    $sth = $self->{DBH}->prepare( q{ 
	SELECT NR,ID,COUNTBM FROM BENCHRUNS WHERE ID=?; 
    }) or die "Can't prepare statement: $DBI::errstr";
    $rc = $sth->execute($id) or die "Can't execute statement: $DBI::errstr";
    if(($nr,$rid,$countbm) = $sth->fetchrow()) {
	$inserted=1;
    } else {
	$inserted=0;
	$countbm=0;
    }
    $sth->finish();
    #update or insert new entry
    $countbm++;
    if($inserted) {
	$sth = $self->{DBH}->prepare( q{ 
	    UPDATE BENCHRUNS SET COUNTBM=? WHERE NR=?; 
	}) or die "Can't prepare statement: $DBI::errstr";
	$rc = $sth->execute($countbm,$nr) or die "Can't execute statement: $DBI::errstr";
    } else {
	$nr=$nextid=$self->get_next_id("BENCHRUNS");
	$sth = $self->{DBH}->prepare( q{ 
	    INSERT INTO BENCHRUNS (NR,ID,COUNTBM) VALUES(?,?,?); 
	}) or die "Can't prepare statement: $DBI::errstr";
	$rc = $sth->execute($nextid,$id,$countbm) or die "Can't execute statement: $DBI::errstr";
	$nextid++;
	$self->set_next_id("BENCHRUNS",$nextid);
    }
    return($nr);
}

sub insert_in_benchmarks {
    my($self) = shift;
    my($id,$fn,$st_size,$st_mtime)=@_;
    my($nr,$nextid,$sth);
    $nr=$nextid=$self->get_next_id();
    $sth = $self->{DBH}->prepare( q{ 
	INSERT INTO BENCHMARKS (NR,ID,FN,FNSIZE,MTIME) VALUES(?,?,?,?,?); 
    }) or die "Can't prepare statement: $DBI::errstr";
    my $rc = $sth->execute($nr,$id,$fn,$st_size,$st_mtime) or die "Can't execute statement: $DBI::errstr";
    $sth->finish();
    $nextid++;
    $self->set_next_id($nextid);
    return($nr);
}


1;


