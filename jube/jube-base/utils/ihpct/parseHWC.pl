#!/usr/bin/perl
use XML::Parser;
use Data::Dumper;

use strict;
use warnings;

my %DataNames;

my @IPDEFList;
my @IPList;

our $currentlabel = "";

my $p = XML::Parser->new(Style => 'Tree');

my $toParseFileNumber = @ARGV;

#print "number of files to parse $toParseFileNumber\n";

my $fileCnt;
for($fileCnt=0; $fileCnt < $toParseFileNumber; $fileCnt++)
{

    my $tree = $p->parsefile($ARGV[$fileCnt]);
    
    $tree = ${$tree}[1];

    my $cnt;

#print "number of vector elements: $#{$tree}\n";
    for($cnt=0; $cnt<$#{$tree}; $cnt++)
    {
#    print "\t$cnt: ${$tree}[$cnt]\n";
    }
    
    for($cnt=1; $cnt < $#{$tree}; $cnt+=2)
    {
#    print "check $cnt -th element\n"; 
	if(${$tree}[$cnt] eq "ipdef") 
	{
#	print "\tfound ipdef element\n";
	    push(@IPDEFList, ${$tree}[$cnt+1]);
	}
	if(${$tree}[$cnt] eq "ip") 
	{
#	print "\tfound ip element\n";
	    push(@IPList, ${$tree}[$cnt+1]);
	}
	
    }
    
    my $IPDEFSize = @IPDEFList;
    my $IPSize    = @IPList;
#print "number of ipdef elements: $IPDEFSize\n";
#print "number of ip elements: $IPDEFSize\n";
    
    my $ipdefCnt;
    for($ipdefCnt=0; $ipdefCnt<$IPDEFSize; $ipdefCnt++)
    {
#    print Dumper($IPDEFList[$ipdefCnt]);
	my $elCnt;
	my $elNum = $#{$IPDEFList[$ipdefCnt]} + 1;
	
#    print "parse IPDEF list no. $ipdefCnt with $elNum elements\n";
	
	for($elCnt=0; $elCnt<$elNum; $elCnt++)
	{
#	print "checking element $elCnt\n";
	    if(${$IPDEFList[$ipdefCnt]}[$elCnt] eq "datadef") 
	    {
#	    print "\tfound datadef element\n";
		my $id    = ${${$IPDEFList[$ipdefCnt]}[$elCnt+1]}[0]{'id'};
		my $label = ${${$IPDEFList[$ipdefCnt]}[$elCnt+1]}[0]{'label'};
#	    print "\t\tid: $id, label: $label\n";
		$DataNames{$id} = $label;
	    }
	}
    }
    
#print Dumper(%DataNames);
    
    my $ipCnt;
    for($ipCnt=0; $ipCnt<$IPSize; $ipCnt++)
    {
#    print Dumper($IPDEFList[$ipdefCnt]);
	my $elCnt;
	my $elNum = $#{$IPList[$ipCnt]} + 1;
	
	my $IPLabel = ${$IPList[$ipCnt]}[0]{'label'};
	
#    print "parse IP list no. $ipdefCnt with $elNum elements, labeled $IPLabel\n";
	
	for($elCnt=0; $elCnt<$elNum; $elCnt++)
	{
#	print "checking element $elCnt\n";
	    if(${$IPList[$ipCnt]}[$elCnt] eq "d") 
	    {
#	    print "\tfound d element\n";
		my $id = ${${$IPList[$ipCnt]}[$elCnt+1]}[0]{'id'};
		my $v  = ${${$IPList[$ipCnt]}[$elCnt+1]}[0]{'v'};
		print "IHPCT: libHPM: in section $IPLabel: $DataNames{$id}: $v\n";
#	    $DataNames{'$id'} = $label;
	    }
	}
    }
}
    
