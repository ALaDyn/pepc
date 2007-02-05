#!/usr/local/bin/perl 

BEGIN {
    unshift(@INC,'/var/local/www/kfa/zam/pepc/cgi-bin');
}

use CGI;
use CGI::Carp 'fatalsToBrowser'; # show perl error messages on browser
use Mail::Send;
use FZJweb;

my $marker="<font color=\"red\">*</font>";
my $cnt,$cntreq,$emailok;
my @data;

$query = new CGI;
#print $query->header;
#print $query->start_html(-title=>'Download PEPC',
#			 -BGCOLOR=>'white');

my $page = new FZJweb(
		      -oe => "Zentralinstitut f&uuml;r Angewandte Mathematik",
		      -oe_short => "ZAM",
		      -author => "Pepc",
		      -author_email => "pepc.zam\@fz-juelich.de",
		      -date => "30.11.2005",
		      -page_title  => "COMPLEX ATOMISTIC MODELLING AND SIMULATION",
		      -menu_file   => "/var/local/www/fzj/zam/pepc/cgi-bin/menu",
		      -menu_pos    => ["PEPC","Download & Install"]
		      );

print $page->header;


print    $query->h1($query->em('PEPC'),'download page');


print    "<P><B>Registration Information</B>";

#print    "<P><font color=red>!!! The download page is under construction !!!</font></B>";

print    "<P>(fields required are marked with $marker)\n";  
print    $query->start_form(POST);
print    "<table border=\"0\" width=\"100%\">";

print    "<tr><td width=\"110\">Full Name$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'fullname', -size => 80),"</td></tr>\n";

print    "<tr><td width=\"110\">Company Name:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'companyname', -size => 70),"</td></tr>\n";

print    "<tr><td width=\"110\">Email$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'email', -size => 80),"</td></tr>\n";

print    "<tr><td width=\"110\">Telephone:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'telephone', -size => 40),"</td></tr>\n";

print    "<tr><td width=\"110\">Address$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'address1', -size => 70),"</td></tr>\n";

print    "<tr><td width=\"110\">Address 2nd:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'address2', -size => 60),"</td></tr>\n";

print    "<tr><td width=\"110\">Address 3rd:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'address3', -size => 70),"</td></tr>\n";

print    "<tr><td width=\"110\">City$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'city', -size => 80),"</td></tr>\n";

print    "<tr><td width=\"110\">Zipcode$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'zip', -size => 20),"</td></tr>\n";

print    "<tr><td width=\"110\">State/Province$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'state', -size => 80),"</td></tr>\n";

print    "<tr><td width=\"110\">Country$marker:</td>\n<td width=\"275\">";
print     $query->textfield(-name => 'country', -size => 80),"</td></tr>\n";

print    "</table>";
print    $query->p;
print    $query->submit('Get Access');
print    $query->end_form;
print    $query->hr;

$data[0]= $query->param('fullname');
$data[1]= $query->param('companyname');
$data[2]= $query->param('email');
$data[3]= $query->param('telephone');
$data[4]= $query->param('address1');
$data[5]= $query->param('address2');
$data[6]= $query->param('address3');
$data[7]= $query->param('city');
$data[8]= $query->param('state');
$data[9]= $query->param('zip');
$data[10]=$query->param('country');

$cnt=0;
for($i=0;$i<11;$i++) {
    $cnt++ if($data[$i]);
}

if($cnt>0) {
    # form filled
    print $query->p;
    print $query->b('Access:');
    print $query->br;
    $cntreq=0;
    for $i (0,2,4,7,8,9,10) {
	if($data[$i]) {
	    $cntreq++;
	} else {
	    print "&#160;<font color=\"red\"> no information in field \#",$i+1,"</font><br>\n";
	}
    }
    # check email
    if($data[2]=~/^(.+)\@(\w+\.)?(.+\.)?([^\.]+)\.([^\.]+)$/) {
#	print "&#160;<font color=\"blue\"> email ok ($1,$2,$3,$4,$5)</font><br>\n";
	$emailok=1;
    } else {
	print "&#160;<font color=\"red\"> email NOT ok</font> ($data[2])<br>\n";
	$emailok=0;
    }

    if(($cntreq==7) && ($emailok==1)) {

#	print "from: ",$query->remote_host();
#	print $query->p;
	my $mypath="/var/local/www/kfa/zam/pepc/tmp";
	
	if (-f "$mypath/cnt.dat") {
	    open(IN,"$mypath/cnt.dat");
	    $cnt=<IN>;
	    chomp($cnt);
	    close(IN);
	} else {
	    $cnt=0;
	}
	$cnt++;
	open(OUT,"> $mypath/cnt.dat");
	print OUT "$cnt\n";
	close(OUT);

	my $nr=$cnt%10;
	my $myname=sprintf("%04d",$nr);
	my $tmpurl="http://www.kfa-juelich.de/zam/pepc/tmp/tmpfull$myname";
	my $tmppath="$mypath/tmpfull$myname";

	
	if (&update_sub_dir($tmppath)) {
	    print $query->p;
#	    print $query->b("$cnt, $nr");
#	    print $query->br;
	    print $query->b('Version 1.1 (January 2006):');
	    
	    print $query->ul(
			     $query->li("Source tar file",
				   $query->br,
					$query->a({href=>"$tmpurl/pepc1.1.tar.gz"},
						  'pepc1.1.tar.gz')
					)
			     );
	    print $query->b('Version 1.4 beta (last update: 26/11/2006):');
	    print $query->ul(
			     $query->li("Source tar file",
				   $query->br,
					$query->a({href=>"$tmpurl/pepc_1.4.tar.gz"},
						  'pepc_1.4.tar.gz')
					)
			     );

	    print $query->b('Demos');
	    print $query->ul(
			     $query->li("Demo tar file",
				   $query->br,
					$query->a({href=>"$tmpurl/pepc_demos1.0.tar.gz"},
						  'pepc_demos1.0.tar.gz')
					)
			     );
		
	
	    print $query->p, $query->i('Use', $query->b('shift+left mouse button'), 
				       'for download (Mozilla only)!'), $query->p;

	    print $query->hr;

	    $printlog.="="x60;$printlog.="\n";
	    $printlog =localtime()."\n";
	    $printlog.="-"x60;$printlog.="\n";
	    $printlog.="from host:      ".$query->remote_host()."\n";
	    $printlog.="cnt:            ".$cnt."\n";
	    $printlog.="nr:             ".$nr."\n";
	    $printlog.="tmpdir:         ".$tmppath."\n";
	    $printlog.="Full Name:      ".$data[0]."\n";
	    $printlog.="Company Name:   ".$data[1]."\n";
	    $printlog.="Email:          ".$data[2]."\n";
	    $printlog.="Telephone:      ".$data[3]."\n";
	    $printlog.="Address1:       ".$data[4]."\n";
	    $printlog.="Address2:       ".$data[5]."\n";
	    $printlog.="Address3:       ".$data[6]."\n";
	    $printlog.="City:           ".$data[7]."\n";
	    $printlog.="State/Province: ".$data[8]."\n";
	    $printlog.="Zip-/Postcode:  ".$data[9]."\n";
	    $printlog.="Country:        ".$data[10]."\n";
	    $printlog.="="x60;$printlog.="\n\n";

	    open(OUT,">>$mypath/access.log");
	    print OUT $printlog;
	    close(OUT);

	    open(OUT,">$mypath/access_$cnt.log");
	    print OUT $printlog;
	    close(OUT);

	    system("/var/local/www/kfa/zam/pepc/cgi-bin/smail.pl $mypath/access_$cnt.log $cnt")

	}
    }
}

print $page->footer;

#print $query->a({href=>'http://www.kfa-juelich.de/zam/pepc'},'Back to PEPC main page');
#print $query->end_html;
    
    
sub update_sub_dir {
    my($tmppath)=@_;

    $softwarepath="/var/local/www/kfa/zam/pepc/software";
    if (! (-d $tmppath)) {
	mkdir($tmppath,0755);
#	print "<BR>--mkdir $tmppath ($?)<BR>";
    } else {
	unlink("$tmppath/pepc1.1.tar.gz");
	unlink("$tmppath/pepc_1.4.tar.gz");
	unlink("$tmppath/pepc_demos1.0.tar.gz");
    }
    symlink("$softwarepath/pepc_1.4.tar.gz","$tmppath/pepc_1.4.tar.gz");
    symlink("$softwarepath/pepc_demos1.0.tar.gz","$tmppath/pepc_demos1.0.tar.gz");
   
    return(1);
}
