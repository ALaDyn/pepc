#!/usr/local/bin/perl 

BEGIN {
    unshift(@INC,'/var/local/www/kfa/zam/pepc/cgi-bin');
}

use Mail::Send;
my $filename=@ARGV[0];
my $cnt=@ARGV[1];

$data="";
open(IN,$filename) or die "file not found $filename\n";
while(<IN>) {
    $data.=$_;
}
close(IN);

# Mail to pepc.zam@fz-juelich.de
$mail=Mail::Send->new(
		      Subject => "Access to PEPC software #$cnt",
		      To      => "pepc.zam\@fz-juelich.de"
				  );
$mailhandle=$mail->open();
print $mailhandle "\n";
print $mailhandle $data;
$mailhandle->close();

