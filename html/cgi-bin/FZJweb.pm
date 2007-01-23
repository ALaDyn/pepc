#
#  FZJweb.pm       (Th.Eickermann, 30.10.2004)
#
#  a module for creating Web-Pages that conform to the
#  FZJ-style. The following public methods are provided:
#
#  $page = new FZJweb(...); # creates a FZJweb Page object
#                           # see usage below for a description of parameters  
#  $page->header;           # returns a string with the HTML-code for
#                           # the page header and the navigation menu
#  $page->footer;           # returns a string with the HTML code for
#                           # the page footer
#
#  The variabel part of the page is defined via a parameter-hash of new.
#  The following parameters are required:
#
#  To appear in the footer:
#   -oe               full name of the OE
#   -oe_short         abbriviation of the OE name
#   -author           full name of responsible author
#   -author_email     email-address of author
#   -date             date of page-creation/update
#
#  To appear in the header:
#   -page_title       title of page: is printed above menu on left side
#   -menu_file        a file that describes the contents of the menu on left side (see below)
#   -menu_pos         description of this page position in the menu (see below)
#
#  Optional parameters (default values in []):
#   -OE               OE name to appear on top of page
#                     [ <upper-case of oe> ]
#   -oe_home          URL of oe-Homepage
#                     "http://www.fz-juelich.de/<lower-case of oe_short>"
#   -uplinks          a list of up-links appearing at the top of the menu on the left
#                     [ <oe_short>-Startseite ]
#   -html_title       title of page (usually printed in browser window-header)
#                     [ <oe_short> - <page_title> ]
#   -description      text appearing as comment in HTML-code
#                     [ <page_title> ]
#   -keywords         text appearing as comment in HTML-code
#                     [ <page_title> ]
#   -navigation_text  text printed above page content
#                     [ current position in menu ]
#   -full_url         URL of current page (not the print-version) ,is only used in print-version
#                     [ 'N/A' ]
#   -print            if non-zero: create a print-version of the page
#   -print_url        URL of print-version of this file
#                     [ no print-button is created ]
#   -year             year to be used in copyright notice
#                     [ year from <date> ]
#   -english          URL of english version of page
#                     [ no english page ]
#   -german           URL of german version of page
#                     [ no german page ]
#   example:
#
#   $page = new FZJweb(-oe => "Zentralinstitut f&uuml;r Angewandte Mathematik",
#                      -oe_short => "ZAM",
#                      -author => "Thomas Eickermann",
#                      -author_email => "th.eickermann@fz-juelich.de",
#                      -date => "30.10.2004",
#
#   ChangeLog:
#   05/04/2005 Th.Eickermann: bug-fixes: menu-file parser (last line problems),
#                                        allow syntax 'FZJweb::new' AND 'FZJweb->new'
#   11/10/2004 Th.Eickermann: added options -english/-german, fixed a missing ${date}
#

use strict;
use IO::File;

package FZJweb;
use Carp;

# zam121
#my $DEFAULT_COMMON_BASE = "https://zam121.zam.kfa-juelich.de/FZJweb/common";
my $DEFAULT_COMMON_BASE = "http://www.fz-juelich.de/common";
my $PICTURE_PATH = "pictures";

# www.fz-juelich.de
#my $DEFAULT_COMMON_BASE = "http://www.fz-juelich.de/common";
#my $PICTURE_PATH = "pictures";

#my $DEFAULT_COMMON_BASE = "common";

my $template_header;
my $template_footer;
my $template_uplink;
my $template_menu1;
my $template_menu2;
my $template_menu2_start;
my $template_menu2_stop;
my $template_pfad;
my $template_print;
my $template_print_header;
my $template_print_footer;


# parameters for the optional language button
my $language_html = '<a href="${language_url}"><img src="${common_base}/icons/${language_icon}" alt="english" name="speech" width="22" height="17" border="0"></a>';
my $no_language   = '';

my $german_icon  = "icn_deut_0.gif";
my $english_icon = "icn_engl_0.gif";

my $german_text   = "deutsche Version";
my $english_text  = "English Version";


sub usage {
    croak("FZJweb::new called with invalid parameters");
}
sub new {
    my $class = "FZJweb";

    $class = shift unless( $_[0] =~ /^-/ );

    my $page = { @_ };

    usage() unless defined $page->{-oe};
    usage() unless defined $page->{-oe_short};
    usage() unless defined $page->{-author};
    usage() unless defined $page->{-author_email};
    usage() unless defined $page->{-date};

    usage() unless defined $page->{-page_title};
    usage() unless defined $page->{-menu_file};
    usage() unless defined $page->{-menu_pos};

    # set default values of optional parameters where required
    #
    if( ! defined $page->{-OE} ) {
	$page->{-OE} = "$page->{-oe} ($page->{-oe_short})";
	$page->{-OE} =~ tr/a-z/A-Z/;
	$page->{-OE} =~ s/UML;/uml;/g;
    }

    if( ! defined $page->{-oe_home} ) {
	my $oe =  $page->{-oe_short};
	$oe =~ tr/A-Z/a-z/;
	$page->{-oe_home} = "http://www.fz-juelich.de/$oe";
    }

    if( ! defined $page->{-uplinks} ) {
	$page->{-uplinks} = [ [ "$page->{-oe_short}-Startseite", $page->{-oe_home}] ];
    }

    if( ! defined $page->{-html_title} ) {
	$page->{-html_title} = "$page->{-oe_short} - $page->{-page_title}";
    }
    
    if( ! defined $page->{-description} ) {
	$page->{-description} = $page->{-page_title};
    }

    if( ! defined $page->{-keywords} ) {
	$page->{-keywords} = $page->{-page_title};
    }

    if( ! defined $page->{-navigation_text} ) {
	$page->{-navigation_text} = join(" > ", @{$page->{-menu_pos}});
    }

    if( ! defined $page->{-print} ) {
	$page->{-print} = 0;
    }

    if( ! defined $page->{-year} ) {
	my @d = split(/\./, $page->{-date});
	$page->{-year} = $d[-1];
    }

    if( ! defined $page->{-full_url} ) {
	$page->{-full_url} = "N/A";
    }

    if( ! defined $page->{-common_base} ) {
	$page->{-common_base} = $DEFAULT_COMMON_BASE;
    }


    if( defined $page->{-english} ) {
	$page->{-language_html} = $language_html;
	$page->{-language_url}  = $page->{-english};
	$page->{-language_icon} = $english_icon;
	$page->{-language_text} = $english_text;
    } elsif( defined $page->{-german} ) {
	$page->{-language_html} = $language_html;
	$page->{-language_url}  = $page->{-german};
	$page->{-language_icon} = $german_icon;
	$page->{-language_text} = $german_text;
    } else {
	$page->{-language_html} = $no_language;
	$page->{-language_text} = $no_language;
 	$page->{-language_url}  = "";   }

    bless $page, "FZJweb";

    $page->_make_menu;

    return $page;
}

sub header {
    my $page = shift;

    if( ! $page->{-print} ) {
	return $page->_replace_vars( $template_header);
    } else {
	return $page->_replace_vars( $template_print_header);
    }
}

sub footer {
    my $page = shift;

    if( ! $page->{-print} ) {
	return $page->_replace_vars( $template_footer);
    } else {
	return $page->_replace_vars( $template_print_footer);
    }
}

sub _replace_vars {
    my($page,$text) = @_;

    $text =~ s/\${html_title}/$page->{-html_title}/g;
    $text =~ s/\${page_title}/$page->{-page_title}/g;
    $text =~ s/\${description}/$page->{-description}/g;
    $text =~ s/\${navigation_text}/$page->{-navigation_text}/g;
    $text =~ s/\${author}/$page->{-author}/g;
    $text =~ s/\${author_email}/$page->{-author_email}/g;
    $text =~ s/\${keywords}/$page->{-keywords}/g;
    $text =~ s/\${date}/$page->{-date}/g;
    $text =~ s/\${oe}/$page->{-oe}/g;
    $text =~ s/\${oe_short}/$page->{-oe_short}/g;
    $text =~ s/\${oe_home}/$page->{-oe_home}/g;
    $text =~ s/\${OE}/$page->{-OE}/g;
    $text =~ s/\${year}/$page->{-year}/g;
    $text =~ s/\${full_url}/$page->{-full_url}/g;
##    $text =~ s/\${}/$page->{-}/g;
    $text =~ s/\${menu_html}/$page->{-menu_html}/;

    if( defined $page->{-print_url} ) {
	$text =~ s/\${print_html}/$template_print/g;
	$text =~ s/\${print_url}/$page->{-print_url}/g;
    } else {
	$text =~ s/\${print_html}//g;
    }
    my $uplinks_html = "";
    if( defined $page->{-uplinks} ) {
	foreach my $uplink ( @{$page->{-uplinks}} ) {
	    my $add = $template_uplink;
	    $add =~ s/\${uplink_url}/$uplink->[1]/;
	    $add =~ s/\${uplink_text}/$uplink->[0]/;
	    $uplinks_html .= $add;
	}
    }
    $text =~ s/\${uplinks_html}/$uplinks_html/g;
    $text =~ s/\${language_html}/$page->{-language_html}/g;
    $text =~ s/\${language_url}/$page->{-language_url}/g;
    $text =~ s/\${language_icon}/$page->{-language_icon}/g;
    $text =~ s/\${language_text}/$page->{-language_text}/g;

    $text =~ s/\${common_base}/$page->{-common_base}/g;
    $text =~ s/\${pictures}/$PICTURE_PATH/g;
    return $text;
}

#  setup the menu
#
sub _make_menu {
    my $page = shift;

    #  read and parse menu-description
    #
    if( ! defined $page->{-menu_file} ) {
	$page->{-menu_html} = "";
	return;
    }
    my $in = new IO::File;
    $in->open("<$page->{-menu_file}") || die "Could not open '$page->{-menu_file}'";

    my $mlines = [];
    my $last_nspace = 0;
    my $ok;

    do {
	$ok = $page->_read_add($in, $mlines, 0, $last_nspace, undef);
    } while( $ok );

    $in->close();

    $page->{-menu} = $mlines;

    #  create the HTML-code for the menu
    #
    my $menu = "";

    foreach my $mline1 ( @{$page->{-menu}} ) {
	my($text1, $mfile1, $level1, $sub1) = @$mline1;

## aimed at beeing generic - but isn't so far
##	my $here = ( ($page->{-menu_pos}[0] eq $text1) && ($level == $level1) ) ? '_here' : '';
	my $here = ( ($page->{-menu_pos}[0] eq $text1) ) ? '_here' : '';

	my $add = $template_menu1;
	$add =~ s/\$menu_text/$text1/;
	$add =~ s/\$menu_html/$mfile1/;
	$add =~ s/\$here/$here/;

	$menu .= $add;
	
## aimed at beeing generic - but isn't so far
##	if( defined $sub1 && (($here ne '') || ($level>0) && ($text1 eq $parents->[$level1][0])) ) {
	if( defined $sub1 && ($here ne '') ) {
	    $menu .= $page->_add_sub(1, $sub1);
	}
    }
    $page->{-menu_html} = $menu;
}

sub _add_sub {
    my $page = shift;
    my($level, $mline1) = @_;
    my($text1, $mfile1, $level1, $sub1) = @$mline1;
    
    my $add = $template_menu2_start;

    foreach my $mline1 ( @$mline1 ) {
	my($text1, $mfile1, $level1, $sub1) = @$mline1;
	
	my $here = ( (defined $page->{-menu_pos}[$level1]) &&
		     ($page->{-menu_pos}[$level1] eq $text1) &&
		     ($level == $level1) ) ? '_here' : '';
	
	$add .= $template_menu2;
	$add =~ s/\$menu2_text/$text1/;
	$add =~ s/\$menu2_html/$mfile1/;
	$add =~ s/\$here/$here/;
    }
    $add .= $template_menu2_stop;
    
    return $add;
}

sub _read_add {
    my $page = shift;
    my( $in, $mlines, $level, $last_nspace, $top_lines) = @_;
    
    while( my $line = <$in> ) {
	
	next unless($line =~ /(\s*)\"([^"]+)\"\s+(\S+)/); #"])) # make emacs happy
	    
	my($space, $text, $mfile) = ($1, $2, $3);

        my $nspace = length($space);
        if( $nspace == $last_nspace ) {
            push( @$mlines, [ $text, $mfile, $level, undef ] );
        } elsif( $nspace > $last_nspace ) {
            my $new_mlines = [ [ $text, $mfile, $level+1, undef ] ];
            my $status = $page->_read_add($in, $new_mlines, $level+1, $nspace, $mlines);
            $mlines->[$#{$mlines}-$status][3] = $new_mlines;
            #push( @$mlines, [ $text, $mfile, $new_mlines ] );
        } else {
            push(@$top_lines, [ $text, $mfile, $level-1, undef ]);
            return 1;
        }
    }
    return 0;
}

################################################################################
#
# start of template texts
#
$template_header = <<'EOT';
Content-Type: text/html; charset=ISO-8859-1


<html>
    <head>
        <title>${html_title}</title>
        
    <!-- Meta-Tags -->
    <!-- Contentory Version 2.1 f-->
    <META HTTP-EQUIV="PRAGMA" CONTENT="no-cache">
    <META HTTP-EQUIV="CACHE-CONTROL" CONTENT="PRIVATE">
    <META HTTP-EQUIV="EXPIRES" CONTENT="0">
    <META HTTP-EQUIV="CONTENT-TYPE" CONTENT="text/html; charset=iso-8859-15">
    <META HTTP-EQUIV="IMAGETOOLBAR" CONTENT="no">
    <META NAME="AUTHOR" CONTENT="${author}">
    <META name="EMAIL" CONTENT="${author_email}">
    <META NAME="KEYWORDS" CONTENT="${keywords}">
    <META NAME="GENERATED" CONTENT="${date}">
    <META NAME="DESCRIPTION" CONTENT="${description}">
    <META NAME="oe" CONTENT="${oe_short}">
    <META NAME="MSSMARTTAGSPREVENTPARSING" CONTENT="TRUE">
    <!-- /Meta-Tags -->
        
        <link rel="stylesheet" href="${common_base}/css/fzj.css" type="text/css">
    </head>
<body bgcolor="#FFFFFF" text="#000000" link="#3366FF" vlink="#663366" alink="#663366" topmargin=0 marginheight=0 marginwidth=0 leftmargin=0>
<!--  Portal -->
<!-- --------------------HEADLINE PORTAL PLUGIN--------------------	-->
<table width="777" border="0" cellspacing="0" cellpadding="0">
  <tr>
    <td valign="top" rowspan="2" width="188">
      <a href="http://www.fz-juelich.de/"><img src="${common_base}/icons/fzj_logo_sm.gif" alt="FZ-Juelich" border="0"></a>
    </td>
    <td width="1" background="${common_base}/icons/portal_line.gif">
      <img src="${common_base}/icons/spacer.gif" width="1" height="19">
    </td>
    
    <td width="144" align="left" valign="top">
      <table width="144" border="0" cellpadding="0" cellspacing="0">
        <tr>
          <td width="8" height="12"><img src="${common_base}/icons/spacer.gif"></td>
          <td width="136"><a href="http://www.fz-juelich.de/portal/forschungszentrum"><img name="portal1" height="19" src="${common_base}/icons/portal_zentrum_0.gif" alt="Das Forschungzentrum" border="0"></a></td>
        </tr>
        
      </table>
    </td>
    
    <td width="1" background="${common_base}/icons/portal_line.gif">
      <img src="${common_base}/icons/spacer.gif" width="1" height="12">
    </td>
    
    
    <td width="144" align="left" valign="top">
      <table width="144" border="0" cellpadding="0" cellspacing="0">
        <tr>
          <td width="8" height="12"><img src="${common_base}/icons/spacer.gif"></td>
          <td width="136"><a href="http://www.fz-juelich.de/portal/forschung"><img name="portal2" height="19" src="${common_base}/icons/portal_forschung_0.gif" alt="Forschung" border="0"></a></td>
        </tr>
        
      </table>
    </td>
    
    <td width="1" background="${common_base}/icons/portal_line.gif">
      <img src="${common_base}/icons/spacer.gif" width="1" height="12">
    </td>
    
    
    <td width="144" align="left" valign="top">
      <table width="144" border="0" cellpadding="0" cellspacing="0">
        <tr>
          <td width="8" height="12"><img src="${common_base}/icons/spacer.gif"></td>
          <td width="136"><a href="http://www.fz-juelich.de/portal/wissen"><img name="portal3" height="19" src="${common_base}/icons/portal_wissen_0.gif" alt="Wissen" border="0"></a></td>
        </tr>
        
      </table>
    </td>
    
    <td width="1" background="${common_base}/icons/portal_line.gif">
      <img src="${common_base}/icons/spacer.gif" width="1" height="12">
    </td>
    
    
    <td width="144" align="left" valign="top">
      <table width="144" border="0" cellpadding="0" cellspacing="0">
        <tr>
          <td width="8" height="12"><img src="${common_base}/icons/spacer.gif"></td>
          <td width="136"><a href="http://www.fz-juelich.de/portal/angebote"><img name="portal4" height="19" src="${common_base}/icons/portal_angebote_0.gif" alt="Angebote" border="0"></a></td>
        </tr>
        
      </table>
    </td>
    
    
    <td valign="TOP" width="9">
 <!--     <a href="/zam/index.php?path=sicherheit%2Fkapitel1&index=319&portal=1"><img src="${common_base}/icons/" width="9" height="9" alt="Portal auf-/zuklappen" title="Portal auf-/zuklappen" border="0"></a> -->
    </td>
  </tr>
  <!-- PORTALBALKEN -->
   <tr>
     <td class="tabellecyan" align="left" colspan="7" background="${common_base}/icons/portal_stripe.jpg"><font face="Arial, Helvetica, sans-serif" size="-1" class="white">&nbsp;</font></td>
     <td class="tabellecyan" valign="top" align="right">

${print_html}     

     <a href="http://www.fz-juelich.de/portal/search"><img name="search" src="${common_base}/icons/icn_search_0.gif" alt="Search" width="19" height="17" border="0"></a>
     
    ${language_html}
     </td>
     <td valign="top" width="9">&nbsp;</td>
   </tr>
   <!-- /PORTALBALKEN -->
</table>
<!-- --------------------HEADLINE PORTAL PLUGIN--------------------	-->
<!-- /Portal -->

<!-- CONTENT TABLE 1-3-1 -->
<table width="768" border="0" cellspacing="0" cellpadding="0">
	<!-- ABSTAND ROW-->
	<tr><td colspan="3" valign="top" height="3"><img src="${common_base}/icons/spacer.gif" width="1" height="3"></td></tr>
	<!-- /ABSTAND ROW-->
	
	<!-- DACHINSTITUTSNAME ROW-->
	<tr><td colspan="3" class="tabelledunkelblau" valign="top" height="38"><font face="Arial, Helvetica, sans-serif" size="+1" class="headlinewhite"><img src="${common_base}/icons/spacer.gif" width="10" height="10"><br><img src="${common_base}/icons/spacer.gif" width="10" height="10">
        ${OE}
    </font></td></tr>
	<!-- /DACHINSTITUTSNAME ROW-->

	<!-- INSTITUTSNAME ROW-->
    <tr><td colspan="3" valign="top" height="38"> <font face="Arial, Helvetica, sans-serif" size="+1" class="headline"> <img src="${common_base}/icons/spacer.gif" width="10" height="10"><br><img src="${common_base}/icons/spacer.gif" width="10" height="10">
    ${page_title}
    </font></td></tr>
	<!-- /INSTITUTSNAME ROW-->

	<!-- SEPERATION -->
	<tr><td colspan="3" height="1"><img src="${common_base}/icons/dots_hor_768_db.gif" width="768" height="1"></td></tr>
	<!-- /SEPERATION -->
	
	<!-- CONTENT ROW -->
	<tr>
    	<td width="188" align="left" valign="top">
            <!-- NAVIGATION -->
            <table width="188" cellpadding="0" border="0" cellspacing="0" class="tabellenav">
<tr>
    <td colspan="2"><img src="${common_base}/icons/spacer.gif" width="10" height="4"></td>
</tr>
<tr>
    <td colspan="2"><img src="${common_base}/icons/spacer.gif" width="10" height="3"><br>
    <img src="${common_base}/icons/dots_hor_188_db.gif" width="188" height="1"></td>
</tr>
${uplinks_html}
${menu_html}
</table>

            <!-- /NAVIGATION -->
    		<br>
        </td>
		
        <td width="10"><img src="${common_base}/icons/spacer.gif" width="10" height="10"></td>
		
        <td width="570" align="left" valign="top">
            <!-- PFAD --> <table><tr height=5><td height=5><font size="1">${navigation_text}</font></td></tr><tr height=2><td height=2><img src="${common_base}/icons/dots_hor_480_db.gif"></td></tr></table> <!-- /PFAD -->
            <!-- CONTENT BLOCK TEXT-->

<!-- Inhalt Anfang -->
EOT

$template_print = << 'EOT';
<a href="${print_url}"><img name="print" src="${common_base}/icons/icn_print_0.gif" alt="Print Version" width="21" height="17" border="0"></a>
EOT

$template_footer = << 'EOT';
<!-- Inhalt Ende -->

<!-- /CONTENT BLOCK TEXT-->
        </td>
    </tr>
    <!-- /CONTENT ROW -->
</table>
<!-- /CONTENT TABLE 1-3-1 -->

<!-- INSTITUT FOOTER TABLE -->
<table width="768" border="0" cellspacing="0" cellpadding="0">
	<tr>
        <td width="20"><img src="${common_base}/${pictures}/spacer.gif" width="20" height="15"></td>
        <td width="168"><img src="${common_base}/${pictures}/spacer.gif" width="168" height="10"></td>
        <td width="5" align="right"><img src="${common_base}/${pictures}/spacer.gif" width="5" height="10"></td>
        <td width="140"><img src="${common_base}/${pictures}/spacer.gif" width="140" height="10"></td>
        <td width="5" align="right"><img src="${common_base}/${pictures}/spacer.gif" width="5" height="10"></td>
        <td width="140"><img src="${common_base}/${pictures}/spacer.gif" width="140" height="10"></td>
        <td width="5" align="right"><img src="${common_base}/${pictures}/spacer.gif" width="5" height="10"></td>
        <td width="140"><img src="${common_base}/${pictures}/spacer.gif" width="140" height="10"></td>
        <td width="5" align="right"><img src="${common_base}/${pictures}/spacer.gif" width="5" height="10"></td>
        <td width="140"><img src="${common_base}/${pictures}/spacer.gif" width="140" height="10"></td>
	</tr>
	<tr>
        <td align="left" valign="top"><img src="${common_base}/${pictures}/copyright.gif" width="13" height="12"></td>
        <td align="left" valign="top"> <font face="Arial, Helvetica, sans-serif" size="-1">
    	    <a href="mailto:fzj@fz-juelich.de" class="footerlink">Forschungszentrum J&uuml;lich<br>D-52425 J&uuml;lich </a></font>
    	</td>
        <td valign="top"><img src="${common_base}/${pictures}/dots_vert_43.gif" width="1" height="43"></td>
        <!-- INSTITUT NAME & URL -->
        <td align="left" valign="top"> <font face="Arial, Helvetica, sans-serif" size="-1"> 
            <a href="${oe_home}" class="footerlink">${oe}</a></font><br>
            </font>
        </td>
        <!-- /INSTITUT NAME & URL -->
        <td valign="top"><img src="${common_base}/${pictures}/dots_vert_43.gif" width="1" height="43"></td>
        <!-- FORSCHUNGSBEREICH NAME & URL-->      
        <td align="left" valign="top"> 
            <font size="-1" face="Arial, Helvetica, sans-serif"> 
            <a href="http://www.fz-juelich.de/portal/forschung/information" class="footerlink"><img src="${common_base}/${pictures}/k_info_sm.gif" width="35" height="30" border="0" align="left"> Information</a></font>
        </td>
        <!-- /FORSCHUNGSBEREICH NAME & URL-->      
        <td valign="top"><img src="${common_base}/${pictures}/dots_vert_43.gif" width="1" height="43"></td>
        <td align="left" valign="top"> 
            <font face="Arial, Helvetica, sans-serif" size="-1"> 
                <a href="${language_url}" class="footerlink">${language_text}</a><br>
                <a href="http://www.fz-juelich.de/portal/impressum" class="footerlink">Impressum</a><br>
            </font>
        </td>
        <td valign="top"><img src="${common_base}/${pictures}/dots_vert_43.gif" width="1" height="43"></td>
        <td align="left" valign="top">
            <font face="Arial, Helvetica, sans-serif" size="-1">${date} <br>
            <a href="mailto:${author_email}" class="footerlink"><img src="${common_base}/${pictures}/go_mail_db.gif" width="9" height="9" border="0"> 
            ${author}</a></font>
        </td>
	</tr>
</table>
<!-- /INSTITUT FOOTER TABLE -->

</body>
</html>
<!-- cached -->
EOT

$template_uplink = << 'EOT';
<tr>
    <td>&nbsp;</td>
    <td align="left" valign="top"><font face="Arial, Helvetica, sans-serif" size="-1"> 
    <img src="${common_base}/icons/spacer.gif" width="10" height="3"><br>
    <img src="${common_base}/icons/go_up_db.gif" width="9" height="9" border="0"> <a href="${uplink_url}" class="institutslink">${uplink_text}</a></font></td>
</tr>
<tr>
    <td colspan="2"><img src="${common_base}/icons/spacer.gif" width="10" height="3"><br>
    <img src="${common_base}/icons/dots_hor_188_db.gif" width="188" height="1"></td>
</tr>
EOT

$template_menu1 = << 'EOT';
<tr>
    <td width="10"><img src="${common_base}/icons/spacer.gif" width="10" height="1"></td>
    <td align="left" valign="top"><font face="Arial, Helvetica, sans-serif" size="-1"> 
    <img src="${common_base}/icons/spacer.gif" width="160" height="3"><br>
    <img src="${common_base}/icons/go$here_db.gif" width="9" height="9" border="0">
    <a href="$menu_html" class="institutslink">$menu_text</a></font></td>
</tr>
EOT

$template_menu2 = << 'EOT';
<tr>
    <td width="10"><img src="${common_base}/icons/spacer.gif" width="10" height="1"></td>
    <td align="left" valign="top"><font face="Arial, Helvetica, sans-serif" size="-1"> 
    <img src="${common_base}/icons/spacer.gif" width="160" height="3"><br>
    <img src="${common_base}/icons/go$here_db.gif" width="9" height="9" border="0"> <a href="$menu2_html" class="institutslink">$menu2_text</a></font></td>
</tr>
EOT

$template_menu2_start = << 'EOT';
<tr>
    <td width="10"><img src="${common_base}/icons/spacer.gif" width="10" height="1"></td>
    <td align="left" valign="top"><font face="Arial, Helvetica, sans-serif" size="-1"> 
        <table width="178" cellpadding="0" border="0" cellspacing="0" class="tabellenav">
EOT

$template_menu2_stop = << 'EOT';
        </table>
</font></td></tr>
EOT

$template_pfad = << 'EOT';
 <font size="-2">></font> <a class="textlink" href="$pfad_file"><font size="-2">$pfad</font></a>
EOT

$template_print_header = << 'EOT';
<html>
    <head>
        <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-15">
        <title></title>
        <link rel="stylesheet" href="${common_base}/css/fzj.css" type="text/css">
    </head>
    <body bgcolor="#ffffff" text="#000000" link="#3366FF" vlink="#663366" alink="#3366FF">
        <font face="Arial, Helvetica, sans-serif" size=-2>
        Forschungszentrum J&uuml;lich Online - ${date}<br>

${oe} (${oe_short})<br>URL:&nbsp;
<a href="${full_url}" target="_blank">${full_url}</a> </font>
<br><br>
<!-- Inhalt Anfang -->
EOT

$template_print_footer = << 'EOT';
<!-- Inhalt Ende -->
<p><br clear=all>
<font face="Arial, Helvetica, sans-serif" size=-2>&copy; Forschungszentrum J&uuml;lich Online ${year}<br>
${oe}<br>
Alle Rechte vorbehalten<br></font> 
    </body>
</html>
EOT

1;
