#!/usr/local/bin/perl

@Ffiles=`ls *.f90`;

$TRACETOOL="vampir";
#$TRACETOOL="hpm";

$outdir="./kojak";

$i=0;

open(INSTRC,"> ./$outdir/instrc.new");

foreach $file (@Ffiles) {
    chomp($file);
#    print "$file\n";
    $prefix=$file;
    $prefix=~s/\.F//s;
    $contents=`cat $file`;
    $FFILE{$prefix}=1;
    ($changed,$changestr,$newcontents)=&mpichange($contents);
    if($changed) {
	printf("%02d: file %20s changed: %s\n",$i++,$file,$changestr);
#	system("mv $file $file.save");
	open(NEWFILE,"> ./$outdir/$file");
	print NEWFILE $newcontents;
	close(NEWFILE);
    }
}

&create_add_files_f90();

&write_configfile();

close(INSTRC);


sub mpichange {
    my($contents)=@_;
    $changed=0;
    $changestr="";

    $VTLOG="";
    $newcontents=&inst_stuff($contents);
    $changestr.=" [VAMPIR: $VTLOG]";
    $changed=1;

    return($changed,$changestr,$newcontents);
}





######################################################
# instrumentiert eine Source-Datei
######################################################
sub inst_stuff {
    my($stuff)=@_;

    my($before,$name,$subr,$insubr,$anf_line,$end_line,$local_sub,%contents);

    $before="";
    
    while($stuff=~m/^\s*(PROGRAM|SUBROUTINE|FUNCTION|INTEGER FUNCTION|REAL FUNCTION) ([A-Za-z0-9_]+).*$/mi) {
	$before.=$`;
	$anf_line=$&;
	$stuff=$';
	$name=$2;

#	print "WF: $name\n";
	
	if  ( ($stuff=~m/^\s*END\s*(PROGRAM|SUBROUTINE|FUNCTION) $name$/mi) || 
	      ($stuff=~m/^\s*END\s*FUNCTION\s*$/mi) || 
	      ($stuff=~m/^\s*END\s*PROGRAM\s*$/mi) ) {
	    $insubr=$`;
	    $end_line=$&;
	    $stuff=$';
	    if($stuff=~/^\n/s) {
		$stuff=$';	    
		$end_line.=$&;
	    }

	    if ($insubr=~/^\s*(PROGRAM|SUBROUTINE|FUNCTION|INTEGER FUNCTION|REAL FUNCTION)/mi) {
		# lokale SUBROUTINE
                ####################################################
		my($lbefore,$lname,$lanf_line,$lend_line,$linsubr);
		$lbefore="";
		%contents=();
		while($insubr=~m/^\s*(SUBROUTINE|FUNCTION|INTEGER FUNCTION|REAL FUNCTION) ([A-Za-z0-9_]+).*$/mi) {
		    $lanf_line=$&;
		    $lbefore.=$`;
		    $lname=$2;
		    print "WF: insubr >$lname<\n";
		    $insubr=$';
		    if($insubr=~m/^\s*END (PROGRAM|SUBROUTINE|FUNCTION) $lname/mi) {
			$linsubr=$`;
			$lend_line=$&;
			$insubr=$';
			$insubr=~/^[^\n]*?\n/s;
			$insubr=$';	    
			$lend_line.=$&;
			$insubr.="C WFINC $lname\n";
			$contents{$lname}=$lanf_line . $linsubr . $lend_line;
		    } else {
			print "\nWARNING: no END SUBROUTINE/FUNCTION for $name\n";
		    }
		}
		$lbefore.=$insubr;
		$insubr=$lbefore;    
		$local_sub=1;
                ####################################################
	    } else {
		$local_sub=0;
	    }
	    $insubr=$anf_line . $insubr . $end_line;
	    $num_subr++;

	    printf (INSTRC "%7d %-31s %25s\n", $num_subr, $name, "$instname off");
	    
	    $VTLOG.="$name($local_sub) ";

	    # Instrumentierung der Subroutine
	    $insubr=inst_subr($insubr,$name,$num_subr);

	    
	    if ($local_sub) {
		foreach $lname (keys(%contents)) {
		    $num_subr++;
		    printf (INSTRC "%7d %-31s %25s\n", $num_subr, $lname, "$instname");
		    print "\t\t(L $lname ";
		    $contents{$lname}=inst_subr($contents{$lname},$lname,$num_subr);
		    print ")\n";
		    $insubr.=s/C WFINC $lname/$contents{$lname}/;
		    
		}
		$insubr =~ s/^C WFINC (.*)$/$contents{$1}/gmi;
	    }

	    $before.=$insubr;
	    
	} else {
	    print "\nWARNING: no END PROGRAM|SUBROUTINE/FUNCTION for $name\n";
	}
	
    }

    $before.=$stuff;
    $stuff=$before;
    
    return($stuff);
}


######################################################
# instrumentiert eine Subroutine
######################################################
sub inst_subr {
    my($contents,$name,$nr)=@_;
    my(@cont,$i,$program,$contains);


    open(SUB,"> ./Sub/$name.ori.f90");
    print SUB $contents;
    close(SUB);

    if ($contents=~m/^\s*(PROGRAM) ([A-Za-z0-9_]+).*$/mi) {
	$program=1;
    } else {
	$program=0;
    }

    if ($contents=~m/^\s*(CONTAINS)/mi) {
	$contains=1;
    } else {
	$contains=0;
    }

    $VTPLACE=0;
    $use_found=0;
    $VTFUNC{"F_$name"}=1;

    my($lline,$llline);

    @cont=split(/\n/,$contents);

    if($contents=~/^\!VTPLACE/m) {
	$VTPLACE=1;
    }

    if($contents=~/^\s*USE\s/mi) {
	$use_found=1;
    }

    $anf_found=0;
    $i=0;
    while($i<=$#cont) {
	$i++ while($cont[$i]=~/^     [&*]/m);
	if($cont[$i]!~/^\!VTPLACE/) {
	    $i++ while($cont[$i]=~/^[!cC]/m);
	}


	$line=$cont[$i];
#	printf("%03d of %03d: %s\n",$i,$#cont,$line);

#	if ($line=~/(call mpi_[a-z_]+)/) {
#	    print $line;
#	    $line=~s/(call mpi_[a-z_]+)/$1_rep/i;
#	    print $line;
#	}

	if($line=~/^\s*(SUBROUTINE|PROGRAM|FUNCTION|INTEGER FUNCTION|REAL FUNCTION)/mi) {
	    if(!$use_found) {
		$i++ while($cont[$i+1]=~/^     &/m);
		@instlines=&inst("include",$name);
		splice(@cont,$i+1,0,@instlines);
		$i+=scalar @instlines;
	    }
	}

	if($line=~/^\s*USE/mi) {
	    $i++ while($cont[$i+1]=~/^\s*USE\s/mi);
	    @instlines=&inst("include",$name);
	    splice(@cont,$i+1,0,@instlines);
	    $i+=scalar @instlines;
	}
	

        # ANFANG
	if (!$anf_found) {
	    if(!$VTPLACE) {
		$anf_found=&is_statement($line);
		if (!$anf_found) {
		    if($line=~/^\#if/) {
			$j=$i;
			$i++;
			$line=$cont[$i];$anf_found=$anf_found||&is_statement($line);
			while($line!~/^\#endif/) {
			    $i++;
			    $line=$cont[$i];$anf_found=$anf_found||&is_statement($line);
			}
			$i=$j if($anf_found);
		    }
		}
	    }
	    if (!$anf_found) {
		if($line=~/^\!VTPLACE/) {
		    $anf_found=1;
		}
	    }
	    if ($anf_found) {
		$VTLOG.="Start ";
		if ($program) {
		    @instlines=&inst("program_start",$name);
		    splice(@cont,$i,0,@instlines);
		    $i+=scalar @instlines;
		} else {
		    @instlines=&inst("subroutine_start",$name);
		    splice(@cont,$i,0,@instlines);
		    $i+=scalar @instlines;
		}
	    }
	}

        # STOP
	if (($line=~/^\s*STOP/) ) {
	    $VTLOG.="STOP ";
	    @instlines=&inst("stop",$name);
	    splice(@cont,$i,0,@instlines);
	    $i+=scalar @instlines;
	}

        # RETURN
	if ( 
	    ($line=~/^\s*(\d+)?\s*RETURN/mi) 
	    ) {
	    $VTLOG.="RETURN ";
	    @instlines=&inst("return",$name);
	    splice(@cont,$i,0,@instlines);
	    $i+=scalar @instlines;
	}

        # IF RETURN
	if ( 
	    ($line=~/^(\s*IF.*?)RETURN$/mi) 
	     ) {
	    $VTLOG.="IF RETURN ";
	    $cont[$i]="$1 Then";
	    @instlines=&inst("return",$name);
	    push(@instlines,"      RETURN");
	    push(@instlines,"      END IF");
	    splice(@cont,$i+1,0,@instlines);
	    $i+=scalar @instlines;
	}
 
        # CONTAINS
	if ( ($line=~/^\s*CONTAINS/i) ) {
	    $VTLOG.="CONTAINS ";
	    if ($program) {
		@instlines=&inst("return",$name);
		splice(@cont,$i,0,@instlines);
		$i+=scalar @instlines;
	    } else {
		@instlines=&inst("return",$name);
		splice(@cont,$i,0,@instlines);
		$i+=scalar @instlines;
	    }
	}

        # END
	if ( 
	     ($line=~/^\s*END ?(SUBROUTINE|FUNCTION|PROGRAM)?/mi) 
	     and	    
	     ($line!~/^\s*END\s*(IF|DO|WHERE)/mi) 
 	    and (!$contains)
	    ) {
	    if ($program) {
		$VTLOG.="End($i PGM) ";
		@instlines=&inst("program_end",$name);
		splice(@cont,$i,0,@instlines);
		$i+=scalar @instlines;
	    } else {
		$VTLOG.="End($i SUB) ";
		@instlines=&inst("subroutine_end",$name);
		splice(@cont,$i,0,@instlines);
		$i+=scalar @instlines;
	    }
	}

	$i++;
    }
    if (!$anf_found) {
			$VTLOG.="\n\t\tWARNING: ANFANG NICHT GEFUNDEN $name\n";
    }

    $contents=join("\n",@cont);

    $contents.="\n";
   
    
    open(SUB,"> ./Sub/$name.f90");
    print SUB $contents;
    close(SUB);
    return($contents);

}

sub is_statement {
    my($line)=@_;
    $statement=0;
    if ( 
	 ($line!~/^\s*\!/) and
	 ($line!~/^[\*C\#]/) and
	 ($line!~/^\s*INTEGER/i) and
	 ($line!~/^\s*CHARACTER/i) and
	 ($line!~/FORMAT/) and
	 ($line!~/     [&!]/) and
	 ($line!~/^\s*REAL/i) and
	 ($line!~/^\s*WF1\s?\(/i) and
	 ($line!~/^\s*WF2\s?LOGICAL\(/i) and
	 ($line!~/^\s*DWF1\s?\(/i) and
	 ($line!~/^\s*DWF2\s?\(/i) and
	 ($line!~/^\s*USE/i)  and
	 ($line!~/^\s*IMPLICIT/i) and 
	 ($line!~/^\s*LOGICAL/i) and 
	 ($line!~/^\s*COMMON/i) and 
	 ($line!~/^\s*PARAMETER/i)
	 ) { 
	$statement=1 if (
			 ($line=~/=/)    or 
			 ($line=~/ IF /i) or
			 ($line=~/ IF\(/i) or
			 ($line=~/ CALL /i) or
			 ($line=~/ ALLOCATE /i) or
			 ($line=~/ WRITE\(/i) or
			 ($line=~/ OPEN\(/i) or
			 ($line=~/CONTINUE/i) or
			 ($line=~/RETURN/i) or
			 ($line=~/ DO /i) or
			 ($line=~/\s+END /i) 
			 );
    }
    return($statement);
}

sub create_add_files {
    open(VTSUB,"> ./$outdir/VTdefs.F");
    open(VTCOM,"> ./$outdir/VTcommon.h");
    print VTSUB "      SUBROUTINE VTdefs\n";
    print VTSUB "#include \"VT.inc\"\n";
    print VTSUB "        INCLUDE 'VTcommon.h'\n";
    print VTSUB "        INTEGER IE\n";

    print VTCOM "        COMMON /VTdefsC/ ICLASSH\n";
    print VTCOM "     & ,VTNOSCL\n";
    foreach $func (sort (keys(%VTFUNC))) {
	print VTCOM "     & , I$func\n";
    }
    foreach $func (sort (keys(%VTFUNC))) {
	print VTCOM "        INTEGER I$func\n";
    }
    print VTCOM "        INTEGER ICLASSH\n";
    print VTCOM "        INTEGER VTNOSCL\n";

    print VTSUB "        VTNOSCL=VT_NOSCL\n";
    print VTSUB "        call VTCLASSDEF('F',ICLASSH,IE)\n";
    foreach $func (keys(%VTFUNC)) {
	print VTSUB "       call VTFUNCDEF('$func',ICLASSH,I$func,IE)\n";
    }
    print VTSUB "      END\n";
    close(VTSUB);
    close(VTCOM);

}

sub create_add_files_f90 {
    open(VTSUB,"> ./$outdir/VTdefs.f90");
    open(VTCOM,"> ./$outdir/VTcommon.h");
    print VTSUB "      SUBROUTINE VTdefs\n";
    print VTSUB "        INCLUDE \"VT.inc\"\n";
    print VTSUB "        INCLUDE 'VTcommon.h'\n";
    print VTSUB "        INTEGER IE\n";

    print VTCOM "        COMMON /VTdefsC/ ICLASSH, &\n";
    print VTCOM "        VTNOSCL, &\n";
    print VTCOM "       I", join(",&\n\tI",(sort (keys(%VTFUNC)))),"\n";
    foreach $func (sort (keys(%VTFUNC))) {
	print VTCOM "        INTEGER*4 I$func\n";
    }
    print VTCOM "        INTEGER*4 ICLASSH\n";
    print VTCOM "        INTEGER*4 VTNOSCL\n";

    print VTSUB "        VTNOSCL=VT_NOSCL\n";
    print VTSUB "        call VTCLASSDEF('F',ICLASSH,IE)\n";
    foreach $func (keys(%VTFUNC)) {
	print VTSUB "       call VTFUNCDEF('$func',ICLASSH,I$func,IE)\n";
    }
    print VTSUB "      END\n";
    close(VTSUB);
    close(VTCOM);

}

sub write_configfile {
    print "export VT_CONFIG=./$outdir/VTconfig.cfg\n";
    open(VTCONF,"> ./$outdir/VTconfig.cfg");
    foreach $func (sort (keys(%VTFUNC))) {
	print VTCONF "SYMBOL $func OFF\n";
    }
    close(VTCONF);

}


sub inst {
    my($type,$name)=@_;
    my(@lines);
#    print "$type,$name\n";
    if($TRACETOOL eq "vampir") {
	return(&inst_vampir($type,$name));
    }
    else {
	return();
    }
}



sub inst_vampir {
    my($type,$name)=@_;
    my(@lines);
#    print "$type,$name\n";
    if($type eq "include") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!")
	     );
	return(@lines);
    }
    if($type eq "program_start") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!")
	     );
	return(@lines);
    }
    if($type eq "subroutine_start") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!\$POMP INST BEGIN($name)",
		     "!")
	     );
	return(@lines);
    }
    if($type eq "stop") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!")
	     );
	return(@lines);
    }
    if($type eq "return") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!\$POMP INST ALTEND($name)",
		     "!")
	     );
	return(@lines);
    }
    if($type eq "subroutine_end") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!\$POMP INST END($name)",
		     "!")
	     );
	return(@lines);
    }
    if($type eq "program_end") {
	push(@lines,(
		     "!VAMPINST $type",
		     "!")
	     );
	return(@lines);
    }
}

