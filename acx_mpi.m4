dnl @synopsis ACX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl @summary figure out how to compile/link code with MPI
dnl
dnl This macro tries to find out how to compile programs that use MPI
dnl (Message Passing Interface), a standard API for parallel process
dnl communication (see http://www-unix.mcs.anl.gov/mpi/)
dnl
dnl On success, it sets the MPICC, MPICXX, or MPIF77 output variable to
dnl the name of the MPI compiler, depending upon the current language.
dnl (This may just be $CC/$CXX/$F77, but is more often something like
dnl mpicc/mpiCC/mpif77.) It also sets MPILIBS to any libraries that are
dnl needed for linking MPI (e.g. -lmpi, if a special
dnl MPICC/MPICXX/MPIF77 was not found).
dnl
dnl If you want to compile everything with MPI, you should set:
dnl
dnl     CC="$MPICC" #OR# CXX="$MPICXX" #OR# F77="$MPIF77"
dnl     LIBS="$MPILIBS $LIBS"
dnl
dnl NOTE: The above assumes that you will use $CC (or whatever) for
dnl linking as well as for compiling. (This is the default for automake
dnl and most Makefiles.)
dnl
dnl The user can force a particular library/compiler by setting the
dnl MPICC/MPICXX/MPIF77 and/or MPILIBS environment variables.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if an MPI
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_MPI.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2005-09-02
dnl @license GPLWithACException

AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

if test -n "$MPI_DIR"; then
  MPIBIN="$MPI_DIR/bin/"
  MPILIBDIR="-L$MPI_DIR/lib/"
fi

dnl	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
if test -n "$MPIBIN"; then
	AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC, $MPIBIN)
else
	AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)
fi
	acx_mpi_save_CC="$CC"
	MPICC="$MPIBIN$MPICC"
	CC="$MPICC"
	AC_SUBST(MPICC)

dnl	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF90,[MPI Fortran compiler command])
if test -n "$MPIBIN"; then
	AC_CHECK_PROGS(MPIF90, mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpif90c, $F90, $MPIBIN)
else
	AC_CHECK_PROGS(MPIF90, mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpif90c, $F90)
fi
	acx_mpi_save_F90="$F90"
	MPIF90="$MPIBIN$MPIF90"
	F90="$MPIF90"
	AC_SUBST(MPIF90)


if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_TRY_LINK([],[      call MPI_Init], [MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"],,$MPILIBDIR)
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"],,$MPILIBDIR)
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

dnl AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
dnl 	[C++], [CXX="$acx_mpi_save_CXX"],
dnl 	[Fortran], [F90="$acx_mpi_save_F90"])
CC="$acx_mpi_save_CC"
F90="$acx_mpi_save_F90"

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI


