# disable page-wise output (to help redirection)
set pagination off

# make breakpoint pending on future shared library
set breakpoint pending on

# break on a specified function call:
#  identify function via nm <executable>, or a prior gdb run, or from source
#  e.g. break on PMPI_Alltoallv or slsort_keys, then backtrace and show
#  arguments passed during the function call.
#  'continue' as last command ensures unsupervised (cluster) debugging and
#  will simply go on.

# BP #1
break slsort_keys
# without breakpoint number, 'commands' acts on the breakpoint set last
commands
  bt
  info args
  continue
end

# you can also inspect specific variables with every breakpoint, getting closer
# to a 'watch' but more targeted.
# BP #2
break mpi_pepckeys.c:154
commands
  print i
  print (int64_t)(*nin)
  continue
end

# mind you, the source code files get pre-processed during the build process,
# so files should be prepended with pp_
# BP 3
break pp_manipulate_particles.f90:919
commands
  print npold 
  print size(indxl)
  continue
end

# start program execution immediatly after displaying breakpoints
info breakpoints
run
