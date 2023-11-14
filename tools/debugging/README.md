# Tipps/Support for debugging

This directory collects scripts and ideas for debugging PEPC.

**Please extend the list/contents if you find something particularly helpful.**

In general, asking the compiler to enable runtime checks of the Fortran runtime
environment should point to the obvious errors and is a first step. This option
is included for debugging builds of PEPC.

If those fail, use the tools listed below.

### Valgrind

`valgrind`'s memcheck tool may prove helpful when looking for memory issues. A
full leak check may involve the following options
```sh
--tool=memcheck --leak-check=full --track-origins=yes
```
Note that especially MPI libraries take deliberate short-cuts that will show up
as memory issues and clobber output. OpenMPI provides so-called 'suppressions'
file to help with that:
```sh
--suppressions=/usr/share/openmpi/openmpi-valgrind.supp
```
(or the correct path for your installation)

Also, splitting output of `valgrind` into separate files for each spawned
process disentables output a lot. So add switches
```sh
--trace-children=yes --log-file=pepc-v.valgrind.%n.%p.`date -Iminutes`
```
to achieve just that. Here `%n` and `%p` indicate the spawned process number.
While `date -Iminutes` will add a timestamp.

Since PEPC's base64 encoding also introduces errors, we ship a suppressions file
for those
```sh
--suppressions=./tools/debugging/pepc.supp
```
It may also be helpful to generate additional suppresions. To get templates for
those, include `--gen-suppressions=all` to the command line.


### GDB

`gdb` can be batch processed to ease debugging multiple ranks or clusters. An
example 'script' can be found as `slsort_keys.gdb`. See further comments
therein. A (single rank) execution can then look like
```sh
gdb --args bin/<EXECUTABLE> <PARAMS> < slsort_keys.gdb
```
While this is close to 'printf-debugging', it's much more versatile and can make
use of all `gdb` features instead of just printing variables at runtime.


### Sanitizers

GCC includes various sanitisers to check runtime behaviour. Those are much more
light-weight than, e.g. `valgrind` and more automatic than `gdb`. They inlude
address-/leak-/thread-sanitizer. Check online for a description of those. Those
need to be added at compile and link time with some examples included in the
`makefiles.defs` example. Running with the environment variable
```sh
ASAN_OPTIONS=help=1
```
will display options for the address-sanitizer.

<!-- vim: set ts=4 sw=4 tw=80 et :-->
