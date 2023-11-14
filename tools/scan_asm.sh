#!/usr/bin/env bash
# Scan PEPC binary for function calls from a module, then filter the disassembly
# of those. This can (?) be used on any binary (indep. of optimisation options)
# and will print whatever it finds. 'inlining' will remove functions, so be
# careful when looking for specific functions.
# `objdump` is used to disassemble the binary. The other option is compiling w/o
# assembling and looking at generated files. This would need changes to
# makefiles and would NOT see inlining (LTO!) so will provide a different
# picture.
# Optionally, the assembly can be piped through llvm-mca to get a guesstimate on
# instructions/cycles/runtime. Interpret w/ caution.
# 

# Binary to run on
BIN=$1
if [ "X$1X" = "XX" ]; then
   echo -e "$0 <binary> [-m/--mca] [-k/--keep]"
   echo -e "  <binary>  \tpath to binary and its name to check"
   echo -e "            \te.g. \"bin/pepc-<frontend>\""
   echo -e "  -m/--mca  \trun \`llvm-mca\` simulation"
   echo -e "  -k/--keep \tkeep disassembly"
   exit
fi

# Temporaries
TMPDIR=`mktemp -d -t ASM_SCAN.XXXXXXXXXX`
ASM="$TMPDIR/asm"
ASM_="$TMPDIR/asm_"

# Check options
for arg in "$@"
do
   case "$arg" in
      -k) KEEP="1"
         ;;
      --keep) KEEP="1"
         ;;
      -m) MCA="1"
         ;;
      --mca) MCA="1"
         ;;
   esac
done

# Module name(s) too look for
# Hint: use `nm` to look at the binary and see how names are formed, e.g.
MODULES="module_interaction_specific module_mirror_boxes module_math_tools module_pepc" # 
MODULES="module_interaction_specific module_tree"
MODULES="module_interaction_specific"

echo "Scanning binary $BIN"
for MODULE in $MODULES
do
   echo "Scanning for module functions ($MODULE)"
   # Disassemble, scanning for all functions
   for FUNC in `objdump -S $BIN  | grep "<.*${MODULE}.*>:" | awk '{print $2}'`
   do
      echo -e "\n==============\n\nScanning for symbol $FUNC"
      FUNC_NAME=`echo $FUNC | sed -e 's/<//' -e 's/>://'`

      echo -e "Counting function's asm lines : $FUNC_NAME"
      # Disassemble again, restricting output only to a single function
      objdump -d $BIN | awk -v RS= "/^[[:xdigit:]]+ <${FUNC_NAME}>/" > "$ASM"
      # Those lines will look like
      #   42dbd1:       c5 31 58 8d 50 ff ff    ASM....
      # Cutting address bits
      cut -c 33- $ASM > $ASM_

      # Count instructions (non-empty lines)
      ASMLINES=`grep   "\S" $ASM_ | wc -l`
      VASMLINES=`grep  "\S" $ASM_ | grep ^v | wc -l`
      JASMLINES=`grep  "\S" $ASM_ | grep ^j | wc -l`
      SDASMLINES=`grep "\S" $ASM_ | grep ^.*sd | wc -l`
      PDASMLINES=`grep "\S" $ASM_ | grep ^.*pd | wc -l`
      CASMLINES=`grep  "\S" $ASM_ | grep ^call | wc -l`
      CALLS=`grep      "\S" $ASM_ | grep ^call | sort | uniq -c`

      echo -en "\nASM count:\n  $ASMLINES lines, "
      echo -en "out of which vector $VASMLINES, jumps $JASMLINES, "
      echo -en "calls $CASMLINES, packed double $PDASMLINES, and "
      echo -en "single double $SDASMLINES\n"

      if `test -n "${MCA}"`; then
         # Now pipe through llvm-mca to get a guesstimate of the number of cycles
         # WE IGNORE ANY ERORRS HERE, LIKELY RUBBISH
         echo -e "\nMCA simulation:"
         llvm-mca --iterations=1 $ASM_ 2>/dev/null | head -n 9
      fi

      if [ "$CASMLINES" -ne 0 ]; then
         echo -e "\nFunction calls:"
         #NO=`grep "$CALLS" $ASM_ | wc -l`
         echo -e "$CALLS"
      fi

      if `test -n "$KEEP"`; then
         cp $ASM_ "${TMPDIR}/${FUNC_NAME}.asm"
      fi
   done
done

if `test -n "$KEEP"`; then
   rm $ASM $ASM_
   echo -e "\nSeparate ASM files can be found in $TMPDIR"
else
   rm -rf $TMPDIR
fi

#  vim: set ts=4 sw=3 tw=0 et :
