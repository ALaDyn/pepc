#!/usr/bin/env bash
# Short script to generate a zip file OF THE CURRENT WORKING TREE for Zenodo.
# MAKE SURE TO RENAME THE FILE BEFORE UPLOADING.

# We use .zip over other packagers since Zenodo then provides a peak in the file

# create temporary directory
TMP=`mktemp -d`
# exclude any build data or binaries, as well as JSC benchmarks
# need to mangle directories for zip to ignore them
find bin >> ${TMP}/excludelist
find bin -type d -printf "%p/\n" >> ${TMP}/excludelist
find lib >> ${TMP}/excludelist
find lib -type d -printf "%p/\n" >> ${TMP}/excludelist
find build >> ${TMP}/excludelist
find build -type d -printf "%p/\n" >> ${TMP}/excludelist
find benchmark >> ${TMP}/excludelist
find benchmark -type d -printf "%p/\n" >> ${TMP}/excludelist
# zip it up
zip --symlinks -r ${TMP}/tmp.zip * --exclude @${TMP}/excludelist \
   && mv ${TMP}/tmp.zip ./zenodo.zip
# remove temporaries
cp ${TMP}/excludelist .
rm -f ${TMP}/excludelist
rmdir ${TMP}

# vim: set ts=3 sw=3 tw=80 expandtab :
