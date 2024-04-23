#!/usr/bin/env bash
# Short script to generate a zip file OF THE CURRENT WORKING TREE for Zenodo.
# MAKE SURE TO RENAME THE FILE BEFORE UPLOADING.

# We use .zip over other packagers since Zenodo then provides a peak in the file

# create temporary directory
TMP=`mktemp -d`
# collect files under version control from git, exclude anything in 'benchmark'
git ls-tree -r master --name-only | grep -vE "^benchmark/" > ${TMP}/includelist
# zip it up
zip --symlinks -r ${TMP}/tmp.zip * --include @${TMP}/includelist \
   && mv ${TMP}/tmp.zip ./zenodo.zip
# create what might be the correct filename
cp ./zenodo.zip ./PEPC-`git tag --sort=taggerdate | tail -n 1`.zip
# remove temporaries
rm -f ${TMP}/includelist
rmdir ${TMP}

# vim: set ts=3 sw=3 tw=80 expandtab :
