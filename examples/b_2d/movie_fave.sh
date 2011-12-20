../extract_density.py
tar cvf ../$1_ave.tar fave.*
for f in fave*png ; do echo $f; convert -quality 100 $f `basename $f png`jpg; done 
mencoder "mf://fave*.jpg" -mf fps=10 -o ../$1_ave.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 
