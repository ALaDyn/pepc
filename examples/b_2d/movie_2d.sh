#../movie_sheath.py
for f in *png ; do echo $f; convert -quality 100 $f `basename $f png`jpg; done 
mencoder "mf://*.jpg" -mf fps=10 -o ../$1.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 
