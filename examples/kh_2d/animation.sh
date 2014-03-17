#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Produces a video from fields/heatmap_XXX_??????.bin.png"
    echo "Usage: ${0} ne|ni|nefromleft|nifromleft|potential|ex|ey|vex|vey|vix|viy"
    exit 1
fi

QUALITY=6
FILEOUT="./heatmap_${1}.bin"

rm -f animation.list*

FFMPEGFILES="concat:"
FFMPEGOUT="./heatmap_"
for datafield in "$@"
do
  FILENAME="animation.list.${datafield}"
  ls -1 fields/heatmap_${datafield}_??????.bin.png | tee animation.list | sed "s/\(.*\)/file '\1'/" > ${FILENAME}
  FFMPEGFILES="${FFMPEGFILES}${FILENAME}|"
  FFMPEGOUT="${FFMPEGOUT}_${datafield}"
done

# see http://www.mplayerhq.hu/DOCS/HTML/de/menc-feat-enc-libavcodec.html

for datafield in "$@"
do
    case $QUALITY in
	0) # even higher: changed bitrate from 800k (approx default) to 3000kb/s
	    type mencoder || { echo >&2 "I require mencoder but it's not installed.  Aborting. Use JuropaGPFS nodes."; exit 1; }
	    mencoder mf://@animation.list -mf fps=16:type=png -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=3000:mbd=2:mv0:trell:v4mv:cbp:last_pred=3:predia=2:dia=2:vmax_b_frames=2:vb_strategy=1:precmp=2:cmp=2:subcmp=2:preme=2:qns=2 -oac copy -o "$FILEOUT.avi"
	    ;;
	1) # highest
	    type mencoder || { echo >&2 "I require mencoder but it's not installed.  Aborting. Use JuropaGPFS nodes."; exit 1; }
	    mencoder mf://@animation.list -mf fps=16:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:mv0:trell:v4mv:cbp:last_pred=3:predia=2:dia=2:vmax_b_frames=2:vb_strategy=1:precmp=2:cmp=2:subcmp=2:preme=2:qns=2 -oac copy -o "$FILEOUT.avi"
	    ;;
	2)
	    type mencoder || { echo >&2 "I require mencoder but it's not installed.  Aborting. Use JuropaGPFS nodes."; exit 1; }
	    mencoder mf://@animation.list -mf fps=16:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:v4mv:last_pred=2:dia=-1:vmax_b_frames=2:vb_strategy=1:cmp=3:subcmp=3:precmp=0:vqcomp=0.6:turbo -oac copy -o "$FILEOUT.avi"
	    ;;
	3)
	    type mencoder || { echo >&2 "I require mencoder but it's not installed.  Aborting. Use JuropaGPFS nodes."; exit 1; }
	    mencoder mf://@animation.list -mf fps=16:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell:v4mv:turbo -oac copy -o "$FILEOUT.avi"
	    ;;
	4) #@ lowest
	    type mencoder || { echo >&2 "I require mencoder but it's not installed.  Aborting. Use JuropaGPFS nodes."; exit 1; }
	    mencoder mf://@animation.list -mf fps=16:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:turbo -oac copy -o "$FILEOUT.avi"
	    ;;
	5) # lossless h.264
	    type ffmpeg || { echo >&2 "I require ffmpeg but it's not installed.  Aborting. Use your workstation computer."; exit 1; }
	    ffmpeg -f concat -i "${FFMPEGFILES}" -r 25 -c:v libx264 -preset veryslow -qp 0 "${FFMPEGOUT}.mkv"
           #compression# = { ultrafast | veryslow }
	    ;;
	6) # high-quality h.264
	    type ffmpeg || { echo >&2 "I require ffmpeg but it's not installed.  Aborting. Use your workstation computer."; exit 1; }
	    ffmpeg -f concat -i "${FFMPEGFILES}" -r 25 -c:v libx264 "${FFMPEGOUT}.mkv"
            #compression# = { ultrafast | veryslow }
	    ;;
	7) # high-quality, half size h.264
	    type ffmpeg || { echo >&2 "I require ffmpeg but it's not installed.  Aborting. Use your workstation computer."; exit 1; }
	    ffmpeg -f concat -i "${FFMPEGFILES}" -r 25 -c:v libx264 -preset veryslow -crf 18 -vf scale=iw/2:-1 -sws_flags bilinear "${FFMPEGOUT}.mkv"
            #compression# = { ultrafast | veryslow }
	    ;;
    esac
done
