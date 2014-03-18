#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Produces a video from fields/TYPE_XXX_??????.bin.png"
    echo "Usage:"
    echo "   ${0} TYPE [FIELD [FIELD [FIELD [..]]]]"
    echo "        - FIELD=ne|ni|nefromleft|nifromleft|potential|ex|ey|vex|vey|vix|viy"
    echo "        - TYPE=heatmap|streamplot"
    echo""
    echo "If more than one field is supplied, in addition to the"
    echo "separate videos, another video is generated that shows"
    echo "all individual movies next to each other."
    exit 1
fi

QUALITY=6

OVFILES_FFMPEG=""
OVDATA_FFMPEG=""
IDX="0"
FILEMERGED=""
TYPE=${1}

shift # get rid of first argument

for datafield in "$@"
do
    echo "################### processing ${datafield} ####################"
    FILENAME_FFMPEG="animation.ffmpeg"
    FILENAME_MPLAYER="animation.mencoder"
    FILEOUT="./${TYPE}_${datafield}"
    FILEMERGED="${FILEMERGED}_${datafield}"
    OVFILES_FFMPEG="${OVFILES_FFMPEG} -i ${FILEOUT}.mkv"
    if [ "${IDX}" -eq 0 ]; then
      OVDATA_FFMPEG="[0:v]pad=iw*$#:ih"
    else
      PREVIDX=$((${IDX} - 1))
      OVDATA_FFMPEG="${OVDATA_FFMPEG}[vid${PREVIDX}];[vid${PREVIDX}][${IDX}:v]overlay=w*${IDX}"
    fi
    IDX=$((${IDX} + 1))
    ls -1 fields/${TYPE}_${datafield}_??????.bin.png | tee ${FILENAME_MPLAYER} | sed "s/\(.*\)/file '\1'/" > ${FILENAME_FFMPEG}

    case $QUALITY in
        5) # lossless h.264
            type ffmpeg || { echo >&2 "I require ffmpeg but it's not installed.  Aborting. Use your workstation computer."; exit 1; }
            ffmpeg -f concat -i "${FILENAME_FFMPE}" -r 25 -c:v libx264 -preset veryslow -qp 0 "${FILEOUT}.mkv"
            #compression# = { ultrafast | veryslow }
            ;;
        6) # high-quality h.264
            type ffmpeg || { echo >&2 "I require ffmpeg but it's not installed.  Aborting. Use your workstation computer."; exit 1; }
            ffmpeg -f concat -i "${FILENAME_FFMPEG}" -r 25 -c:v libx264 "${FILEOUT}.mkv"
            #compression# = { ultrafast | veryslow }
            ;;
        7) # high-quality, half size h.264
            type ffmpeg || { echo >&2 "I require ffmpeg but it's not installed.  Aborting. Use your workstation computer."; exit 1; }
            ffmpeg -f concat -i "${FILENAME_FFMPEG}" -r 25 -c:v libx264 -preset veryslow -crf 18 -vf scale=iw/2:-1 -sws_flags bilinear "${FILEOUT}.mkv"
            #compression# = { ultrafast | veryslow }
            ;;
    esac
done

if [ "$#" -gt 2 ]; then
  echo "################### merging videos ####################"
  echo ffmpeg ${OVFILES_FFMPEG} -filter_complex "${OVDATA_FFMPEG}" -r 25 -c:v libx264 "heatmap${FILEMERGED}.mkv"
  ffmpeg ${OVFILES_FFMPEG} -filter_complex "${OVDATA_FFMPEG}" -r 25 -c:v libx264 "heatmap${FILEMERGED}.mkv"
fi

