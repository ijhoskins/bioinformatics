#!/bin/bash

usage() {
cat << EOF  
Usage: ./extract_transcript_sequences [-h] -b REGION_BED

Gets transcript lengths from a region BED file.

-h	Help
-b	Full path to region BED file

EOF
}

while getopts ":b:" o; do
	case "${o}" in
		b)
			REGION_BED=${OPTARG}
			;;
		*)
			usage
			exit
			;;
	esac
done
shift $((OPTIND-1))

echo "Started $0"

CURR_DIR="$PWD"

touch transcript_lengths.txt

awk -v OFS="\t" '{
if(FNR==1){LAST_TRX=$1; TRX_LEN=$3-$2; next}

if($1==LAST_TRX){TRX_LEN+=$3-$2} else {print LAST_TRX, TRX_LEN; LAST_TRX=$1; TRX_LEN=$3-$2}}

END {print LAST_TRX, TRX_LEN}
' "$REGION_BED" >> transcript_lengths.txt

echo "Completed $0"
