#!/bin/bash

usage() {
cat << EOF  
Usage: ./extract_transcript_sequences [-h] -i TRANSCRIPT_IDs -t GENCODE_GTF -g GENOME_FASTA

Extracts transcript sequences given a list of Ensembl transcript IDs. 

Bedtools must be installed.

-h	Help
-i	Full path to Ensembl transcript ID file; specify IDs without version suffix, one per line
-t	Full path to Ensembl or Gencode GTF file
-g	Full path to genome FASTA

EOF
}

while getopts ":i:t:g:" o; do
	case "${o}" in
		i)
			TRX_ID_FILE=${OPTARG}
			;;
		t)
			GTF=${OPTARG}
			;;
		g)
			GENOME=${OPTARG}
			;;
		*)
			usage
			exit
			;;
	esac
done
shift $((OPTIND-1))

echo "Started $0"

# TODO: split sequence into multiple lines so there is no limit to the transcript 
# that can be extracted
MAX_TRX_LEN=$(getconf ARG_MAX)
echo Maximum transcript length able to extract is "$MAX_TRX_LEN"

LOGFILE="$PWD/extract_transcript_sequences_stderr.log"
touch "$LOGFILE"

CURR_DIR="$PWD"
TEMP_DIR=$(mktemp -d -t temp_XXXXXXXXXX)

cd "$TEMP_DIR"

# Extract records that match the Ensembl IDs
echo "Filtering GTF. This could take some time..."
fgrep -f "$TRX_ID_FILE" "$GTF" | awk '{if($3=="exon") {print}}' - > filt.gtf

# Ensure the exons are sorted from 5' to 3' 
more filt.gtf | awk '{if($7=="+") {print}}' | sort -k4n > filt_pos.gtf
more filt.gtf | awk '{if($7=="-") {print}}' | sort -k4nr > filt_neg.gtf
cat filt_pos.gtf filt_neg.gtf > filt_sorted.gtf

while read TRX_ID
do
	fgrep "$TRX_ID" filt_sorted.gtf > trx_id.gtf
	
	if fgrep 'transcript_type "protein_coding"' trx_id.gtf > /dev/null
	then
		touch "$TRX_ID".fa
		echo \>"$TRX_ID" >> "$TRX_ID".fa
		bedtools getfasta -s -fi "$GENOME" -bed trx_id.gtf | fgrep -v \> | tr -d '\n' >> "$TRX_ID".fa
		echo >> "$TRX_ID".fa
		mv "$TRX_ID".fa "$CURR_DIR"
	else
		echo "$TRX_ID was not found in the GTF or is not protein coding and was filtered out" >> "$LOGFILE"
	fi
done < "$TRX_ID_FILE"
	
cd "$CURR_DIR"
rm -rf "$TEMP_DIR"

echo "Completed $0"
