#!/bin/bash

usage() {
cat << EOF  
Usage: ./get_region_bed.sh [-h] -i TRANSCRIPT_IDs -t ENSEMBL_GTF -p filter_gtf.py

Gets positions of 5' splice sites for each transcript.

Bedtools and pybedtools must be installed.

-h	Help
-i	Full path to Ensembl transcript ID file; specify IDs without version suffix, one per line
-t	Full path to Ensembl GTF file
-p Python script filter_gtf.py

EOF
}

while getopts ":i:t:p:" o; do
	case "${o}" in
		i)
			TRX_ID_FILE=${OPTARG}
			;;
		t)
			GTF=${OPTARG}
			;;
		p)
			PYTHON_SCRIPT=${OPTARG}
			;;
		*)
			usage
			exit
			;;
	esac
done
shift $((OPTIND-1))

echo "Started $0"

LOGFILE="$PWD/get_splice_site_stderr.log"
touch "$LOGFILE"

CURR_DIR="$PWD"
TEMP_DIR=$(mktemp -d -t temp_XXXXXXXXXX)

cd "$TEMP_DIR"

# Extract records that match the Ensembl IDs
python3 "$PYTHON_SCRIPT" -g "$GTF" -i "$TRX_ID_FILE" -o "$TEMP_DIR" -e "filt.gtf"

GTF_BASENAME="${GTF##*/}"
FILT_GTF="${GTF_BASENAME%.*}.filt.gtf"
awk '{if($3=="exon") {print}}' "$FILT_GTF" > filt_exon.gtf

# Ensure the exons are sorted from 5' to 3' 
more filt_exon.gtf | awk '{if($7=="+") {print}}' | sort -k4n > filt_pos.gtf
more filt_exon.gtf | awk '{if($7=="-") {print}}' | sort -k4nr > filt_neg.gtf
cat filt_pos.gtf filt_neg.gtf > filt_sorted.gtf

while read TRX_ID
do
	fgrep "$TRX_ID" filt_sorted.gtf > trx_id.gtf
	
	if fgrep 'transcript_biotype "protein_coding"' trx_id.gtf > /dev/null
	then
		awk -v OFS="\t" -v trx_id="$TRX_ID" '
		BEGIN {TRX_LEN=0; EXON=1}
	
		{EXON_LEN=$5-$4+1; a[EXON]=TRX_LEN+EXON_LEN; TRX_LEN+=EXON_LEN; EXON+=1}
	
		END {
		if(length(a)==1){print trx_id, "NA"} else {
		
		sep=";"; ss=a[1]; for ( i=2; i<length(a); i++) {ss = ss sep a[i]}; 
		print trx_id, ss}
		
		}' trx_id.gtf > "$TRX_ID"_splice_site.txt

		mv "$TRX_ID"_splice_site.txt "$CURR_DIR"
	else
		echo "$TRX_ID was not found in the GTF or is not protein coding and was filtered out" >> "$LOGFILE"
	fi
done < "$TRX_ID_FILE"
	
cd "$CURR_DIR"
rm -rf "$TEMP_DIR"

echo "Completed $0"

