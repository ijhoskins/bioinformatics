#!/bin/bash

usage() {
cat << EOF  
Usage: ./get_region_bed.sh [-h] -i TRANSCRIPT_IDs -t GENCODE_GTF
Gets UTR and CDS regions as a BED file. Bedtools must be installed.

-h	Help
-i	Full path to Ensembl transcript ID file; specify IDs without version suffix, one per line
-t	Full path to Gencode GTF file

EOF
}

while getopts ":i:t:" o; do
	case "${o}" in
		i)
			TRX_ID_FILE=${OPTARG}
			;;
		t)
			GTF=${OPTARG}
			;;
		*)
			usage
			exit
			;;
	esac
done
shift $((OPTIND-1))

echo "Started $0"

LOGFILE="$PWD/region_bed_stderr.log"
touch "$LOGFILE"

CURR_DIR="$PWD"
TEMP_DIR=$(mktemp -d -t temp_XXXXXXXXXX)

cd "$TEMP_DIR"

# Extract records that match the Ensembl IDs
echo "Filtering GTF. This could take some time..."
fgrep -f "$TRX_ID_FILE" "$GTF" | awk '{if($3=="UTR" || $3=="CDS") {print}}' - > filt.gtf

# Ensure the exons are sorted from 5' to 3' 
more filt.gtf | awk '{if($7=="+") {print}}' | sort -k4n > filt_pos.gtf
more filt.gtf | awk '{if($7=="-") {print}}' | sort -k4nr > filt_neg.gtf
cat filt_pos.gtf filt_neg.gtf > filt_sorted.gtf

while read TRX_ID
do
	fgrep "$TRX_ID" filt_sorted.gtf > trx_id.gtf
	
	if fgrep 'transcript_type "protein_coding"' trx_id.gtf > /dev/null
	then
		awk -v OFS="\t" -v trx_id="$TRX_ID" '
		BEGIN {FIVEPRIME_LEN=0; THREEPRIME_LEN=0; FIVEPRIME_SEEN="False"}
	
		{if($3=="UTR" && FIVEPRIME_SEEN=="False"){FIVEPRIME_LEN+=$5-$4+1}; 
		if($3=="CDS"){CDS_LEN+=$5-$4+1; FIVEPRIME_SEEN="True"};
		if($3=="UTR" && FIVEPRIME_SEEN=="True"){THREEPRIME_LEN+=$5-$4+1}}
	
		END {
		if(FIVEPRIME_LEN==0){print trx_id, 0, CDS_LEN, "CDS", 0, "+"};
	
		if(FIVEPRIME_LEN>0){print trx_id, 0, FIVEPRIME_LEN, "UTR5", 0, "+";
		print trx_id, FIVEPRIME_LEN, FIVEPRIME_LEN+CDS_LEN, "CDS", 0, "+"};
	
		if(THREEPRIME_LEN>0){print trx_id, FIVEPRIME_LEN+CDS_LEN,FIVEPRIME_LEN+CDS_LEN+THREEPRIME_LEN, "UTR3", 0, "+"}
	
		}' trx_id.gtf > "$TRX_ID"_region.bed

		mv "$TRX_ID"_region.bed "$CURR_DIR"
	else
		echo "$TRX_ID was not found in the GTF or is not protein coding and was filtered out" >> "$LOGFILE"
	fi
done < "$TRX_ID_FILE"
	
cd "$CURR_DIR"
rm -rf "$TEMP_DIR"

echo "Completed $0"
