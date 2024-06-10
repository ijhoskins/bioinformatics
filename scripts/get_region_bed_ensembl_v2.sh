#!/bin/bash

usage() {
cat << EOF  
Usage: ./get_region_bed.sh [-h] -i TRANSCRIPT_IDs -t ENSEMBL_GTF -p filter_gtf.py

Gets UTR and CDS regions as a BED file.

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

LOGFILE="$PWD/region_bed_stderr.log"
touch "$LOGFILE"

CURR_DIR="$PWD"
TEMP_DIR=$(mktemp -d -t temp_XXXXXXXXXX)

cd "$TEMP_DIR"

# Extract records that match the Ensembl IDs
python3 "$PYTHON_SCRIPT" -g "$GTF" -i "$TRX_ID_FILE" -o "$TEMP_DIR" -e "filt.gtf"

GTF_BASENAME="${GTF##*/}"
FILT_GTF="${GTF_BASENAME%.*}.filt.gtf"
awk '{if($3=="five_prime_utr" || $3=="three_prime_utr" || $3=="CDS") {print}}' "$FILT_GTF" > filt_exon.gtf

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
		BEGIN {FIVEPRIME_LEN=0; THREEPRIME_LEN=0}
	
		{if($3=="five_prime_utr"){FIVEPRIME_LEN+=$5-$4+1}; 
		if($3=="CDS"){CDS_LEN+=$5-$4+1};
		if($3=="three_prime_utr"){THREEPRIME_LEN+=$5-$4+1}}
	
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

