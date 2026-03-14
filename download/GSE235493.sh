#!/bin/bash

# Download script for GSE235493 data
# paper: https://doi.org/10.1016/j.neuron.2025.04.025
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235493
# code: https://doi.org/10.5281/zenodo.15243470

set -e

source "$(dirname "$0")/../functions/utils.sh"

DATA_DIR="../../data/BrainOmicsData/raw/GSE235493"

log_message "Starting GSE235493 data download..."

DOWNLOAD_LIST="
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235493&format=file|GSE235493_RAW.tar|0
"
# other files
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235493&format=file&file=GSE235493%5FOrganotypicCHD8Patchseq%2Ecsv%2Egz|GSE235493_OrganotypicCHD8Patchseq.csv.gz|0
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235493&format=file&file=GSE235493%5FAcuteRhesusPatchseq%2Ecsv%2Egz|GSE235493_AcuteRhesusPatchseq.csv.gz|0


# Perform batch download
batch_download "$DOWNLOAD_LIST" "$DATA_DIR" 5


# Extract the tar file (idempotent)
if [ ! -f "$DATA_DIR/.GSE235493_extracted" ]; then
	log_message "Extracting $DATA_DIR/GSE235493_RAW.tar ..."
	tar -xvf "$DATA_DIR/GSE235493_RAW.tar" -C "$DATA_DIR"
	# create a marker so we don't re-extract next runs
	touch "$DATA_DIR/.GSE235493_extracted"
else
	log_message "Archive already extracted (marker exists), skipping tar extraction."
fi

# Extract the csv files safely (only if .gz exists and uncompressed file missing)
for gzfile in GSE235493_OrganotypicCHD8Patchseq.csv.gz GSE235493_AcuteRhesusPatchseq.csv.gz; do
	if [ -f "$DATA_DIR/$gzfile" ]; then
		csvfile="${gzfile%.gz}"
		if [ ! -f "$DATA_DIR/$csvfile" ]; then
			log_message "Decompressing $gzfile -> $csvfile"
			gunzip -c "$DATA_DIR/$gzfile" > "$DATA_DIR/$csvfile"
			# keep the original .gz to preserve source unless you prefer removal; comment next line to keep .gz
			rm -f "$DATA_DIR/$gzfile"
		else
			log_message "$csvfile already exists, skipping decompression."
		fi
	else
		log_message "$gzfile not present in $DATA_DIR, skipping."
	fi
done

# Organize extracted files into per-sample folders (idempotent)
# Folder name is derived from filename prefix: first_two_fields joined by '_', e.g. GSM7503706_pcd87
if [ ! -f "$DATA_DIR/.GSE235493_organized" ]; then
	log_message "Organizing files into sample folders under $DATA_DIR ..."
	# Only consider regular files directly under DATA_DIR (avoid recursing into subdirs)
		find "$DATA_DIR" -maxdepth 1 -type f -print0 | while IFS= read -r -d '' f; do
		base=$(basename "$f")
		# skip archive, marker files and the two CSV metadata files we keep in DATA_DIR
		case "$base" in
				.*|GSE235493_RAW.tar|.GSE235493_extracted|.GSE235493_organized|GSE235493_OrganotypicCHD8Patchseq.csv|GSE235493_AcuteRhesusPatchseq.csv)
				continue
				;;
		esac

		# derive sample prefix: use first two underscore-separated fields if present
		IFS='_' read -r p1 p2 rest <<< "$base"
		if [ -n "$p2" ]; then
			prefix="${p1}_${p2}"
		else
			prefix="$p1"
		fi

			destdir="$DATA_DIR/$prefix"
			if [ -f "$destdir" ]; then
				log_message "Cannot organize $base: destination $destdir already exists as a file. Skipping."
				continue
			fi

			mkdir -p "$destdir"

		if [ -e "$destdir/$base" ]; then
			log_message "Already moved: $base -> $prefix/"
		else
			log_message "Moving $base -> $prefix/"
			mv -v "$f" "$destdir/"
		fi
	done

	# create organization marker
	touch "$DATA_DIR/.GSE235493_organized"
else
	log_message "Files already organized (marker exists), skipping organization step."
fi


cleanup_temp_files "$DATA_DIR"

log_message "GSE235493 data download and organization completed!" --message-type success
