# Add reconciled source identities without replacing historical source donor IDs.
add_canonical_donor_identity <- function(
  meta, crosswalk_file = "data/canonical_donor_crosswalk.tsv"
) {
  required <- c("Dataset", "Global_Donor_ID")
  if (!all(required %in% names(meta))) {
    stop("Canonical donor identities require Dataset and Global_Donor_ID")
  }
  crosswalk <- read.delim(crosswalk_file, sep = "\t", quote = "",
                         stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(required, "Canonical_Donor_ID", "Identity_Status")
  if (!all(required %in% names(crosswalk)) ||
      anyNA(crosswalk[, required]) ||
      any(!nzchar(as.matrix(crosswalk[, required])))) {
    stop("Canonical donor crosswalk has missing required identities")
  }
  key <- paste(crosswalk$Dataset, crosswalk$Global_Donor_ID, sep = "\r")
  if (anyDuplicated(key)) stop("Canonical donor crosswalk has duplicate source keys")
  donor <- as.character(meta$Global_Donor_ID)
  index <- match(paste(meta$Dataset, donor, sep = "\r"), key)
  matched <- !is.na(index)
  canonical <- donor
  canonical[matched] <- crosswalk$Canonical_Donor_ID[index[matched]]
  status <- rep("source_specific_identity_retained", length(donor))
  status[is.na(donor) | donor %in% c("", "NA", "Unknown", "unknown")] <- "unknown"
  status[matched] <- crosswalk$Identity_Status[index[matched]]
  meta$Canonical_Donor_ID <- canonical
  meta$Donor_Identity_Resolution <- status
  meta
}
