# required metadata columns
.req_cols_meta <- c("sample.id",
                    "file.name",
                    "chain")

#' Read a set of mixcr samples
#'
#' @description Reads all samples listed in a metadata table containing
#' 'sample.id', 'file.name' and 'chain' columns:
#' 'sample.id' is a unique identifier of a sample,
#' 'file.name' contains full sample path,
#' 'chain' should specify TRA,TRB,IGH,IGK,IGL or their combination
#'
#' @param metadata a metadata table
#' @param ... extra parameters for sample reader
#'
#' @return sample bundle containing metadata and a list of sample data tables
#'
#' @export
read_mixcr_dataset <- function(metadata, ...) {
  found_cols <- intersect(.req_cols_meta, colnames(metadata))

  if (!setequal(.req_cols_meta, found_cols)) {
    missing_cols <- paste0(setdiff(found_cols, .req_cols_meta), collapse = ",")
    stop(paste0("The following required columns in metadata are missing: ",
                missing_cols))
  }

  if (length(metadata$sample.id) != length(unique(metadata$sample.id))) {
    stop("Duplicate sample.id not allowed")
  }

  samples <- metadata$file.name %>%
    as.list %>%
    lapply(function(x) read_mixcr_sample(x, ...))

  names(samples) <- metadata$sample.id

  list("metadata" = metadata,
       "samples" = samples)
}

# Auxilliary gz-friendly fread function
.fread_gz <- function(filename) {
  if (endsWith(filename, "gz")) {
    if ("command -v gzcat" %>%
        system(ignore.stdout = T, ignore.stderr = T) == 0) {
      filename <- paste("gzcat", filename)
    } else if ("command -v gzcat" %>%
               system(ignore.stdout = T, ignore.stderr = T) == 0) {
      filename <- paste("zcat", filename)
    } else {
      stop("Neither gzcat or zcat is present - cannot ungzip")
    }
  }
  fread(filename)
}

# required columns
.req_cols_mixcr <- c("refPoints",
                     "cloneId",
                     "cloneCount",
                     "allVHitsWithScore",
                     "allDHitsWithScore",
                     "allJHitsWithScore",
                     "allCHitsWithScore")

# temp columns that can be removed in further analysis
.tmp_cols_mixcr <- c("refPoints")

#' Read mixcr sample
#'
#' @description Reads mixcr clonotype table into a data table
#'
#' @param filename full path to sample
#' @param dropExtraColumns drop all columns except ones needed for this package
#'
#' @return sample data table
#'
#' @export
read_mixcr_sample <- function(filename, dropExtraColumns = F) {
  data <- .fread_gz(filename)

  # check missing cols
  found_cols <- intersect(.req_cols_mixcr, colnames(data))

  if (!setequal(.req_cols_mixcr, found_cols)) {
    missing_cols <- paste0(setdiff(found_cols, .req_cols_mixcr), collapse = ",")
    stop(paste0("The following required columns in sample ",
                sample_id,
                " are missing: ",
                missing_cols))
  }

  # remove columns unused by package if requested
  if (dropExtraColumns) {
    data <- data %>%
      select(one_of(.req_cols_mixcr))
  }

  # most likely v/d/j
  data$v         <- str_split_fixed(data$allVHitsWithScore, fixed("*"), Inf)[,1]
  data$d         <- str_split_fixed(data$allDHitsWithScore, fixed("*"), Inf)[,1]
  data$j         <- str_split_fixed(data$allJHitsWithScore, fixed("*"), Inf)[,1]

  # vdj rearrangement markup
  ref_points_tmp <- str_split_fixed(data$refPoints, fixed(":"), Inf)
  data$vDel      <- as.integer(ref_points_tmp[,11])
  data$vdIns     <- as.integer(ref_points_tmp[,13]) - as.integer(ref_points_tmp[,12])
  data$dDel5     <- as.integer(ref_points_tmp[,14])
  data$dDel3     <- as.integer(ref_points_tmp[,15])
  data$djIns     <- as.integer(ref_points_tmp[,17]) - as.integer(ref_points_tmp[,16])
  data$vjIns     <- as.integer(ref_points_tmp[,17]) - as.integer(ref_points_tmp[,12])
  data$jDel      <- as.integer(ref_points_tmp[,18])

  # non-negative insertions
  data$vdIns     <- ifelse(data$vdIns < 0, 0, data$vdIns)
  data$djIns     <- ifelse(data$djIns < 0, 0, data$djIns)
  data$vjIns     <- ifelse(data$vjIns < 0, 0, data$vjIns)

  if (dropExtraColumns) {
    data <- data %>%
      select(-one_of(.tmp_cols_mixcr))
  }

  # replace blank segment names
  data %>%
    mutate(v = ifelse(v == "", "undef", v),
           d = ifelse(d == "", "undef", d),
           j = ifelse(j == "", "undef", j))
}
