# all possible chains
.chain.list <- c("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")

# list of chains from corresponding metadata column / all chains if missing
.get_chains <- function(metadata) {
  ifelse(is.null(metadata) | is.na(metadata),
         .chain.list,
         metadata$chain %>% strsplit(","))
}

# Compute 1D segment usage by type
.compute_segment_usage <- function(sample, type) {
  if (type[1] == "V") {
    sample$segment <- sample$allVHitsWithScore
  } else if (type[1] == "D") {
    sample$segment <- sample$allDHitsWithScore
  } else if (type[1] == "J") {
    sample$segment <- sample$allJHitsWithScore
  } else if (type[1] == "C") {
    sample$segment <- sample$allCHitsWithScore
  } else {
    stop(paste0("Unknown segment type ", type[1]))
  }

  sample$segment <- with(sample,
                         ifelse(segment == "",
                                "undef*00(1)",
                                segment))

  sample %>%
    # remove all extra data
    select(cloneId, cloneCount, segment) %>%
    # unnest
    separate_rows(segment, sep = ",") %>%
    # format values
    mutate(weight = str_split_fixed(segment, "[()]", n = 3)[,2] %>%
             as.numeric,
           segment = str_split_fixed(segment, fixed("*"), n = 2)[,1]) %>%
    # normalize weight
    group_by(cloneId) %>%
    mutate(weight = weight / sum(weight)) %>%
    # compute weighted read/clonotype count/frequency
    group_by(segment, segment.type = type) %>%
    summarise(count.clonotypes = sum(weight),
              count.reads = sum(weight * cloneCount)) %>%
    ungroup %>%
    mutate(freq.clonotypes = count.clonotypes / sum(count.clonotypes),
           freq.reads = count.reads / sum(count.reads))
}

#' Compute 1D segment usage
#'
#' @description Computes the distribution of segment usage frequencies
#' for a list of segment types: 'V', 'D', 'J' or 'C'. By default, all
#' segments are analyzed. In case metadata is provided, applicable segment
#' types are deduced from the list of chains present in data
#'
#' @param sample clonotype table, imported using 'read_mixcr_dataset' or
#' 'read_mixcr_sample' functions
#' @param metadata metadata, a named list with chain identifiers
#' stored under the 'chain' entry - a comma-separated list of chains such as
#' 'TRA', 'TRB', etc is expected. If set to NA will assume all chains are present
#'
#' @return segment usage table with unweighted 'freq.clonotypes' and
#' weighted 'freq.reads' frequencies, unnormalized 'count.*' values
#' are also provided
#'
#' @export
compute_segment_usage <- function(sample,
                                  metadata = NA) {
  chains <- .get_chains(metadata)

  types <- c("V", "J")
  if (any(c("TRB", "TRD", "IGH") %in% chains)) {
    types <- c(types, "D")
  }
  if (any(c("IGH", "IGK", "IGL") %in% chains)) {
    types <- c(types, "C")
  }

  types %>%
    as.list %>%
    lapply(function(x) .compute_segment_usage(sample, x)) %>%
    rbindlist
}

# Compute 2D segment usage by type
.compute_segment_usage2 <- function(sample, type1, type2) {
  type = paste0(type1[1], type2[1])

  if (type1[1] == "V") {
    sample$segment1 <- sample$allVHitsWithScore
  } else if (type1[1] == "D") {
    sample$segment1 <- sample$allDHitsWithScore
  } else if (type1[1] == "J") {
    sample$segment1 <- sample$allJHitsWithScore
  } else if (type1[1] == "C") {
    sample$segment1 <- sample$allCHitsWithScore
  } else {
    stop(paste0("Unknown segment type ", type1[1]))
  }

  if (type2[1] == "V") {
    sample$segment2 <- sample$allVHitsWithScore
  } else if (type2[1] == "D") {
    sample$segment2 <- sample$allDHitsWithScore
  } else if (type2[1] == "J") {
    sample$segment2 <- sample$allJHitsWithScore
  } else if (type2[1] == "C") {
    sample$segment2 <- sample$allCHitsWithScore
  } else {
    stop(paste0("Unknown segment type ", type2[1]))
  }

  sample$segment1 <- with(sample,
                          ifelse(segment1 == "",
                                 "undef*00(1)",
                                 segment1))
  sample$segment2 <- with(sample,
                          ifelse(segment2 == "",
                                 "undef*00(1)",
                                 segment2))

  sample %>%
    # remove all extra data
    select(cloneId, cloneCount, segment1, segment2) %>%
    # unnest
    separate_rows(segment1, sep = ",") %>%
    separate_rows(segment2, sep = ",") %>%
    # format values
    mutate(segment1.weight = str_split_fixed(segment1, "[()]", n = 3)[,2] %>%
             as.numeric,
           segment1 = str_split_fixed(segment1, fixed("*"), n = 2)[,1],
           segment2.weight = str_split_fixed(segment2, "[()]", n = 3)[,2] %>%
             as.numeric,
           segment2 = str_split_fixed(segment2, fixed("*"), n = 2)[,1]) %>%
    # normalize weight
    group_by(cloneId) %>%
    mutate(weight = segment1.weight * segment2.weight /
             sum(segment1.weight * segment2.weight)) %>%
    # compute weighted read/clonotype count/frequency
    group_by(segment1, segment2, segment2.type = type) %>%
    summarise(count.clonotypes = sum(weight),
              count.reads = sum(weight * cloneCount)) %>%
    ungroup %>%
    mutate(freq.clonotypes = count.clonotypes / sum(count.clonotypes),
           freq.reads = count.reads / sum(count.reads))
}

#' Compute 2D segment usage
#'
#' @description Computes the distribution of segment usage frequencies
#' for a list of segment type pairs: 'VD', 'DJ' or 'VJ'. By default, all
#' segments are analyzed. In case metadata is provided, applicable segment
#' types are deduced from the list of chains present in data
#'
#' @param sample clonotype table, imported using 'read_mixcr_dataset' or
#' 'read_mixcr_sample' functions
#' @param metadata metadata, a named list with chain identifiers
#' stored under the 'chain' entry - a comma-separated list of chains such as
#' 'TRA', 'TRB', etc is expected. If set to NA will assume all chains are present
#'
#' @return segment pair usage table with unweighted 'freq.clonotypes' and
#' weighted 'freq.reads' frequencies, unnormalized 'count.*' values
#' are also provided
#'
#' @export
compute_segment_usage2 <- function(sample,
                                   metadata = NA) {
  chains <- .get_chains(metadata)

  types <- c()
  if (any(c("TRA", "TRG", "IGL", "IGK") %in% chains)) {
    types <- "VJ"
  }
  if (any(c("TRB", "TRD", "IGH") %in% chains)) {
    types <- c(types, "VD", "DJ")
  }

  types %>%
    strsplit("") %>%
    lapply(function(x) .compute_segment_usage2(sample, x[1], x[2])) %>%
    rbindlist
}

# Compute insertion size histogram
.compute_insertions <- function(sample, type) {
  if (type[1] == "VD") {
    sample$insertions <- sample$vdIns
  } else if (type[1] == "DJ") {
    sample$insertions <- sample$djIns
  } else if (type[1] == "VJ") {
    sample$insertions <- sample$vjIns
  } else {
    stop(paste0("Unknown insertions type ", type[1]))
  }

  sample %>%
    group_by(insertions, ins.type = type[1]) %>%
    summarise(count.clonotypes = n(),
              count.reads = sum(cloneCount)) %>%
    ungroup %>%
    mutate(freq.clonotypes = count.clonotypes / sum(count.clonotypes),
           freq.reads = count.reads / sum(count.reads))
}

#' Compute insertion size histogram
#'
#' @description Computes the distribution of insert sizes in various
#' CDR3 regions: 'VD', 'DJ' or 'VJ'. In case metadata is provided, applicable
#' insert regions/types are deduced from the list of chains present in data
#'
#' @param sample clonotype table, imported using 'read_mixcr_dataset' or
#' 'read_mixcr_sample' functions
#' @param metadata metadata, a named list with chain identifiers
#' stored under the 'chain' entry - a comma-separated list of chains such as
#' 'TRA', 'TRB', etc is expected. If set to NA will assume all chains are present
#'
#' @return insert size table with unweighted 'freq.clonotypes' and
#' weighted 'freq.reads' frequencies, unnormalized 'count.*' values
#' are also provided
#'
#' @export
compute_insertions <- function(sample,
                               metadata = NA) {
  chains <- .get_chains(metadata)

  types <- c()
  if (any(c("TRA", "TRG", "IGL", "IGK") %in% chains)) {
    types <- "VJ"
  }
  if (any(c("TRB", "TRD", "IGH") %in% chains)) {
    types <- c(types, "VD", "DJ")
  }

  types %>%
    as.list %>%
    lapply(function(x) .compute_insertions(sample, x)) %>%
    rbindlist
}

# Compute mean deletion sizes
.compute_deletions_m <- function(sample, type) {
  if (type[1] == "V3") {
    sample$deletions <- sample$vDel
    sample$segment   <- sample$v
  } else if (type[1] == "D5") {
    sample$deletions <- sample$dDel5
    sample$segment   <- sample$d
  } else if (type[1] == "D3") {
    sample$deletions <- sample$dDel3
    sample$segment   <- sample$d
  } else if (type[1] == "J5") {
    sample$deletions <- sample$jDel
    sample$segment   <- sample$j
  } else {
    stop(paste0("Unknown deletions type ", type[1]))
  }

  sample %>%
    group_by(segment, del.type = type[1]) %>%
    summarise(mean.clonotypes = mean(deletions),
              mean.reads = sum(deletions * cloneCount) / sum(cloneCount),
              count.clonotypes = n(),
              count.reads = sum(cloneCount))
}

#' Compute mean deletion sizes by segment
#'
#' @description Computes the means of segment trimming sizes for
#' various segments: 'V3', 'D5', 'D3' or 'J5'.  In case metadata is provided,
#' segments applicable to computing mean deletion/trimming size are deduced
#' from the list of chains present in data. Note that values for trimmings
#' are negative, however, positive values are recorded in case of P-segments
#' - polyndromic duplications.
#'
#' @param sample clonotype table, imported using 'read_mixcr_dataset' or
#' 'read_mixcr_sample' functions
#' @param metadata metadata, a named list with chain identifiers
#' stored under the 'chain' entry - a comma-separated list of chains such as
#' 'TRA', 'TRB', etc is expected. If set to NA will assume all chains are present
#'
#' @return deletion size by segment table with unweighted 'mean.clonotypes' and
#' weighted 'mean.reads' means, total number of clonotypes/reads for a given
#' segment is stored in 'count.*' columns
#'
#' @export
compute_deletions_m <- function(sample,
                              metadata = NA) {
  chains <- .get_chains(metadata)

  types <- c("V3", "J5")
  if (any(c("TRB", "TRD", "IGH") %in% chains)) {
    types <- c(types, "D5", "D3")
  }

  types %>%
    as.list %>%
    lapply(function(x) .compute_deletions_m(sample, x)) %>%
    rbindlist
}

# Compute deletion size by segment histogram
.compute_deletions <- function(sample, type) {
  if (type[1] == "V3") {
    sample$deletions <- sample$vDel
    sample$segment   <- sample$v
  } else if (type[1] == "D5") {
    sample$deletions <- sample$dDel5
    sample$segment   <- sample$d
  } else if (type[1] == "D3") {
    sample$deletions <- sample$dDel3
    sample$segment   <- sample$d
  } else if (type[1] == "J5") {
    sample$deletions <- sample$jDel
    sample$segment   <- sample$j
  } else {
    stop(paste0("Unknown deletions type ", type[1]))
  }

  sample %>%
    group_by(segment, deletions, del.type = type[1]) %>%
    summarise(count.clonotypes = n(),
              count.reads = sum(cloneCount)) %>%
    ungroup %>%
    mutate(freq.clonotypes = count.clonotypes / sum(count.clonotypes),
           freq.reads = count.reads / sum(count.reads))
}

#' Compute the distribution of deletion sizes by segment
#'
#' @description Computes the histogram of segment trimming sizes for
#' various segments: 'V3', 'D5', 'D3' or 'J5'.  In case metadata is provided,
#' segments applicable to computing mean deletion/trimming size are deduced
#' from the list of chains present in data. Note that values for trimmings
#' are negative, however, positive values are recorded in case of P-segments
#' - polyndromic duplications.
#'
#' @param sample clonotype table, imported using 'read_mixcr_dataset' or
#' 'read_mixcr_sample' functions
#' @param metadata metadata, a named list with chain identifiers
#' stored under the 'chain' entry - a comma-separated list of chains such as
#' 'TRA', 'TRB', etc is expected. If set to NA will assume all chains are present
#'
#' @return deletion size by segment 2D histogram/table with unweighted 'freq.clonotypes' and
#' weighted 'freq.reads' frequencies, unnormalized 'count.*' values
#' are also provided
#'
#' @export
compute_deletions <- function(sample,
                              metadata = NA) {
  chains <- .get_chains(metadata)

  types <- c("V3", "J5")
  if (any(c("TRB", "TRD", "IGH") %in% chains)) {
    types <- c(types, "D5", "D3")
  }

  types %>%
    as.list %>%
    lapply(function(x) .compute_deletions(sample, x)) %>%
    rbindlist
}
