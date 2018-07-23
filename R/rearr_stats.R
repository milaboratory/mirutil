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
compute_segment_usage <- function(sample,
                                  metadata = NA,
                                  ...) {
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
compute_segment_usage2 <- function(sample,
                                   metadata = NA,
                                   ...) {
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
compute_insertions <- function(sample,
                               metadata = NA,
                               ...) {
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
compute_deletions_m <- function(sample,
                                metadata = NA,
                                ...) {
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
compute_deletions <- function(sample,
                              metadata = NA,
                              ...) {
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


#
.bases <- c("A", "T", "G", "C")
.mock_count <- rep(1/4, 4)
.df_base <- data.table(expand.grid(nt.1 = .bases, nt.2 = .bases),
                       count = 0) %>%
  mutate(nt.1 = as.character(nt.1),
         nt.2 = as.character(nt.2)) %>%
  as.data.table

# Compute insertion probabilities from a table of insert sequences
# always returns 4x4=16 row table
.insert_prop_table <- function(bases,
                               reverse) {
  nbases <- length(bases)

  df.lst = list(.df_base)

  if (reverse) {
    df.lst <- c(df.lst,
                list(data.table(nt.1  = bases[1],
                           nt.2  = .bases,
                           count = .mock_count)))
  } else {
    df.lst <- c(df.lst,
                list(data.table(nt.1  = .bases,
                                nt.2  = bases[1],
                                count = .mock_count)))
  }

  if (nbases > 1) {
    if (reverse) {
      df.lst <- c(df.lst,
                  list(data.table(nt.1 = bases[2:nbases],
                                  nt.2 = bases[1:(nbases-1)],
                                  count = 1)))
    } else {
      df.lst <- c(df.lst,
                  list(data.table(nt.1 = bases[1:(nbases-1)],
                                  nt.2 = bases[2:nbases],
                                  count = 1)))
    }
  }

  df.lst %>%
    rbindlist %>%
    group_by(nt.1, nt.2) %>%
    summarise(count = sum(count))
}

# Compute insertion base probability histogram
.compute_insert_prop <- function(sample, type,
                                 insert.prob.inner.cores) {
  reverse <- F

  if (type[1] == "nVJ") {
    sample <- sample %>%
      filter(vEnd >= 0, jStart >= 0, vEnd < jStart) %>%
      mutate(insert = substr(nSeqCDR3, vEnd + 1, jStart)) %>%
      select(insert, cloneCount)
  } else if (type[1] == "nVD") {
    sample <- sample %>%
      filter(vEnd >= 0, dStart >= 0, vEnd < dStart) %>%
      mutate(insert = substr(nSeqCDR3, vEnd + 1, dStart)) %>%
      select(insert, cloneCount)
  } else if (type[1] == "nDJ") {
    reverse <- T
    sample <- sample %>%
      filter(dEnd >= 0, jStart >= 0, dEnd < jStart) %>%
      mutate(insert = substr(nSeqCDR3, dEnd + 1, jStart)) %>%
      select(insert, cloneCount)
  } else {
    stop(paste0("Unknown insert base type ", type[1]))
  }

  # just in case remove inserts with non-ATGC chars
  sample <- sample %>%
    filter(grepl("^[ATGC]+$", insert))

  # get frequencies
  with(sample,
       getDiNtFreq(insert, cloneCount, reverse)) %>%
    mutate(ins.profile.type = type[1],
           freq.clonotypes = count.clonotypes / sum(count.clonotypes),
           freq.reads = count.reads / sum(count.reads))

}


# # compute insert prob tables, all have len of 16
# res <- sample %>%
#   .$insert %>%
#   strsplit("") %>%
#   mclapply(function(bases) .insert_prop_table(bases, reverse),
#            mc.cores = insert.prob.inner.cores) %>%
#   rbindlist
#
# # append cloneCount by hand!
# res$cloneCount <- rep(sample$cloneCount, each = 16)
#
# res %>%
#   group_by(nt.1, nt.2) %>%
#   summarise(count.clonotypes = sum(count),
#             count.reads = sum(cloneCount * count)) %>%
#   ungroup %>%
#   mutate(ins.profile.type = type[1],
#          freq.clonotypes = count.clonotypes / sum(count.clonotypes),
#          freq.reads = count.reads / sum(count.reads))

compute_insertion_profile <- function(sample,
                                      metadata = NA,
                                      insert.prob.inner.cores = 1,
                                      ...) {
  chains <- .get_chains(metadata)

  types <- c()
  if (any(c("TRA", "TRG", "IGL", "IGK") %in% chains)) {
    types <- "nVJ"
  }
  if (any(c("TRB", "TRD", "IGH") %in% chains)) {
    types <- c(types, "nVD", "nDJ")
  }

  types %>%
    as.list %>%
    lapply(function(x) .compute_insert_prop(sample,
                                            x,
                                            insert.prob.inner.cores)) %>%
    rbindlist
}
