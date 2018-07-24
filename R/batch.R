#' Apply an analysis function to all samples in a dataset
#'
#' @description Applies an analysis function to each sample and binds together
#' resulting statistics
#'
#' @param dataset a set of samples and metadata produced by 'read_mixcr_dataset'
#' function
#' @param fun a function that takes a sample and metadata list as first and
#' second arguments
#' @param cores number of cores for parallelization
#'
#' @return a summary data table of analysis results
batch_analysis <- function(dataset,
                           fun,
                           cores = 4) {
  dataset$metadata %>%
    .to_rowlist %>%
    mclapply(function(x) {
      #print(x$sample.id);
      fun(dataset$samples[[x$sample.id]], metadata = x) %>%
        mutate(sample.id = x$sample.id, chain = x$chain)
    }, mc.cores = cores) %>%
    rbindlist
}

# transform data frame to a list of named lists by row
.to_rowlist <- function(df) {
  df %>%
    split(1:nrow(df)) %>%
    lapply(as.list)
}

#' Computes rearrangement statistics for a dataset
#'
#' @description Produces a set of data frames holding various rearrangement statistics:
#' segment usage, insert size histograms, segment trimming histograms
#'
#' @param dataset a set of samples and metadata produced by 'read_mixcr_dataset'
#' function
#' @param cores number of cores for parallelization
#'
#' @return a list holding 'segment.usage', 'insert.size' and 'deletion.size' tables
#'
#' @export
compute_rearr_stats <- function(dataset,
                                cores = 4) {
  list("segment.usage"  = batch_analysis(dataset,
                                         compute_segment_usage,
                                         cores),
       "insert.size"    = batch_analysis(dataset,
                                         compute_insertions,
                                         cores),
       "deletion.size"  = batch_analysis(dataset,
                                         compute_deletions,
                                         cores),
       "insert.profile" = batch_analysis(dataset,
                                         compute_insertion_profile,
                                         cores)
       )
}

#' Computes distances between samples in a dataset based on pre-computed
#' rearrangement statistics
#'
#' @description Produces a data frame with pairwise distances between samples computed
#' for segment usage, insert size and segment trimming size histograms. A square
#' root of Jensen-Shannon Divergence is used as a metric. Note that no distances
#' are computed for samples having distinct antigen receptor chain lists
#'
#' @param stats_bundle a set of rearrangement statistics produced by
#' 'compute_rearr_stats' function
#' @param value.types specifies whether to use weighted histograms - 'reads' or
#' to count each clonotype once - 'clonotypes'; can specify both
#' @param add.pseudocounts if true will re-compute histograms by adding a pseudocount
#' of 1 to each histogram, including missing bins, and re-normalize them to new total
#' @param cores number of cores for parallelization
#'
#' @return a data frame holding pairwise distances between samples for each
#' rearrangement statistic
#'
#' @export
compute_rearr_stat_dist <- function(stats_bundle,
                                    value.types = c("reads", "clonotypes"),
                                    add.pseudocounts = F,
                                    cores = 2) {
  list(
    stats_bundle$segment.usage %>%
      segment_usage_dist(value.types, add.pseudocounts, cores, filter.by.chain = T),
    stats_bundle$insert.size %>%
      insert_size_dist(value.types, add.pseudocounts, cores, filter.by.chain = T),
    stats_bundle$deletion.size %>%
      deletion_size_dist(value.types, add.pseudocounts, cores, filter.by.chain = T),
    stats_bundle$insert.profile %>%
      insert_profile_dist(value.types, add.pseudocounts, cores, filter.by.chain = T)
  ) %>%
    rbindlist
}

#' Performs a multidimensional scaling placing samples in a 2D plane
#' according to distance between distributions of rearrangement statistics
#'
#' @description Produces a data frame with 2D coordinates for each sample
#' computed for segment usage, insert size and segment trimming size statistics.
#' Note that this routine will produce separate placements for each antigen
#' receptor chain present in dataset
#'
#' @param dists_bundle a set of distances produced by
#' 'compute_rearr_stat_dist' function
#' @param value.types specifies whether to use weighted histograms - 'reads' or
#' to count each clonotype once - 'clonotypes'; can specify both
#' @param add.pseudocounts if true will re-compute histograms by adding a pseudocount
#' of 1 to each histogram, including missing bins, and re-normalize them to new total
#' @param cores number of cores for parallelization
#'
#' @return a data frame holding sample coordinates for each
#' rearrangement statistic and chain
#'
#' @export
compute_rearr_stat_mds <- function(dists_bundle,
                                   cores = 2) {
  dists_bundle %>%
    select(chain, statistic, type, value.type) %>%
    unique %>%
    .to_rowlist %>%
    mclapply(function(x) {
      dists_bundle %>%
        filter(chain == x$chain,
               statistic == x$statistic,
               type == x$type,
               value.type == x$value.type) %>%
        dists_to_mds %>%
        mutate(chain = x$chain,
               statistic = x$statistic,
               type = x$type,
               value.type = x$value.type)
    }, mc.cores = cores) %>%
    rbindlist
    #group_by(chain, statistic, type, value.type) %>%
    #do(dists_to_mds(.))
}

.merge_metadata <- function(data, metadata) {
  if ("sample.id" %in% colnames(data)) {
    return(data %>% merge(metadata))
  } else if ("sample.id.1" %in% colnames(data)) {
    metadata.1 <- metadata
    colnames(metadata.1) <- paste0(colnames(metadata.1), ".1")
    metadata.2 <- metadata
    colnames(metadata.2) <- paste0(colnames(metadata.2), ".2")
    return(data %>%
             merge(metadata.1) %>%
             merge(metadata.2))
  } else {
    stop("No sample id in data frame")
  }
}

#' @export
add_metadata <- function(data, metadata) {
  if (is.data.frame(data) | is.data.table(data)) {
    return(.merge_metadata(data, metadata))
  } else if (is.list(data)) {
    for (nn in names(data)) {
      data[[nn]] <- data[[nn]] %>% .merge_metadata(metadata)
    }
    return(data)
  } else {
    stop("Unknown data format")
  }
}
