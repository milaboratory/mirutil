#' Apply an analysis function to all samples in a dataset
#'
#' @description Applies an analysis function to each sample and binds together
#' resulting statistics
#'
#' @param dataset a metadata table
#' @param fun a function that takes a sample and metadata list as first and
#' second arguments
#' @param cores number of cores for parallelization
#'
#' @return a summary data table of analysis results
batch_analysis <- function(dataset,
                           fun,
                           cores = 2) {
  dataset$metadata %>%
    .to_rowlist %>%
    mclapply(function(x)
      fun(dataset$samples[[x$sample.id]], metadata = x) %>%
        mutate(sample.id = x$sample.id, chain = x$chain),
      mc.cores = cores) %>%
    rbindlist
}

.to_rowlist <- function(df) {
  df %>%
    split(1:nrow(df)) %>%
    lapply(as.list)
}

compute_rearr_stats <- function(dataset,
                                cores = 2) {
  list("segment.usage" = batch_analysis(dataset,
                                        compute_segment_usage,
                                        cores),
       "insert.size"   = batch_analysis(dataset,
                                        compute_insertions,
                                        cores),
       "deletion.size" = batch_analysis(dataset,
                                        compute_deletions,
                                        cores))
}

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
      deletion_size_dist(value.types, add.pseudocounts, cores, filter.by.chain = T)
  ) %>%
    rbindlist
}

compute_rearr_stat_mds <- function(dists_bundle) {
  dists_bundle %>%
    group_by(chain, statistic, type, value.type) %>%
    do(dists_to_mds(.))
}
