## Various Jensen-Shannon Divergence implementations

# frequencies / zero-safe
.jsd_freq_safe <- function(p, q) {
  p <- ifelse(is.na(p), 0, p)
  q <- ifelse(is.na(q), 0, q)
  m <- 0.5 * (p + q)
  0.5 * sum(ifelse(p == 0, 0, p * log(p / m)) +
              ifelse(q == 0, 0, q * log(q / m)))
}

# frequencies / zero-unsafe
.jsd_freq_unsafe <- function(p, q) {
  m <- 0.5 * (p + q)
  0.5 * sum(p * log(p / m) + q * log(q / m))
}

# count / with pseudocounts
.jsd_count_safe <- function(x, y) {
  x <- ifelse(is.na(x), 0, x) + 1
  y <- ifelse(is.na(y), 0, y) + 1
  .jsd_freq_unsafe(x / sum(x), y / sum(y))
}

# distance based on count / with pseudocounts
.jsd_count_dist <- function(p, q) {
  sqrt(.jsd_count_safe(p, q) / log(2))
}

# distance based on frequencies
.jsd_freq_dist <- function(p, q) {
  sqrt(.jsd_freq_safe(p, q) / log(2))
}

# compute a distance between two histograms corresponding to different samples
.histogram_dist_inner <- function(data,
                                  sample.id.1,
                                  sample.id.2,
                                  .type,
                                  add.pseudocounts,
                                  filter.by.chain) {
  # select distance metric to use
  if (add.pseudocounts) {
    distfun <- .jsd_count_dist
  } else {
    distfun <- .jsd_freq_dist
  }

  # filter sample 1
  d1 <- data %>%
    mutate(value.1 = value) %>%
    filter(sample.id == sample.id.1,
           type == .type)

  if (filter.by.chain) {
    chain.1 <- (d1 %>% .$chain)[1]
  } else {
    chain.1 <- "any"
  }

  d1 <- d1 %>%
    select(variable, value.1)

  if (nrow(d1) == 0) {
    return(data.frame())
  }

  # sample for sample2
  d2 <- data %>%
    mutate(value.2 = value) %>%
    filter(sample.id == sample.id.2,
           type == .type)

  if (filter.by.chain) {
    chain.2 <- (d2 %>% .$chain)[1]
  } else {
    chain.2 <- "any"
  }

  d2 <- d2 %>%
    select(variable, value.2)

  if (nrow(d2) == 0 | chain.1 != chain.2) { # chain match check here
    return(data.frame())
  }

  # compute distance
  res <- d1 %>%
    merge(d2,
          all = T, allow.cartesian = T, by = "variable") %>%
    summarise(d = distfun(value.1, value.2)) %>%
    mutate(sample.id.1 = sample.id.1,
           sample.id.2 = sample.id.2,
           type = .type)

  if (filter.by.chain) {
    res$chain <- chain.1
  }

  res
}

#' Compute the pairwise Jensen-Shannon distance between samples for a set of histograms
#'
#' @description Computes pairwise set of Jensen-Shannon distance provided a set of
#' histograms. Several histogram groups can be provided having different 'type' value.
#' If specified, will also only compare samples with the same antigen receptor chain
#'
#' @param data a set of histograms, should contain 'sample.id', statistic type - 'type',
#' bin identifier - 'variable' and frequency/count - 'value' columns
#' @param cores number of cores for parallelization
#' @param add.pseudocounts if True, will add pseudocounts and re-normalize histograms,
#' in this case 'value' should contain count, not frequency
#' @param filter.by.chain if True, will only compare samples with the same chain,
#' requires additional 'chain' column in 'data'
#'
#' @return a pairwise distance table, containing sample ids in 'sample.id.1' and
#' 'sample.id.2' columns, statistic type in 'type' column and distance in 'd' column;
#' if 'filter.by.chain' is set to true will also include 'chain' column
histogram_dist <- function(data,
                           cores,
                           add.pseudocounts,
                           filter.by.chain) {
  if (filter.by.chain & !("chain" %in% colnames(data))) {
    stop("Cannot filter by chain with chain column missing")
  }

  with(data,
       expand.grid(
         sample.id.1 = unique(sample.id),
         sample.id.2 = unique(sample.id),
         type = unique(type),
         stringsAsFactors = F)
       ) %>%
    filter(sample.id.1 > sample.id.2) %>%
    mutate(x = paste(sample.id.1, sample.id.2, type)) %>%
    .$x %>%
    strsplit(" ") %>%
    mclapply(function(x) .histogram_dist_inner(data, x[1], x[2], x[3],
                                               add.pseudocounts,
                                               filter.by.chain),
             mc.cores = cores) %>%
    rbindlist
}

# selects between reads (weighted) and clonotypes (unweighted),
# count values (add.pseudocounts = T) and frequency values (add.pseudocounts = F)
.select_value <- function(data,
                          value.type,
                          add.pseudocounts) {
  if (value.type == "reads") {
    if (add.pseudocounts) {
      data$value <- data$count.reads
    } else {
      data$value <- data$freq.reads
    }
  } else if (value.type == "clonotypes") {
    if (add.pseudocounts) {
      data$value <- data$count.clonotypes
    } else {
      data$value <- data$freq.clonotypes
    }
  } else {
    stop(paste0("Unknown value type ", value.type))
  }
  data
}

# a generic function for computing distances for a given histogram
.dist_value <- function(data,
                        value.types,
                        add.pseudocounts,
                        cores,
                        filter.by.chain,
                        name) {
  value.types %>%
    as.list %>%
    lapply(function(value.type) {
      data %>%
        .select_value(value.type, add.pseudocounts) %>%
        histogram_dist(cores, add.pseudocounts, filter.by.chain) %>%
        mutate(value.type = value.type)
    }) %>%
    rbindlist %>%
    mutate(statistic = name)
}

#' Compute the pairwise Jensen-Shannon distance between samples for segment usage
#'
#' @description Computes pairwise set of Jensen-Shannon distance provided a set of
#' segment usage histograms. If specified, will also only compare samples with the
#' same antigen receptor chain
#'
#' @param data a set of histograms produced by 'compute_segment_usage' function
#' @param cores number of cores for parallelization
#' @param add.pseudocounts if True, will add pseudocounts and re-normalize histograms,
#' in this case 'value' should contain count, not frequency
#' @param filter.by.chain if True, will only compare samples with the same chain,
#' requires additional 'chain' column in 'data'
#'
#' @return a pairwise distance table, containing sample ids in 'sample.id.1' and
#' 'sample.id.2' columns, statistic type in 'type' column and distance in 'd' column;
#' if 'filter.by.chain' is set to true will also include 'chain' column
segment_usage_dist <- function(data,
                               value.types = c("reads", "clonotypes"),
                               add.pseudocounts = F,
                               cores = 2,
                               filter.by.chain = F) {
  data$variable <- data$segment
  data$type <- data$segment.type

  .dist_value(data,
              value.types,
              add.pseudocounts,
              cores,
              filter.by.chain,
              "segment.usage")
}

#' Compute the pairwise Jensen-Shannon distance between samples for paired segment usage
#'
#' @description Computes pairwise set of Jensen-Shannon distance provided a set of
#' segment pair usage histograms. If specified, will also only compare samples with the
#' same antigen receptor chain
#'
#' @param data a set of histograms produced by 'compute_segment_usage2' function
#' @param cores number of cores for parallelization
#' @param add.pseudocounts if True, will add pseudocounts and re-normalize histograms,
#' in this case 'value' should contain count, not frequency
#' @param filter.by.chain if True, will only compare samples with the same chain,
#' requires additional 'chain' column in 'data'
#'
#' @return a pairwise distance table, containing sample ids in 'sample.id.1' and
#' 'sample.id.2' columns, statistic type in 'type' column and distance in 'd' column;
#' if 'filter.by.chain' is set to true will also include 'chain' column
segment2_usage_dist <- function(data,
                               value.types = c("reads", "clonotypes"),
                               add.pseudocounts = F,
                               cores = 2,
                               filter.by.chain = F) {
  data$variable <- paste(data$segment.1, data$segment.2)
  data$type <- data$segment2.type

  .dist_value(data,
              value.types,
              add.pseudocounts,
              cores,
              filter.by.chain,
              "segment2.usage")
}

#' Compute the pairwise Jensen-Shannon distance between samples for insert
#' size distribution
#'
#' @description Computes pairwise set of Jensen-Shannon distance provided a set of
#' insert size distribution histograms. If specified, will also only compare samples
#' with the same antigen receptor chain
#'
#' @param data a set of histograms produced by 'compute_insertions' function
#' @param cores number of cores for parallelization
#' @param add.pseudocounts if True, will add pseudocounts and re-normalize histograms,
#' in this case 'value' should contain count, not frequency
#' @param filter.by.chain if True, will only compare samples with the same chain,
#' requires additional 'chain' column in 'data'
#'
#' @return a pairwise distance table, containing sample ids in 'sample.id.1' and
#' 'sample.id.2' columns, statistic type in 'type' column and distance in 'd' column;
#' if 'filter.by.chain' is set to true will also include 'chain' column
insert_size_dist <- function(data,
                             value.types = c("reads", "clonotypes"),
                             add.pseudocounts = F,
                             cores = 2,
                             filter.by.chain = F) {
  data$variable <- data$insertions
  data$type <- data$ins.type

  .dist_value(data,
              value.types,
              add.pseudocounts,
              cores,
              filter.by.chain,
              "insert.size")
}

#' Compute the pairwise Jensen-Shannon distance between samples for segment trimming
#' size distribution
#'
#' @description Computes pairwise set of Jensen-Shannon distance provided a set of
#' per-segment trimming size distribution histograms. If specified, will also only
#' compare samples with the same antigen receptor chain
#'
#' @param data a set of histograms produced by 'compute_deletions' function
#' @param cores number of cores for parallelization
#' @param add.pseudocounts if True, will add pseudocounts and re-normalize histograms,
#' in this case 'value' should contain count, not frequency
#' @param filter.by.chain if True, will only compare samples with the same chain,
#' requires additional 'chain' column in 'data'
#'
#' @return a pairwise distance table, containing sample ids in 'sample.id.1' and
#' 'sample.id.2' columns, statistic type in 'type' column and distance in 'd' column;
#' if 'filter.by.chain' is set to true will also include 'chain' column
deletion_size_dist <- function(data,
                               value.types = c("reads", "clonotypes"),
                               add.pseudocounts = F,
                               cores = 2,
                               filter.by.chain = F) {
  data$variable <- paste(data$segment, data$deletions)
  data$type <- data$del.type

  .dist_value(data,
              value.types,
              add.pseudocounts,
              cores,
              filter.by.chain,
              "deletion.size")
}

#' Convert pairwise distance tables to dist object
#'
#' @description Converts between distance tables and R dist object
#'
#' @param data distance table produced by '*_dist' functions
#'
#' @return an R dist object
dists_to_rdist <- function(dists) {
  mat <- dists %>%
    rbind(dists %>%
            mutate(tmp = sample.id.1,
                   sample.id.1 = sample.id.2,
                   sample.id.2 = tmp) %>%
            select(-tmp)) %>%
    dcast(sample.id.1 ~ sample.id.2, value.var = 'd', fill = 0)
  rownames(mat) <- mat$sample.id.1
  mat$sample.id.1 <- NULL
  mat %>% as.matrix %>% as.dist
}

#' Multidimensional scaling for pairwise sample distances
#'
#' @description Performs a multidimensional scaling for pairwise sample distance
#' tables produced by '*_dist' functions
#'
#' @param data distance table produced by '*_dist' functions
#'
#' @return a data table holding sample ids and their 2D coordinates
dists_to_mds <- function(dists) {
  if (nrow(dists) == 1) {
    return(data.frame(sample.id = c(dists$sample.id.1, dists$sample.id.2),
                      mds.x = c(0, 0),
                      mds.y = c(1, 1)))
  }

  xy <- dists %>%
    dists_to_rdist %>%
    isoMDS(maxit = 500) %>%
    .$points

  data.frame(sample.id = rownames(xy),
             mds.x = xy[,1],
             mds.y = xy[,2])
}
