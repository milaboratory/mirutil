.jsd_freq_safe <- function(p, q) {
  p <- ifelse(is.na(p), 0, p)
  q <- ifelse(is.na(q), 0, q)
  m <- 0.5 * (p + q)
  0.5 * sum(ifelse(p == 0, 0, p * log(p / m)) +
              ifelse(q == 0, 0, q * log(q / m)))
}

.jsd_freq_unsafe <- function(p, q) {
  m <- 0.5 * (p + q)
  0.5 * sum(p * log(p / m) + q * log(q / m))
}

.jsd_count_safe <- function(x, y) {
  x <- ifelse(is.na(x), 0, x) + 1
  y <- ifelse(is.na(y), 0, y) + 1
  .jsd_freq_unsafe(x / sum(x), y / sum(y))
}

.jsd_count_dist <- function(p, q) {
  sqrt(.jsd_count_safe(p, q) / log(2))
}

.jsd_freq_dist <- function(p, q) {
  sqrt(.jsd_freq_safe(p, q) / log(2))
}

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

histogram_dist <- function(data,
                           cores,
                           add.pseudocounts,
                           filter.by.chain) {
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
