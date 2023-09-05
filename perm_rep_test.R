PERMANOVA_repeat_measures <- function(
    D, permute_within, blocks = NULL, block_data, permutations=999,
    metadata_order = c(names(permute_within), names(block_data)),
    na.rm=F) {
  
  # Make sure D is a dist object
  if (class(D) != "dist") {
    stop("D must be a dist object")
  }
  
  # Default to free permutations if blocks is not given
  if (!missing(block_data) && is.null(blocks)) {
    stop("blocks must be given if block_data is present")
  } else if (is.null(blocks)) {
    blocks <- rep(1, nrow(permute_within))
    block_data <- as.data.frame(matrix(0, nrow=1, ncol=0))
  } else if (length(unique(blocks)) == 1) {
    warning("blocks only contains one unique value")
  }
  
  # Ensure no metadata overlap between permute_within and block_data
  if (length(intersect(names(permute_within), names(block_data))) > 0) {
    stop("metadata is repeated across permute_within and block_data")
  }
  
  # Ensure that metadata_order only contains stuff in permute_within and block_data
  if(length(setdiff(metadata_order, union(names(permute_within), names(block_data)))) > 0) {
    stop("metadata_order contains metadata not in permute_within and block_data")
  }
  
  # Ensure that the data in permute_within matches that in dist
  ord <- rownames(as.matrix(D))
  if (length(ord) != nrow(permute_within) || length(blocks) != length(ord)) {
    stop("blocks, permute_within, and D are not the same size")
  }
  if (is.null(rownames(permute_within))) {
    warning("permute_within has no rownames - can't verify sample orders")
  } else if (!all(ord == rownames(permute_within))) {
    stop("rownames do not match between permute_within and D")
  }
  
  # Ensure matching between blocks and block_data
  if (any(is.na(blocks))) {
    stop("NAs are not allowed in blocks")
  }
  if (is.factor(blocks)) {
    if (any(!(levels(blocks) %in% rownames(block_data)))) {
      stop("not all block levels are contained in block_data")
    }
    # Match blocks with block_data and discard level information
    block_data <- block_data[match(levels(blocks), rownames(block_data)), , drop=F]
    blocks <- as.numeric(blocks)
  } else if (is.numeric(blocks)) {
    if (blocks < 1 || max(blocks) > nrow(block_data)) {
      stop("Numeric blocks has indices out of range")
    }
  } else if (is.character(blocks)) {
    if (is.null(rownames(block_data)) || !all(blocks %in% rownames(block_data))) {
      stop("blocks does not match the rownames of block_data")
    }
    # Transform to numeric
    blocks <- match(blocks, rownames(block_data))
  } else {
    stop("blocks must be a numeric, factor, or character vector")
  }
  
  # Error out on NA metadata rather than allowing adonis to error out with
  # a totally nonsensical error message
  na.removed <- 0
  if (any(is.na(permute_within)) || any(is.na(block_data))) {
    if (na.rm) {
      n_prerm <- length(blocks)
      
      # Remove NAs in block_data
      hasna <- (rowSums(is.na(block_data)) > 0) | (sapply(split(rowSums(is.na(permute_within)) > 0, blocks), mean) == 1)
      block_data <- block_data[!hasna,, drop=F]
      keep <- !hasna[blocks]
      blocks <- cumsum(!hasna)[blocks]
      
      blocks <- blocks[keep]
      permute_within <- permute_within[keep,, drop=F]
      D <- as.matrix(D)[keep, keep]
      # block_data is not subset, as the rows with NAs are no longer referenced in blocks
      
      # Remove NAs in permute_within
      keep <- rowSums(is.na(permute_within)) == 0
      blocks <- blocks[keep]
      permute_within <- permute_within[keep,, drop=F]
      D <- as.dist(D[keep, keep])
      
      if (length(blocks) < ncol(permute_within) + ncol(block_data)) {
        stop(sprintf("After omitting samples with NAs, the number of samples (%d) is less than the number of metadata (%d)",
                     length(blocks), ncol(permute_within) + ncol(block_data)))
      } else if (length(blocks) < n_prerm * 0.5) {
        warning(sprintf("Removed %d samples with NA metadata", n_prerm - length(blocks)))
      }
      na.removed <- n_prerm - length(blocks)
    } else {
      stop("Some metadata is NA! adonis does not support any NA in the metadata")
    }
  }
  
  # Warn on some suspicious input
  persample <- apply(permute_within, 1, function(x)is.factor(x) && !any(duplicated(x)))
  if (any(persample)) {
    warning(sprintf("%s in permute_within has one DOF per sample.", colnames(permute_within)[which(persample)[1]]))
  }
  if (length(unique(blocks)) < nrow(block_data)) {
    warning("Not all blocks have a sample associated with them. Block permutations will still be performed over the full set of blocks - if this is not desired, subset block_data to only the blocks which appear in the data.")
  }
  if (!any(duplicated(blocks))) {
    warning("blocks contains no duplicated elements")
  }
  
  library(vegan)
  library(permute)
  
  # Test statistic from non-permuted data
  mtdat <- cbind(permute_within, block_data[blocks,,drop=F])
  ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
  R2 <- ad$aov.tab$R2
  names(R2) <- rownames(ad$aov.tab)
  
  # Permutations
  nullsamples <- matrix(NA, nrow=length(R2), ncol=permutations)
  for (i in seq_len(permutations)) {
    within.i <- shuffle(nrow(permute_within), control=how(blocks=blocks))
    block.i <- sample(seq_len(nrow(block_data)))
    mtdat <- cbind(
      permute_within[within.i,,drop=F],
      block_data[block.i,,drop=F][blocks,,drop=F])
    perm.ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
    
    nullsamples[,i] <- perm.ad$aov.tab$R2
  }
  
  # For residuals, test the other direction (i.e. p-value of all covariates)
  n <- length(R2)
  R2[n-1] <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]
  
  # P value calculation similar to adonis's
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  
  P[n] <- NA    # No p-values for "Total"
  ad$aov.tab$`Pr(>F)` <- P
  
  if (na.rm) {
    ad$na.removed <- na.removed
  }
  
  return (ad)
}
