
##########################################################################
################ Functions from R Package PRECISION 0.2.0 ################
################ https://github.com/rebeccahch/precision ################
##########################################################################
############################# Not to export ##############################
##########################################################################



##################### functions for reducing signals #####################

"estimate.biological.effect" <- function(uhdata) {
    biological.effect <- uhdata
    return(biological.effect)
  }


"med.sum.pbset" <- function(data, pbset.id = NULL,
                            num.per.unipbset = 10) {
  stopifnot(length(unique(table(rownames(data)))) == 1)

  if(length(unique(rownames(data))) == length(rownames(data))){
    cat("Already probe-set level\n")
    data.ps <- data
  } else {
    if(is.null(pbset.id)) pbset.id <- unique(rownames(data))

    data <- data[rownames(data) %in% pbset.id, ]
    data.ps <- apply(data, 2, function(x) tapply(x, rep(1:length(unique(rownames(data))),
                                                        each = num.per.unipbset), median))
    rownames(data.ps) <- pbset.id
  }
  return(data.ps)
}


"limma.pbset" <- function(data, group.id,
                          group.id.level = c("E", "V"),
                          pbset.id = NULL){

  stopifnot(length(unique(rownames(data))) == nrow(data))
  stopifnot(is.character(group.id))
  stopifnot(group.id %in% group.id.level)

  if(is.null(pbset.id)) pbset.id <- rownames(data)

  limma.level <- factor(group.id, levels = group.id.level)
  design.mat <- model.matrix(~0 + limma.level)
  colnames(design.mat) <- group.id.level
  cont.mat <- limma::makeContrasts(contrasts = paste0(group.id.level[2], "-", group.id.level[1]),
                                   levels = design.mat)
  fit.temp <- limma::lmFit(data, design.mat)
  contr.temp <- limma::contrasts.fit(fit.temp, cont.mat)
  eb.temp <- limma::eBayes(contr.temp)
  final.temp.1 <- limma::topTable(eb.temp,number = nrow(data))
  final.temp <- final.temp.1[match(rownames(data), rownames(final.temp.1)),]

  g1.mean <- apply(data[, group.id == group.id.level[1]], 1, mean)
  g2.mean <- apply(data[, group.id == group.id.level[2]], 1, mean)
  g1.sd <- apply(data[, group.id == group.id.level[1]], 1, sd)
  g2.sd <- apply(data[, group.id == group.id.level[2]], 1, sd)

  format.t <- data.frame(final.temp, g1.mean, g2.mean, g1.sd, g2.sd)
  name.a <- ncol(final.temp) + 1
  names(format.t)[name.a:ncol(format.t)] <- c("g1.mean","g2.mean", "g1.sd", "g2.sd")

  return(format.t)
}


"reduce.signal" <- function(biological.effect,
                                 group.id,
                                 group.id.level = c("E", "V"),
                                 reduce.multiplier = 1/2,
                                 pbset.id = NULL){
  stopifnot(nrow(biological.effect) != length(unique(rownames(biological.effect)))) # probe level
  stopifnot(length(unique(table(rownames(biological.effect)))) == 1)
  if(is.null(pbset.id)) pbset.id <- unique(rownames(biological.effect))

  n.p.u <- unique(table(rownames(biological.effect)))
  biological.effect.psl <- med.sum.pbset(biological.effect, num.per.unipbset = n.p.u)
  s.e.limma <- limma.pbset(data = biological.effect.psl,
                           group.id = group.id,
                           group.id.level = group.id.level,
                           pbset.id = pbset.id)
  de.ind <- s.e.limma$P.Value < 0.01

  sample.g1 <- rowMeans(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[1]])
  sample.g2 <- rowMeans(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[2]])

  half.signal <- (sample.g1 - sample.g2)*reduce.multiplier

  reduced.biological.effect.de <- cbind(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[1]],
                             biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[2]] + half.signal)

  # combine and colnames, rownames back to original order
  temp <- biological.effect

  temp[rownames(temp) %in% pbset.id[de.ind], ] <- reduced.biological.effect.de[, colnames(biological.effect)]
  redhalf.biological.effect.pl.p10 <- temp; rm(temp)

  return(redhalf.biological.effect.pl.p10)
}





"amplify.signal" <- function(biological.effect,
                                 group.id,
                                 group.id.level = c("E", "V"),
                                 amplify.multiplier = 2.4,
                                 pbset.id = NULL){
  stopifnot(nrow(biological.effect) != length(unique(rownames(biological.effect)))) # probe level
  stopifnot(length(unique(table(rownames(biological.effect)))) == 1)
  if(is.null(pbset.id)) pbset.id <- unique(rownames(biological.effect))

  n.p.u <- unique(table(rownames(biological.effect)))
  biological.effect.psl <- med.sum.pbset(biological.effect, num.per.unipbset = n.p.u)
  s.e.limma <- limma.pbset(data = biological.effect.psl,
                           group.id = group.id,
                           group.id.level = group.id.level,
                           pbset.id = pbset.id)
  de.ind <- s.e.limma$P.Value < 0.01

  sample.g1 <- rowMeans(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[1]])
  sample.g2 <- rowMeans(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                group.id == group.id.level[2]])
  
  amplify_ind = (sample.g2>sample.g1) + 1
  
  half.signal <- abs(sample.g1 - sample.g2)*amplify.multiplier

  reduced.biological.effect.de <- cbind(biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[1]] + half.signal*(amplify_ind==1),
                             biological.effect[rownames(biological.effect) %in% pbset.id[de.ind],
                                  group.id == group.id.level[2]] + half.signal*(amplify_ind==2))

  # combine and colnames, rownames back to original order
  temp <- biological.effect

  temp[rownames(temp) %in% pbset.id[de.ind], ] <- reduced.biological.effect.de[, colnames(biological.effect)]
  redhalf.biological.effect.pl.p10 <- temp; rm(temp)

  return(redhalf.biological.effect.pl.p10)
}



###################### functions for normalizaiton ######################

"med.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))

  if(!is.null(train)){
    # median normalization training
    ref.dis <- median(train)
    temp <- apply(train, 2, median) - ref.dis
    shifts.train <- matrix(rep(temp, each = nrow(train)), ncol = ncol(train))
    train.mn <- train - shifts.train
  } else train.mn <- NULL

  if(is.null(test)) {
    test.fmn <- NULL
  } else{
    temp <- apply(test, 2, median) - ref.dis
    shifts.test <- matrix(rep(temp, each = nrow(test)), ncol = ncol(test))
    test.fmn <- test - shifts.test
  }

  return(list("train.mn" = train.mn,
              "test.fmn" = test.fmn,
              "ref.dis" = ref.dis))
}


"quant.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))

  if(!is.null(train)){
    # quantile normalization training
    train.qn <- preprocessCore::normalize.quantiles(train, copy = TRUE)
    dimnames(train.qn) <- dimnames(train)

    ref.dis <- as.numeric(sort(train.qn[, 1]))
  } else train.qn <- NULL

  if(is.null(test)){
    test.fqn <- NULL
  } else{
    ## frozen quantile normalize test
    test.fqn <- apply(test, 2, function(x){ord <- rank(x); ref.dis[ord]})
    dimnames(test.fqn) <- dimnames(test)
  }

  return(list("train.qn" = train.qn,
              "test.fqn" = test.fqn,
              "ref.dis" = ref.dis))
}


"vs.norm" <- function(train = NULL, test = NULL, ref.dis = NULL){
  if(!is.null(train) && !is.null(test)) stopifnot(nrow(train) == nrow(test))
  if(is.null(train)) stopifnot(!is.null(ref.dis))

  if(!is.null(train)){
    # vsn training
    train0 <- 2^train
    ref.dis <- vsn::vsn2(train0, verbose = FALSE)
    train.vsn <- log2(exp(as.matrix(ref.dis)))
  } else train.vsn <- NULL

  if(is.null(test)) {
    test.fvsn <- NULL
  } else {
    test0 <- 2^test
    test.fvsn0 <- vsn::vsn2(test0, ref.dis, verbose = FALSE)
    test.fvsn <- log2(exp(as.matrix(test.fvsn0)))
  }

  return(list("train.vsn" = train.vsn,
              "test.fvsn" = test.fvsn,
              "ref.dis" = ref.dis))
}
