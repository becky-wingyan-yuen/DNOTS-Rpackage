
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



###################### functions for rehybridization ######################

"estimate.handling.effect" <- function(uhdata, nuhdata){

  stopifnot(rownames(uhdata) == rownames(nuhdata))
  stopifnot(dim(uhdata) == dim(nuhdata))

  ## the list of 7 box 3 arrays and 2 bad arrays that require removal + median imputation
  ## from the rest of the arrays in the same slide
  rm.list <- c("JB4160V.b3","JB4387V.b3","JB4388V.b3",
               "JB4650V.b3","JB4757V.b3","JB4833V.b3",
               "GL3793V.b3","JB4933E","JB4952V")

  temp.handling.effect <- NULL
  for(i in 1:ncol(nuhdata)){
    temp.name <- substr(colnames(nuhdata)[i], 1, 7)
    temp.data <- uhdata[, colnames(uhdata) == temp.name]
    temp.diff <- nuhdata[, i] - temp.data
    temp.handling.effect <- cbind(temp.handling.effect, temp.diff)
  }
  colnames(temp.handling.effect) <- colnames(nuhdata)[1:ncol(nuhdata)]

  handling.effect <- NULL
  for(i in 1:24){
    begin.n <- (i-1)*8 + 1
    end.n <- i*8
    temp.data <- temp.handling.effect[, begin.n:end.n]
    temp.name <- colnames(temp.data)
    indi.vec <- temp.name %in% rm.list
    md.array <- apply(temp.data[, !(temp.name %in% rm.list)], 1, median)
    for(j in 1:8){
      if(indi.vec[j]){
        handling.effect <- cbind(handling.effect, md.array)
      } else{
        handling.effect <- cbind(handling.effect, temp.data[, j])
      }
    }
  }

  colnames(handling.effect) <- paste(substr(colnames(temp.handling.effect), 1, 7),
                             seq(1, ncol(nuhdata)), sep = ".")

  return(handling.effect)
}


"stratification.design" <- function(seed, num.array,
                                    batch.id){

  stopifnot(length(unlist(batch.id)) == num.array)
  set.seed(seed)

  g1.sample <- g2.sample <- NULL
  for(j in 1:length(batch.id)){
    sample.number <- batch.id[[j]]
    g1 <- sample(sample.number, size = length(sample.number)/2)
    g2 <- sample.number[!sample.number %in% g1]
    g1.sample <- c(g1.sample, g1)
    g2.sample <- c(g2.sample, g2)
  }
  g1.sample <- sample(g1.sample)
  g2.sample <- sample(g2.sample)

  ind <- c(g1.sample, g2.sample)
  return(ind)
}


"confounding.design" <- function(seed, num.array,
                                 degree = "complete",
                                 rev.order = FALSE){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(num.array))
  stopifnot(degree %in% c("complete", "partial"))
  stopifnot(num.array %% 2 == 0)

  set.seed(seed)
  if(degree == "complete"){
    g1 <- sample(1:(num.array/2))
    g2 <- sample((num.array/2 + 1):num.array)
    if(!rev.order){
      ind <- c(g1, g2)
    } else{
      ind <- c(g2, g1)
    }
  } else{ # partial
    swapsize <- ceiling(num.array/2/10)
    temp1.ind <- sample(1:(num.array/2), size = swapsize)
    temp2.ind <- sample((num.array/2+1):num.array, size = swapsize)
    g1 <- sample(1:(num.array/2))
    g2 <- sample(((num.array/2) + 1):num.array)
    g1[g1 %in% temp1.ind] <- temp2.ind
    g2[g2 %in% temp2.ind] <- temp1.ind
    rm(temp1.ind, temp2.ind)

    if(!rev.order){
      ind <- c(g1, g2)
    } else{
      ind <- c(g2, g1)
    }
  }
  return(ind)
}

"rehybridize" <- function (biological.effect,
                                handling.effect,
                                group.id,
                                group.id.level = c("E", "V"),
                                array.to.sample.assign,
                                icombat = FALSE,
                                isva = FALSE,
                                iruv = FALSE,
                                biological.effect.ctrl = NULL,
                                handling.effect.ctrl = NULL) {
  stopifnot(dim(biological.effect) == dim(handling.effect))
  stopifnot(rownames(biological.effect) == rownames(handling.effect))
  stopifnot(group.id %in% group.id.level)
  stopifnot(sum(icombat, isva, iruv) < 2)
  if(iruv) stopifnot(!is.null(biological.effect.ctrl) & !is.null(handling.effect.ctrl))

  halfcut <- length(array.to.sample.assign)/2
  out <- cbind(biological.effect[, group.id == group.id.level[1]] +
                 handling.effect[, array.to.sample.assign[1:halfcut]],
               biological.effect[, group.id == group.id.level[2]] +
                 handling.effect[, array.to.sample.assign[(halfcut + 1):length(array.to.sample.assign)]])

  if(icombat){ # ComBat
    cat("ComBat adjusting\n")
    mod.tr = model.matrix(~rep(1, ncol(out)))
    batch = cut(array.to.sample.assign, 8*(0:(ncol(out)/8))) # adjust by array slide
    table(batch, substr(colnames(out), 7, 7))

    combat_dat = sva::ComBat(dat = out,
                        batch = batch,
                        mod = mod.tr,
                        par.prior = TRUE)
    out <- combat_dat
  }

  if(isva){ # sva
    mod0 = model.matrix(~1, data = data.frame(colnames(out)))
    mod = model.matrix(~rep(c(0, 1), each = halfcut)) # current "out" is not in original order
    n.sv = sva::num.sv(dat = out, mod = mod, method = "leek")
    sva_dat = sva::sva(out, mod, mod0, n.sv=n.sv)
    out <- out[, colnames(biological.effect)]
    return(list(trainData = out, trainMod = mod, trainSV = sva_dat))
    stop()
  }

  if(iruv){ # ruv
    halfcut <- length(array.to.sample.assign)/2
    ctrl <- cbind(biological.effect.ctrl[, group.id == group.id.level[1]] +
                    handling.effect.ctrl[, array.to.sample.assign[1:halfcut]],
                  biological.effect.ctrl[, group.id == group.id.level[2]] +
                    handling.effect.ctrl[, array.to.sample.assign[(halfcut + 1):length(array.to.sample.assign)]])
    combine_dat <- rbind(out, ctrl)
    ctrl.ind <- rownames(combine_dat) %in% unique(rownames(ctrl))

    cat("RUV4 normalize \n")
    Y <- t(combine_dat)
    X <- rep(c(0, 1), each = halfcut) # current "combine_dat" is not in original order
    temp <- ruv::getK(Y = Y, X = matrix(X), ctl = ctrl.ind)
    xx <- ruv::RUV4(Y = Y, X = matrix(X), ctl = ctrl.ind, k = temp$k)
    ruv_dat <- Y - xx$W %*% solve(t(xx$W) %*% xx$W, t(xx$W) %*% Y) # no shrinkage
    rm(temp, xx, X, Y)
    ruv_dat <- t(ruv_dat)[!ctrl.ind, ]

    out <- ruv_dat
  }

  out <- out[, colnames(biological.effect)]
  return(out)
}

