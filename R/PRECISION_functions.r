
##########################################################################
################ Functions from R Package PRECISION 0.2.0 ################
################ https://github.com/rebeccahch/precision ################
##########################################################################
############################# Not to export ##############################
##########################################################################


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
