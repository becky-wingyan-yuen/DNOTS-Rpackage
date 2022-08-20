
######################################################################
###################### Functions for clustering ######################
######################################################################


cluster_pam = function(dat, true_cluster, distance = c("euclidean", "spearman", "pearson"), ncluster = 2, normalize = FALSE, 
	norm_method = c("quantile", "median", "vsn"), combat = FALSE, batch = NULL, eliminate = TRUE, eliminate_lv = 6 ){
	
	if( is.vector(true_cluster) == FALSE | is.numeric(true_cluster) == FALSE ){
		stop("Please format true cluster labels as a numeric vector.")
	}
	
	if( is.data.frame(dat) == FALSE && is.matrix(dat) == FALSE ){
		stop("Please format data as a data frame or matrix.")
	}else if( is.data.frame(dat) ){
		dat = as.matrix(dat)
	}
	
	#### normalization ####
	
	if( normalize == FALSE ){
		
		uhdata_averaged = dat
		
	}else{
		
		if( norm_method == "quantile" ){
			uhdata_median_norm = quant.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.qn
		}else if( norm_method == "median" ){
			uhdata_median_norm = med.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.mn
		}else if( norm_method == "vsn" ){
			uhdata_median_norm = vs.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.vsn
		}
		
	}
	
	############# ComBat to eliminate batch effects #############
	
	if( combat == TRUE && is.null(batch) == FALSE ){
		uhdata_averaged_combat = sva::ComBat(uhdata_averaged, batch)
		uhdata_averaged = uhdata_averaged_combat
	}else if( combat == TRUE && is.null(batch) ){
		stop("Please provide batch information for Combat.")
	}
	
	###### eliminate probes with >95% expressions <6 ######
	
	if( eliminate == TRUE ){
		
		proportion = rowSums(uhdata_averaged<eliminate_lv)/dim(dat)[2]*100
		uhdata_eliminated_probe = rownames(uhdata_averaged)[which(proportion>95)]
		
		if( length(uhdata_eliminated_probe) > 0 ){
			uhdata_averaged = uhdata_averaged[-which(rownames(uhdata_averaged)%in%uhdata_eliminated_probe),]
		}
		
	}
	
	###### clustering ######
	
	if( distance == "euclidean" ){
		uhdata_pam = cluster::pam(t(uhdata_averaged), k = ncluster, diss = FALSE)
	}else if( distance == "spearman" ){
		uhdata_dcorr = 1 - cor(uhdata_averaged, method = "spearman")
		uhdata_pam = cluster::pam(uhdata_dcorr, k = ncluster, diss = TRUE)
	}else if( distance == "pearson" ){
		uhdata_dcorr = 1 - cor(uhdata_averaged, method = "pearson")
		uhdata_pam = cluster::pam(uhdata_dcorr, k = ncluster, diss = TRUE)
	}
	
	###### cluster evaluation ######
	
	uni_est_cluster = as.numeric(uhdata_pam$clustering)
	
	indexes = fpc::cluster.stats(d = dist(t(uhdata_averaged)), clustering = true_cluster, 
		alt.clustering = uni_est_cluster)
	
	###### extract results ######
	
	result = list()
	result$clustering = uni_est_cluster
	result$cluster = list()
	
	for( i in 1:ncluster ){
		result$cluster[[i]] = names(uhdata_pam$clustering)[which(uhdata_pam$clustering==i)]
	}
	
	result$ARI = indexes$corrected.rand
	
	return(result)
	
}





cluster_Skmeans = function(dat, true_cluster, ncluster = 2, normalize = FALSE, 
	norm_method = c("quantile", "median", "vsn"), combat = FALSE, batch = NULL, eliminate = TRUE, eliminate_lv = 6 ){
	
	if( is.vector(true_cluster) == FALSE | is.numeric(true_cluster) == FALSE ){
		stop("Please format true cluster labels as a numeric vector.")
	}
	
	if( is.data.frame(dat) == FALSE && is.matrix(dat) == FALSE ){
		stop("Please format data as a data frame or matrix.")
	}else if( is.data.frame(dat) ){
		dat = as.matrix(dat)
	}
	
	#### normalization ####
	
	if( normalize == FALSE ){
		
		uhdata_averaged = dat
		
	}else{
		
		if( norm_method == "quantile" ){
			uhdata_median_norm = quant.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.qn
		}else if( norm_method == "median" ){
			uhdata_median_norm = med.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.mn
		}else if( norm_method == "vsn" ){
			uhdata_median_norm = vs.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.vsn
		}
		
	}
	
	############# ComBat to eliminate batch effects #############
	
	if( combat == TRUE && is.null(batch) == FALSE ){
		uhdata_averaged_combat = sva::ComBat(uhdata_averaged, batch)
		uhdata_averaged = uhdata_averaged_combat
	}else if( combat == TRUE && is.null(batch) ){
		stop("Please provide batch information for Combat.")
	}
	
	###### eliminate probes with >95% expressions <6 ######
	
	if( eliminate == TRUE ){
		
		proportion = rowSums(uhdata_averaged<eliminate_lv)/dim(dat)[2]*100
		uhdata_eliminated_probe = rownames(uhdata_averaged)[which(proportion>95)]
		
		if( length(uhdata_eliminated_probe) > 0 ){
			uhdata_averaged = uhdata_averaged[-which(rownames(uhdata_averaged)%in%uhdata_eliminated_probe),]
		}
		
	}
	
	###### clustering ######
	
	km.perm = sparcl::KMeansSparseCluster.permute(t(uhdata_averaged), K = ncluster, nperms=25, nvals = 10, silent = TRUE)
	uhdata_skmeans = sparcl::KMeansSparseCluster(t(uhdata_averaged), K = ncluster, wbounds = km.perm$bestw, silent = TRUE)
	
	###### cluster evaluation ######
	
	uni_est_cluster = as.numeric(uhdata_skmeans[[1]]$Cs)
	
	indexes = fpc::cluster.stats(d = dist(t(uhdata_averaged)), clustering = true_cluster, 
		alt.clustering = uni_est_cluster)
	
	###### extract results ######
	
	result = list()
	result$clustering = uni_est_cluster
	result$cluster = list()
	
	for( i in 1:ncluster ){
		result$cluster[[i]] = names(uhdata_skmeans[[1]]$Cs)[which(uhdata_skmeans[[1]]$Cs==i)]
	}
	
	result$ARI = indexes$corrected.rand
	result$tuning_par = km.perm$bestw
	
	return(result)
	
}





cluster_other = function(dat, true_cluster, clust_method = c("SOM", "kmeans", "MNM"), ncluster = 2, normalize = FALSE, 
	norm_method = c("quantile", "median", "vsn"), combat = FALSE, batch = NULL, eliminate = TRUE, eliminate_lv = 6 ){
	
	if( is.vector(true_cluster) == FALSE | is.numeric(true_cluster) == FALSE ){
		stop("Please format true cluster labels as a numeric vector.")
	}
	
	if( is.data.frame(dat) == FALSE && is.matrix(dat) == FALSE ){
		stop("Please format data as a data frame or matrix.")
	}else if( is.data.frame(dat) ){
		dat = as.matrix(dat)
	}
	
	#### normalization ####
	
	if( normalize == FALSE ){
		
		uhdata_averaged = dat
		
	}else{
		
		if( norm_method == "quantile" ){
			uhdata_median_norm = quant.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.qn
		}else if( norm_method == "median" ){
			uhdata_median_norm = med.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.mn
		}else if( norm_method == "vsn" ){
			uhdata_median_norm = vs.norm(train = dat)
			uhdata_averaged = uhdata_median_norm$train.vsn
		}
		
	}
	
	############# ComBat to eliminate batch effects #############
	
	if( combat == TRUE && is.null(batch) == FALSE ){
		uhdata_averaged_combat = sva::ComBat(uhdata_averaged, batch)
		uhdata_averaged = uhdata_averaged_combat
	}else if( combat == TRUE && is.null(batch) ){
		stop("Please provide batch information for Combat.")
	}
	
	###### eliminate probes with >95% expressions <6 ######
	
	if( eliminate == TRUE ){
		
		proportion = rowSums(uhdata_averaged<eliminate_lv)/dim(dat)[2]*100
		uhdata_eliminated_probe = rownames(uhdata_averaged)[which(proportion>95)]
		
		if( length(uhdata_eliminated_probe) > 0 ){
			uhdata_averaged = uhdata_averaged[-which(rownames(uhdata_averaged)%in%uhdata_eliminated_probe),]
		}
		
	}
	
	###### clustering ######
	
	if( clust_method == "SOM" ){
		uhdata_som = som::som(t(uhdata_averaged), xdim = ncluster, ydim = 1)
		uni_est_cluster = as.numeric(uhdata_som$visual$x+1)
	}else if( clust_method == "kmeans" ){
		uhdata_kmeans = kmeans(t(uhdata_averaged), centers = ncluster)
		uni_est_cluster = as.numeric(uhdata_kmeans$cluster)
	}else if( clust_method == "MNM" ){
		uhdata_mclust = mclust::Mclust(t(uhdata_averaged), G = 1:ncluster, verbose = FALSE)
		uni_est_cluster = as.numeric(uhdata_mclust$classification)
	}
	
	###### cluster evaluation ######
	
	indexes = fpc::cluster.stats(d = dist(t(uhdata_averaged)), clustering = true_cluster, 
		alt.clustering = uni_est_cluster)
	
	###### extract results ######
	
	result = list()
	result$clustering = uni_est_cluster
	result$cluster = list()
	
	if( clust_method == "SOM" ){
		for( i in 1:ncluster ){
			result$cluster[[i]] = names(uhdata_som$data[,1])[which(uhdata_som$visual$x+1==i)]
		}
	}else if( clust_method == "kmeans" ){
		for( i in 1:ncluster ){
			result$cluster[[i]] = names(uhdata_kmeans$cluster)[which(uhdata_kmeans$cluster==i)]
		}
	}else if( clust_method == "MNM" ){
		for( i in 1:ncluster ){
			result$cluster[[i]] = names(uhdata_mclust$classification)[which(uhdata_mclust$classification==i)]
		}
	}
	
	result$ARI = indexes$corrected.rand
	
	return(result)
	
}

