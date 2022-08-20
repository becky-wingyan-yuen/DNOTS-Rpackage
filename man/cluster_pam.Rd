\name{cluster_pam}
\alias{cluster_pam}
\title{
Clustering by Partition Around Medoids (PAM) 
}
\description{
The function gives clustering results from a data frame or matrix by Partition Around Medoids (PAM) method.
}
\usage{
cluster_pam(dat, true_cluster, distance = c("euclidean", "spearman", "pearson"), 
			ncluster = 2, normalize = FALSE, 
			norm_method = c("quantile", "median", "vsn"), 
			combat = FALSE, batch = NULL, eliminate = TRUE, 
			eliminate_lv = 6 )
}
\arguments{
  \item{dat}{
  data matrix or data frame, each row corresponds to a probe and each column corresponds to a sample.
}
  \item{true_cluster}{
  numeric vector indicating the true cluster labels corresponding to the samples, prespecified by the user.
}
  \item{distance}{
  character string specifying the method to be used for calculating dissimilarities between observations in \code{pam()}. The available options are "euclidean", "spearman" and "pearson".
}
  \item{ncluster}{
  an integer indicating the number of clusters. The default value is 2.
}
  \item{normalize}{
  logical indicator. If \code{normalize = TRUE}, the data will be normalized by the chosen normalization method. The default value is \code{normalize = FALSE}.
}
  \item{norm_method}{
  character string specifying the normalization method to be used. The available options are "quantile", "median" and "vsn".
}
  \item{combat}{
  logical indicator. If \code{combat = TRUE}, \code{ComBat()} will be applied after normalization. The default value is \code{combat = FALSE}.
}
  \item{batch}{
  numeric vector indicating the batch labels for \code{ComBat()}, prespecified by the user.
}
  \item{eliminate}{
  logical indicator. If \code{eliminate = TRUE}, probes with more than 95\% expressions smaller than \code{eliminate_lv}. The default value is \code{eliminate = TRUE}.
}
  \item{eliminate_lv}{
  default value is 6.
}
}
\value{
List with the following elements:
	\item{clustering}{numeric vector indicating the estimated cluster labels corresponding to the samples.}
	\item{cluster}{\code{nclust} character string vectors storing the names of samples in each cluster.}
	\item{ARI}{corrected Rand index evaluating the accuracy of clustering results.}
}
\references{
	Yilin Wu*, Becky Wing-Yan Yuen*, Yingying Wei, Li-Xuan Qin^. On data normalization for tumor subtyping with microRNA data, Manuscript.
}
\author{
	Becky Wing-Yan Yuen
}
\examples{
data("example_data", package = "DNOTS")

uni_true_cluster = rep(NA,dim(example_data)[2])
uni_true_cluster[which(substr(colnames(example_data),7,7)=="E")] = 1
uni_true_cluster[which(substr(colnames(example_data),7,7)=="V")] = 2

results = cluster_pam(dat = example_data, true_cluster = uni_true_cluster, 
						distance = "euclidean")
						
str(results)
}