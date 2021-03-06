% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hclust_tree.R
\name{hclust_tree}
\alias{hclust_tree}
\title{Build the hierarchical clustering tree.}
\usage{
hclust_tree(pinmat, mat_fdr, mat_dist, hcmethod = "average")
}
\arguments{
\item{pinmat}{The incidence table generated by \code{calc_pinmat}.}

\item{mat_fdr}{The FDR matrix generated by \code{fisher_fdr}}

\item{mat_dist}{The dissmilarity based on Fisher's test p-values for 
hierarchical clustering.}

\item{hcmethod}{Default: average}
}
\value{
A hclust objects with new items added.
}
\description{
Hierarchical clustering with Fisher's test p-values as distance matrix.
Also add feature coverage information for each node in the tree.
}
