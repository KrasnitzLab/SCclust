% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_clones.R
\name{find_clones}
\alias{find_clones}
\title{Identify nodes in a hierarchical tree which qualify as clones.
Identify 'hard' clones first, then expand them to 'soft' clones. Expansion may 
result in clone mergers.
Based on hierarchical clustering, identify the hard/soft clones.}
\usage{
find_clones(hc, fdrthresh = -2, sharemin = 0.85, nshare = 3, bymax = T,
  climbfromsize = 2, climbtoshare = 3)
}
\arguments{
\item{hc}{An hclust object with additional items generated by \code{hclust_tree}.}

\item{fdrthresh}{maximal allowed value for log10(FDR) for any pair of leaves 
in a clone node. Default: -2.}

\item{sharemin}{A feature is considered 'widely shared' if present in sharemin 
fraction of leaves in a node.Default: 0.90.}

\item{nshare}{Minimal number of 'widely shared' features in a hard clone. 
Default: 3.}

\item{bymax}{Logical. If TRUE (default), use maximal, and otherwise mean, 
FDR for the node as a criterion for a hard clone.}

\item{climbfromsize}{An integer: minimal size of a hard clone allowed to be 
expanded}

\item{climbtoshare}{An integer: expand the clone as long as the number of 
widely shared features is at least this value}
}
\value{
An hclust object, with hard/soft clones indicated
}
\description{
Identify nodes in a hierarchical tree which qualify as clones.
Identify 'hard' clones first, then expand them to 'soft' clones. Expansion may 
result in clone mergers.
Based on hierarchical clustering, identify the hard/soft clones.
}
