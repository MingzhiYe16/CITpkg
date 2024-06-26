% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/C_CIT.R
\name{cit.bp.m.v2}
\alias{cit.bp.m.v2}
\title{Do CIT on binary output with multiple mediators and version 2}
\usage{
cit.bp.m.v2(
  L,
  G,
  T,
  C = NULL,
  maxit = 10000,
  n.perm = 0,
  perm.index = NULL,
  rseed = NULL,
  v2 = TRUE
)
}
\arguments{
\item{L}{A numeric vector}

\item{G}{A numeric vector}

\item{T}{A numeric vector}

\item{C}{A numeric vector}

\item{maxit}{A integer}

\item{n.perm}{A integer}

\item{perm.index}{a dataframe}

\item{rseed}{A number for seed.}

\item{v2}{A bool for whether use v2}
}
\value{
A data frame of pvalue
}
\description{
This function do CIT on binary output with multiple mediators and version 2
}
