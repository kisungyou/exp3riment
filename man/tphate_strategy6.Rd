% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tphate_strategy6.R
\name{tphate_strategy6}
\alias{tphate_strategy6}
\title{Transition PHATE : Strategy 6}
\usage{
tphate_strategy6(data, ndim = 2, nbdk = 5, alpha = 2, scaler = 0.5)
}
\description{
set 1 for same-cluster affinities
otherwise, use log-like decay but maximum of those as Jessie suggested
}
