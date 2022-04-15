README

asotran.RData are the preprocessed data:
(3-day moving average, transformed to unit Frechet by the average ecdf value)

prcp.RData contains the original downloaded data.  

tpdm.R contains functions which are called by run2.R, can be sourced

run2.R is the file which should recreate the analysis.  After loading the libraries
(around line 60) and loading asotran.RData file and result.RData, I think you can skip until about line 100.  The estimated TPDM is the R object sigma.aso.