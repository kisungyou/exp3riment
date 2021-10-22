#' load 'Digits' data of multiple features
#' 
#' UCI >> Multiple Features Data Set
#' 
#' @usage data(Digits)
#' 
#' @format a list\describe{
#' \item{fou}{76 Fourier coefficients of the character shapes.}
#' \item{fac}{216 profile correlations.}
#' \item{kar}{64 Karhunen-Loève coefficients.}
#' \item{pix}{240 pixel averages in 2 x 3 windows.}
#' \item{zer}{47 Zernike moments.}
#' \item{mor}{6 morphological features.}
#' \item{label}{length-\eqn{2000} vector of class labels.}
#' }
#' 
#' 
"Digits"



# data generation ---------------------------------------------------------
# dat_zer = as.matrix(read.table("~/Desktop/mfeats/mfeat-zer", quote="\"", comment.char=""))
# dat_pix = as.matrix(read.table("~/Desktop/mfeats/mfeat-pix", quote="\"", comment.char=""))
# dat_mor = as.matrix(read.table("~/Desktop/mfeats/mfeat-mor", quote="\"", comment.char=""))
# dat_kar = as.matrix(read.table("~/Desktop/mfeats/mfeat-kar", quote="\"", comment.char=""))
# dat_fou = as.matrix(read.table("~/Desktop/mfeats/mfeat-fou", quote="\"", comment.char=""))
# dat_fac = as.matrix(read.table("~/Desktop/mfeats/mfeat-fac", quote="\"", comment.char=""))
# label   = rep(0:9, each=200)
# 
# 
# Digits = list()
# Digits$fou = dat_fou
# Digits$fac = dat_fac
# Digits$kar = dat_kar
# Digits$pix = dat_pix
# Digits$label = label

