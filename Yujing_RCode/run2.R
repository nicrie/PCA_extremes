## Prepare a data file. Just for precipitation data.


## 1. Download data from ftp.
##    The downloaded data consists of files each containing data for a
##    specific station.
##    Relevant files in the compressed file including the XX XX , where
##    longitude and latitude of each station, plus the actual working period
##    is recorded.
##    Used bash script to remove variables other than PRCP (precipitation) from each file.
##    Read in data for stations with a record length more than a specific
##    level, and combine all of PRCP data into a matrix, with each column
##    corresponding to each station, and each row corresponding to each
##    date. The records from each station have been trimmed to share the
##    same length. For entries where there is no value, NA is put in.

##   (Put code files for this step into a separate folder.)


## 2. Data preprocessing.
##    Three steps:
##       i. Moving average each time series.
##       ii. ECDF smoothing.
##       iii. Probability integral transformation.
##    i. Moving average each time series.
##       Goal is to align asynchronicous extreme precipitation at different
##    locations but are due to the same extreme event. Ignoring this fact
##    may lead to underestimation of the extremal dependence.
##    Function used: ma
##    ii. ECDF smoothing.
##        Moving average will lead to clusters in the time series, which in
##    turn will affect the tail of the PIT transformed distribution. E.g.,
##    one large precipitation may contribute multiple times to the estimate
##    of the extremal dependence.
##    Function used: ECDF.smoothing
##    iii. Probability integral transformation.
##         A routine step in multivariate extreme analysis. Regular
##    multivariate extreme analysis, at least when studying the tail
##    dependence, requires each margin to have the same distribution, which
##    is almost impossible for real data. So a transformation has to be done
##    to make them so. In OUR case, marginal series are transformed to be
##    unit Frechet distributed.
##    Function used: to.alpha.2


## 3. Estimation of Tail Pairwise Dependence Matrix (TPDM).
##    See the paper.
##    Function used: rw.Sigma.
##                   if the matrix is not positive definite yet: nearPD from
##    package (Matrix)


## 4. Principal Component Analysis of TPDM.
##    Function used: eigen


##   Following is the complete sequence of operations analysing extreme
##   prcipitation over continental US. Corresponding to the above 2 -- 4.
##   Two code files.

x1 = c(6,  10,  3, NaN, 5,  9, 11,  1, 1, 10)
x2 = c(12,  7,  1,  9,  4,  6,  NaN,  1, 11, 11)
x3 = c(1,  8,  1,  3, 12,  9, 10, 2,  7,  5)
x4 = c( 8,  NaN,  1,  3, 12,  10, 10, 2,  7,  5)
x = cbind(x1, x2, x3, x4) 

k = 5
th = 5
u = 0.80

## Compute TPDM
P <- ncol(x)                      # Number of stations
M <- nrow(x)                      # Number of obs
Sigma <- matrix(0, P, P)
for ( i in 1 : P) {
  if (i %% 5 == 0) print(i)
  for ( j in 1 : P) {
    r <- sqrt(x[, i] ^ 2 + x[, j] ^ 2)
    w1 <- x[, i] / r
    w2 <- x[, j] / r
    th <- quantile(r, u, na.rm = TRUE)
    id <- decls(r, th, 5)
    Sigma[i, j] <- sum(w1[id] * w2[id], na.rm = TRUE) / (length(id)) * 2
  }
}
Sigma



## Data loading.
library(Matrix)                         # nearPD
library(ggplot2)                        # plotting
library(ggpubr)                         # Align several plots from ggplot2
                                        # into a single panel
library(data.table)                     # To read data file

load("desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/prcp.RData")

saveRDS(prcp$long, "desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/long.RData")
saveRDS(prcp$lat, "desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/lat.RData")
saveRDS(prcp$obs, "desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/obs.RData")
saveRDS(prcp$time, "desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/time.RData")

## Precipitation data, 1950-1-1 to 2016-12-31
## 1 object:
##         prcp: a list with the following elements:
##                time: an array of dates of observation
##                 obs: observed precipitation (10th of cm)
##           lat, long: latitude and longitude of the stations

source("desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/tpdm.R")
## Codes for computing Tail Pairwise Dependence Matrix (TPDM)


## For Hurricane season data (August, September, October), Calculate TPDM,
## its eigen decomposition, and other quantities like principal components

aso.id <- take.month(prcp$time, c(8, 9, 10))
data.aso <- prcp$obs[aso.id, ]          # Data corresponding to the
                                        # Hurricane season
aso.ave <- ma(data.aso, 3)              # 3-day moving average
aso.smo <- ECDF.smoothing(aso.ave, 3)   # ECDF smoothing
aso.tran <- to.alpha.2(aso.smo)         # Transfer the margins to have the
                                        # same regular variation index: 2

sigma.aso <- rw.Sigma(aso.tran)         # Compute the TPDM

## Compute the positive definite matrix that is closest to sigma.aso, in
## case it is not PD.
## From package Matrix.
## sigma.aso.p <- nearPD(sigma.aso)
## sigma.aso <- sigma.aso.p[[1]]

eigen.aso <- eigen(sigma.aso)           # Eigen-decomposition of TPDM

U.aso <- -eigen.aso$vectors             # Basis vectors for PCA of TPDM

aso.tran.inv <- invTrans(aso.tran)      # Inverse transformation.
V.aso <- pc(U.aso, aso.tran.inv)        # Each column is a principal component.



##########################################################################################################
##########################################################################################################
##########################################################################################################
#source(tpdm.R)
#load libraries at line 60-ish
#load data sets (asotran.RData, prcp.RData, result.RData)

aso.id <- take.month(prcp$time, c(8, 9, 10))
aso.date <- prcp$time[aso.id]
N <- dim(V.aso)[1]


## To plot Basis vectors and Reconstructions.
long.aso <- prcp$long
lat.aso <- prcp$lat

## To plot time series of PCs with El-Nino and La Niña and ENSO-medium
## periods added.
enso <- fread(file = "desktop/fisica_mates/cinque/tfg-PCA_extremes/Yujing_RCode/enso.txt", select = c(1, 10), header = TRUE)

aso.date <- prcp$time[aso.id]
aso.date <- as.Date(aso.date)
N <- length(aso.date)

ENSO <- rep(NA, 67)                     # 67: From 1951 - 2016
for ( i in 1 : 67) {
    if (enso$ASO[i] <= -0.5) ENSO[i] <- 1
    else if (enso$ASO[i] <= 0.5) ENSO[i] <- 2
    else ENSO[i] <- 3
}
ENSO <- factor(ENSO, labels = c("low", "medium", "high"))
enso$ENSO <- ENSO

## Prepare year break for the time series of PC plot.
aso.yr <- take.year(prcp$time[aso.id])

yr.st <- rep(0, 67)
aso.yr.u <- unique(aso.yr)
for ( YR in 1 : 67) {
    yr.st[YR] <- which(aso.yr == aso.yr.u[YR])[1]
}

yr.bk <- rep(0, 3)
yri <- 1
for ( YR in c(1960, 1980, 2000)) {
    yr.bk[yri] <- which(aso.yr == YR)[1]
    yri <- yri + 1
}

mon.st <- rep(0, 67)
aso.mon <- take.month1(prcp$time[aso.id])
j <- 1
for ( i in 1 : 67) {
    mon.st[j] <- which((aso.yr == aso.yr.u[i]) & (aso.mon %in% c(8, 9, 10)))[1]
    j <- j + 1

}
enso$start <- mon.st
enso$end <- c(mon.st[-1]-1, length(aso.yr))





plot.list <- vector("list", 12)
j <- 1
## Legend
d <- scale_color_gradientn(colours = c("blue", "white", "red"),
                           limits = c(-0.11, 0.1), breaks = c(-0.11, 0,
                                                              0.1),
                           labels = c("0.64", "log2", "0.74"), name =
                                                                   "Basis\nSize")

plot.locVec2(U.aso[, 3], long.aso, lat.aso, limit = c(-0.11, 0.1), size = 2) + 
                            guides(fill = guide_legend(title = NULL, keywidth = 0.5,
                             keyheight = 1, reverse = TRUE)) + d

## Plot the first 6 basis and first 6 PCs
for ( i in 1 : 6) {
    ## Plot of the basis
    u.aso <- U.aso[, i]
    pu <- plot.locVec2(U.aso[, i], long.aso, lat.aso, limit = c(-0.11, 0.1), size = 2) +
        guides(fill = guide_legend(title = NULL, keywidth = 0.5,
                                   keyheight = 1, reverse = TRUE))+ d
    pu <- pu + ggtitle(bquote(e[.(i)]))
    plot.list[[j]] <- pu + theme(legend.title = element_text(size = 10))


    j <- j + 1
    pc.aso <- data.frame(date = aso.date, pc = V.aso[, i], order = c(1 : N))
    ppc <- ggplot(pc.aso, aes(x = order, y = pc))
    ppc <- ppc + scale_x_continuous(breaks = c(yr.bk), labels = c(1960, 1980,
                                                               2000))
    ppc <- ppc + geom_line(size = 0.2)
    yrng <- range(pc.aso$pc)
    ## Add ENSO periods
    ppc <- ppc + geom_rect(aes(NULL, NULL, xmin = start, xmax = end,
                               fill = ENSO), ymin = yrng[1], ymax = yrng[2],
                           data = enso)
    ## Add cross at the date of Hurricane Floyd
    data0 <- data.frame(x = 4555, y = V.aso[4555, i])
    ppc <- ppc + scale_fill_manual(values = alpha(c("blue", "grey", "red"), .3))  +
        guides(fill = guide_legend(title = "ONI", keywidth = 0.5, keyheight = 1, reverse = TRUE)) +
        labs(x = "Year", y = element_blank()) + geom_point(data = data0, aes(x = x, y = y), shape = 4, colour = "red", size = 3) + ggtitle(bquote(v[t*","*.(i)]))
 plot.list[[j]] <- ppc
    j <- j + 1
}

rec = reconstruct(V.aso, U.aso, 1140, 4555)
max(rec)

## Add quantile regression to PC1
library(quantreg)
pc1 <- data.frame(pc = V.aso[, 1], t = c(1:N))
pc1.qr <- rq(pc ~ t, 0.95, pc1)
pc1.qr.pred <- predict(pc1.qr)
pc1.qr.pred <- data.frame(pc = pc1.qr.pred, t = c(1: N))
pc1.plot <- plot.list[[2]] +
    geom_line( aes(x = t, y = pc), data = pc1.qr.pred, linetype = "dashed")
pc1.plot
plot.list[[2]] <- pc1.plot

##ADDED DC Aug 13:  Bootstrap p-value accounting for autocorrelation
#build a block bootstrap sample
blockSize <- 4
numBlocks <- ceiling(N/blockSize)
iter <- 1000
tvalVector <- numeric(iter)
for(i in seq(1, iter))
{
	blockMtx <- matrix(nrow = blockSize, ncol = numBlocks)
	blockMtx[1,] <- sample(1:N, size = numBlocks, replace = T)
	blockMtx[2,] <- min(blockMtx[1,] + 1, N)
	blockMtx[3,] <- min(blockMtx[1,] + 2, N)
	blockMtx[4,] <- min(blockMtx[1,] + 3, N)
	blockVec <- as.numeric(blockMtx)
	pc1Boot <- data.frame(pc = V.aso[blockVec, 1], t = c(1:N))
	bootTest <- rq(pc ~ t, 0.95, pc1Boot)
	tvalVector[i] <- summary(bootTest)$coefficients[2,3]
}
tvalVector <- c(tvalVector, summary(pc1.qr)$coefficients[2,3])



## Multiple plots.
## Basis vectors.
a <- ggarrange(plotlist = plot.list[seq(1, 11, 2)], ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom", font.label = list(size = 10, color = "black", face = "bold", family = NULL))
#pdf("pc.pdf", width = 4)

## Principal component time serieses.
b <- ggarrange(plotlist = plot.list[seq(2, 12, 2)], ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom", font.label = list(size = 10, color = "black", face = "bold", family = NULL))

pdf("Desktop/eigenvectors.pdf", width = 4)
a
dev.off()

pdf("Desktop/pc.pdf", width = 4)
b
dev.off()

## Reconstruction of the observations on Sep.16,1999, Hurricane Floyd
plot.list <- vector("list", 6)
index <- c(-1, 1140, 2, 6, 10, 20)      # -1: original obs, others: number
                                        # of basis vectors used in reconstruction
limit <- c(-10, 80)                     # Limit of the observations

id <- 4555


for ( i in 1 : 6) {
    k <- index[i]
    if (k != -1) {
    x1 <- (reconstruct(V.aso, U.aso, k, id))
    print(max(x1, na.rm = TRUE))


    plot.list[[i]] <- plot.locVec(x1, long.aso, lat.aso, limit, Size = 2, axis = FALSE)

    } else {
        plot.list[[i]] <- plot.locVec(aso.tran[id, ], long.aso, lat.aso, limit, Size = 2, axis = FALSE)
    }

}
plot.list[[1]] <- plot.list[[1]] + ggtitle("Transformed Data")
plot.list[[2]] <- plot.list[[2]] + ggtitle("All Eigenvectors")
for ( i in 3 : 6) {
    plot.list[[i]] <- plot.list[[i]] + ggtitle(paste(index[i], "Eigenvectors"))
}


d <- ggarrange(plotlist = plot.list, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")
pdf("Desktop/reconstruction.pdf")
d
dev.off()

## PC2-vs-PC6
score26 <- data.frame(date = aso.date, scorei = V.aso[, 2],
                             scorej = V.aso[, 6], order = c(1 : N))
score26$Year <- aso.yr
score26 <- merge(score26, enso, by = "Year")

pic.s <- ggplot(score26, aes(x = scorei, y = scorej, color = ENSO)) +
    geom_point(size = 2)
myColors <- c( "#619CFF", "#00BA38", "#F8766D")
names(myColors) <- levels(enso$ENSO)
colScale <- scale_color_manual(name = "ONI", values = myColors)
pic.s <- pic.s + colScale
pic.s <- pic.s +  guides(color = guide_legend(reverse = TRUE)) + labs(x = "PC 2", y = "PC 6")


range0 <- c(-150, -50)                  # Range for the lower box, score 2
range1 <- c(-150, -50)                  # Range for the upper box, score 2
range2 <- c(-100, -25)                   # Range for the lower box, score 6
range3 <- c(25, 100)                    # Range for the upper box, score 6


enso0 <- score26$ENSO
bid <- boxid(V.aso, 2, 6, range0, range1, range2, range3)
ensop <- enso.p(enso0, bid)
boxPlot <- add.box(pic.s, range0, range1, range2, range3)

pdf("boxplot.pdf")
boxPlot
dev.off()

