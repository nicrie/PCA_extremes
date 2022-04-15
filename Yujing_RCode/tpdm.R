## Data preprocessing.


take.month <- function(time, month) {
    ## Obtain the index of observations that in the corresponding months.
    ## Input:
    ##   time: an array of dates
    ##  month: a vector of months
    ## Output:
    ##   a vector of integers, index of observations that in the corresponding months.
    library(stringr)
    ymd <- str_split(time, "-")
    Month <- as.integer(sapply(ymd, function(x) x[2]))
    which(Month %in% month)
}


## Moving average.
## To align days of observations that is related to the same
## extreme event but happen on consecutive but different days.
ma.one <- function(x, k = 3) {
    ## k-day moving average for one station
    res <- x
    n <- length(x)
    k2 <- floor(k / 2)
    K <- k2 * 2
    x <- c(rep(NA, k2), x)
    for( i in 1 : n) {
        if (is.na(x[i + k2])) {
            res[i] <- NA
        } else {
            id <- c(i : (i + K))
            res[i] <- mean(x[id], na.rm = TRUE)
        }
    }
    res
}

ma <- function(X, k = 3) {
    ## k-day moving average for all stations
    apply(X, 2, ma.one)
}

## ECDF smoothing.
## To remove the temporal depdence induced by the moving-average process.
ECDF.smoothing.one <- function(x, k = 3) {
    ## ECDF smoothing for one station
    n <- length(x)
    dat <- vector("list", (k + 1))
    ecdf.list <- vector("list", k)
    each <- matrix(0, n, k)
    for ( i in 1 : k) {
        id <- seq(i, n, by = k)
        dat[[i]] <- sort(as.vector(x[id]))
    }
    ni <- sapply(dat, length)

    for ( i in 1 : k) {
        ecdf.list[[i]] <- approxfun(dat[[i]],
                                    y = seq(1, ni[i]) / (ni[i] + 1),
                                    yleft = 0, yright = 1 - 1e-9, ties="ordered")
        each[, i] <- ecdf.list[[i]](x)
    }
    apply(each, 1, mean)
}

ECDF.smoothing <- function(X, k = 3) {
    ## ECDF smoothing for all stations.
    apply(X, 2, ECDF.smoothing.one, k = k)
}

## Probability integral transformation
## To make all of the margins follow a Frechet distribution with alpha = 2.

to.alpha.2 <- function(X) {
    sqrt(1 / (-log(X)))
}

## Compute TPDM

decls <- function(x, th, k) {
    ## Ordinary decluster.
    id.big <- which(x > th)
    id.dif <- diff(id.big)
    tick <- which(id.dif >= k)
    start <- id.big[c(1, tick + 1)]              # Where a new cluster begins
    end <- c(id.big[tick], last(id.big, 1))
    n <- length(start)
    id.res <- rep(0, n)
    for ( i in 1 : n) {
        temp <- x[start[i] : end[i]]
        id.res[i] <- which(temp == max(temp, na.rm = TRUE))[1] + start[i] - 1
    }
    id.res
}

rw.Sigma <- function(X, u = 0.80) {
    ## Compute TPDM
    P <- ncol(X)                      # Number of stations
    M <- nrow(X)                      # Number of obs
    Sigma <- matrix(0, P, P)
    for ( i in 1 : P) {
        if (i %% 5 == 0) print(i)
        for ( j in 1 : P) {
            r <- sqrt(X[, i] ^ 2 + X[, j] ^ 2)
            w1 <- X[, i] / r
            w2 <- X[, j] / r
            th <- quantile(r, u, na.rm = TRUE)
            id <- decls(r, th, 5)
            Sigma[i, j] <- sum(w1[id] * w2[id], na.rm = TRUE) / (length(id)) * 2
        }
    }
    Sigma
}


##some functions
##applies the transformation t
trans <- function(x)
{
    ##because it takes an exponential, this function flakes out if x is too big
    ##hence for big values of x, we return x
    v <- log(1 + exp(x))
    id <- which(x < -20)
    v[!is.finite(v)] <- x[!is.finite(v)]
    v[id] <- exp(x[id])
    return(v)
}

##applies the inverse transformation t^{-1}
invTrans <- function(v)
{
    ##same trickeration for big values of v
    ##still returns -Inf if v is machine zero
    x <- log(exp(v) - 1)
    x[!is.finite(x) & v > 1 & !is.na(x)] <- v[!is.finite(x) & v > 1 &
                                                  !is.na(x)]

    return(x)
}

#does matrix multiplication
AMult <- function(A, v)
{
    id <- which(!is.na(v))
    x <- invTrans(v)
    if (length(x) > 1) A <- A[, id]
    x <- x[id]

    if(length(x) == 1) Ax <- A * x
    else  Ax <- A %*% x
    return(list(trans(Ax), which(is.na(v))))
}



## Calculate Principal Component
pc.one <- function(i, U, invX) {
    id.na <- is.na(invX[i, ])
   t(U[!id.na, ]) %*% invX[i, !id.na]
}

pc <- function(U, invX) {
    res <- invX
    n <- nrow(invX)
    for ( i in 1 : n) {
        res[i, ] <- pc.one(i, U, invX)
    }
    res
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

                                        # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



plot.locVec <- function(vec, long, lat,
                        limit = range(vec, na.rm = TRUE), Size = 0.5, axis = FALSE) {
    frame <- data.frame(vec = vec, long = long, lat = lat)
    all_states <- map_data("state")

    p <- ggplot()
    p <- p + geom_polygon(data = all_states, aes(x = long, y = lat,
                                                 group = group))
    p <- p + geom_point(data = frame,
                        aes(x = long, y = lat, color = vec), size = Size)

    p <- p + scale_color_gradientn(colours = rev(rainbow(7))[-c(1, 2)],
                                   limits = limit,
                                   name = "Transformed \n Precipitation") +
        theme(legend.title=element_text(size=rel(0.8)),
              legend.text=element_text(size=rel(0.7)))
    if (axis == FALSE) {
        p <- p + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank())
    }
    p

}



plot.locVec2 <- function(vec, long, lat,
                         limit, size = 1, axis = FALSE) {
    frame <- data.frame(vec = vec, long = long, lat = lat)
    all_states <- map_data("state")

    p <- ggplot()
    p <- p + geom_polygon(data = all_states, aes(x = long, y = lat,
                                                 group = group))
    p <- p + geom_point(data = frame,
                        aes(x = long, y = lat, color = vec), size = size)

    ## p <- p + scale_color_gradient2(limits  =  limit, low = "blue",
    ##                                high = "red", na.value = "yellow")  +
    ##     theme(legend.title=element_text(size=rel(0.8)),
    ##           legend.text=element_text(size=rel(0.7)))
    if (axis == FALSE) {
        p <- p + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank())
    }
    p

}


take.month1 <- function(time){
    ## For each observation, output the month it is from
    library(stringr)
    ymd <- str_split(time, "-")
    year <- as.integer(sapply(ymd, function(x) x[2]))
    year
}


take.year <- function(time) {
    ## For each observation, output the year it is from
    library(stringr)
    ymd <- str_split(time, "-")
    year <- as.integer(sapply(ymd, function(x) x[1]))
    year
}


reconstruct <- function(V, U, k, day) {
    trans(U[, 1 : k] %*% V[day, 1 : k])
}


plot.trans <- function(x) {
    l <- log(2) / (0.74 - log(2))
    if (x > log(2)) {
        res <- l * (x - log(2)) + log(2)
    } else {
        res <- x
    }
    return(res)
}

plot.trans.inv <- function(x) {
    l <- log(2) / (0.74 - log(2))
    if (x > log(2)) {
        res <- (x - log(2)) / l + log(2)
    } else {
        res <- x
    }
    return(res)
}

ranid <- function(V, k, range) {
    a <- V.aso[, k] > min(range)
    b <- V.aso[, k] < max(range)
    which(a & b)
}


boxid <- function(V, k1, k2, r0, r1, r2, r3) {
    a <- ranid(V, k1, r0)
    b <- ranid(V, k1, r1)

    c <- ranid(V, k2, r2)
    d <- ranid(V, k2, r3)

    list(upper = intersect(b, d), lower = intersect(a, c))
}

enso.p <- function(enso0, boxid) {
    up <- boxid$upper
    lo <- boxid$lower

    eup <- enso0[up]
    elo <- enso0[lo]
    xup <- sum(eup == "low")
    nup <- length(up)

    xlo <- sum(elo == "low")
    nlo <- length(lo)
    prop.test(c(xup, xlo), c(nup, nlo), correct = FALSE)
}


add.box <- function(pic, r0, r1, r2, r3) {
    x1 <- c(r0[1], r1[1])
    x2 <- c(r0[2], r1[2])
    y1 <- c(r2[1], r3[1])
    y2 <- c(r2[2], r3[2])
    d <- data.frame(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
    pic.s + geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1,
    ymax=y2, fill = c('a', 'b')), alpha = 0, inherit.aes=FALSE, color="black") + guides(fill = FALSE)
}
