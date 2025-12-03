library(simest)

## Test smooth.pen.reg(*, agcv=TRUE) ---------------

## -- Using a subset of R's sunspots data  ------
osV <- abbreviate(gsub("[^[:alnum:]]", '', sub("\\(.*", '', osVersion)), 12)

if(!dev.interactive(TRUE)) pdf(paste0("smooth.pen_sunsp__", osV, ".pdf"), width = 9, height=5)

str(ssp <- window(sunspot.m2014, start = 1900))
plot(ssp)

yss <- as.numeric(ssp)
xss <- as.numeric(time(ssp))

fm1 <- smooth.pen.reg(xss, yss, lambda = 0.1)
plot(fm1)

fm1s <- smooth.pen.reg(xss, sqrt(yss), lambda = 0.1)
plot(fm1s) ## residual plots a really good

## Compare with Stahel's version of "started log"
##
## This is from example(u.log, package = "sfsmisc") [MM: ~/R/Pkgs/sfsmisc/man/u.log.Rd ]

  ## by MM: "modularized" by providing a threshold-computer function separately:
  logst_thrWS <- function(x, mult = 1) {
      lq <- quantile(x[x > 0], probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
      if (lq[1] == lq[2])
          lq[1] <- lq[2]/2
      lq[1]^(1 + mult)/lq[2]^mult
  }
  logst0 <- function(x, calib = x, threshold = thrFUN(calib, mult=mult),
                     thrFUN = logst_thrWS, mult = 1, base = 10) {
      ## Using pmax.int() instead of logical indexing -- NA's work automatically - even faster
      xm <- pmax.int(threshold, x)
      res <- log(xm, base) + (x - xm)/(threshold * log(base))
      attr(res, "threshold") <- threshold
      attr(res, "base") <- base
      res
  }

fm1l <- smooth.pen.reg(xss, logst0(yss), lambda = 0.1)
plot(fm1l)

##------- Use sqrt() transformed yss for a bit --------

set.seed(101)
(t1 <- system.time(
fm1sG <- smooth.pen.reg(xss, sqrt(yss), lambda = 0.1, agcv=TRUE, fit = TRUE)
))
(t05 <- system.time(
fm05sG <- smooth.pen.reg(xss, sqrt(yss), lambda = 0.05, agcv=TRUE, fit = TRUE)
))
(t2 <- system.time(
fm2sG <- smooth.pen.reg(xss, sqrt(yss), lambda = 0.2, agcv=TRUE, fit = TRUE)
))
(t4 <- system.time(
fm4sG <- smooth.pen.reg(xss, sqrt(yss), lambda = 0.4, agcv=TRUE, fit = TRUE)
))

str(fm4sG); fm4sG
str(fm2sG); fm2sG
str(fm1sG); fm1sG
str(fm05sG);fm05sG

p1 <- predict(fm1sG)
## should give the same but uses different code
p1.<- predict(fm1sG, newdata = xss)
all.equal(p1, p1., tolerance = 0) # is TRUE, too [Fedora 42 Lnx, x86_64]
stopifnot(all.equal(p1, p1.))

pmat <- cbind(p01 = predict(fm05sG), p1, p2 = predict(fm2sG), p4 = predict(fm4sG))
cbind(xss, rt.y = sqrt(yss), pmat)
plot(xss, sqrt(yss), cex = 3/4, col = adjustcolor(1, 1/3))
title("smooth.pen.reg(xss, sqrt(yss), lambda = L)   for L = 0.05, 0.1, 0.2, 0.4")
matlines(xss, pmat, lwd = 3, col = adjustcolor(2:5, 1/2))
legend("topleft", paste("lambda =", c(0.05, 0.1, 0.2, 0.4)), col=2:5, lty=1:4, lwd=3, bty="n")
