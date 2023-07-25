n <- 500
p <- 10

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- as.matrix(exp(X %*% rep(1, p)) + rnorm(n, mean = 0, sd = 0.2))

SDR1 <- cumuSIR(X = X, Y = Y)
SDR2 <- cumuSAVE(X = X, Y = Y)

SIM1 <- SIMR(X = X, Y = Y)
SIM2 <- SIMR(X = X, Y = Y, initial = SDR1$basis[, 1])
SIM3 <- SIMR(X = X, Y = Y, initial = rep(1, p))
SIM4 <- SIMR(X = X, Y = Y, initial = rep(1, p), kernel = "K2_Biweight")
SIM5 <- SIMR(X = X, Y = Y, initial = rep(1, p), bandwidth = 1)
SIM6 <- SIMR(X = X, Y = Y, initial = rep(1, p), method = "nlminb")
SIM7 <- SIMR(X = X, Y = Y, initial = rep(1, p), method = "nmk")

MIM1 <- MIMR(X = X, Y = Y, n.index = 1)
MIM2 <- MIMR(X = X, Y = Y, n.index = 2)
MIM3 <- MIMR(X = X, Y = Y, n.index = 2, initial = SDR1$basis[, 1:2])
MIM4 <- MIMR(X = X, Y = Y, n.index = 2, initial = SDR1$basis[, 1:2], bandwidth = 1)
MIM5 <- MIMR(X = X, Y = Y, n.index = 2, initial = SDR1$basis[, 1:2], method = "nmk")

MDR1 <- CVMDR(X = X, Y = Y)
MDR2 <- CVMDR(X = X, Y = Y, initial = SDR1$basis)
MDR3 <- CVMDR(X = X, Y = Y, initial = rep(1, p))
MDR4 <- CVMDR(X = X, Y = Y, initial = rep(1, p), method = "nlminb")

SID1 <- SIDRuniY(X = X, Y = Y)
SID2 <- SIDRuniY(X = X, Y = Y, initial = SDR1$basis[, 1])
SID3 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p))
SID4 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p), kernel = "K2_Biweight")
SID5 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p), bandwidth = 1)
SID6 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p), method = "nlminb")
SID7 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p), method = "nmk")
SID8 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p), dist.mode = "sample")
SID9 <- SIDRuniY(X = X, Y = Y, initial = rep(1, p), dist.mode = "quantile")
SID10 <- SIDRmultiY(X = X, Y = Y, initial = rep(1, p))
SID11 <- SIDRmultiY(X = X, Y = Y, initial = rep(1, p), dist.mode = "sample")

MID1 <- MIDRuniY(X = X, Y = Y, n.index = 1)
MID2 <- MIDRuniY(X = X, Y = Y, n.index = 2)
MID3 <- MIDRmultiY(X = X, Y = Y, n.index = 2)

SDR1 <- CVSDRuniY(X = X, Y = Y)
SDR2 <- CVSDRuniY(X = X, Y = Y, initial = rep(1, p))
SDR3 <- CVSDRmultiY(X = X, Y = Y, initial = rep(1, p))



