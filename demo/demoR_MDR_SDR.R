n <- 1000
p <- 5

X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- as.matrix(sin(X %*% rep(1, p))+rnorm(n, mean = 0, sd = 0.2))

SDR1 <- cumuSIR(X = X, Y = Y)
SDR2 <- cumuSAVE(X = X, Y = Y)
SDR3 <- CVSDR(X = X, Y = Y)
SDR4 <- CVSDR(X = X, Y = Y, initial = rep(1, p))
MDR1 <- CVMDR(X = X, Y = Y)
MDR2 <- CVMDR(X = X, Y = Y, initial = rep(1, p))
SID1 <- SIDR(X = X, Y = Y)
SID2 <- SIDR(X = X, Y = Y, initial = SDR1$basis[, 1])
SID3 <- SIDR(X = X, Y = Y, initial = rep(1, p))
SID4 <- SIDR(X = X, Y = Y, initial = rep(1, p), kernel = "K2_Biweight")
SID5 <- SIDR(X = X, Y = Y, initial = rep(1, p), bandwidth = 1)
SIM1 <- SIMR(X = X, Y = Y)
SIM2 <- SIMR(X = X, Y = Y, initial = SDR1$basis[, 1])
SIM3 <- SIMR(X = X, Y = Y, initial = rep(1, p))
SIM4 <- SIMR(X = X, Y = Y, initial = rep(1, p), kernel = "K2_Biweight")
SIM5 <- SIMR(X = X, Y = Y, initial = rep(1, p), bandwidth = 1)
MID1 <- MIDR(X = X, Y = Y, n.index = 1)
MID2 <- MIDR(X = X, Y = Y, n.index = 2)
MIM1 <- MIMR(X = X, Y = Y, n.index = 1)
MIM2 <- MIMR(X = X, Y = Y, n.index = 2)

SDR3 <- CVSDR(X = X, Y = Y)
SDR4 <- CVSDR(X = X, Y = Y, initial = rep(1, p))

system.time(
  SDR3 <- SIDR(X = X, Y = Y)
)

system.time(
  SDR5 <- SIDR(X = X, Y = Y, initial = rep(1, p))
)
