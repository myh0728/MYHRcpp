n <- 100
p <- 5

Xi <- matrix(rnorm(n * p), n, p)
Ci <- 1 + rexp(n)
Zi <- rexp(n)
beta0 <- rep(1, p)
tau <- 3
test.data <- simRecur.GFPI(Xi = Xi, Zi = Zi, Ci = Ci,
                           beta0 = beta0, lambda0 = function(t){0.5},
                           end.study = tau)

shuf <- sample(1:length(test.data$id), length(test.data$id), replace = FALSE)
test.recurReg.proData <- recurReg.proData(id = test.data$id[shuf],
                                          t.stop = test.data$t.stop[shuf],
                                          covariate = test.data[shuf, c("X1", "X2")],
                                          sort = TRUE)

################################################################################

data(readmission, package = "frailtypack")

test.data <- data.frame(id = readmission$id,
                        t.start = readmission$t.start,
                        t.stop = readmission$t.stop,
                        is.event = readmission$event,
                        chemo = readmission$chemo == "Treated",
                        gender = readmission$sex == "Male",
                        charlson12 = readmission$charlson == "1-2",
                        charlson3 = readmission$charlson == "3",
                        dukesC = readmission$dukes == "C",
                        dukesD = readmission$dukes == "D")

##### Combining pseudo-partial likelihood scores estimation

Rij <- recurReg.proData(id = test.data$id,
                        t.start = test.data$t.start,
                        t.stop = test.data$t.stop,
                        is.event = test.data$is.event,
                        covariate = test.data[c("chemo", "gender",
                                                "charlson12", "charlson3",
                                                "dukesC", "dukesD")])

# certain options of w.t: c("unit", "Gehan", "cumbase", "S0t")
test.PPL <- CWPPL(Rij = Rij) # default w.t = list("unit)
test.PPL <- CWPPL(Rij = Rij, w.t = list("Gehan"))
test.PPL <- CWPPL(Rij = Rij, w.t = list(function(t){1 / (1 + t)}))
test.PPL <- CWPPL(Rij = Rij, w.t = list(Rij$t.event ^ 0.5))
test.PPL <- CWPPL(Rij = Rij, w.t = list("unit", function(t){1 / (1 + t)})) # default method = "EL"
test.PPL <- CWPPL(Rij = Rij, w.t = list("cumbase", Rij$t.event ^ 0.5), method = "GMM")
test.PPL <- CWPPL(Rij = Rij, w.t = list("unit", function(t){1 / (1 + t)}, "cumbase"))

inferCWPPL.boot(esti = test.PPL, method = "naive", clevel = 0.95,
                n.boot = 10, seed = 123)
inferCWPPL.boot(esti = test.PPL, method = "RWB.gamma42", clevel = 0.95,
                n.boot = 10, seed = 123)
inferCWPPL.boot(esti = test.PPL, method = "RWB.exp", clevel = 0.95,
                n.boot = 10, seed = 123)
inferCWPPL.boot(esti = test.PPL, method = "ELcali", clevel = 0.95,
                n.boot = 10, seed = 123)
# default method = "ELcali"
# default confidence level: clevel = 0.95
# default n.boot = 500
# default seed = NULL

inferCWPPL.SEperturb(esti = test.PPL, method = "RWB.gamma42",
                     n.boot = 10, seed = 123)
inferCWPPL.SEperturb(esti = test.PPL, method = "RWB.exp",
                     n.boot = 10, seed = 123)
inferCWPPL.SEperturb(esti = test.PPL, method = "naive",
                     n.boot = 10, seed = 123)
# default method = "RWB.gamma42"
# default n.boot = 500
# default seed = NULL




