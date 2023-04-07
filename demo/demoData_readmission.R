data(readmission, package = "frailtypack")

test.data <- data.frame(id = readmission$id,
                        t.start = readmission$t.start,
                        t.stop = readmission$t.stop,
                        is.event = readmission$event,
                        chemo = readmission$chemo=="Treated",
                        gender = readmission$sex=="Male",
                        charlson12 = readmission$charlson=="1-2",
                        charlson3 = readmission$charlson=="3",
                        dukesC = readmission$dukes=="C",
                        dukesD = readmission$dukes=="D")

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
test.PPL <- CWPPL(Rij = Rij, w.t = list(function(t){1/(1+t)}))
test.PPL <- CWPPL(Rij = Rij, w.t = list(Rij$t.event^0.5))
test.PPL <- CWPPL(Rij = Rij, w.t = list("unit", function(t){1/(1+t)})) # default method = "EL"
test.PPL <- CWPPL(Rij = Rij, w.t = list("cumbase", Rij$t.event^0.5), method = "GMM")
test.PPL <- CWPPL(Rij = Rij, w.t = list("unit", function(t){1/(1+t)}, "cumbase"))

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



