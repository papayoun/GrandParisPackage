source("vignettes/betaCode_unknownSigma.R")
source("vignettes/debugSelectedParticles.R")

particleSet <- cbind(Step1 = Res$Particles[, 1], Step2 = Res$Particles[,2],
                     Selec1 = Sel)
