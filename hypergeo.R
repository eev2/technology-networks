folder = 'classes' # 'classes'

if (!require(BiasedUrn)) {
    install.packages("BiasedUrn", repos = "http://www.stats.bris.ac.uk/R/")
    library(BiasedUrn)
}

OBS = as.matrix(read.csv(paste0(folder, "_OBSCnt.csv"), header = FALSE))
nsec = scan(paste0(folder, "_nsec.csv"))
N = sum(nsec)
nyr = ncol(OBS)

fn = function(lgod, OBS) -log(dMFNCHypergeo(OBS, nsec, sum(OBS),
  exp(lgod), precision = 1E-7))

wghty = array(0,dim(OBS))
mxlly = hylly = numeric(nyr)
fity = list()

for (i in 1:nyr) {
#  fity[[i]] = fit = nlm(fn, rep(0,8), OBS = OBS[,i], iterlim = 25)
  fity[[i]] = fit = optim(rep(0,8), fn, OBS = OBS[,i], method = "L-BFGS-B",
                          control = list(maxit = 25))
  mxlly[i] = -fit$value
  wghty[,i] = exp(fit$par - fit$par[8])
  hylly[i] = -fn(rep(0,8), OBS[,i])
}

LRS = 2*(mxlly - hylly)

png(paste0(folder, "_variety.png"), height = 300)
par(mar = c(1.5,1.5,.1,.1), tcl = NA, mgp = c(0,.5,0))
plot(1980:2011, LRS, xlab = "", ylab = "")
abline(h = qchisq(.95,7), col = 2)
dev.off()

save.image(file = paste0(folder, "_hypergeo.RData"))
