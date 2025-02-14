library(ggplot2)
library(ggrepel, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R") 
library(tidyverse)
library(viridisLite)
#library(viridis)
library(ggsignif, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")
library(reshape2, lib.loc = "/dss/dsshome1/lxc0E/di67kah/R")

com <- c("W1x")
header <- c("RUN","NPOP1","NPOP2","NPOP3","NPOP4","NANC1","NANC2","NANC3","TMIG2S","MIG01","MIG10",
            "MIG12","MIG21","MIG23","MIG32","MIGA0A1","MIGA1A0","MIGA1A2","MIGA2A1","MIGA2A3","MIGA3A2", 
            "TMIG2E","TDIV3","TDIV2","TDIV1","MaxEstLhood","MaxObsLhood")

setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/CCW_CCE_COW_COE_folded/fastsimcoal2/bestruns")
standardmodel <- read.delim(paste("./",com[1],"/",com[1],".bestlhoods",sep=""), header=TRUE, sep="\t")
RUN <- 0
standardmodel <- cbind(RUN,standardmodel)
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2023__MigratoryGenomics/sim_final/CCW_CCE_COW_COE_folded/fastsimcoal2/bestruns/modelfit/W1x_finite")
modelbs <- read.delim(paste(com[1],"_finite.bestlhoods.all",sep=""), header=FALSE, sep="\t")
colnames(modelbs) <- header
model <- rbind(standardmodel,modelbs)
model <- model[,c(-26, -27)]

# bias correction

bc <-  function(t0, tt){2*t0-colMeans(tt)}
ci <- function(tt){quantile(tt, c(0.025, (1 - 0.025)))}
t0 = model[1,][2:(length(model))]
tt = modelbs[2:(length(modelbs)-2)]

for (i in 1:length(t0)) {
  print(bc(t0[i],tt[i]))
  print(ci(tt[,i]))  
}

# boxplot

pdf(file="Model1parabs_boxplot.pdf")
modelSPA$RUN <- c(TRUE, rep(FALSE, nrow(modelSPA) - 1L))  # tag first row as interesting
df.2 <- melt(modelSPA[1:(length(modelSPA)-1)])  # convert df to long format
ggplot(subset(df.2, !RUN), aes(x=variable, y=value)) + 
  geom_boxplot() + scale_y_log10() +
  geom_point(data=subset(df.2, RUN), aes(x=variable, y=value), color="red", size=2)
dev.off()

