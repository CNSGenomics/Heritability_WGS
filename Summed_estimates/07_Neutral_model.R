

##### This code calculates the standard error when merging the estimates from different bins
# It requires:
# - The raw estimates from the GCTA --reml option (a *.hsq file), ie the V(Gx), V(e) and V(p) without the V(Gx)/Vp rows
# - The variance-covariance matrix of the variance components at the end of the .log file.
# - The summed estimates within each MAF bin (summed across all LD bins)

library(data.table)
library(tidyverse)

#Load the raw variance components estimates for height and BMI
raw_bmi <- fread("Raw_est_BMI_Quartile48pcs", h=F)
raw_height<- fread("Raw_est_height_Quartile48pcs", h=F)

#Load the variance-covariance matrix of the variance estimates
varcov_bmi <- fread("Varcovar_BMI_Quartile48pcs", h=F, select=c(1:17))
varcov_height <- fread("Varcovar_height_Quartile48pcs", h=F, select=c(1:17))

#Load the summed scaled estimates (summed across all LD bins within each MAF bin)
var_main_evol <- fread("07_Neutral_model_summed_setimates.txt", h=T)

seheightcumquant <- data.frame(matrix(ncol=1,nrow = 4))
ngrm <- 16 # Total number of GRMs
nld <- 4 # Number of LD bins

v <- sum(raw_height$V2[1:ngrm])
for (i in 2:4){
  varU <- sum(varcov_height[c(1:(nld*i)),c(1:(nld*i)),with=F])
  varV <- sum(varcov_height[c(1:ngrm),c(1:ngrm),with=F])
  u <- sum(raw_height[1:(nld*i),2])
  covuv <- sum(varcov_height[c(1:ngrm),c(1:(nld*i)),with=F])
  seheightcumquant[i,1] <- sqrt((u/v)^2 * (varU/u^2 + varV/v^2 - 2*covuv/(u*v) ))
}

v <- sum(raw_bmi$V2[1:ngrm])

seBMIcumquant <- data.frame(matrix(ncol=1,nrow = 7))
for (i in 2:4){
  varU <- sum(varcov_bmi[c(1:(nld*i)),c(1:(nld*i)),with=F])
  varV <- sum(varcov_bmi[c(1:ngrm),c(1:ngrm),with=F])
  u <- sum(raw_bmi[1:(nld*i),2])
  covuv <- sum(varcov_bmi[c(1:ngrm),c(1:(nld*i)),with=F])
  seBMIcumquant[i,1] <- sqrt((u/v)^2 * (varU/u^2 + varV/v^2 - 2*covuv/(u*v) ))
}


var_main_evol$se <- c(seBMIcumquant[2:4,], seheightcumquant[2:4,])
var_main_evol$se[c(3,6)] <- 0.00001
supp_plot_neutral_model <- ggplot(var_main_evol %>% mutate(cv=cumsum(n)/sum(n)), 
  aes(x = Binnum, y = cv, colour = Trait)) + 
  geom_smooth(method='lm', se=FALSE) +
  geom_point(size=3, aes(shape= Trait)) +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color="black", linetype=3) +
  scale_y_continuous(limits =c(0, 1.01))+
  geom_errorbar(aes(ymin=cv-se, ymax=cv), width=.005) +
  labs(x ="MAF", y = expression("Cumulative contribution to genetic variance" %+-% SE)) + 
  scale_x_continuous(limits =c(0, 0.5)) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="black"),
        legend.position = c(0.7, 0.3)) 

