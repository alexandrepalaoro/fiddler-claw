library(ape)
library(phytools)
library(geiger)
library(tidyverse)
library(bbmle)
library(mvMORPH)
library(vioplot)
library(ellipse)
library(phylolm)
library(ratematrix)
library(viridis)
library(scales)
library(vioplot)
library(performance)
library(corncob)
library(lattice)
library(latticeExtra)
library(RColorBrewer)
library(sp)
library(pander)
library(scatterplot3d)

#### Along with all the packages, we are also an analysis provided in the following paper:

# Slater, G. J., & Friscia, A. R. (2019). Hierarchy in adaptive radiation: a case study using the Carnivora (Mammalia). Evolution, 73(3), 524-539. 

#### I will explain more about it below.

source("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/code/RateMatrixFunctions.R")

## Loading data

Setting up the code and loading the phylogeny and data. Since we are using the datasheet that contains the measurements from all individuals, we are loading it and summarizing it by species. That way we can run the comparative methods.

#### Loading tree ####

phy = read.nexus("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/Fiddler_complete.nex")

plot(phy); axisPhylo()

#### Loading raw data ####

uca = read.csv("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/def_linear_01-08-25.csv" , h = T, sep = ';')

#### Summarizing data ####
uca$prop.claw = uca$claw_size / uca$carapace
uca$prop.polex = uca$polex / uca$carapace
uca$prop.manus = uca$manus / uca$carapace
uca$prop.inlever = uca$in.lever / uca$carapace
uca$prop.out1 = uca$out.lever1 / uca$carapace
uca$prop.out2 = uca$out.lever2 / uca$carapace

uca.short = uca %>% group_by(sp) %>%
  summarise(
    mean.carapace=mean( carapace , na.rm=T),
    mean.claw = mean(claw_size, na.rm=T),
    mean.tip=mean(tip_adv, na.rm=T),
    mean.tub=mean(tub_adv, na.rm=T),
    mean.polex=mean(polex, na.rm=T),
    mean.manus=mean(manus, na.rm=T),
    mean.out.tub = mean(out.lever1, na.rm=T),
    mean.dactyl = mean(out.lever2, na.rm = T),
    mean.inl = mean(in.lever, na.rm = T),
    sd.carapace = sd(carapace, na.rm=T),
    sd.claw = sd(prop.claw, na.rm=T),
    sd.tip = sd(tip_adv, na.rm=T),
    sd.tub = sd(tub_adv, na.rm=T),
    sd.polex = sd(prop.polex, na.rm=T),
    sd.manus = sd(prop.manus, na.rm=T),
    sd.out.tub = sd(prop.out1, na.rm=T),
    sd.dactyl = sd(prop.out2, na.rm = T),
    prop.inlever = mean( prop.inlever ),
    prop.out1 = mean( prop.out1 , na.rm = T ),
    prop.out2 = mean( prop.out2 , na.rm = T),
    prop.manus = mean( prop.manus , na.rm = T ),
    prop.polex = mean( prop.polex , na.rm = T ),
    prop.claw = mean( prop.claw , na.rm = T ),
    
    
  )

uca.short = as.data.frame(uca.short)

rownames(uca.short) = uca.short$sp

name.check(phy, uca.short)

#### Loading mating strategies ####

uca.mate = read.csv("data/def_mating_short.csv" , h = T, sep = ';')
uca.mate = uca.mate %>% arrange (sp)

uca.short$mating = uca.mate$mating


#### PGLS ####

# Our first series of tests regard claw evolutionary allometry and changes in the biomechanical relationships 
# within components of the claw. To do this, we performed four PGLS regressions you can see below.


### Claw evolutionary allometry

# I saved the results in RDS because I bootstraped the model. It saves time in the long run...
claw1_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw1_BM.rds")
claw1_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw1_OUf.rds")
claw1_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw1_OUr.rds")


### The models are below

#claw1_BM = phylolm(log(mean.claw) ~ mating * log(mean.carapace), boot = 1000, 
#            data=uca.short, phy=phy, model='BM', measurement_error=T)

#claw1_OUf = phylolm(log(mean.claw) ~ mating * log(mean.carapace), boot = 1000, 
#                data=uca.short, phy=phy, model='OUfixedRoot', measurement_error=T)

#claw1_OUr =  phylolm(log(mean.claw) ~ mating * log(mean.carapace), boot = 1000, 
#                data=uca.short, phy=phy, model='OUrandomRoot', measurement_error=T,
#                upper.bound = 300)

###

compare_performance(claw1_BM, claw1_OUf, claw1_OUr)

## BM was chosen because it has the lowest AIC value

plot(claw1_BM)

summary(claw1_BM)

plot( log(mean.claw) ~ log(mean.carapace), data = uca.short,
      cex = 2 , pch = 21, bg = c("black","white")[as.numeric(as.factor(uca.short$mating))])

curve( coef(claw1_BM)[1] + (coef(claw1_BM)[3] * x), add = T,
       col = 'red', lty = 2, lwd = 4)

### Claw lever versus claw strength

###### CLAW DACTYL LENGTH VERSUS CLAW MANUS LENGTH - PROPORTIONAL #####

claw5_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw5_BM.rds")
claw5_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw5_OUf.rds")
claw5_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw5_OUr.rds")

### The models are below
#claw5_BM = phylolm(prop.out2 ~ prop.manus * mating, data = uca.short, boot = 1000,
#                   phy = phy, model = "BM", measurement_error = T)

#claw5_OUf = phylolm(prop.out2 ~ prop.manus * mating, data = uca.short, boot = 1000,
#                    phy = phy, model = "OUfixedRoot", measurement_error = T)

#claw5_OUr = phylolm(prop.out2 ~ prop.manus * mating, data = uca.short, boot = 1000,
#                    phy = phy, model = "OUrandomRoot", measurement_error = T,
#                    upper.bound = 20)  

###

compare_performance(claw5_BM, claw5_OUf, claw5_OUr)

## Once again, BM was the model with the lowest AIC

summary(claw5_BM)

#Filled dots = burrow mating system; white dots = surface mating system

plot( prop.out2 ~ prop.manus, data = uca.short,
      cex = 2 , pch = 21, bg = c("black","white")[as.numeric(as.factor(uca.short$mating))])

#### Mechanical advantage versus relative claw muscle ####

# Given that not all species have tubercles on their dactyl (the movable finger), we performed one test of the 
# mechanical advantage calculated at the tip, and another test where mechanical advantage was calculated at the 
# tubercule. The second test excluded 9 species our sample. 

#### Mechanical advantage at the tip ####

claw6_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw6_BM.rds")
claw6_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw6_OUf.rds")
claw6_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw6_OUr.rds")

### The models are below
#claw6_BM = phylolm(mean.tip ~ prop.manus * mating, data = uca.short, boot = 1000,
#                   phy = phy, model = "BM", measurement_error = T)

#claw6_OUf = phylolm(mean.tip ~ prop.manus * mating, data = uca.short, boot = 1000,
#                    phy = phy, model = "OUfixedRoot", measurement_error = T)

#claw6_OUr = phylolm(mean.tip ~ prop.manus * mating, data = uca.short, boot = 1000,
#                    phy = phy, model = "OUrandomRoot", measurement_error = T,
#                    upper.bound = 10^3)  

compare_performance(claw6_BM, claw6_OUf, claw6_OUr)

summary(claw6_BM)

plot( mean.tip ~ prop.manus, data = uca.short,
      cex = 2 , pch = 21, 
      bg = c("black","white")[as.numeric(as.factor(uca.short$mating))])


#### Mechanical advantage at the tubercle ####

claw7_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw7_BM.rds")
claw7_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw7_OUf.rds")
claw7_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw7_OUr.rds")

### The models are below
#claw7_BM = phylolm(mean.tub ~ prop.manus * mating, data = uca.short, boot = 1000,
#                   phy = phy, model = "BM", measurement_error = T)

#claw7_OUf = phylolm(mean.tub ~ prop.manus * mating, data = uca.short,  boot = 1000,
#                    phy = phy, model = "OUfixedRoot", measurement_error = T)

#claw7_OUr = phylolm(mean.tub ~ prop.manus * mating, data = uca.short, boot = 10,
#                    phy = phy, model = "OUrandomRoot", measurement_error = T)  
###

compare_performance(claw7_BM, claw7_OUf, claw7_OUr)

## Even though the OU fixed root and BM are tied, BM is a much simpler model

summary(claw7_BM)

plot( mean.tub ~ prop.manus, data = uca.short,
      cex = 2 , pch = 21, bg = c("black","white")[as.numeric(as.factor(uca.short$mating))])


#### PHYLO REGRESSION GRAPHS ####

## Evolutionary allometry

plot( log(mean.claw) ~ log(mean.carapace), data = uca.short,
      cex = 2.5 , pch = 21, bg = c("black","grey")[as.numeric(as.factor(uca.short$mating))],
      las = 1, cex.axis = 1.5, cex.lab = 1.5,
      bty = 'l',
      ylab = "Mean claw length (log)",
      xlab = "Mean carapace width (log)")

curve( coef(claw1_BM)[1] + (coef(claw1_BM)[3] * x), add = T,
       col = 'black', lty = 2, lwd = 4)

## Dactyl and manus

plot( prop.out2 ~ prop.manus, data = uca.short,
      cex = 2.5 , pch = 21, bg = c("black","grey")[as.numeric(as.factor(uca.short$mating))],
      las = 1, cex.axis = 1.5, cex.lab = 1.5,
      bty = 'l',
      ylab = "Proportional dactyl length",
      xlab = "Proportional manus length")


### Mech Adv at the tip

plot( mean.tip ~ prop.manus, data = uca.short,
      cex = 2.5 , pch = 21, bg = c("black","grey")[as.numeric(as.factor(uca.short$mating))],
      las = 1, cex.axis = 1.5, cex.lab = 1.5,
      bty = 'l',
      ylab = "Mean mechanical advantage (tip)",
      xlab = "Proportional manus length")

### Mech Adv at the tub

plot( mean.tub ~ prop.manus, data = uca.short,
      cex = 2.5 , pch = 21, bg = c("black","grey")[as.numeric(as.factor(uca.short$mating))],
      las = 1, cex.axis = 1.5, cex.lab = 1.5,
      bty = 'l',
      ylab = "Mean mechanical advatange (tubercle)",
      xlab = "Proportional manus length")

#### Intraspecific variation ####

# Since we had measures of intraspecific variation (standard variations), we decided to run PGLS models on 
# the SDs as well. Our goal was to test if the surface mating system had more variation than the burrow 
# mating system. To do this, we used the SD of one variable (e.g., carapace width, claw length) as the 
# dependent variable and the mating system as the independent factor. We ran the same three evolutionary models 
# without measurement error this time (because we are looking at the measurement error itself). The analyses are 
# as shown below.

## Separating the groups to make plotting easier later on

burrow.mean = uca.short[uca.short$mating == 'burrow',]
surface.mean = uca.short[uca.short$mating == 'surface',]

### Carapace width ####

#### Carapace SD ####

sd1_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd1_BM.rds")
sd1_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd1_OUf.rds")
sd1_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd1_OUr.rds")

#sd1_BM = phylolm(sd.carapace ~ mating , boot = 1000, 
#            data=uca.short, phy=phy, model='BM')

#sd1_OUf = phylolm(sd.carapace ~ mating, boot = 1000, 
#                data=uca.short, phy=phy, model='OUfixedRoot')

#sd1_OUr =  phylolm(sd.carapace ~ mating, boot = 1000, 
#                data=uca.short, phy=phy, model='OUrandomRoot',
#                upper.bound = 300)

compare_performance(sd1_BM, sd1_OUf, sd1_OUr)

summary(sd1_OUf)


#### Claw length ####

sd2_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd2_BM.rds")
sd2_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd2_OUf.rds")
sd2_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd2_OUr.rds")

#sd2_BM = phylolm(sd.claw ~ mating , boot = 1000, 
#                 data=uca.short, phy=phy, model='BM')

#sd2_OUf = phylolm(sd.claw ~ mating, boot = 1000, 
#                  data=uca.short, phy=phy, model='OUfixedRoot',
#                  upper.bound = 10)

#sd2_OUr =  phylolm(sd.claw ~ mating, boot = 1000, 
#                   data=uca.short, phy=phy, model='OUrandomRoot',
#                   upper.bound = 20)


compare_performance(sd2_BM, sd2_OUf, sd2_OUr) 

summary(sd2_OUf)


#### Dactyl length ####

sd3_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd3_BM.rds")
sd3_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd3_OUf.rds")
sd3_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd3_OUr.rds")

#sd3_BM = phylolm(sd.dactyl ~ mating , boot = 1000, 
#                 data=uca.short, phy=phy, model='BM')

#sd3_OUf = phylolm(sd.dactyl ~ mating, boot = 1000, 
#                  data=uca.short, phy=phy, model='OUfixedRoot',
#                  upper.bound = 10)

#sd3_OUr =  phylolm(sd.dactyl ~ mating, boot = 1000, 
#                   data=uca.short, phy=phy, model='OUrandomRoot',
#                   upper.bound = 20)


compare_performance(sd3_BM, sd3_OUf, sd3_OUr) 

summary(sd3_OUf)

#### Manus length ####

sd4_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd4_BM.rds")
sd4_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd4_OUf.rds")
sd4_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd4_OUr.rds")

#sd4_BM = phylolm(sd.manus ~ mating , boot = 1000, 
#                 data=uca.short, phy=phy, model='BM')

#sd4_OUf = phylolm(sd.manus ~ mating, boot = 1000, 
#                  data=uca.short, phy=phy, model='OUfixedRoot')

#sd4_OUr =  phylolm(sd.manus ~ mating, boot = 1000, 
#                   data=uca.short, phy=phy, model='OUrandomRoot')

compare_performance(sd4_BM, sd4_OUf, sd4_OUr) 

summary(sd4_OUf)


#### Mechanical advantage at the tip ####

sd5_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd5_BM.rds")
sd5_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd5_OUf.rds")
sd5_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd5_OUr.rds")

#sd5_BM = phylolm(sd.tip ~ mating , boot = 1000, 
#                 data=uca.short, phy=phy, model='BM')

#sd5_OUf = phylolm(sd.tip ~ mating, boot = 1000, 
#                  data=uca.short, phy=phy, model='OUfixedRoot')

#sd5_OUr =  phylolm(sd.tip ~ mating, boot = 1000, 
#                   data=uca.short, phy=phy, model='OUrandomRoot')

compare_performance(sd5_BM, sd5_OUf, sd5_OUr) 

summary(sd5_OUf)


#### Mechanical advantage at the tubercle ####

sd6_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd6_BM.rds")
sd6_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd6_OUf.rds")
sd6_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/sd6_OUr.rds")

#sd6_BM = phylolm(sd.tub ~ mating , boot = 1000, 
#                 data=uca.short, phy=phy, model='BM')

#sd6_OUf = phylolm(sd.tub ~ mating, boot = 1000, 
#                  data=uca.short, phy=phy, model='OUfixedRoot',
#                  upper.bound = 20)

#sd6_OUr =  phylolm(sd.tub ~ mating, boot = 10, 
#                   data=uca.short, phy=phy, model='OUrandomRoot',
#                   upper.bound = 1000) ## it's not estimating well. not going to use it


compare_performance(sd6_BM, sd6_OUf, sd6_OUr) 

summary(sd6_OUf)

#### Intraspecific variation graphs ####

```{r, fig.width = 8, fig.height = 12}
par(bty = 'l', mfrow = c(3,2), cex.lab = 1.5, cex.axis = 1.3,
    mar = c(5.1,7.1,4.1,2.1), mgp = c(4,1,0))

plot(density(burrow.mean$sd.carapace, na.rm = T), 
     main = "", las = 1, bty = 'l',
     xlab = "Standard deviation of carapace width (cm)",
     xlim = c(-0.2, 0.8))
polygon( density(burrow.mean$sd.carapace, na.rm = T), col = alpha( "black", 0.8), border = "black" )
polygon( density(surface.mean$sd.carapace, na.rm = T),
         col = alpha('grey', 0.5) , border = alpha('grey', 0.7) )
abline(v = coef(sd1_OUf)[1], col = 'red', lwd = 2)
abline(v = coef(sd1_OUf)[1] + coef(sd1_OUf)[2], col = 'grey', lwd = 3, lty = 2)
legend("topleft", legend = "A",bty = 'n', cex = 2)

plot(density(burrow.mean$sd.claw, na.rm = T), 
     main = "", las = 1, bty = 'l',
     xlab = "Standard deviation of proportional claw length",
     ylim = c( 0, 5 ))
polygon( density(burrow.mean$sd.claw, na.rm = T), col = "black", border = 'black' )
polygon( density(surface.mean$sd.claw, na.rm = T),
         col = alpha('grey', 0.5) , border = alpha('grey', 0.5) )
abline(v = coef(sd2_OUf)[1], col = 'red', lwd = 2)
legend("topleft", legend = "B",bty = 'n', cex = 2)


plot(density(burrow.mean$sd.dactyl, na.rm = T), 
     main = "", las = 1, bty = 'l',
     xlab = "Standard deviation of proportional dactyl length",
     ylim = c(0,6),
     xlim = c(-0.1,0.5))
polygon( density(burrow.mean$sd.dactyl, na.rm = T), col = "black", border = 'black' )
polygon( density(surface.mean$sd.dactyl, na.rm = T),
         col = alpha('grey', 0.5) , border = alpha('grey', 0.5) )
abline(v = coef(sd3_OUf)[1], col = 'red', lwd = 2)
legend("topleft", legend = "C",bty = 'n', cex = 2)

plot(density(burrow.mean$sd.manus, na.rm = T), 
     main = "", las = 1, bty = 'l',
     xlab = "Standard deviation of proportional manus length",
     xlim = c(-0.05,0.2))
polygon( density(burrow.mean$sd.manus, na.rm = T) , col = "black", border = 'black' )
polygon( density(surface.mean$sd.manus, na.rm = T) , 
         col = alpha('grey', 0.5) , border = alpha('grey', 0.5) )
abline( v = coef( sd4_OUf )[1] , col = 'red' , lwd = 2)
legend("topleft", legend = "D",bty = 'n', cex = 2)

plot(density(burrow.mean$sd.tip, na.rm = T), 
     main = "", las = 1, bty = 'l',
     xlab = "Standard deviation of mechanical 
     advantage at the tip",
     xlim = c( -0.02, 0.1))
polygon( density(burrow.mean$sd.tip, na.rm = T), col = "black", border = 'black' )
polygon( density(surface.mean$sd.tip, na.rm = T),
         col = alpha('grey', 0.5) , border = alpha('grey', 0.5) )
abline( v = coef(sd5_OUf)[1] , col = 'red' , lwd = 2)
legend("topleft", legend = "E",bty = 'n', cex = 2)


plot(density(burrow.mean$sd.tub, na.rm = T), 
     main = "", las = 1, bty = 'l',
     xlab = "Standard deviation of mechanical 
     advantage at the tubercle",
     xlim = c(-0.1, 0.3) )
polygon( density(burrow.mean$sd.tub, na.rm = T), col = "black", border = 'black' )
polygon( density(surface.mean$sd.tub, na.rm = T),
         col = alpha('grey', 0.5) , border = alpha('grey', 0.5) )

abline(v = coef(sd6_OUf)[1], col = 'red', lwd = 2)
abline(v = coef(sd6_OUf)[1] + coef(sd6_OUf)[2], col = 'grey', lwd = 3, lty = 2)
legend("topleft", legend = "F",bty = 'n', cex = 2)
 
dev.off() #to remove par configurations

#### Ancestral state reconstruction ####

# We calculated the ancestral state of both the mating system and proportional claw length. 
# For the mating system, we ran the "Equal rates" model and the "All rates differ" model and compared the 
# best fitted model using AIC. For the continuous variable (proportional claw length), we also used ace() 
# but simulated trait evolution using a Brownian motion model of evolution for each. However, to improve 
# the amount of information in our tests, we performed the ancestral state reconstructions with all the species 
# in the phylogeny. For the species we could not determine the mating system, we calculated the probability of 
# them being in each state using the "make.simmap()" function of the phytools package. Further, these are simply 
# visual aids and not tests. Thus, we did not perform any comparison between different evolutionary models 
# for the continuous variable.


#### ANC RECON #####

uca.full = read.csv("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/medidas_full.csv" , h = T, sep = ';')

uca.mate.full = read.csv("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/mating_short_all_spp.csv", h = T, sep = ';')

phy2 = read.nexus("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/combined_50_burn_MCC_2.nex")

uca.mate.full = uca.mate.full %>% arrange (sp)

rownames(uca.mate.full) = uca.mate.full$sp

uca.sim = uca.mate.full[,-c(1)]

sim_ard = make.simmap(tree = phy2, x = as.matrix(uca.sim), model = "ARD", nsim = 1000)

sim_er = make.simmap(tree = phy2, x = as.matrix(uca.sim), model = "ER", nsim = 1000)

### Model comparison

AIC.er = (-2 * sim_er[[1]]$logL) + (2 * 1)

AIC.ard = (-2 * sim_ard[[1]]$logL) + (2 * 2)

AIC.er
AIC.ard

AIC.er - AIC.ard

### The ER model either had the lowest AIC.

## Check the results

sum_sim = describe.simmap(sim_er)

## We now run ace for continuous variables. It does this by running a BM model of evolution on our data. We will run this on the claw to carapace ratio and the mechanical advantage.


uca.short2 = uca.full %>% group_by(sp) %>%
  summarise(
    mean.carapace = mean( carapace, na.rm = T),
    mean.claw = mean( claw_size, na.rm=T)
  )

uca.short2 = as.data.frame(uca.short2)

rownames(uca.short2) = uca.short2$sp

name.check(phy, uca.short2)

uca.short2$claw.ratio = uca.short2$mean.claw/uca.short2$mean.carapace
uca.short2$claw.ratio

claw.r = setNames(uca.short2$claw.ratio, rownames(uca.short2)) 

contC = contMap(phy2, claw.r, plot=FALSE)

ace.ratio = ace(claw.r, phy2, type="continuous")


### Preparing for plotting

ace_for_plot = sum_sim$ace[-c(45:90),]

col = viridis(200,alpha=1,option="D",direction=-1)
cols = setNames(c("black","grey"), c("burrow", "surface")) 

## Now, plotting both trees combined. We merged them on photoshop later.

plot(contC, lwd=5, ftype="i", xlim=c(0,65), leg.txt='Claw to carapace ratio', res = 300,
     offset = 1.6)

nodelabels(node = 1:phy2$Nnode+Ntip(phy), pie = ace_for_plot, # ACE do mating
           piecol = cols, cex= .85)

tiplabels(pie = sum_sim$tips, offset = 1.5,
          piecol = cols, cex=0.65)


text("Mating system",x=44,y=-1,cex=1,font=1, pos=1)
text(c('Burrow','Surface'), x=42,y=c(-3,-4),cex=.9,font=1,pos=4)
points(x=41.5,y=-3,pch=21,bg='black',cex=2.5)
points(x=41.5,y=-4,pch=21,cex=2.5, bg = 'grey')

# One additional test and plot regarding the proportional claw length. Does the mating system influence 
# proportional claw length?
  
uca.short$prop.claw = uca.short$mean.claw / uca.short$mean.carapace


prop_BM = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/prop_BM.rds")
prop_OUf = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/prop_OUf.rds")
prop_OUr = readRDS("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/prop_OUr.rds")

### These are the models
#prop_BM = phylolm(prop.claw ~ mating, data = uca.short, phy = phy, model = "BM",
#                  measurement_error = T, boot = 1000)
#prop_OUf = phylolm(prop.claw ~ mating, data = uca.short, phy = phy, model = "OUfixedRoot",
#                  upper.bound = 300, measurement_error = T, boot = 1000)
#prop_OUr = phylolm(prop.claw ~ mating, data = uca.short, phy = phy, model = "OUrandomRoot",
#                  upper.bound = 300, measurement_error = T, boot = 1000)
###


compare_performance(prop_BM, prop_OUf, prop_OUr)

summary(prop_BM)

par(bty = 'l',  cex.lab = 1.5, cex.axis = 1.3)
vioplot(prop.claw ~ as.factor(mating), data = uca.short,
        las = 1,  names = c("Burrow", "Surface"),
        ylab = "Claw length divided by carapace width",
        xlab = 'Mating system',
        col = c('black', 'grey'),
        rectCol = 'white', lineCol = 'grey30',
        colMed = 'grey30', cex = 2,
        areaEqual = T)
dev.off() # to remove par() configurations

## Ratematrix

# We used the ratematrix package to calculate the rates of evolution of claw traits in each mating systems. 
# The package also calculated the evolutionary covariance and evolutionary correlation of the traits in each 
# mating system. After all of this was calculated, we used the functions implemented in Slater & Friscia (2019) 
# to test if the evolutionary covariation of the traits differed in each mating system. We performed three analyses: 
# one for claw evolutionary allometry, another for the claw lever and muscle, and another for claw mechanical 
# advantage and relative manus length (the same models we ran PGLS regressions for). 

# For all models, we used a custom flat prior that encompassed the distribution of our data. 

#### CLAW MANUS AND DACTYL ####

### PREPARING THE DATA FOR RATEMATRIX ###

mat = setNames(uca.short$mating, rownames(uca.short)) 

resp.data = uca.short[c("prop.out2","prop.manus")]

resp.scale = scale(resp.data,scale = F)

colsb = setNames(c("black","gray"),levels(as.factor(mat)))

phy.sim = make.simmap(phy, mat, model="ER")
class(phy.sim)

#### BUILDING A FLAT PRIOR ####

par.mu = rbind( c(-2,2), c(-2,2) )
par.sd = rbind( c(0,2), c(0,2) )
prior = makePrior(r=2, p=2, den.mu = "unif", par.mu = par.mu, den.sd = "unif", par.sd=par.sd)
sample.prior = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)

# After building the prior, the next step is running the model to get the mcmc chain for the analyses
# We are running four independent chains, each with a different prior sampling. After we run each chain, 
# we analyze the run, making sure it converged and run the next one.

##### MODELS - DO NOT RUN #####
### DO NOT RUN UNLESS YOU WANT NEW MCMC CHAINS ###
### THE MCMC CHAINS GENERATED ARE LOADED AFTER THIS PART ###

## First MCMC chain

#handle = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                        v = 50, prior = prior, start = sample.prior
#)

#logAnalyzer(handle, burn = 0.3, thin = 1)

#mcmc = readMCMC(handle,burn = 0.3, thin = 1)

### Reading the mcmc chains to make sure everything is in order and the roots have been sampled appropriately

#plotRatematrix(mcmc, set.xlim = c(-0.02,0.3))
#plotRootValue(mcmc)

## Second MCMC chain

#sample.prior2 = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle2 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                         v = 50, prior = prior, start = sample.prior2
#)

#logAnalyzer(handle2, burn = 0.3, thin = 1)

#mcmc2 = readMCMC(handle2,burn = 0.3, thin = 1)

#plotRatematrix(mcmc2)
#plotRootValue(mcmc2)


## Third MCMC chain

#sample.prior3 = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle3 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                         v = 50,  prior = prior, start = sample.prior3
#)

#logAnalyzer(handle3, burn = 0.3, thin = 1)

#mcmc3 = readMCMC(handle3,burn = 0.3, thin = 1)

#plotRatematrix(mcmc3)
#plotRootValue(mcmc3)

## Fourth MCMC chain

#sample.prior4= samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle4 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                         v = 50, prior = prior, start = sample.prior4
#)

#logAnalyzer(handle4, burn = 0.3, thin = 1)
#
#mcmc4 = readMCMC(handle4,burn = 0.3, thin = 1)

#plotRatematrix(mcmc4)
#plotRootValue(mcmc4)

# Now, saving all mcmc chains in the same file so we don't have to run it every time.

#save(mcmc, mcmc2, mcmc3, mcmc4, file = "data/mcmc_dactyl.Rdata")

##### LOAD MCMC CHAINS GENERATED! ####

load(file = "D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/mcmc_dactyl.Rdata")


# Measuring the convergence of the four chains
Rfactor = checkConvergence(mcmc, mcmc2, mcmc3, mcmc4)
Rfactor

# Merging them all in the same object
merged_mcmc = mergePosterior(mcmc, mcmc2, mcmc3, mcmc4)

# Plotting the results. 

# grey = surface mating 
# black = burrow mating

plotRatematrix(merged_mcmc, set.xlim = c(-0.003,0.03), alphaOff = 0.7, alphaEll = 0.7,
               colors=c("grey","black")[c(2,1)], n.lines = 100, point.wd = 1,
               set.leg = c("Proportional dactyl length", "Proportional manus length"))

# This graph requires some explaining. The plots on the main diagonal show the rate of evolution of one trait in each mating system. The spread on the surface mating (grey) suggests that there is more uncertainty on the estimation of the rate of evolution of that regime when compared to the burrow mating system. The plot on the off-diagonal (the one that is a little bit transparent) show the evolutionary covariance of the traits in each mating system. The ellipsis shows the 95% confidence interval of each bivariate distribution. The size, direction, and angle of the ellipsis shows how these two traits are covarying.    

plotRootValue(merged_mcmc)

# This second graph shows how the posterior distribution of the values retrieved from the root of the tree.


### Now, we test how much the two mating systems overlap with one another

# This first test shows the probability of overlap between the posterior distributions. With the small values denoting that there is very low probability of the two posteriors to overlap.
testRatematrix(chain = merged_mcmc, par = "all")

# This test shows the overlap of the correlation between covariance matrices. 
testRatematrix(chain = merged_mcmc, par = "correlation")

# To know the correlations, we need to first extract these correlations and then plot them. 

corr_dact_manus = extractCorrelation(merged_mcmc)

hist(x = corr_dact_manus[,1], xlim = c(-1,1), main = "", col = "black"
     , border = "white", breaks = 20, freq = FALSE)
hist(x = corr_dact_manus[,2], xlim = c(-1,1), main = "", col = "grey"
     , border = "white", breaks = 20, freq = FALSE,add=T)

boxplot(x = corr_dact_manus)

## None of these, however, are formal statistical tests.
## That is why we use the same framework in Slater & Friscia (2019) to test if the covariance matrices differ between the mating systems. The analysis is based on Ovaskainen et al. (2008) work, but basically it calculates pairwise differences between the covariance matrices.
## Below, we provide the same results.

complete.dactyl.post <- merged_mcmc$matrix

### Ovaskainen Test ###
## If you remove the # from the code, it will save the results in a txt file

#zz <- file("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/dactyl_OvaskainenTest_all.txt", open = "wt")
#sink(zz)
ovd.test <- pairwise.Ovaskainen(complete.dactyl.post)
lapply(ovd.test, quantile, c(0.025, 0.5, 0.975))
#sink()


#### CLAW LENGTH AND CARAPACE WIDTH ####

# Now we will do the same as before, but for the evolutionary allometry

### PREPARING THE DATA FOR RATEMATRIX ###

mat = setNames(uca.short$mating, rownames(uca.short)) 

resp.data = uca.short[c("mean.claw","mean.carapace")]

resp.scale = scale(resp.data,scale = F)

colsb = setNames(c("black","gray"),levels(as.factor(mat)))

phy.sim = make.simmap(phy, mat, model="ER")
class(phy.sim)

#### BUILDING A FLAT PRIOR ####

par.mu <- rbind( c(-5,5), c(-5,5) )
par.sd <- rbind( c(0,3), c(0,3) )
prior <- makePrior(r=2, p=2, den.mu = "unif", par.mu = par.mu, den.sd = "unif", par.sd = par.sd)
sample.prior = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)


##### MODELS - DO NOT RUN #####
### DO NOT RUN UNLESS YOU WANT NEW MCMC CHAINS ###
### THE MCMC CHAINS GENERATED ARE LOADED AFTER THIS PART ###

## First MCMC chain

#handle = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 10000000, 
#                        v = 50, prior = prior, start = sample.prior
#)

#logAnalyzer(handle, burn = 0.3, thin = 1)

#mcmc = readMCMC(handle,burn = 0.3, thin = 1)

#plotRatematrix(mcmc, set.xlim = c(-0.02,0.3))
#plotRootValue(mcmc)

## Second MCMC chain

#sample.prior2 = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle2 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 10000000, 
#                         v = 50, prior = prior, start = sample.prior2
#)

#logAnalyzer(handle2, burn = 0.3, thin = 1)

#mcmc2 = readMCMC(handle2,burn = 0.3, thin = 1)

#plotRatematrix(mcmc2)
#plotRootValue(mcmc2)

## Third MCMC chain

#sample.prior3 = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle3 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 10000000, 
#                         v = 50,  prior = prior, start = sample.prior3
#)

#logAnalyzer(handle3, burn = 0.3, thin = 1)

#mcmc3 = readMCMC(handle3,burn = 0.3, thin = 1)

#plotRatematrix(mcmc3)
#plotRootValue(mcmc3)

## Fourth MCMC chain

#sample.prior4= samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle4 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 10000000, 
#                         v = 50, prior = prior, start = sample.prior4
#)

#logAnalyzer(handle4, burn = 0.3, thin = 1)

#mcmc4 = readMCMC(handle4,burn = 0.3, thin = 1)

#plotRatematrix(mcmc4)
#plotRootValue(mcmc4)

#save(mcmc, mcmc2, mcmc3, mcmc4, file = "data/mcmc_claw.Rdata")

##### LOAD MCMC CHAINS GENERATED! ####

load(file = "D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/mcmc_claw.Rdata")

# Measuring the convergence of the four chains
Rfactor = checkConvergence(mcmc, mcmc2, mcmc3, mcmc4)
Rfactor

# Merging them all in the same objetc
merged_mcmc = mergePosterior(mcmc, mcmc2, mcmc3, mcmc4)

# Plotting the results. 

# grey = surface mating 
# black = burrow mating

plotRatematrix(merged_mcmc, set.xlim = c(0,0.6), alphaOff = 0.7, alphaEll = 0.7,
               colors=c("grey","black")[c(2,1)], n.lines = 100, point.wd = 1,
               set.leg = c("Mean claw length (scaled)", "Mean carapace width (scaled)"))

# This graph requires some explaining. The plots on the main diagonal show the rate of evolution of one trait in each mating system. The spread on the surface mating (grey) suggests that there is more uncertainty on the estimation of the rate of evolution of that regime when compared to the burrow mating system. The plot on the off-diagonal (the one that is a little bit transparent) show the evolutionary covariance of the traits in each mating system. The ellipsis shows the 95% confidence interval of each bivariate distribution. The size, direction, and angle of the ellipsis shows how these two traits are covarying.    

plotRootValue(merged_mcmc)

# This second graph shows how the posterior distribution of the values retrieved from the root of the tree.

### Now, we test how much the two mating systems overlap with one another

# This first test shows the probability of overlap between the posterior distributions. With the small values denoting that there is very low probability of the two posteriors to overlap.
testRatematrix(chain = merged_mcmc, par = "all")

# This test shows the overlap of the correlation between covariance matrices. 
testRatematrix(chain = merged_mcmc, par = "correlation")

# To know the correlations, we need to first extract these correlations and then plot them. 

corr_claw_carapace = extractCorrelation(merged_mcmc)

hist(x = corr_claw_carapace[,1], xlim = c(-1,1), main = "", col = "black"
     , border = "white", breaks = 20, freq = FALSE)
hist(x = corr_claw_carapace[,2], xlim = c(-1,1), main = "", col = "grey"
     , border = "white", breaks = 20, freq = FALSE,add=T)

boxplot(x = corr_claw_carapace)

## None of these, however, are formal statistical tests.
## That is why we use the same framework in Slater & Friscia (2019) to test if the covariance matrices differ between the mating systems. The analysis is based on Ovaskainen et al. (2008) work, but basically it calculates pairwise differences between the covariance matrices.
## Below, we provide the same results.

complete.claw.post <- merged_mcmc$matrix

### Ovaskainen Test ###
## If you remove the # from the code, it will save the results in a txt file

#zz <- file("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/claw_OvaskainenTest_all.txt", open = "wt")
#sink(zz)
ovd.test <- pairwise.Ovaskainen(complete.claw.post)
lapply(ovd.test, quantile, c(0.025, 0.5, 0.975))
#sink()

```

#### MECHANICAL ADVANTAGE AND RELATIVE MANUS LENGTH #####
### PREPARING THE DATA FOR RATEMATRIX ###

mat = setNames(uca.short$mating, rownames(uca.short)) 

resp.data = uca.short[c("mean.tip","prop.manus")]

resp.scale = scale(resp.data,scale = F)

colsb = setNames(c("black","gray"),levels(as.factor(mat)))

phy.sim = make.simmap(phy, mat, model="ER")
class(phy.sim)

#### BUILDING A FLAT PRIOR ####

par.mu = rbind( c(-1.5,1.5), c(-1.5,1.5) )
par.sd = rbind( c(0,1), c(0,1) )
prior = makePrior(r=2, p=2, den.mu = "unif", par.mu=par.mu,den.sd = "unif", par.sd=par.sd)
sample.prior = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)


##### MODELS - DO NOT RUN #####
### DO NOT RUN UNLESS YOU WANT NEW MCMC CHAINS ###
### THE MCMC CHAINS GENERATED ARE LOADED AFTER THIS PART ###

## First MCMC chain

#handle = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                        v = 50, prior = prior, start = sample.prior
#)

#logAnalyzer(handle, burn = 0.3, thin = 1)

#mcmc = readMCMC(handle,burn = 0.3, thin = 1)

#plotRatematrix(mcmc)
#plotRootValue(mcmc)

## Second MCMC chain

#sample.prior2 = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle2 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                         v = 50, prior = prior, start = sample.prior2
#)

#logAnalyzer(handle2, burn = 0.3, thin = 1)

#mcmc2 = readMCMC(handle2,burn = 0.3, thin = 1)

#plotRatematrix(mcmc2)
#plotRootValue(mcmc2)

## Third MCMC chain

#sample.prior3 = samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle3 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                         v = 50,  prior = prior, start = sample.prior3
#)

#logAnalyzer(handle3, burn = 0.3, thin = 1)

#mcmc3 = readMCMC(handle3,burn = 0.3, thin = 1)

#plotRatematrix(mcmc3)
#plotRootValue(mcmc3)

## Fourth MCMC chain

#sample.prior4= samplePrior(n=1, prior = prior, sample.sd = T, rebuild.R = F)
#handle4 = ratematrixMCMC(data = resp.scale, phy = phy.sim, dir = '/result_mcmc',gen = 1e+07, 
#                         v = 50, prior = prior, start = sample.prior4
#)

#logAnalyzer(handle4, burn = 0.3, thin = 1)

#mcmc4 = readMCMC(handle4,burn = 0.3, thin = 1)

#plotRatematrix(mcmc4)
#plotRootValue(mcmc4)

#save(mcmc, mcmc2, mcmc3, mcmc4, file = "data/mcmc_ratio.Rdata")

##### LOAD MCMC CHAINS GENERATED! ####

load(file = "D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/mcmc_ratio.Rdata")

# Measuring the convergence of the four chains
Rfactor = checkConvergence(mcmc, mcmc2, mcmc3, mcmc4)
Rfactor

# Merging them all in the same objetc
merged_mcmc = mergePosterior(mcmc, mcmc2, mcmc3, mcmc4)

# Plotting the results. 

# grey = surface mating 
# black = burrow mating

plotRatematrix(merged_mcmc, set.xlim = c(-0.00015,0.003), alphaOff = 0.7, alphaEll = 0.7,
               colors=c("grey","black")[c(2,1)], n.lines = 100, point.wd = 1,
               set.leg = c("Mean mechanical advantage", "Mean proportional manus length"))

# This graph requires some explaining. The plots on the main diagonal show the rate of evolution of one trait in each mating system. The spread on the surface mating (grey) suggests that there is more uncertainty on the estimation of the rate of evolution of that regime when compared to the burrow mating system. The plot on the off-diagonal (the one that is a little bit transparent) show the evolutionary covariance of the traits in each mating system. The ellipsis shows the 95% confidence interval of each bivariate distribution. The size, direction, and angle of the ellipsis shows how these two traits are covarying.    

plotRootValue(merged_mcmc)

# This second graph shows how the posterior distribution of the values retrieved from the root of the tree.

### Now, we test how much the two mating systems overlap with one another

# This first test shows the probability of overlap between the posterior distributions. With the small values denoting that there is very low probability of the two posteriors to overlap.
testRatematrix(chain = merged_mcmc, par = "all")

# This test shows the overlap of the correlation between covariance matrices. 
testRatematrix(chain = merged_mcmc, par = "correlation")

# To know the correlations, we need to first extract these correlations and then plot them. 

corr_tip_ratio = extractCorrelation(merged_mcmc)

hist(x = corr_tip_ratio[,1], xlim = c(-1,1), main = "", col = "black"
     , border = "white", breaks = 20, freq = FALSE)
hist(x = corr_tip_ratio[,2], xlim = c(-1,1), main = "", col = "grey"
     , border = "white", breaks = 20, freq = FALSE,add=T)

boxplot(x = corr_tip_ratio)

## None of these, however, are formal statistical tests.
## That is why we use the same framework in Slater & Friscia (2019) to test if the covariance matrices differ between the mating systems. The analysis is based on Ovaskainen et al. (2008) work, but basically it calculates pairwise differences between the covariance matrices.
## Below, we provide the same results.

complete.ratio.post <- merged_mcmc$matrix

### Ovaskainen Test ###
## If you remove the # from the code, it will save the results in a txt file

#zz <- file("D:/OneDrive/Orientados/Mestrado - Jonatas/Mestrado/new_paper/claw.evol/data/ratio_OvaskainenTest_all.txt", open = "wt")
#sink(zz)
ovd.test = pairwise.Ovaskainen(complete.ratio.post)
lapply(ovd.test, quantile, c(0.025, 0.5, 0.975))
#sink()


### EVOLUTIONARY CORRELATION GRAPHS

par(bty = 'l', cex.lab = 1.5, cex.axis = 1.3,
    mar = c(5.1,7.1,4.1,2.1), mgp = c(4,1,0))


hist(x = corr_claw_carapace[,1], xlim = c(-1,1), main = "", col = "black"
     , border = "white", breaks = 20, freq = FALSE,
     ylim = c(0,10),
     xlab = "Evolutionary correlation between carapace width
     and claw length"
)
abline( v = median(corr_claw_carapace[,1]), lwd = 2, col = 'red', lty = 2)

hist(x = corr_claw_carapace[,2], xlim = c(-1,1), main = "", col = alpha("grey", 0.8)
     , border = "white", breaks = 20, freq = FALSE,add=T)
abline( v = median(corr_claw_carapace[,2]), lwd = 2, col = 'red', lty = 2)

legend(x = -1, y = 8, 
       legend = c("Surface", "Burrow"),
       pt.bg = c("grey", "black"), 
       pch = 22,
       bty = "n",
       pt.cex = 2,
       title = "Mating System",
       title.cex = 1.5,
       cex = 1.2)
legend("topleft", legend = "A", bty = 'n', cex = 3)



hist(x = corr_dact_manus[,1], xlim = c(-1,1), main = "", col = alpha("black", 1),
     border = "white", breaks = 20, freq = FALSE,
     xlab = "Evolutionary correlation between proportional
     manus and porportional dactyl length")
abline( v = median(corr_dact_manus[,1]), lwd = 2, col = 'red', lty = 2)

hist(x = corr_dact_manus[,2], xlim = c(-1,1), main = "", col = alpha("grey",0.8)
     , border = "white", breaks = 20, freq = FALSE,add=T)
abline( v = median(corr_dact_manus[,2]), lwd = 2, col = 'red', lty = 2)

legend( "topleft", legend = "B", bty = 'n', cex = 3)




hist(x = corr_tip_ratio[,1], main = "", col = "black"
     , border = "white", breaks = 20, freq = FALSE,
     xlim = c(-1,1), ylim = c(0,2.5),
     xlab = "Evolutionary correlation between MA (tip) 
     and proportional manus length")
abline( v = median(corr_tip_ratio[,1]), lwd = 2, col = 'red', lty = 2)

hist(x = corr_tip_ratio[,2], main = "", col = alpha("grey", 0.8)
     , border = "white", breaks = 20, freq = FALSE,add=T)
abline( v = median(corr_tip_ratio[,2]), lwd = 2, col = 'red', lty = 2)


legend("topleft", legend = "C", bty = 'n', cex = 3)

dev.off() # Run this command if you want to remove the par() configuration

## Covariance matrices with intraspecific variation

# One potential source of noise in any evolutionary analyses is intraspecific variation. However, the ratematrix approach we used cannot account for variation in its current formulation. Thus, we used a different package to do that. As we explain thoroughly in our Supplementary Files, we will use the "mvMORPH" package, which allows us to calculate evolutionary covariance matrices with intraspecific variation. We will run two models, one with intraspecific variation and another without intraspecific variation and compare them. Similar results should suggest that intraspecific variation does not influence our results much.

# First, we need to build the data and the standard error matrices. 
# We will use only claw width and carapace length for this.

# Getting the data and the sd we calculated previously
data_mvBM = uca.short[,c(2,3)]
error = uca.short[,c(11,12)]

# Now, transforming it into standard error
error_sq = (error/sqrt(uca.mate$n))^2

# We have some NAs in our sample because we only had one or two individuals of some species. Thus, we added an error compartible with the rest of our sample
error_sq[is.na(error_sq)] <- 10^-3

# Turn the data into matrix
data_matrix = as.matrix(data_mvBM)
error_matrix = as.matrix(error_sq)


# Now, we need the mating systems AND an ancestral reconstruction tree to make this work
mat = setNames(uca.short$mating, rownames(uca.short)) 
sim = make.simmap(phy, mat, model="ER", nsim=100)

# Running the model with intraspecific variation first
t1 = mvBM(sim[[1]], data = data_matrix, error = error_matrix,
          model = "BMM")

# Now, the same model but without intraspecific variation
t1_n = mvBM(sim[[1]], data = data_matrix, model = "BMM")

### This is what we obtained ###
## WITH INTRASPECIFIC VARIATION ##
# Estimated rate matrix 
# ______________________ 
#  , , burrow
#
#                mean.carapace   mean.claw
# mean.carapace    0.01057907   0.02201515
# mean.claw        0.02201515   0.05179114
#
# , , surface
#
#                mean.carapace    mean.claw
# mean.carapace    0.05757604    0.1180860
# mean.claw        0.11808598    0.2464731



## WITHOUT INTRASPECIFIC VARIATION
# Estimated rate matrix 
# ______________________ 
# , , burrow
#
#                mean.carapace   mean.claw
# mean.carapace    0.01024554   0.02138044
# mean.claw        0.02138044   0.05369690
#
# , , surface
#
#                mean.carapace   mean.claw
# mean.carapace    0.06070039   0.1174544
# mean.claw        0.11745444   0.2502149


# Now, for the sake of the argument, let's check the evolutionary correlations between carapace and claw. First, we will compare the burrow mating system with and without intraspecific variation.

# BURROW
cov2cor(t1$sigma[,,1]) 
cov2cor(t1_n$sigma[,,1])

## WITH ##
#                mean.carapace   mean.claw
# mean.carapace     1.0000000    0.9368431
# mean.claw         0.9368431    1.0000000

## WITHOUT ##
#                mean.carapace   mean.claw
# mean.carapace     1.0000000    0.9097224
# mean.claw         0.9097224    1.0000000

# SURFACE
cov2cor(t1$sigma[,,2]) 
cov2cor(t1_n$sigma[,,2])

## WITH ##
#                mean.carapace   mean.claw
# mean.carapace     1.0000000    0.9961823
# mean.claw         0.9961823    1.0000000

## WITHOUT ##
#                mean.carapace   mean.claw
# mean.carapace     1.0000000    0.9529771
# mean.claw         0.9529771    1.0000000


# As can be seen, the difference in those models is minimal. 
# For more information and interpretation, see our Supplementary Files.


```

## 3D visualization of claw morphospace

# Now, as a visual aid of how each mating systme is exploring the morphospace, we will plot the claw length, manus length, and dactyl length both in raw format and proportionally to body size. We will also calculate a bounding box around the two mating systems to visualize how each mating system explores the morphospace. There is no analysis, just visual representation.
# We would also like to disclose that we used AI to help with to code for the bounding boxes.

## First plot, proportional sizes

s3d = scatterplot3d(
  x = uca.short$prop.claw , 
  y = uca.short$prop.manus , 
  z = uca.short$prop.out2,
  
  pch = 21,
  bg = c("black" , "grey")[as.numeric( as.factor( uca.short$mating ) )],
  grid = F,
  cex.symbols = 2,
  type = 'h',
  angle = 80,
  lty.hplot = 2,
  xlab = "Mean proportional claw length",
  ylab = 'Mean proportional manus length',
  zlab = "Mean proportional dactyl length",
)

add_bounding_boxes = function(data, group_var) {
  mating_systems = unique(data[[group_var]])
  
  for (system in mating_systems) {
    system_data <- data[data[[group_var]] == system, ]
    
    # Calculate min and max for each dimension
    x_range = range(system_data$prop.claw)
    y_range = range(system_data$prop.manus)
    z_range = range(system_data$prop.out2)
    
    # Define the 8 corners of the bounding box
    corners = matrix(c(
      x_range[1], y_range[1], z_range[1], # bottom front left
      x_range[2], y_range[1], z_range[1], # bottom front right
      x_range[2], y_range[2], z_range[1], # bottom back right
      x_range[1], y_range[2], z_range[1], # bottom back left
      x_range[1], y_range[1], z_range[2], # top front left
      x_range[2], y_range[1], z_range[2], # top front right
      x_range[2], y_range[2], z_range[2], # top back right
      x_range[1], y_range[2], z_range[2]  # top back left
    ), ncol = 3, byrow = TRUE)
    
    # Define which corners to connect to form the 12 edges of the box
    edges = list(
      c(1,2), c(2,3), c(3,4), c(4,1), # bottom face
      c(5,6), c(6,7), c(7,8), c(8,5), # top face
      c(1,5), c(2,6), c(3,7), c(4,8)  # vertical edges
    )
    
    # Draw all edges of the bounding box
    for (edge in edges) {
      s3d$points3d(
        corners[edge, ],
        type = "l",
        col = ifelse(system == "surface", "grey", "black"),
        lty = ifelse(system == "surface", 1, 2),
        lwd = 2
      )
    }
  }
}

# Add bounding boxes to the plot
add_bounding_boxes(uca.short, "mating")

## Now the second plot with raw values.

s3d = scatterplot3d(
  x = uca.short$mean.claw , 
  y = uca.short$mean.manus , 
  z = uca.short$mean.dactyl,
  
  pch = 21,
  bg = c("black" , "grey")[as.numeric( as.factor( uca.short$mating ) )],
  grid = F,
  cex.symbols = 2,
  type = 'h',
  angle = 80,
  lty.hplot = 2,
  xlab = "Mean claw length (cm)",
  ylab = 'Mean manus length (cm)',
  zlab = "Mean  dactyl length (cm)",
)



add_bounding_boxes = function(data, group_var) {
  mating_systems = unique(data[[group_var]])
  
  for (system in mating_systems) {
    system_data = data[data[[group_var]] == system, ]
    
    # Calculate min and max for each dimension
    x_range = range(system_data$mean.claw)
    y_range = range(system_data$mean.manus)
    z_range = range(system_data$mean.dactyl)
    
    # Define the 8 corners of the bounding box
    corners = matrix(c(
      x_range[1], y_range[1], z_range[1], # bottom front left
      x_range[2], y_range[1], z_range[1], # bottom front right
      x_range[2], y_range[2], z_range[1], # bottom back right
      x_range[1], y_range[2], z_range[1], # bottom back left
      x_range[1], y_range[1], z_range[2], # top front left
      x_range[2], y_range[1], z_range[2], # top front right
      x_range[2], y_range[2], z_range[2], # top back right
      x_range[1], y_range[2], z_range[2]  # top back left
    ), ncol = 3, byrow = TRUE)
    
    # Define which corners to connect to form the 12 edges of the box
    edges = list(
      c(1,2), c(2,3), c(3,4), c(4,1), # bottom face
      c(5,6), c(6,7), c(7,8), c(8,5), # top face
      c(1,5), c(2,6), c(3,7), c(4,8)  # vertical edges
    )
    
    # Draw all edges of the bounding box
    for (edge in edges) {
      s3d$points3d(
        corners[edge, ],
        type = "l",
        col = ifelse(system == "surface", "grey", "black"),
        lty = ifelse(system == "surface", 1, 2),
        lwd = 2
      )
    }
  }
}

# Add bounding boxes to the plot
add_bounding_boxes(uca.short, "mating")

# DONE :D