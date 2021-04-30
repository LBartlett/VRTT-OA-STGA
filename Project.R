# Script by Ethan Hackmeyer with a lot of help from Lewis Bartlett
# ejh16320@uga.edu , lewis.bartlett@uga.edu
# Oxalic Acid Shop Towel study

#Reading in the 2020 Shearer Trial Data

UGAdata.2020 <- read.csv(file = 'ShearerTrialData.csv', header = T, stringsAsFactors = F)

#Creating a column for the change in percent mite increase

UGAdata.2020$DeltaPMI <- NA

#Populating DeltaPMI by taking the PMI at the end minus the PMI at the beginning of the experiment

for(E in 1:NROW(UGAdata.2020)){
  UGAdata.2020$DeltaPMI[E] <- (UGAdata.2020$MPB2[which(UGAdata.2020$Colony == UGAdata.2020$Colony[E])]
                               -
                                UGAdata.2020$MPB1[which(UGAdata.2020$Colony == UGAdata.2020$Colony[E])]
                               )
}

shapiro.test(UGAdata.2020$DeltaPMI)
hist(UGAdata.2020$DeltaPMI)

boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$TreatmentCode)
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$OxalicDose)
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$GlycerinTowel)

Model1 <- glm(
  DeltaPMI ~ Yard + TreatmentCode,
  family = 'gaussian',
  data = UGAdata.2020
)
anova(Model1, test = 'F')
summary(Model1)

hist(resid(Model1))
qqnorm(resid(Model1))
qqline(resid(Model1), col = "blue1", lwd = 2)

Model2 <- glm(
  DeltaPMI ~ Yard + OxalicDose - 1,
  family = 'gaussian',
  data = UGAdata.2020
)
anova(Model2, test = 'F')
summary(Model2)

hist(resid(Model2))
qqnorm(resid(Model2))
qqline(resid(Model2), col = "blue1", lwd = 2)

library(afex)

FullMiteMod <- mixed(DeltaPMI ~ OxalicDose + (1|Yard),
                     data = UGAdata.2020)

nice(FullMiteMod)

MirrorMod <- lmer(DeltaPMI ~ OxalicDose + (1|Yard),
                  data = UGAdata.2020)
summary(MirrorMod)


Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

ColRef <- data.frame(Treatment = levels(as.factor(UGAdata.2020$Treatment)), Col =  c('blue4','pink2','orange3','red2'))

boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$TreatmentCode,main = 'Treatment Influence on Mite Numbers', ylab = expression(paste('Percent Mite Change',sep='')),
        xlab = 'Treatment Type', border = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50), 
        cex.axis = 1.3, cex.lab = 1.5, outline = T, lwd = 1.2, outcol = 'white',
        boxlty = 1, whisklty = 1, staplelty = 1, boxwex = 0.7, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50))
stripchart(UGAdata.2020$DeltaPMI ~ UGAdata.2020$TreatmentCode,
           col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 60),
           vertical = T, add = T, pch = 16, cex = 1.2, 
           method = 'jitter', lwd = 2)

boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$OxalicDose, main = 'Oxalic Dosage', ylab = expression(paste('Percent Mite Change',sep='')),
        xlab = 'Treatment', border = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50), 
        cex.axis = 1.3, cex.lab = 1.5, outline = T, lwd = 1.2, outcol = 'white',
        boxlty = 1, whisklty = 1, staplelty = 1, boxwex = 0.7, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50))
stripchart(UGAdata.2020$DeltaPMI ~ UGAdata.2020$OxalicDose,
           col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 60),
           vertical = T, add = T, pch = 16, cex = 1.2, 
           method = 'jitter', lwd = 2)
          
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$GlycerinTowel, main = 'Glycerin Towel Presence', ylab = expression(paste('Percent Mite Change',sep='')),
        xlab = 'Towel Present', border = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50), 
        cex.axis = 1.3, cex.lab = 1.5, outline = T, lwd = 1.2, outcol = 'white',
        boxlty = 1, whisklty = 1, staplelty = 1, boxwex = 0.7, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50))
stripchart(UGAdata.2020$DeltaPMI ~ UGAdata.2020$GlycerinTowel,
           col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 60),
           vertical = T, add = T, pch = 16, cex = 1.2, 
           method = 'jitter', lwd = 2)


## Just wanted to see if the presence of the glycerin towel alone with no OA looks any different from the glycerin towel with OA         
ZeroOxalic <- UGAdata.2020[which(UGAdata.2020$OxalicDose ==0),]
hist(ZeroOxalic$DeltaPMI)
boxplot(ZeroOxalic$DeltaPMI ~ ZeroOxalic$GlycerinTowel)
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$GlycerinTowel)
## From the comparison of the boxplots above,it seems like the presence of the glycerin towel alone has slightly different results from the OA soaked glycerin towel

GlycerinModel <- glm(
  DeltaPMI ~ Yard + TreatmentCode,
  family = 'gaussian',
  data = ZeroOxalic
)
anova(GlycerinModel, test = 'F')
summary(GlycerinModel)

hist(resid(GlycerinModel))
qqnorm(resid(GlycerinModel))
qqline(resid(GlycerinModel), col = "blue1", lwd = 2)

# Results are kinda interesting, p value was .3074 but I'm not sure if I did it correctly