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

#Constructing a linear model and creating an ANOVA

Model1 <- glm(
  DeltaPMI ~ Yard + TreatmentCode,
  family = 'gaussian',
  data = UGAdata.2020
)
anova(Model1, test = 'F')

boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$TreatmentCode)
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$OxalicDose)
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$GlycerinTowel)

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
