# Analysis Script for Oxalic Acid 'Slow-Release' Trials
# Script by Ethan J Hackmeyer & Lewis J Bartlett
# ejh16320@uga.edu , lewis.bartlett@uga.edu
# Maintained on GitHUb: https://github.com/LBartlett/VRTT-OA-STGA.git

### Load Package(s)
library(afex)

## Load in Data ('Shearer' was collaborating beekeeper)
UGAdata.2020 <- read.csv(file = 'ShearerTrialData.csv', header = T, stringsAsFactors = F)

# Calculate a change in per-capita mite loads, 'Delta PMI' (change in percent mite intensity from pre- to post- treatment)
for(E in 1:NROW(UGAdata.2020)){
  UGAdata.2020$DeltaPMI[E] <- (UGAdata.2020$MPB2[which(UGAdata.2020$Colony == UGAdata.2020$Colony[E])]
                               -
                                 UGAdata.2020$MPB1[which(UGAdata.2020$Colony == UGAdata.2020$Colony[E])]
  )
}

# factorise treatments
UGAdata.2020$TreatmentCode <- factor(UGAdata.2020$TreatmentCode, levels = rev(unique(UGAdata.2020$TreatmentCode)))

# Check data frame is ready for analysis
head(UGAdata.2020)

# Quick transparency function for plotting

Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

# Colour scheme
ColRef <- data.frame(Treatment = levels(as.factor(UGAdata.2020$Treatment)), Col =  c('blue4','pink3','orange3','red2'))

# Set margings to accomodate plots

par(mar = c(5,5,4,2))

################
#### Analysis here

### Break down by hypotheses: each section has a graph [boxplot] and matching analysis

## Start by treatment code [as factors]

# Plot
boxplot(UGAdata.2020$DeltaPMI ~ UGAdata.2020$TreatmentCode, main = NA, ylab = expression(paste(Delta, ' Percent Mite Intensity',sep='')),
        xlab = 'Treatment Type', border = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50), 
        cex.axis = 1.3, cex.lab = 1.5, outline = T, lwd = 1.75, outcol = 'white',
        boxlty = 1, whisklty = 1, staplelty = 1, boxwex = 0.7, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 35))
stripchart(UGAdata.2020$DeltaPMI ~ UGAdata.2020$TreatmentCode,
           col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 50),
           vertical = T, add = T, pch = 16, cex = 1.2, 
           method = 'jitter', lwd = 2)
abline(h = 0, lty = 2)

# Analysis

TreatMiteMod <- mixed(DeltaPMI ~ TreatmentCode + (1|Yard),
                     data = UGAdata.2020)

nice(TreatMiteMod)

# No evidence of any differences across the four treatments (true control, no oxalic acid control, both oxalic doses) on change in mite levels

## Strip out the 'CTRL' level, just work on oxalic dose
UGAdata.2020.NC <- UGAdata.2020[which(UGAdata.2020$TreatmentCode != 'CTRL'),]

boxplot(UGAdata.2020.NC$DeltaPMI ~ UGAdata.2020.NC$OxalicDose, main = NA, ylab = expression(paste(Delta, ' Percent Mite Intensity',sep='')),
        xlab = 'Oxalic Acid Dose (g)', border = sapply(X = as.character(ColRef$Col[2:4]), FUN = Transpa, percent = 50), 
        cex.axis = 1.3, cex.lab = 1.5, outline = T, lwd = 1.2, outcol = 'white',
        boxlty = 1, whisklty = 1, staplelty = 1, boxwex = 0.7, 
        col = sapply(X = as.character(ColRef$Col[2:4]), FUN = Transpa, percent = 50))
stripchart(UGAdata.2020.NC$DeltaPMI ~ UGAdata.2020.NC$OxalicDose,
           col = sapply(X = as.character(ColRef$Col[2:4]), FUN = Transpa, percent = 60),
           vertical = T, add = T, pch = 16, cex = 1.2, 
           method = 'jitter', lwd = 2)
abline(h = 0, lty = 2)

DoseMiteMod <- mixed(DeltaPMI ~ OxalicDose + (1|Yard),
                      data = UGAdata.2020.NC)

nice(DoseMiteMod)

# No evidence of effectiveness of oxalic acid delivered via these methods across these doses in mite control

## Check effect of glycerin towel by comparing CTRL and 0OA
UGAdata.2020.NOA <- UGAdata.2020[which(UGAdata.2020$OxalicDose == 0),]

boxplot(UGAdata.2020.NOA$DeltaPMI ~ UGAdata.2020.NOA$GlycerinTowel, main = NA, ylab = expression(paste(Delta, ' Percent Mite Intensity',sep='')),
        xlab = 'Glycerine Towel', border = sapply(X = as.character(ColRef$Col[1:2]), FUN = Transpa, percent = 50), 
        cex.axis = 1.3, cex.lab = 1.5, outline = T, lwd = 1.2, outcol = 'white',
        boxlty = 1, whisklty = 1, staplelty = 1, boxwex = 0.7, 
        col = sapply(X = as.character(ColRef$Col[1:2]), FUN = Transpa, percent = 50))
stripchart(UGAdata.2020.NOA$DeltaPMI ~ UGAdata.2020.NOA$GlycerinTowel,
           col = sapply(X = as.character(ColRef$Col[1:2]), FUN = Transpa, percent = 60),
           vertical = T, add = T, pch = 16, cex = 1.2, 
           method = 'jitter', lwd = 2)
abline(h = 0, lty = 2)

GTMiteMod <- mixed(DeltaPMI ~ GlycerinTowel + (1|Yard),
                     data = UGAdata.2020.NOA)

nice(GTMiteMod)


# The glycerin towel itself didn't seem to have an effect, perhaps as expected.


### Note of the graphs above, the second two graphs are just subsets of the first.

# Get whether any treatments showed any change in PMI at all
TreatMiteModNI <- mixed(DeltaPMI ~ TreatmentCode - 1 + (1|Yard),
                      data = UGAdata.2020)

summary(TreatMiteModNI)

# look at coefficient compared to standard error; 95%CI is +- 2*STE 
# only evidence for change =/= 0 is in 'true-control' control group, which significantly increased (> 0)
# note that 0OA (shop towel control) doesn't align with this.

# POI :: 'coefficient is actually just mean: see this comparison)
mean(UGAdata.2020$DeltaPMI[which(UGAdata.2020$TreatmentCode == 'CTRL')])








