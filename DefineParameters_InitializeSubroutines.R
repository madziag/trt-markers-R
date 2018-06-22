# SET DIRECTORY. IMPORT DATA. INITIALIZE PACKAGES
setwd("C:/Users/Eigenaar/Desktop/R")

d <- read.table("tbs.txt")

install.packages("pROC")
install.packages("MASS")
install.packages("rms")
install.packages("splines")

library(MASS); library(pROC); library(rms); library(splines)


# CATEGORICAL FACTORS -> AS FACTORS
d$meth.conc   <- as.factor(d$meth.conc)
d$parity.grp  <- as.factor(d$parity.grp)
d$twinb.pres  <- as.factor(d$twinb.pres)
d$chorio      <- as.factor(d$chorio)
d$cv.grp2     <- as.factor(d$cv.grp2)
d$inlab       <- as.factor(d$inlab)
d$rupmem.rand <- as.factor(d$rupmem.rand)

# TREATMENT MUST BE 0 OR 1
d$group.x     <- ifelse(d$group.x == 1, 0, 1)


############################################
### DEFINE DATA AND CONTROL PARAMETERS #####
############################################

marker <- d[,c("age", "efw.a", "efw.b", "gest.agec2", "g.age.del1",
               "dia.gest", "diab.t", "prev.cs", "prob.preecl", "meth.conc",
               "twinb.pres", "chorio", "cv.grp2", "inlab", "rupmem.rand", "parity.grp")]
trt        <- d$group.x
outcome    <- d$outcome.28d.any

# CONTROL PARAMETERS
trtA       <- "Planned Cesearan"   # NAME OF GROUP 0
trtB       <- "Planned Vaginal"    # NAME OF GROUP 1
threshold  <- 0                    # TREATMENT THRESHOLD
nrsplits   <- 5                 # NUMBER OF SPLITS FOR CROSSVALIDATION
propindevpset <- 0.5               # PROPORTION DESIRED IN DEVELOPMENT SET
nrquantiles <- 10                  # NUMBER OF QUNATILES FOR SUBGROUP ANALYSIS
splinedf    <- 4
rule        <- "p"
sls         <- .05
k.aic       <- 2
nrbstrapstovalidate  <- 100
percentageinvalidate <- 0.40
validate   <- "y"
plotsize   <- 0.4
yaxislimits <- c(0, 0.2)  # DESIRED YLIM(MIN = 0)(MARKER VS OBSERVED RISK/RD-COMPLETE DATASET)



###########################
### SUBROUTINES ###########
###########################

# TREATMENT AND OUTCOME
check1 <- function(x){
  OKAY <- "More than two categories"
  if (length(table(x)) == 2){
    OKAY <- "2 categories and the codes are good"
  }
  if (names(table(x)[1]) != "0" | names(table(x)[2]) != "1"){
    OKAY= "2 categories but the codes are not good"
  }
  return(OKAY)
}


# ADD COLUMN FOR CADIT VALUE 
cadit <- ifelse((trt == 0 & outcome == 0)|(trt == 1 & outcome == 1), 0, 
                ifelse((trt == 0 & outcome == 1)|(trt == 1 & outcome == 0), 1, NA))

###############################
### MISCELLANEOUS FUNCTIONS ###
###############################
# CALCULATES MIDPOINTS OF QUANTILES FOR PLOTS
midpoints <- function(x, dp = 5){
  lower <- as.numeric(sub("\\((.+),.*", "\\1", x))
  upper <- as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", x))
  return(round(lower + (upper - lower)/2, dp))
}

# CLEANS UP GRAPHS
cleanup <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), axis.line = element_line(color = "black"),
                 panel.border = element_rect(colour = "black", fill = NA, size = 1))

#################################
### ANALYSIS ON COMPLETE DATA ###
#################################
# OVERALL EFFECT COMPARISON
OverallTrtEffectF <- function(trtxx, outcomexx){
  yy <- table(trt, outcome)
  xx <- prop.table(table(trt, outcome), 1)
  RD <- xx[2,2] - xx[1,2]
  RDlow   <- RD - 1.96 * sqrt(xx[1,2] * (1 - xx[1,2])/sum(yy[1,1] + yy[1,2]) + xx[2,2] * (1 - xx[2,2])/sum(yy[2,1] + yy[2,2]))
  RDupp   <- RD + 1.96 * sqrt(xx[1,2] * (1 - xx[1,2])/sum(yy[1,1] + yy[1,2]) + xx[2,2] * (1 - xx[2,2])/sum(yy[2,1] + yy[2,2]))
  RDstats <- c(RD = RD, Lower = RDlow, Upper = RDupp)
  RR <- xx[2,2] / xx[1,2]
  RRlow   <- log(RR) - 1.96 * sqrt((1 - xx[1,2])/(xx[1,2] * sum(yy[1,1] + yy[1,2])) + (1 - xx[2,2])/(xx[2,2] * sum(yy[2,1] + yy[2,2])))
  RRupp   <- log(RR) + 1.96 * sqrt((1 - xx[1,2])/(xx[1,2] * sum(yy[1,1] + yy[1,2])) + (1 - xx[2,2])/(xx[2,2] * sum(yy[2,1] + yy[2,2])))
  RRstats <- c(RR = RR, Lower = exp(RRlow),  Upper = exp(RRupp))
  chiSqrTest <- chisq.test(outcome, trt)
  fisherTest <- fisher.test(outcome, trt)
  summaryGLM <- summary(glm(outcome ~ trt, family = binomial))
  return(list(RD = RDstats, RR = RRstats, ChiSqrTest = chiSqrTest, 
              FisherTest = fisherTest, SummaryGLM = summaryGLM))
}   

#UNIVARIATE ANALYSIS
# TESTS STATISTICAL INTERACTION IN THE ENTIRE DATA SET BY BIOMARKER
InteractionsF <- function(trtxx, outcomexx, biomarkersxx) {
  res <- data.frame(matrix(ncol = 9, nrow = 0))
  for (i in 1:ncol(biomarkersxx)) {
    fo  <- paste("yyy ~ trt*", paste(names(biomarkersxx)[i], collapse = "+"))
    ddc <- summary(glm(as.formula(fo), data = data.frame(trt = trtxx, yyy = outcomexx, biomarkersxx), family = binomial)) 
    res <- rbind(res, c(ddc$coefficients[4,4], as.numeric(ddc$coefficients[,1]), as.numeric(ddc$coefficients[,2])))
  }
  res <- data.frame(biomarker = names(biomarkersxx), res = res)
  names(res) <- c("biomarker", "interaction p-value", "estimate.intercept", 
                  "estimate.trt", "estimate.BM", "estimate.trt*BM", "SD.intercept", 
                  "SD.trt", "SD.BM", "SD.trt*BM")
  return(res)
}

# PLOTS STATISTICAL INTERACTIONS BY BIOMARKER
ORplotF <- function(resxx){
  layout(matrix(c(1,2), nrow = 1, ncol = 2), widths = c(1.5,0.5))
  par(mar = c(5, 4, 4, 2) + 0.1)
  plot(resxx[,6], resxx[,10], 
       xlab = "log(OR ratio)", ylab = "standard error", main = "Complete Dataset",
       pch = 16, cex = 1.25, col = 1:length(resxx[,6]), 
       xlim = c(), ylim = c(min(0, min(resxx[,10])), max(0, max(resxx[,10]))))
  abline(a = 0, b = 0.5, lty = 2)
  abline(a = 0, b = -0.5, lty = 2)
  par(mar = c(5, 0, 4, 0) + 0.1)
  plot(0:1, 0:1, col = NULL, axes = F, xlab = "", ylab = "")
  legend("center", legend = (resxx[,1]), pch = 16, col = 1:length(resxx[,6]),
         bty = "n", cex = 0.75)
  par(mar = c(5, 4, 4, 2) + 0.1)
  layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T), widths = c(1))
}

# PLOT OF MOST SIGNIFICANT INTERACTION AND CROSSVALIDATED FIT STATS
MostSignificantBMF <- function(resxx, biomarkersxx, outcomexx, trtxx, nbootxx, yaxislimitsxx) {
  BMminPvalue <- biomarkersxx[,which(as.numeric(as.character(resxx[,2])) == min(as.numeric(as.character(resxx[,2]))))[1]]
  orderbypos  <- order(BMminPvalue)
  BMminPvalue <- BMminPvalue[orderbypos]
  outcomexx   <- outcomexx[orderbypos]
  trtxx       <- trtxx[orderbypos]
  hh0         <- glm(outcomexx ~ trtxx * BMminPvalue, data = data.frame(outcomexx, trtxx, BMminPvalue), family = binomial)
  hh0ext      <- lrm(outcomexx ~ trtxx * BMminPvalue, data = data.frame(outcomexx, trtxx, BMminPvalue), x = T, y = T)
  predstats   <- validate(hh0ext, B = nbootxx)
  hh0a        <- predict.glm(hh0, newdata = data.frame(BMminPvalue, trtxx = rep(0,length(BMminPvalue))), se.fit=TRUE, type = "response")
  hh0b        <- predict.glm(hh0, newdata = data.frame(BMminPvalue, trtxx = rep(1,length(BMminPvalue))), se.fit=TRUE, type = "response")
  yaxislimitsN <- c(-yaxislimitsxx[2], yaxislimitsxx[2])
  plot(BMminPvalue, (hh0a$fit - hh0b$fit), 
       ylab = paste("Risk Difference(",trtA, "-", trtB,")"), type = "l", lwd = 3, col = 1, lty = 1, 
       xlab = as.character(resxx[which(as.numeric(as.character(resxx[,2])) == min(as.numeric(as.character(resxx[,2]))))[1],1]),
       main = "Complete Dataset", ylim = yaxislimitsN)
  abline(h = 0)
  sefit = sqrt(hh0a$se.fit^2 + hh0b$se.fit^2)
  lines(BMminPvalue, (hh0a$fit - hh0b$fit) + 2 * sefit, lwd = 2.5, lty = 2, col = 1)
  lines(BMminPvalue, (hh0a$fit - hh0b$fit) - 2 * sefit, lwd = 2.5, lty = 2, col = 1)
  return(predstats)
}

ChecktablesBinaryF <- function(selectedbiomarker, biomarkersxx, trtxx, outcomexx){
  groupbms <- biomarkersxx[,selectedbiomarker]
  groupbmsxx <- na.omit(groupbms)
  from <- as.numeric(unique(groupbmsxx)[1])
  x <- unique(groupbmsxx)
  to <- as.numeric(x[length(x)])
  for (i in from:to){
    trtxx <- trt[groupbmsxx == i]
    outcomexx <- outcome[groupbmsxx == i]
    if (length(table(trtxx, outcomexx)) == 4){
      CHECKTABLE <- "good"
    }else{
      CHECKTABLE <- "bad"
    }
  }
  return (CHECKTABLE)
}

ChecktablesFactorsF <- function(selectedbiomarker, biomarkersxx, trtxx, outcomexx){
  groupbms <- biomarkersxx[,selectedbiomarker]
  groupbmsxx <- na.omit(groupbms)
  from <- as.numeric(levels(groupbmsxx)[1])
  x <- levels(groupbmsxx)
  to <- as.numeric(x[length(x)])
  for (i in from:to){
    trt <- trtxx[groupbmsxx == i]
    outcome <- outcomexx[groupbmsxx == i]
    if (length(table(trt, outcome)) == 4){
      CHECKTABLE <- "good"
    }else{
      CHECKTABLE <- "bad"
    }
  }
  return (CHECKTABLE)
}

# BIOMARKER VS OBSERVED RISK DIFFERENCE: CONTINUOUS VARIABLES
RDSingleMarkerContF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx,
                                nrquantilesxx, splinedfxx, yaxislimitsxx, thresholdxx){
  groupbms     <- biomarkersxx[,selectedbiomarker]
  groupbms     <- as.numeric(unlist(groupbms))
  tempdf       <- na.omit(data.frame(trtxx, outcomexx, groupbms))
  cutpoints    <- unique(c(-1, quantile(tempdf$groupbms, probs = c((1:nrquantilesxx)/nrquantilesxx))))
  subpops      <- cut(tempdf$groupbms, cutpoints)
  subpops      <- midpoints(subpops)
  hh2          <- table(tempdf$trtxx, tempdf$outcomexx, subpops)[,2,] /apply(table(tempdf$trtxx, tempdf$outcomexx, subpops), c(1,3), sum)
  diffs        <- hh2[1,] - hh2[2,]
  cutpoints    <- as.numeric(cutpoints[-1])
  cutpoints    <- cutpoints[is.na(diffs) == FALSE]
  diffs        <- diffs[is.na(diffs) == FALSE]
  yaxislimitsN <- c(-yaxislimitsxx[2], yaxislimitsxx[2])
  hh3          <- predict(lm(diffs ~ ns(cutpoints,df= splinedfxx)), se.fit = TRUE)
  par(mfrow = c(1,1))
  par(mar = c(5, 4, 4, 2) + 0.1)
  plot(cutpoints, diffs, col = NULL, 
       xlab = paste("subgroups:", selectedbiomarker), ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
       ylim =  yaxislimitsN, main = "Risk Difference: Continuous Markers")
  abline(h = thresholdxx)
  points(cutpoints,diffs)
  lines(cutpoints, hh3$fit, col = 1,lty = 1,lwd = 2)
  lines(cutpoints, hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
  lines(cutpoints, hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
}

# BIOMARKER VS OBSERVED RISK: CONTINUOUS VARIABLES IN SEPARATE TREATMENT AND NON TREATMENT GROUPS
RiskIndividualTrtsContF <- 
  function(selectedbiomarkerxx, trtxx, outcomexx, biomarkersxx, nrquantilesxx, splinedfxx, 
           yaxislimitsxx, thresholdxx) {
    groupbms   <- biomarkersxx[,selectedbiomarkerxx]
    groupbms   <- as.numeric(unlist(groupbms))
    tempdf     <- na.omit(data.frame(trtxx, outcomexx, groupbms))
    # TRT0
    tempdf0    <- subset(tempdf, tempdf$trtxx == 0)
    cutpoints0 <- unique(c(-1, quantile(tempdf0$groupbms, probs = c((1:nrquantilesxx) / nrquantilesxx))))
    subpops0   <- cut(tempdf0$groupbms, cutpoints0)
    subpops0   <- midpoints(subpops0)
    props0     <- table(tempdf0$trtxx, tempdf0$outcomexx, subpops0)[,2,] / apply(table(tempdf0$trtxx, tempdf0$outcomexx, subpops0),c(1, 3),sum)
    cutpoints0 <- as.numeric(cutpoints0[-1])
    cutpoints0 <- cutpoints0[is.na(props0) == FALSE]
    props0     <- props0[is.na(props0) == FALSE]
    hh0        <- predict(lm(props0 ~ ns(cutpoints0, df = splinedfxx)), se.fit = TRUE)
    # TRT1
    tempdf1    <- subset(tempdf, tempdf$trtxx == 1)
    cutpoints1 <- unique(c(-1, quantile(tempdf1$groupbms, probs = c((1:nrquantilesxx) / nrquantilesxx))))
    subpops1   <- cut(tempdf1$groupbms, cutpoints1)
    subpops1   <- midpoints(subpops1)
    props1     <- table(tempdf1$trtxx, tempdf1$outcomexx, subpops1)[,2,] / apply(table(tempdf1$trtxx, tempdf1$outcomexx, subpops1),c(1, 3),sum)
    cutpoints1 <- as.numeric(cutpoints1[-1])
    cutpoints1 <- cutpoints1[is.na(props1) == FALSE]
    props1     <- props1[is.na(props1) == FALSE]
    hh1        <- predict(lm(props1 ~ ns(cutpoints1, df = splinedfxx)),se.fit = TRUE)
    par(mfrow = c(1,1))
    par(mar = c(5, 4, 4, 2) + 0.1)
    plot(cutpoints0, props0, col = NULL, xlab = paste("subgroups:", selectedbiomarkerxx),
         ylab = "Risk of Outcome", main = "Risk: Continuous Markers - Individual Treatment Groups", 
         ylim = yaxislimitsxx)
    abline(h = thresholdxx)
    points(cutpoints0, props0)
    points(cutpoints1, props1, col = 2)
    lines(cutpoints0, hh0$fit, col = 1, lty = 1, lwd = 2)
    lines(cutpoints1, hh1$fit, col = 2, lty = 1, lwd = 2)
    legend("topright", inset = .05, cex = 1, title = "Legend", c(trtA,trtB),
           horiz = FALSE, lty = c(1,1), lwd = c(2,2), col = c("black","red"), bg ="grey96", text.font = 3)
  }

# BIOMARKER VS OBSERVED RISK DIFFERENCE: DICHOTOMOUS VARIABLES
RDSingleMarkerDichF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx, 
                                yaxislimitsxx){
  
  ChecktablesBinary<- ChecktablesBinaryF(selectedbiomarker = selectedbiomarker, 
                                         biomarkersxx =marker, trtxx = trt, 
                                         outcomexx = outcomexx)
  
  if(ChecktablesBinary == "good"){  
    groupbms     <- biomarkersxx[,selectedbiomarker]
    groupbms     <- as.numeric(unlist(groupbms))
    tempdf       <- na.omit(data.frame(trtxx, outcomexx, groupbms))
    
    # GROUP 0 (BMS = 0)
    group0 <- subset(tempdf, tempdf$groupbms == 0)
    t0     <- table(group0$trtxx, group0$outcomexx)
    pt0    <- prop.table(table(group0$trtxx, group0$outcomexx), 1)
    RDGroup0  <- pt0[2,2] - pt0[1,2]
    RD0low    <- RDGroup0 - 1.96 * sqrt(pt0[1,2] * (1 - pt0[1,2])/sum(t0[1,1] + t0[1,2]) + pt0[2,2] * (1 - pt0[2,2])/sum(t0[2,1] + t0[2,2]))
    RD0upp    <- RDGroup0 + 1.96 * sqrt(pt0[1,2] * (1 - pt0[1,2])/sum(t0[1,1] + t0[1,2]) + pt0[2,2] * (1 - pt0[2,2])/sum(t0[2,1] + t0[2,2]))
    # GROUP 1 (BMS = 1)
    group1 <- subset(tempdf, tempdf$groupbms == 1)
    t1     <- table(group1$trtxx, group1$outcomexx)
    pt1    <- prop.table(table(group1$trtxx, group1$outcomexx), 1)
    RDGroup1  <- pt1[2,2] - pt1[1,2]
    RD1low    <- RDGroup1 - 1.96 * sqrt(pt1[1,2] * (1 - pt1[1,2])/sum(t1[1,1] + t1[1,2]) + pt1[2,2] * (1 - pt1[2,2])/sum(t1[2,1] + t1[2,2]))
    RD1upp    <- RDGroup1 + 1.96 * sqrt(pt1[1,2] * (1 - pt1[1,2])/sum(t1[1,1] + t1[1,2]) + pt1[2,2] * (1 - pt1[2,2])/sum(t1[2,1] + t1[2,2]))
    # ADD TO DATAFRAME   
    tempdf$riskdiff <- ifelse(tempdf$groupbms == 0 & (tempdf$trtxx == 0 | tempdf$trtxx == 1), RDGroup0, 
                              ifelse(tempdf$groupbms == 1 & (tempdf$trtxx == 0 | tempdf$trtxx == 1), RDGroup1, NA))
    tempdf$riskdiffupper <- ifelse(tempdf$groupbms == 0 & (tempdf$trtxx == 0 | tempdf$trtxx == 1), RD0upp, 
                                   ifelse(tempdf$groupbms == 1 & (tempdf$trtxx == 0 | tempdf$trtxx == 1), RD1upp, NA))
    tempdf$riskdifflower <- ifelse(tempdf$groupbms == 0 & (tempdf$trtxx == 0 | tempdf$trtxx == 1), RD0low, 
                                   ifelse(tempdf$groupbms == 1 & (tempdf$trtxx == 0 | tempdf$trtxx == 1), RD1low, NA))
    
    yaxislimitsN <- c(-yaxislimitsxx[2], yaxislimitsxx[2])
    pRD <- ggplot(tempdf, aes(x = factor(tempdf$groupbms), y = tempdf$riskdiff)) 
    pRD +
      labs(title = "Risk Difference: Binary Markers", 
           x = paste("Biomarker:", selectedbiomarker), y = "Risk Difference") +
      scale_x_discrete(labels = c((paste(selectedbiomarker, "= 0")), (paste(selectedbiomarker, "= 1")))) +
      ylim(yaxislimitsN) +
      geom_errorbar(aes(ymin = tempdf$riskdifflower, ymax = tempdf$riskdiffupper), width = 0.05) +
      geom_point(colour = "black", size = 4) + 
      theme(axis.text = element_text(size = 12), title = element_text(size = 14, face = "bold")) +
      cleanup 
  }else{
    print("Cannot calculate Risk Difference - Check your marker")
  }
}

# BIOMARKER VS OBSERVED RISK:DICHOTOMOUS VARIABLES IN SEPARATE TREATMENT AND NON TREATMENT GROUPS
RiskIndividualTrtsDichF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx, yaxislimitsxx){
  ChecktablesBinary<- ChecktablesBinaryF(selectedbiomarker = selectedbiomarker, 
                                         biomarkersxx = marker, trtxx = trt, outcomexx = outcome)
  if(ChecktablesBinary == "good"){  
    groupbms  <- biomarkersxx[,selectedbiomarker]
    groupbms  <- as.numeric(unlist(groupbms))
    tempdf    <- na.omit(data.frame(trtxx, outcomexx, groupbms))
    # GROUP 0 (BMS = 0)
    tempdf0   <- subset(tempdf, tempdf$groupbms == 0)
    t0        <- table(tempdf0$trtxx, tempdf0$outcomexx)
    pt0       <- prop.table(table(tempdf0$trtxx, tempdf0$outcomexx), 1)
    IRB0T0    <- pt0[1,2]
    IRB0T0low <- IRB0T0 - 1.96 * sqrt((pt0[1,2]*(1-pt0[1,2]))/(t0[1,1]+t0[1,2])) 
    IRB0T0low <- ifelse(IRB0T0low < 0, 0, IRB0T0low)
    IRB0T0upp <- IRB0T0 + 1.96 * sqrt((pt0[1,2]*(1-pt0[1,2]))/(t0[1,1]+t0[1,2])) 
    IRB0T1    <- pt0[2,2]
    IRB0T1low <- IRB0T1 - 1.96 * sqrt((pt0[2,2]*(1-pt0[2,2]))/(t0[2,1]+t0[2,2])) 
    IRB0T1low <- ifelse(IRB0T1low < 0, 0, IRB0T1low)
    IRB0T1upp <- IRB0T1 + 1.96 * sqrt((pt0[2,2]*(1-pt0[2,2]))/(t0[2,1]+t0[2,2])) 
    # GROUP 1 (BMS = 1)
    tempdf1   <- subset(tempdf, tempdf$groupbms == 1)
    t1        <- table(tempdf1$trtxx, tempdf1$outcomexx)
    pt1       <- prop.table(table(tempdf1$trtxx, tempdf1$outcomexx), 1)
    pt1       <- prop.table(table(tempdf1$trtxx, tempdf1$outcomexx), 1)
    IRB1T0    <- pt1[1,2]
    IRB1T0low <- IRB1T0 - 1.96 * sqrt((pt1[1,2]*(1-pt1[1,2]))/(t1[1,1]+t1[1,2])) 
    IRB1T0low <- ifelse(IRB1T0low < 0, 0, IRB1T0low)
    IRB1T0upp <- IRB1T0 + 1.96 * sqrt((pt1[1,2]*(1-pt1[1,2]))/(t1[1,1]+t1[1,2])) 
    IRB1T1    <- pt0[2,2]
    IRB1T1low <- IRB1T1 - 1.96 * sqrt((pt1[2,2]*(1-pt1[2,2]))/(t1[2,1]+t1[2,2])) 
    IRB1T1low <- ifelse(IRB1T1low < 0, 0, IRB1T1low)
    IRB1T1upp <- IRB1T1 + 1.96 * sqrt((pt1[2,2]*(1-pt1[2,2]))/(t1[2,1]+t1[2,2])) 
    # ADD TO DATAFRAME
    tempdf$IncidenceRate <- ifelse(tempdf$groupbms == 0 & tempdf$trtxx == 0, IRB0T0, 
                                   ifelse(tempdf$groupbms== 0 & tempdf$trtxx == 1, IRB0T1, 
                                          ifelse(tempdf$groupbms == 1 & tempdf$trtxx == 0, IRB1T0, 
                                                 ifelse(tempdf$groupbms == 1 & tempdf$trtxx == 1, IRB1T1, NA))))
    tempdf$IRLower <- ifelse(tempdf$groupbms == 0 & tempdf$trtxx == 0, IRB0T0low, 
                             ifelse(tempdf$groupbms == 0 & tempdf$trtxx == 1, IRB0T1low, 
                                    ifelse(tempdf$groupbms == 1 & tempdf$trtxx == 0, IRB1T0low, 
                                           ifelse(tempdf$groupbms == 1 & tempdf$trtxx == 1, IRB1T1low, NA))))
    tempdf$IRUpper <- ifelse(tempdf$groupbms == 0 & tempdf$trtxx == 0, IRB0T0upp, 
                             ifelse(tempdf$groupbms == 0 & tempdf$trtxx == 1, IRB0T1upp, 
                                    ifelse(tempdf$groupbms == 1 & tempdf$trtxx == 0, IRB1T0upp, 
                                           ifelse(tempdf$groupbms == 1 & tempdf$trtxx == 1, IRB1T1upp, NA))))
    
    limits <- aes(ymax = tempdf$IRUpper, ymin = tempdf$IRLower)
    yaxislimitsN <- c(- yaxislimitsxx[2], yaxislimitsxx[2])
    dodge <- position_dodge(width = 0.3)  
    p <- ggplot(tempdf, aes(x = as.factor(groupbms), y = IncidenceRate, group = as.factor(trtxx), color = as.factor(trtxx)))
    p + 
      labs(title = "Risk: Binary Markers - Individual Treatment Groups", 
           x = paste("Biomarker:", selectedbiomarker), y = "Risk") +   
      scale_x_discrete(labels = c((paste(selectedbiomarker, "= 0")), (paste(selectedbiomarker, "= 1")))) +
      ylim(yaxislimitsxx) +
      geom_errorbar(limits, width = 0.1, position = dodge) + 
      geom_point(size = 4, position = dodge) +
      scale_color_manual(values = c("black", "red"), name  = "Treatment", breaks = c(0,1), labels = c(trtA, trtB)) +
      theme(axis.text = element_text(size = 12), title = element_text(size = 12, face = "bold")) + cleanup 
  }else{
    print("Cannot calculate Risk Difference - Check your marker")
  }
}


RDSingleMarkerFactorsF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx, 
                                   yaxislimitsxx){
  
  ChecktablesFactors <- ChecktablesFactorsF(selectedbiomarker = selectedbiomarker, 
                                            biomarkersxx = marker, trtxx = trt, 
                                            outcomexx = outcome)
  
  if(ChecktablesFactors == "good"){
    groupbms     <- biomarkersxx[,selectedbiomarker]
    tempdf       <- na.omit(data.frame(trt, outcome, groupbms))
    RDGroup <- c() ;RDlow <- c(); RDupp <- c()
    from <- as.numeric(levels(groupbmsxx)[1])
    x <- levels(groupbmsxx)
    to <- as.numeric(x[length(x)])
    for (i in from:to){
      group <- subset(tempdf, tempdf$groupbms == i)
      t     <- table(group$trt, group$outcome)
      pt    <- prop.table(table(group$trt, group$outcome), 1)
      RDGroup[i] <- pt[2,2] - pt[1,2]
      RDlow[i]   <- RDGroup[i] - 1.96 * sqrt(pt[1,2] * (1 - pt[1,2])/sum(t[1,1] + t[1,2]) + pt[2,2] * (1 - pt[2,2])/sum(t[2,1] + t[2,2]))
      RDupp[i]   <- RDGroup[i] + 1.96 * sqrt(pt[1,2] * (1 - pt[1,2])/sum(t[1,1] + t[1,2]) + pt[2,2] * (1 - pt[2,2])/sum(t[2,1] + t[2,2]))
    }
    
    RDTable = cbind(RDlow, RDGroup, RDupp) 
    valuesRDL <- c(RDTable[,1])
    tempdf$riskdifflower <- valuesRDL[tempdf$groupbms]  
    valuesRD <- c(RDTable[,2])
    tempdf$riskdiff <- valuesRD[tempdf$groupbms]      
    valuesRDU <- c(RDTable[,3])
    tempdf$riskdiffupper <- valuesRDU[tempdf$groupbms]  
    
    yaxislimitsN <- c(-yaxislimits[2], yaxislimits[2])
    pRD <- ggplot(tempdf, aes(x = factor(tempdf$groupbms), y = tempdf$riskdiff)) 
    pRD +
      labs(title = "Risk Difference: Dichotomous Markers", 
           x = paste("Biomarker:", selectedbiomarker), y = "Risk Difference") +
      scale_x_discrete(labels = c(levels(groupbms))) +
      ylim(yaxislimitsN) +
      geom_errorbar(aes(ymin = tempdf$riskdifflower, ymax = tempdf$riskdiffupper), width = 0.05) +
      geom_point(colour = "black", size = 4) + 
      theme(axis.text = element_text(size = 12), title = element_text(size = 14, face = "bold")) +
      cleanup  
  }else{
    print("Cannot calculate Risk Difference - Check your marker levels")
  }
}


RiskIndividualTrtsFactorsF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx, 
                                       yaxislimitsxx){
  
  ChecktablesFactors <- ChecktablesFactorsF(selectedbiomarker = selectedbiomarker, 
                                            biomarkersxx = marker, trtxx = trt, 
                                            outcomexx = outcome)
  
  if(ChecktablesFactors == "good"){
    groupbms <- biomarkersxx[,selectedbiomarker]
    tempdf   <- na.omit(data.frame(trt, outcome, groupbms))
    IRT0 <- c(); IRT0low <- c(); IRT0upp <- c()
    IRT1 <- c(); IRT1low <- c(); IRT1upp <- c()
    from <- as.numeric(levels(groupbmsxx)[1])
    x <- levels(groupbmsxx)
    to <- as.numeric(x[length(x)])
    for (i in from:to){
      group <- subset(tempdf, tempdf$groupbms == i)
      t     <- table(group$trt, group$outcome)
      pt    <- prop.table(table(group$trt, group$outcome), 1)
      IRT0[i]    <- pt[1,2]
      IRT0low[i] <- IRT0[i] - 1.96 * sqrt((pt[1,2]*(1-pt[1,2]))/(t[1,1]+t[1,2])) 
      IRT0low[i] <- ifelse(IRT0low[i] < 0, 0, IRT0low[i])
      IRT0upp[i] <- IRT0[i] + 1.96 * sqrt((pt[1,2]*(1-pt[1,2]))/(t[1,1]+t[1,2])) 
      IRT1[i]    <- pt[2,2]
      IRT1low[i] <- IRT1[i] - 1.96 * sqrt((pt[2,2]*(1-pt[2,2]))/(t[2,1]+t[2,2])) 
      IRT1low[i] <- ifelse(IRT1low[i]  < 0, 0, IRT1low[i])
      IRT1upp[i] <- IRT1[i] + 1.96 * sqrt((pt[2,2]*(1-pt[2,2]))/(t[2,1]+t[2,2])) 
    }
    IRTable = cbind(IRT0, IRT0low, IRT0upp, IRT1, IRT1low, IRT1upp) 
    
    tempdf0 <- tempdf[tempdf$trt == 0,]
    tempdf1 <- tempdf[tempdf$trt == 1,]
    
    valuesIRT0low  <- c(IRTable[,2])
    tempdf0$IRlow   <- valuesIRT0low[tempdf0$groupbms]  
    valuesIRT0     <- c(IRTable[,1])
    tempdf0$IR      <- valuesIRT0[tempdf0$groupbms]      
    valuesIRT0upp  <- c(IRTable[,3])
    tempdf0$IRupp   <- valuesIRT0upp[tempdf0$groupbms] 
    valuesIRT1low <- c(IRTable[,5])
    tempdf1$IRlow  <- valuesIRT1low[tempdf1$groupbms]  
    valuesIRT1    <- c(IRTable[,4])
    tempdf1$IR     <- valuesIRT1[tempdf1$groupbms]      
    valuesIRT1upp <- c(IRTable[,6])
    tempdf1$IRupp  <- valuesIRT1upp[tempdf1$groupbms] 
    tempdf <- rbind(tempdf0, tempdf1)
    
    dodge <- position_dodge(width = 0.3)  
    limits <- aes(ymax = tempdf$IRupp, ymin = tempdf$IRlow)
    PI <- ggplot(tempdf, aes(x = as.factor(groupbms), y = IR, group = as.factor(trt), color = as.factor(trt)))
    
    PI +
      labs(title = "Risk: Categorical Markers - Individual Treatment Groups", 
           x = paste("Biomarker:", selectedbiomarker), y = "Risk") +
      scale_x_discrete(labels = c(levels(groupbms))) +
      ylim(yaxislimitsxx) +
      geom_errorbar(limits, width = 0.1, position = dodge) +
      geom_point(size = 4, position = dodge) + 
      scale_color_manual(values = c("black", "red"), name  = "Treatment", breaks = c(0,1), labels = c(trtA, trtB)) +
      theme(axis.text = element_text(size = 12), title = element_text(size = 12, face = "bold")) +
      cleanup  
  }else{
    print("Cannot calculate Risk Difference - Check your marker levels")
  }
}


RDSingleMarkerFactorsF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx, 
                                   yaxislimitsxx){
  
  ChecktablesFactors <- ChecktablesFactorsF(selectedbiomarker = selectedbiomarker, 
                                            biomarkersxx = marker, trtxx = trt, 
                                            outcomexx = outcome)
  
  if(ChecktablesFactors == "good"){
    groupbms     <- biomarkersxx[,selectedbiomarker]
    tempdf       <- na.omit(data.frame(trt, outcome, groupbms))
    RDGroup <- c() ;RDlow <- c(); RDupp <- c()
    
    from <- as.numeric(levels(tempdf$groupbms)[1])
    x <- levels(tempdf$groupbms)
    to <- as.numeric(x[length(x)])
    for (i in from:to){
      group <- subset(tempdf, tempdf$groupbms == i)
      t     <- table(group$trt, group$outcome)
      pt    <- prop.table(table(group$trt, group$outcome), 1)
      RDGroup[i] <- pt[2,2] - pt[1,2]
      RDlow[i]   <- RDGroup[i] - 1.96 * sqrt(pt[1,2] * (1 - pt[1,2])/sum(t[1,1] + t[1,2]) + pt[2,2] * (1 - pt[2,2])/sum(t[2,1] + t[2,2]))
      RDupp[i]   <- RDGroup[i] + 1.96 * sqrt(pt[1,2] * (1 - pt[1,2])/sum(t[1,1] + t[1,2]) + pt[2,2] * (1 - pt[2,2])/sum(t[2,1] + t[2,2]))
    }
    
    RDTable = cbind(RDlow, RDGroup, RDupp) 
    valuesRDL <- c(RDTable[,1])
    tempdf$riskdifflower <- valuesRDL[tempdf$groupbms]  
    valuesRD <- c(RDTable[,2])
    tempdf$riskdiff <- valuesRD[tempdf$groupbms]      
    valuesRDU <- c(RDTable[,3])
    tempdf$riskdiffupper <- valuesRDU[tempdf$groupbms]  
    
    yaxislimitsN <- c(-yaxislimits[2], yaxislimits[2])
    pRD <- ggplot(tempdf, aes(x = factor(tempdf$groupbms), y = tempdf$riskdiff)) 
    pRD +
      labs(title = "Risk Difference: Categorical Markers - Individual Treatment Groups", 
           x = paste("Biomarker:", selectedbiomarker), y = "Risk Difference") +
      scale_x_discrete(labels = c(levels(groupbms))) +
      ylim(yaxislimitsN) +
      geom_errorbar(aes(ymin = tempdf$riskdifflower, ymax = tempdf$riskdiffupper), width = 0.05) +
      geom_point(colour = "black", size = 4) + 
      theme(axis.text = element_text(size = 12), title = element_text(size = 14, face = "bold")) +
      cleanup  
  }else{
    print("Cannot calculate Risk Difference for all factor levels")
  }
}

RiskIndividualTrtsFactorsF <- function(selectedbiomarker, trtxx, outcomexx, biomarkersxx, 
                                       yaxislimitsxx){
  
  ChecktablesFactors <- ChecktablesFactorsF(selectedbiomarker = selectedbiomarker, 
                                            biomarkersxx = marker, trtxx = trt, 
                                            outcomexx = outcome)
  
  if(ChecktablesFactors == "good"){
    groupbms <- biomarkersxx[,selectedbiomarker]
    tempdf   <- na.omit(data.frame(trt, outcome, groupbms))
    IRT0 <- c(); IRT0low <- c(); IRT0upp <- c()
    IRT1 <- c(); IRT1low <- c(); IRT1upp <- c()
    
    from <- as.numeric(levels(tempdf$groupbms)[1])
    x <- levels(tempdf$groupbms)
    to <- as.numeric(x[length(x)])
    for (i in from:to){
      group <- subset(tempdf, tempdf$groupbms == i)
      t     <- table(group$trt, group$outcome)
      pt    <- prop.table(table(group$trt, group$outcome), 1)
      IRT0[i]    <- pt[1,2]
      IRT0low[i] <- IRT0[i] - 1.96 * sqrt((pt[1,2]*(1-pt[1,2]))/(t[1,1]+t[1,2])) 
      IRT0low[i] <- ifelse(IRT0low[i] < 0, 0, IRT0low[i])
      IRT0upp[i] <- IRT0[i] + 1.96 * sqrt((pt[1,2]*(1-pt[1,2]))/(t[1,1]+t[1,2])) 
      IRT1[i]    <- pt[2,2]
      IRT1low[i] <- IRT1[i] - 1.96 * sqrt((pt[2,2]*(1-pt[2,2]))/(t[2,1]+t[2,2])) 
      IRT1low[i] <- ifelse(IRT1low[i]  < 0, 0, IRT1low[i])
      IRT1upp[i] <- IRT1[i] + 1.96 * sqrt((pt[2,2]*(1-pt[2,2]))/(t[2,1]+t[2,2])) 
    }
    IRTable = cbind(IRT0, IRT0low, IRT0upp, IRT1, IRT1low, IRT1upp) 
    
    tempdf0 <- tempdf[tempdf$trt == 0,]
    tempdf1 <- tempdf[tempdf$trt == 1,]
    
    valuesIRT0low  <- c(IRTable[,2])
    tempdf0$IRlow   <- valuesIRT0low[tempdf0$groupbms]  
    valuesIRT0     <- c(IRTable[,1])
    tempdf0$IR      <- valuesIRT0[tempdf0$groupbms]      
    valuesIRT0upp  <- c(IRTable[,3])
    tempdf0$IRupp   <- valuesIRT0upp[tempdf0$groupbms] 
    valuesIRT1low <- c(IRTable[,5])
    tempdf1$IRlow  <- valuesIRT1low[tempdf1$groupbms]  
    valuesIRT1    <- c(IRTable[,4])
    tempdf1$IR     <- valuesIRT1[tempdf1$groupbms]      
    valuesIRT1upp <- c(IRTable[,6])
    tempdf1$IRupp  <- valuesIRT1upp[tempdf1$groupbms] 
    tempdf <- rbind(tempdf0, tempdf1)
    
    dodge <- position_dodge(width = 0.3)  
    limits <- aes(ymax = tempdf$IRupp, ymin = tempdf$IRlow)
    PI <- ggplot(tempdf, aes(x = as.factor(groupbms), y = IR, group = as.factor(trt), color = as.factor(trt)))
    
    PI +
      labs(title = "Risk: Categorical Markers - Individual Treatment Groups", 
           x = paste("Biomarker:", selectedbiomarker), y = "Risk") +
      scale_x_discrete(labels = c(levels(groupbms))) +
      ylim(yaxislimitsxx) +
      geom_errorbar(limits, width = 0.1, position = dodge) +
      geom_point(size = 4, position = dodge) + 
      scale_color_manual(values = c("black", "red"), name  = "Treatment", breaks = c(0,1), labels = c(trtA, trtB)) +
      theme(axis.text = element_text(size = 12), title = element_text(size = 12, face = "bold")) +
      cleanup  
  }else{
    print("Cannot calculate Risk for all factor levels")
  }
}
# BIOMARKER VS OBSERVED RISK DIFFERENCE: ALL VARIABLES 
RDSingleMarkerForAllTypesF <- function(selectedbiomarker){
  SelBiomrk <- marker[,selectedbiomarker]
  if((length(table(SelBiomrk)) == 2) & 
     (names(table(SelBiomrk))[1] == "0") & 
     (names(table(SelBiomrk))[2] == "1") &
     (is.factor(SelBiomrk) == FALSE)){
    RDSingleMarkerDich <- RDSingleMarkerDichF(selectedbiomarker = selectedbiomarker, trtxx = trt, 
                                              outcomexx = outcome, biomarkersxx = marker, 
                                              yaxislimitsxx = yaxislimits)
  } else if (length(table(SelBiomrk)) > 2 & (is.factor(SelBiomrk) == FALSE)){
    RDSingleMarkerCont <- RDSingleMarkerContF(selectedbiomarker = selectedbiomarker, trtxx = trt,
                                              outcomexx = outcome, biomarkersxx = marker,
                                              nrquantilesxx = nrquantiles, splinedfxx = splinedf,
                                              yaxislimitsxx = yaxislimits, thresholdxx = threshold)
  } else if(is.factor(SelBiomrk) == TRUE){
    RDSingleMarkerFactors <- RDSingleMarkerFactorsF(selectedbiomarker = selectedbiomarker, 
                                                    trtxx = trt, outcomexx = outcome, 
                                                    biomarkersxx = marker, 
                                                    yaxislimitsxx = yaxislimits)
  } else{
    print("Marker can only be Continuous, Binary or a Factor")
  }
}
# BIOMARKER VS OBSERVED RISK DIFFERENCE: ALL VARIABLES 
RiskIndTrt1MarkerF <- function(selectedbiomarker){
  SelBiomrk <- marker[,selectedbiomarker]
  if((length(table(SelBiomrk)) == 2) & 
     (names(table(SelBiomrk))[1] == "0") & 
     (names(table(SelBiomrk))[2] == "1") &
     (is.factor(SelBiomrk) == FALSE)){
    RiskIndividualTrtsDich <- RiskIndividualTrtsDichF(selectedbiomarker = selectedbiomarker ,  trtxx = trt,  
                                                      outcomexx = outcome, biomarkersxx=marker,
                                                      yaxislimitsxx = yaxislimits)
  } else if (length(table(SelBiomrk)) > 2 & (is.factor(SelBiomrk) == FALSE)){
    RiskIndividualTrtsCont <- RiskIndividualTrtsContF(selectedbiomarker = selectedbiomarker, trtxx = trt, 
                                                      outcomexx = outcome, biomarkersxx = marker,
                                                      nrquantilesxx = nrquantiles, 
                                                      splinedfxx = splinedf, 
                                                      yaxislimitsxx = yaxislimits, 
                                                      thresholdxx = threshold)
  } else if(is.factor(SelBiomrk) == TRUE){
    RiskIndividualTrtsFactors<- RiskIndividualTrtsFactorsF(selectedbiomarker = selectedbiomarker, trtxx = trt, 
                                                           outcomexx = outcome, biomarkersxx=marker, 
                                                           yaxislimitsxx = yaxislimits)
  } else{
    print("Marker can only be Continuous, Binary or a Factor")
  }
}

############################################
### SPLIT-BUILD-VALIDATE-CROSSVALIDATE #####
############################################
# SPLIT DATA
splitdata <- function(biomarkersxx, trtxx, outcomexx, caditxx, propindevpsetxx) {
  outcome0 <- subset(outcomexx, trtxx == 0); outcome1 <- subset(outcomexx, trtxx == 1)
  trt0 <- trtxx[trtxx == 0]; trt1 <- trtxx[trtxx == 1]
  biomarkers0 <- subset(biomarkersxx, trtxx == 0); biomarkers1 <- subset(biomarkersxx, trtxx == 1)
  cadit0 <- subset(caditxx, trtxx == 0); cadit1 <- subset(caditxx, trtxx == 1)
  help0 = rnorm(nrow(biomarkers0), 0, 1)
  help1 = rnorm(nrow(biomarkers1), 0, 1)
  
  outcomeD0   = subset(outcome0, help0    <= quantile(help0, prob = propindevpset))
  trtD0       = subset(trt0, help0        <= quantile(help0, prob = propindevpset))
  bmsD0       = subset(biomarkers0, help0 <= quantile(help0, prob = propindevpset))
  caditD0     = subset(cadit0, help0      <= quantile(help0, prob = propindevpset))
  outcomeD1   = subset(outcome1, help1    <= quantile(help1, prob = propindevpset))
  trtD1       = subset(trt1, help1        <= quantile(help1, prob = propindevpset))
  bmsD1       = subset(biomarkers1, help1 <= quantile(help1, prob = propindevpset))
  caditD1     = subset(cadit1, help1      <= quantile(help1, prob = propindevpset))
  
  outcomeV0   = subset(outcome0, help0     > quantile(help0, prob = propindevpset))
  trtV0       = subset(trt0, help0         > quantile(help0, prob = propindevpset))
  bmsV0       = subset(biomarkers0, help0  > quantile(help0, prob = propindevpset))
  caditV0     = subset(cadit0, help0       > quantile(help0, prob = propindevpset))
  outcomeV1   = subset(outcome1, help1     > quantile(help1, prob = propindevpset))
  trtV1       = subset(trt1, help1         > quantile(help1, prob = propindevpset))
  bmsV1       = subset(biomarkers1, help1  > quantile(help1, prob = propindevpset))
  caditV1     = subset(cadit1, help1       > quantile(help1, prob = propindevpset))
  
  DevpGroupCadit   <- c(caditD0, caditD1)
  ValGroupbms      <- rbind(bmsV0, bmsV1)
  ValGrouptrt      <- c(trtV0, trtV1)
  ValGroupoutcome  <- c(outcomeV0, outcomeV1)
  ValGroupCadit    <- c(caditV0, caditV1)
  
  return(list(OutcomeD0 = outcomeD0, OutcomeD1 = outcomeD1, 
              TrtD0 = trtD0, TrtD1 = trtD1,  
              BmsD0 = bmsD0, BmsD1 = bmsD1, 
              CaditD01 = DevpGroupCadit, 
              OutcomeV01 = ValGroupoutcome, 
              TrtV01 = ValGrouptrt, 
              BmsV01 = ValGroupbms))
}

# SELECT BIOMARKERS
selectBMSRDorCADIT <- function(BMSD0, BMSD1, OUTCOMED0, OUTCOMED1, CADITD01, 
                               rulexx, slsxx, k.aicxx, validatexx){
  selectedRDGR0 <- NULL; selectedRDGR1 <- NULL; selectedCadit <- NULL 
  BMSD01 <- rbind(BMSD0, BMSD1) # FOR CADIT
  ### RD#1 
  fo <- paste("yyy~", paste(names(BMSD0), collapse = "+"))
  treatment0  <- lrm(as.formula(fo), data = data.frame(yyy = OUTCOMED0, BMSD0), x = TRUE, y = TRUE)
  treatment0a <- fastbw(treatment0, rule = rulexx, sls = slsxx, k.aic = k.aicxx)
  if (length(treatment0a$names.kept) > 0) selectedRDGR0 = treatment0a$names.kept
  if (validatexx == "y") {
    hhx0 <- validate(treatment0, B = nrbstrapstovalidate, bw = TRUE, estimates = FALSE, rule = rulexx, sls = slsxx, k.aic = k.aicxx)
    if (length(names(which(apply(attr(hhx0, "kept"), 2, sum) > percentageinvalidate * nrbstrapstovalidate)))) {
      selectedRDGR0 = names(which(apply(attr(hhx0, "kept"), 2, sum) > percentageinvalidate * nrbstrapstovalidate))
    }
  }
  ### RD#2
  fo <- paste("yyy~", paste(names(BMSD1), collapse = "+"))
  treatment1 <- lrm(as.formula(fo), data = data.frame(yyy = OUTCOMED1, BMSD1), x = TRUE, y = TRUE)
  treatment1a <- fastbw(treatment1, rule = rulexx, sls = slsxx, k.aic = k.aicxx)
  if (length(treatment1a$names.kept) > 0) selectedRDGR1 = treatment1a$names.kept
  if (validatexx == "y") {
    hhx1 <- validate(treatment1, B = nrbstrapstovalidate, bw = TRUE, estimates = FALSE, rule = rulexx, sls = slsxx, k.aic = k.aicxx)
    if (length(names(which(apply(attr(hhx1, "kept"), 2, sum) > percentageinvalidate * nrbstrapstovalidate)))) {
      selectedRDGR1 <- names(which(apply(attr(hhx1, "kept"), 2, sum) > percentageinvalidate * nrbstrapstovalidate))
    }
  }
  
  ### CADIT
  fo <- paste("yyy~", paste(names(BMSD01), collapse = "+"))
  treatmentC  <- lrm(as.formula(fo), data = data.frame(yyy = CADITD01, BMSD01), x = TRUE, y = TRUE)
  treatmentCa <- fastbw(treatmentC, rule = rulexx, sls = slsxx, k.aic = k.aicxx)
  if (length(treatmentCa$names.kept) > 0) selectedCadit = treatmentCa$names.kept
  if (validatexx == "y") {
    hhxC <- validate(treatmentC, B = nrbstrapstovalidate, bw = TRUE, estimates = FALSE, rule = rulexx, sls = slsxx, k.aic = k.aicxx)
    if (length(names(which(apply(attr(hhxC, "kept"), 2, sum) > percentageinvalidate * nrbstrapstovalidate)))){
      selectedCadit <- names(which(apply(attr(hhxC, "kept"), 2, sum) > percentageinvalidate * nrbstrapstovalidate))
    }
  }
  return(list(SelectedGR0 = selectedRDGR0, SelectedGR1 = selectedRDGR1, SelectedCadit = selectedCadit))
}

# BUILD MODEL
buildmodelRDorCadit <- function(SelBMSGR0, SelBMSGR1, SelBMSCadit, OUTCOMED0, OUTCOMED1, CADITD01, BMSD0, BMSD1){
  BMSD01 <- rbind(BMSD0, BMSD1)  
  # RD
  if (length(SelBMSGR0) > 0) {fo = paste("yyy~", paste(SelBMSGR0, collapse = "+"))} else {fo = paste("yyy~1")}
  model0 = glm(as.formula(fo), family = binomial(link = 'logit'), data = data.frame(yyy = OUTCOMED0, BMSD0))
  
  if (length(SelBMSGR1) > 0) {fo = paste("yyy~", paste(SelBMSGR1, collapse = "+"))} else {fo = paste("yyy~1")}
  model1 = glm(as.formula(fo), family = binomial(link = 'logit'), data = data.frame(yyy = OUTCOMED1, BMSD1))
  
  # CADIT
  if (length(SelBMSCadit) > 0) {fo = paste("yyy~", paste(SelBMSCadit, collapse = "+"))} else {fo = paste("yyy~1")}
  modelC = glm(as.formula(fo), family = binomial(link = 'logit'), data = data.frame(yyy = CADITD01, BMSD01))
  
  return(list(modelRD0 = model0, modelRD1 = model1, modelCadit = modelC))
}

# VALIDATE MODEL
validatemodel <- function(MODEL0, MODEL1, MODELCADIT, BMSV01) {
  # RD 
  PredModel0     <- predict.glm(MODEL0, newdata = BMSV01, type = "response")
  PredModel1     <- predict.glm(MODEL1, newdata = BMSV01, type = "response")
  # CADIT
  PredModelCadit <- predict.glm(MODELCADIT, newdata = BMSV01, type = "response")
  CADITDELTA = (PredModelCadit * 2) - 1
  
  return(list(PREDMODEL0 = PredModel0, PREDMODEL1 = PredModel1, PREDMODELCADIT = CADITDELTA))
}

# CROSSVALIDATION
CV <- function(nrsplits) {
  IntermediateResults1  = array(NA, dim = c(nrsplits, length(outcome), 5))
  IntermediateResults2a = array(NA, dim = c(nrsplits, ncol(marker)))
  IntermediateResults2b = array(NA, dim = c(nrsplits, ncol(marker)))
  IntermediateResults2c = array(NA, dim = c(nrsplits, ncol(marker)))
  for (split in 1:nrsplits) {
    splitsdata <- splitdata(biomarkersxx = marker, trtxx = trt, outcomexx = outcome, caditxx = cadit, propindevpsetxx = propindevpset)
    selectsBMSRDorCADIT <- selectBMSRDorCADIT(BMSD0 = splitsdata$BmsD0, BMSD1 = splitsdata$BmsD1, 
                                              OUTCOMED0 = splitsdata$OutcomeD0, OUTCOMED1 = splitsdata$OutcomeD1, 
                                              CADITD01 = splitsdata$CaditD01, rulexx = rule, slsxx = sls, k.aicxx = k.aic, validatexx = validate)
    if (length(selectsBMSRDorCADIT$SelectedGR0) > 0){IntermediateResults2a[split, 1:length(selectsBMSRDorCADIT$SelectedGR0)] = selectsBMSRDorCADIT$SelectedGR0}
    if (length(selectsBMSRDorCADIT$SelectedGR1) > 0){IntermediateResults2b[split, 1:length(selectsBMSRDorCADIT$SelectedGR1)] = selectsBMSRDorCADIT$SelectedGR1}
    if (length(selectsBMSRDorCADIT$SelectedCadit) > 0){IntermediateResults2c[split, 1:length(selectsBMSRDorCADIT$SelectedCadit)] = selectsBMSRDorCADIT$SelectedCadit}
    buildsmodelRDorCadit <- buildmodelRDorCadit(SelBMSGR0 = selectsBMSRDorCADIT$SelectedGR0, 
                                                SelBMSGR1 = selectsBMSRDorCADIT$SelectedGR1, 
                                                SelBMSCadit = selectsBMSRDorCADIT$SelectedCadit,  
                                                OUTCOMED0 = splitsdata$OutcomeD0, OUTCOMED1 = splitsdata$OutcomeD1, 
                                                CADITD01 = splitsdata$CaditD01, 
                                                BMSD0 = splitsdata$BmsD0, BMSD1 = splitsdata$BmsD1)
    validatesmodel <- validatemodel(MODEL0 = buildsmodelRDorCadit$modelRD0, MODEL1 = buildsmodelRDorCadit$modelRD1, 
                                    MODELCADIT = buildsmodelRDorCadit$modelCadit, BMSV01 = splitsdata$BmsV01)
    IntermediateResults1[split, 1:length(splitsdata$TrtV01), 1:5] = cbind(splitsdata$TrtV01, 
                                                                          splitsdata$OutcomeV01, 
                                                                          validatesmodel$PREDMODEL0, 
                                                                          validatesmodel$PREDMODEL1, 
                                                                          validatesmodel$PREDMODELCADIT)  
  } 
  return(list(IntermediateResults1xx  = IntermediateResults1, IntermediateResults2axx = IntermediateResults2a, 
              IntermediateResults2bxx = IntermediateResults2b, IntermediateResults2cxx = IntermediateResults2c))
}

### PLOTS CODE #####
####################
# FREQUENCY OF SELECTED BIOMARKERS AFTER CROSSVALIDATION: RD AND CADIT
funxie1 <- function(xx,hhx) {length(which(hhx %in% xx))}
chosenMarkersPlotRDF <- function(intermediateresults2a, intermediateresults2b, bmsnames, plotsizexx) {
  hh0 <- as.character(intermediateresults2a)
  hh0 <- hh0[is.na(hh0) == F]
  v0  <- sapply(bmsnames, funxie1, hhx = hh0)
  hh1 <- as.character(intermediateresults2b)
  hh1 <- hh1[is.na(hh1) == F]
  v1  <- sapply(bmsnames, funxie1, hhx = hh1)
  layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = T), widths = c(1, plotsizexx, 1))
  par(mar = c(5, 4, 4, 0))
  barplot(height = -v0, names.arg = c(""), horiz = T, xlab = "nr of times selected", xaxt = "n", main = trtA, xlim = c(-max(v0,v1), 0))
  axis(1, at = -c(0:max(v0)), labels = sort(0:max(v0), decreasing = F))
  par(mar = c(5, 0, 4, 0))
  plot(c(1,7), c(1,length(v0)), xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", axes = F)
  text(x = 1,y = c(1:length(v0)), labels = names(v0), cex = 1.00, adj = 0)
  par(mar = c(5, 0, 4, 2))
  barplot(height = v1, names.arg = c(""), horiz = T, xlab = "nr of times selected", main = trtB, xlim = c(0, max(v0, v1)))
  par(mar = c(5, 4, 4, 2) + 0.1)
  layout(matrix(c(1), nrow = 1, ncol = 1, byrow = T),widths = c(1))
  return(list(v0 = v0, v1 = v1))
}

chosenMarkersPlotCaditF <- function(intermediateresults2c, bmsnames, plotsizexx) {
  hhC <- as.character(intermediateresults2c)
  hhC <- hhC[is.na(hhC) == F]
  vC  <- sapply(bmsnames, funxie1, hhx = hhC)
  layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = T), widths = c(1, plotsizexx, 1))
  par(mar = c(5, 4, 4, 0))
  barplot(height = -vC, names.arg = c(""), horiz = T, xlab = "nr of times selected", xaxt = "n", main = "CADIT: All Groups", xlim = c(-max(vC), 0))
  axis(1, at = -c(0:max(vC)), labels = sort(0:max(vC), decreasing = F))
  par(mar = c(5, 0, 4, 0))
  plot(c(1,7), c(1,length(vC)), xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n", axes = F)
  text(x = 1, y = c(1:length(vC)), labels = names(vC), cex = 1.00, adj = 0)
  par(mar = c(5, 0, 4, 2))
  return(vC = vC)
}

####################
# PREDICTED VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED RISK DIFFERENCE)
# OPTION TO CHOOSE PLOT YOU WANT TO SEE
subpopplotQ <- function(splitnr, IntermediateResults1, nrquantilesxx, splinedfxx, thresholdxx) {
  trtValGr     <- IntermediateResults1[splitnr,,1]
  outcomeValGr <- IntermediateResults1[splitnr,,2]
  outcomeValGr <- outcomeValGr[is.na(trtValGr) == F]
  effOfTrt     <- IntermediateResults1[splitnr,,3] - IntermediateResults1[splitnr,,4]
  effOfTrt     <- effOfTrt[is.na(trtValGr) == F]
  trtValGr     <- trtValGr[is.na(trtValGr) == F]
  if (var(effOfTrt, na.rm = T) > 0) {
    cutpoints  <- unique(c(-1, quantile(effOfTrt, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
    subpops    <- cut(effOfTrt, cutpoints)
    hh2        <- table(trtValGr, outcomeValGr, subpops)[,2,] / apply(table(trtValGr, outcomeValGr, subpops), c(1,3), sum)
    diffs      <- hh2[1,] - hh2[2,]
    Ncutpoints <- as.numeric(cutpoints[-1])
    Ncutpoints <- Ncutpoints[is.na(diffs) == FALSE]
    if(length(Ncutpoints)<= 1){
      plot(c(1:10), rep(0,10), 
           col = NULL,
           xlab = "Benefit Function Quantiles:RD", 
           ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
           ylim = c(-1,1),
           main = paste("Test Set: Split", splitnr))
      text(x = c(3.05, 3), y = c(0,-0.1), 
           labels = c("RD MODEL: Only one cutpoint created", " Unable to create subpop plot"), adj = 0)
    }else{
      diffs      <- diffs[is.na(diffs) == FALSE]
      hh3        <- predict(lm(diffs ~ ns(Ncutpoints, df = splinedfxx)), se.fit = TRUE)
      par(mfrow = c(1,1))
      par(mar=c(4, 5, 4, 1) + 0.1)
      par(oma= c(0,0,0,2))
      plot(Ncutpoints, diffs, 
           col = NULL,  
           xlab = "Predicted Benefit Function RD: Quantiles", 
           ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"), 
           ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T), 
                    max(c(diffs, hh3$fit + 2 * hh3$se.fit), na.rm = T)), 
           main = paste("Split", splitnr))
      abline(h = thresholdxx)
      points(Ncutpoints, diffs)
      lines(Ncutpoints, hh3$fit, col = 1, lty = 1, lwd = 2)
      lines(Ncutpoints, hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
      lines(Ncutpoints, hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
      
    }
  }
  if (var(effOfTrt, na.rm = T) == 0) {
    plot(c(1:10), rep(0,10), 
         col = NULL,
         xlab = "Benefit Function Quantiles:RD", 
         ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
         ylim = c(-1,1),
         main = paste("Test Set: Split", splitnr))
    text(x = c(3.05, 3), y = c(0,-0.1), 
         labels = c("RD MODEL: No variation in Risk Difference", " Probably both models are empty."), adj = 0)
  }
}

# PREDICTED CADIT VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED CADIT) 
# OPTION TO CHOOSE PLOT YOU WANT TO SEE
subpopplotQCADIT <- function(splitnr, IntermediateResults1, nrquantilesxx, splinedfxx, thresholdxx) {
  trtValGr     <- IntermediateResults1[splitnr,,1]
  outcomeValGr <- IntermediateResults1[splitnr,,2]
  outcomeValGr <- outcomeValGr[is.na(trtValGr) == F]
  CaditDelta   <- IntermediateResults1[splitnr,,5] 
  CaditDelta   <- CaditDelta[is.na(trtValGr) == F]
  trtValGr     <- trtValGr[is.na(trtValGr) == F]
  if (var(CaditDelta, na.rm = T) > 0) {
    cutpoints  <- unique(c(-1, quantile(CaditDelta, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
    subpops    <- cut(CaditDelta, cutpoints)
    hh2        <- table(trtValGr, outcomeValGr, subpops)[,2,] / apply(table(trtValGr, outcomeValGr, subpops), c(1,3), sum)
    diffs      <- hh2[1,] - hh2[2,]
    Ncutpoints <- as.numeric(cutpoints[-1])
    Ncutpoints <- Ncutpoints[is.na(diffs) == FALSE]; 
    diffs      <- diffs[is.na(diffs) == FALSE]
    if(length(Ncutpoints) <= 1){
      plot(c(1:10), rep(0,10), 
           col = NULL,
           xlab = "Benefit Function Quantiles: CADIT", 
           ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
           ylim = c(-1,1),
           main = paste("Test Set: Split", splitnr))
      text(x = c(3.05, 3), y = c(0,-0.1), 
           labels = c("CADIT MODEL: Only one cutpoint created", " Unable to create subpop plot"), adj = 0)
    }else{
      hh3        <- predict(lm(diffs ~ ns(Ncutpoints, df = splinedfxx)), se.fit = TRUE)
      par(mfrow = c(1,1))
      par(mar=c(4, 5, 4, 1) + 0.1)
      plot(Ncutpoints, diffs, 
           col = NULL,  
           xlab = "Predicted Benefit Function CADIT: Quantiles", 
           ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"), 
           ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T), 
                    max(c(diffs, hh3$fit + 2 * hh3$se.fit), na.rm = T)), 
           main = paste("Split", splitnr))
      abline(h = thresholdxx)
      points(Ncutpoints, diffs)
      lines(Ncutpoints, hh3$fit, col = 1, lty = 1, lwd = 2)
      lines(Ncutpoints, hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
      lines(Ncutpoints, hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
    }
  }
  if (var(CaditDelta, na.rm = T) == 0) {
    plot(c(1:10), rep(0,10), 
         col = NULL,
         xlab = "Benefit Function Quantiles: CADIT", 
         ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
         ylim = c(-1,1),
         main = paste("Test Set: Split", splitnr))
    text(x = c(3.05, 3), y = c(0,-0.1), 
         labels = c("CADIT MODEL: No variation in Risk Difference", " Probably both models are empty."), adj = 0)
  }
}

# PREDICTED VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED RISK DIFFERENCE) 
# SEE ALL PLOTS IN CROSSVALIDATION
subpopplotQAll <- function(nrsplitsxx, IntermediateResults1, nrquantilesxx, splinedfxx, thresholdxx){
  par(mfrow = c(3,3))
  for (splitnr in 1:nrsplitsxx) {
    trtValGr     <- IntermediateResults1[splitnr,,1]
    outcomeValGr <- IntermediateResults1[splitnr,,2]
    outcomeValGr <- outcomeValGr[is.na(trtValGr) == F]
    effOfTrt     <- IntermediateResults1[splitnr,,3] - IntermediateResults1[splitnr,,4]
    effOfTrt     <- effOfTrt[is.na(trtValGr) == F]
    trtValGr     <- trtValGr[is.na(trtValGr) == F]
    if (var(effOfTrt, na.rm = T) > 0) {
      cutpoints  <- unique(c(-1, quantile(effOfTrt, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
      subpops    <- cut(effOfTrt, cutpoints)
      hh2        <- table(trtValGr, outcomeValGr, subpops)[,2,] / apply(table(trtValGr, outcomeValGr, subpops), c(1,3), sum)
      diffs      <- hh2[1,] - hh2[2,]
      Ncutpoints <- as.numeric(cutpoints[-1])
      Ncutpoints <- Ncutpoints[is.na(diffs) == FALSE] 
      if(length(Ncutpoints) <= 1){
        plot(c(1:10), rep(0,10), 
             col = NULL,
             xlab = "Benefit Function Quantiles: CADIT", 
             ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
             ylim = c(-1,1),
             main = paste("Test Set: Split", splitnr))
        text(x = c(1.05, 1), y = c(0.2,-0.2), 
             labels = c("CADIT MODEL: Only one cutpoint created", " Unable to create subpop plot"), adj = 0)
      }else{
        diffs    <- diffs[is.na(diffs) == FALSE]
        hh3      <- predict(lm(diffs ~ ns(Ncutpoints, df = splinedfxx)), se.fit = TRUE)
        plot(Ncutpoints, diffs, 
             col = NULL,  
             xlab = "Predicted Benefit Function RD: Quantiles", 
             ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"), 
             ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T), 
                      max(c(diffs, hh3$fit + 2 * hh3$se.fit), na.rm = T)), 
             main = paste("Split", splitnr))
        abline(h = thresholdxx)
        points(Ncutpoints, diffs)
        lines(Ncutpoints, hh3$fit, col = 1, lty = 1, lwd = 2)
        lines(Ncutpoints, hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
        lines(Ncutpoints, hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
      }
    }
    if (var(effOfTrt, na.rm = T) == 0) {
      plot(c(1:10), rep(0,10), 
           col = NULL,
           xlab = "Benefit Function Quantiles:RD", 
           ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
           ylim = c(-1,1),
           main = paste("Test Set: Split", splitnr))
      text(x = c(1.05, 1), y = c(0.2,-0.2), 
           labels = c("RD MODEL: No variation in Risk Difference", " Probably both models are empty."), adj = 0)
    }
  }
}

# PREDICTED CADIT VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED CADIT) 
# SEE ALL PLOTS IN CROSSVALIDATION
subpopplotQAllCADIT <- function(nrsplitsxx, IntermediateResults1, nrquantilesxx, splinedfxx, thresholdxx) {
  par(mfrow = c(3,3))
  for (splitnr in 1:nrsplitsxx) {
    trtValGr     <- IntermediateResults1[splitnr,,1]
    outcomeValGr <- IntermediateResults1[splitnr,,2]
    outcomeValGr <- outcomeValGr[is.na(trtValGr) == F]
    CaditDelta   <- IntermediateResults1[splitnr,,5] 
    CaditDelta   <- CaditDelta[is.na(trtValGr) == F]
    trtValGr     <- trtValGr[is.na(trtValGr) == F]
    if (var(CaditDelta, na.rm = T) > 0) {
      cutpoints  <- unique(c(-1, quantile(CaditDelta, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
      subpops    <- cut(CaditDelta, cutpoints)
      hh2        <- table(trtValGr, outcomeValGr, subpops)[,2,] / apply(table(trtValGr, outcomeValGr, subpops), c(1,3), sum)
      diffs      <- hh2[1,] - hh2[2,]
      Ncutpoints <- as.numeric(cutpoints[-1])
      Ncutpoints <- Ncutpoints[is.na(diffs) == FALSE] 
      if(length(Ncutpoints) <= 1){
        plot(c(1:10), rep(0,10), 
             col = NULL,
             xlab = "Benefit Function Quantiles: CADIT", 
             ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
             ylim = c(-1,1),
             main = paste("Test Set: Split", splitnr))
        text(x = c(1.05, 1), y = c(0.2,-0.2), 
             labels = c("CADIT MODEL: Only one cutpoint created", " Unable to create subpop plot"), adj = 0)
      }else{
        diffs      <- diffs[is.na(diffs) == FALSE]
        hh3        <- predict(lm(diffs ~ ns(Ncutpoints, df = splinedfxx)), se.fit = TRUE)
        plot(Ncutpoints, diffs, 
             col = NULL,  
             xlab = "Predicted Benefit Function CADIT: Quantiles", 
             ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"), 
             ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T), 
                      max(c(diffs, hh3$fit + 2 * hh3$se.fit), na.rm = T)), 
             main = paste("Split", splitnr))
        abline(h = thresholdxx)
        points(Ncutpoints, diffs)
        lines(Ncutpoints, hh3$fit, col = 1, lty = 1, lwd = 2)
        lines(Ncutpoints, hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
        lines(Ncutpoints, hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
      }
    }
    
    if (var(CaditDelta, na.rm = T) == 0) {
      plot(c(1:10), rep(0,10), 
           col = NULL,
           xlab = "Benefit Function Quantiles: CADIT", 
           ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
           ylim = c(-1,1),
           main = paste("Test Set: Split", splitnr))
      text(x = c(1.05, 1), y = c(0.2,-0.2), 
           labels = c("CADIT MODEL: No variation in Risk Difference", " Probably both models are empty."), adj = 0)
    }
  }
}

# SUBPOPULATION ANALYSIS USING TREATMENT EFFECT: AVERAGE OF ALL PLOTS - RD MODEL
avgsubpopplotRDQ <- function(nrquantilesxx, Intermedresults1xx, splinedfxx, trt2xx, outcome2xx) {
  par(mfrow = c(1,1))
  par(mar = c(4, 3, 4, 2) + 0.1)
  MeanResults1 <- apply(Intermedresults1xx, c(2,3), mean, na.rm = T)
  MeanResults2 <- apply(MeanResults1, 1, mean, na.rm = T)
  MeanResults1 <- MeanResults1[is.na(MeanResults2) == FALSE,]
  
  MeanResults2 <- data.frame(trt2 = trt2xx, outcome2 = outcome2xx, MeanResults1 = MeanResults1[,c(3,4)])
  trt2xx       <- MeanResults2[,1]
  outcome2xx   <- MeanResults2[,2]
  trteffectxx  <- MeanResults2[,3] - MeanResults2[,4]
  if (var(trteffectxx, na.rm = T) > 0) {
    cutpoints <- unique(c(-1, quantile(trteffectxx, probs = c((1:nrquantilesxx) /nrquantilesxx), na.rm = TRUE)))
    subpops   <- cut(trteffectxx, cutpoints)
    hh2       <- table(trt2xx, outcome2xx, subpops)[,2,] / apply(table(trt2xx, outcome2xx, subpops), c(1,3), sum)
    diffs     <- hh2[1,] - hh2[2,]
    cutpoints <- as.numeric(cutpoints[-1])
    cutpoints <- cutpoints[is.na(diffs) == FALSE]
    diffs     <- diffs[is.na(diffs) == FALSE]
    hh3       <- predict(lm(diffs ~ ns(cutpoints, df = splinedfxx)), se.fit = TRUE)
    plot(cutpoints, diffs, 
         col = NULL, 
         xlab = "average predicted benefit function quantiles", 
         ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"), 
         ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T), 
                  max(c(diffs, hh3$fit + 2 * hh3$se.fit), na.rm = T)),
         main = "Averaged over all test sets: RD Model")
    abline(h = 0)
    points(cutpoints,diffs)
    lines(cutpoints,hh3$fit,col = 1,lty = 1,lwd = 2)
    lines(cutpoints,hh3$fit + 2 * hh3$se.fit,col = 1,lty = 2,lwd = 1.5)
    lines(cutpoints,hh3$fit - 2 * hh3$se.fit,col = 1,lty = 2,lwd = 1.5)
  }
  if (var(trteffectxx, na.rm = T) == 0) {
    plot(c(1:10), rep(0,10),
         col = NULL,
         ylim = c(-1,1), 
         xlab = "benefit function quantiles",
         ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
         text(x = c(1.05, 1), y = c(0.2,-0.2),
              labels = c("RD MODEL: No variation in Risk Difference", " Probably both models are empty."), adj = 0))
  }
}

# SUBPOPULATION ANALYSIS USING TREATMENT EFFECT: AVERAGE OF ALL PLOTS - RD MODEL
avgsubpopplotCADITQ <- function(nrquantilesxx, Intermedresults1xx, splinedfxx, trt2xx, outcome2xx) {
  par(mfrow = c(1,1))
  par(mar = c(4, 3, 4, 2) + 0.1)
  MeanResults1 <- apply(Intermedresults1xx, c(2,3), mean, na.rm = T)
  MeanResults2 <- apply(MeanResults1, 1, mean, na.rm = T)
  MeanResults1 <- MeanResults1[is.na(MeanResults2) == FALSE,]
  MeanResults2 <- data.frame(trt2 = trt2xx, outcome2 = outcome2xx, MeanResults1 = MeanResults1[,5])
  trt2xx       <- MeanResults2[,1]
  outcome2xx   <- MeanResults2[,2]
  DeltaCadit   <- MeanResults2[,3]
  if (var(DeltaCadit, na.rm = T) > 0) {
    cutpoints <- unique(c(-1, quantile(DeltaCadit, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
    subpops   <- cut(DeltaCadit, cutpoints)
    hh2       <- table(trt2xx, outcome2xx, subpops)[,2,] / apply(table(trt2xx, outcome2xx, subpops), c(1,3), sum)
    diffs     <- hh2[1,] - hh2[2,]
    cutpoints <- as.numeric(cutpoints[-1])
    cutpoints <- cutpoints[is.na(diffs) == FALSE]
    diffs     <- diffs[is.na(diffs) == FALSE]
    hh3       <- predict(lm(diffs ~ ns(cutpoints, df = splinedfxx)), se.fit = TRUE)
    plot(cutpoints, diffs, 
         col = NULL, 
         xlab = "average predicted benefit function quantiles", 
         ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"), 
         ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T), 
                  max(c(diffs, hh3$fit + 2 * hh3$se.fit),na.rm = T)),
         main = "Averaged over all test sets: CADIT Model")
    abline(h = 0)
    points(cutpoints,diffs)
    lines(cutpoints,hh3$fit,col = 1,lty = 1,lwd = 2)
    lines(cutpoints,hh3$fit + 2 * hh3$se.fit,col = 1,lty = 2,lwd = 1.5)
    lines(cutpoints,hh3$fit - 2 * hh3$se.fit,col = 1,lty = 2,lwd = 1.5)
  }
  if (var(DeltaCadit, na.rm = T) == 0) {
    plot(c(1:10), rep(0,10),
         col = NULL,
         ylim = c(-1,1), 
         xlab = "benefit function quantiles",
         ylab = paste("Observed Risk Difference(",trtA, "-", trtB,")"),
         text(x = c(1.05, 1), y = c(0.2,-0.2),
              labels = c("CADIT MODEL: No variation in Risk Difference", " Probably both models are empty."), adj = 0))
  }
}

# TAIL ORIENTED STEPP PLOT: RD - SEE INDIVIDUAL PLOTS
tailOrientedSTEPPPlotF <- function(splitnr, intermediateResults1xx, nrquantilesxx, splinedfxx){
  trt2xx      <- intermediateResults1xx[splitnr,,1]
  outcome2xx  <- intermediateResults1xx[splitnr,,2]
  outcome2xx  <- outcome2xx[is.na(trt2xx)==F]
  trteffectxx <- intermediateResults1xx[splitnr,,3] - intermediateResults1xx[splitnr,,4]
  trteffectxx <- trteffectxx[is.na(trt2xx)==F]
  trt2xx      <- trt2xx[is.na(trt2xx)==F]
  
  if (var(trteffectxx, na.rm = T) > 0) {
    cutpoints <- unique(c(-1, quantile(trteffectxx, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
    subpops   <- as.numeric(cut(trteffectxx, cutpoints))
    ggggg <- c()
    for (jj in 1:(nrquantilesxx - 1)) {
      hh10  <- table(trt2xx[(subpops <= jj)], outcome2xx[(subpops <= jj)])
      hh11  <- hh10[,2] /apply(hh10, 1, sum)
      ggggg <- rbind(ggggg, c(jj, hh11[1] - hh11[2]))
      hh10  <- table(trt2xx[(subpops > jj)], outcome2xx[(subpops > jj)])
      hh11  <- hh10[,2] / apply(hh10,1,sum)
      ggggg <- rbind(ggggg, c(jj+10, hh11[1]- hh11[2]))
    }
    jj        <- nrquantilesxx
    hh10      <- table(trt2xx, outcome2xx)
    hh11      <- hh10[,2] / apply(hh10,1,sum)
    ggggg     <- rbind(ggggg,c(jj,hh11[1]-hh11[2]))
    rangorde  <- order(ggggg[,1])
    ggggg     <- ggggg[rangorde,]
    weights   <- c(1:nrquantilesxx, (nrquantilesxx - 1):1)
    hh3       <- predict(lm(ggggg[,2] ~ ns(ggggg[,1], df = splinedfxx), weights = weights), se.fit = TRUE)
    maxy      <- max(ggggg[,2], hh3$fit + 2 * hh3$se.fit, na.rm = TRUE)
    miny      <- min(ggggg[,2], hh3$fit - 2 * hh3$se.fit, na.rm = TRUE)
    par(mar = c(4,3,4,2) + 0.1)
    plot(ggggg[,1], ggggg[,2],
         ylim = c(miny, maxy), 
         xlim = c(min(ggggg[,1], na.rm = T), max(ggggg[,1], na.rm = T)),
         xlab = "cumulative predicted benefit function quantiles",
         ylab = "",
         main = "Tail Oriented STEPP", xaxt = "n")
    axis(1, at = 1:(nrquantilesxx + (nrquantilesxx - 1)), labels = weights)
    abline(h = 0, lty = 1)
    abline(h = ggggg[(ggggg[,1] == nrquantilesxx),],lty = 3, lwd = 2, col = 2)
    lines(ggggg[,1], hh3$fit, col = 1, lty = 1,lwd = 2)
    lines(ggggg[,1], hh3$fit + 2 * hh3$se.fit, col = 1,lty = 2, lwd = 1.5)
    lines(ggggg[,1], hh3$fit - 2 * hh3$se.fit, col = 1,lty = 2, lwd = 1.5)
  }
  if (var(trteffectxx, na.rm = T) ==0) {
    plot(c(1:10), rep(0,10), 
         col = NULL,
         ylim = c(-1,1),
         xlab = "predicted benefit function quantiles",
         ylab = "",
         main = "Tail Oriented STEPP")
    text(x = c(1,1), y = c(0,-0.1),labels = c("There is no variation in the risk-differences;", "   probably both prediction models are empty."), adj = 0)
  }
}

# TAIL ORIENTED STEPP PLOT: RD - SEE ALL PLOTS
tailOrientedSTEPPPlotAllF <- function(nrsplitsxx, intermediateResults1xx, nrquantilesxx, splinedfxx){
  par(mfrow = c(3,3))
  for (splitnr in 1:nrsplitsxx) {
    trt2xx      <- intermediateResults1xx[splitnr,,1]
    outcome2xx  <- intermediateResults1xx[splitnr,,2]
    outcome2xx  <- outcome2xx[is.na(trt2xx)==F]
    trteffectxx <- intermediateResults1xx[splitnr,,3] - intermediateResults1xx[splitnr,,4]
    trteffectxx <- trteffectxx[is.na(trt2xx)==F]
    trt2xx      <- trt2xx[is.na(trt2xx)==F]
    
    if (var(trteffectxx, na.rm = T) > 0) {
      cutpoints <- unique(c(-1, quantile(trteffectxx, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
      subpops   <- as.numeric(cut(trteffectxx, cutpoints))
      ggggg <- c()
      for (jj in 1:(nrquantilesxx - 1)) {
        hh10  <- table(trt2xx[(subpops <= jj)], outcome2xx[(subpops <= jj)])
        hh11  <- hh10[,2] /apply(hh10, 1, sum)
        ggggg <- rbind(ggggg, c(jj, hh11[1] - hh11[2]))
        hh10  <- table(trt2xx[(subpops > jj)], outcome2xx[(subpops > jj)])
        hh11  <- hh10[,2] / apply(hh10,1,sum)
        ggggg <- rbind(ggggg, c(jj+10, hh11[1]- hh11[2]))
      }
      jj        <- nrquantilesxx
      hh10      <- table(trt2xx, outcome2xx)
      hh11      <- hh10[,2] / apply(hh10,1,sum)
      ggggg     <- rbind(ggggg,c(jj,hh11[1]-hh11[2]))
      rangorde  <- order(ggggg[,1])
      ggggg     <- ggggg[rangorde,]
      weights   <- c(1:nrquantilesxx, (nrquantilesxx - 1):1)
      hh3       <- predict(lm(ggggg[,2] ~ ns(ggggg[,1], df = splinedfxx), weights = weights), se.fit = TRUE)
      maxy      <- max(ggggg[,2], hh3$fit + 2 * hh3$se.fit, na.rm = TRUE)
      miny      <- min(ggggg[,2], hh3$fit - 2 * hh3$se.fit, na.rm = TRUE)
      par(mar = c(4,3,4,2) + 0.1)
      plot(ggggg[,1], ggggg[,2],
           ylim = c(miny, maxy), 
           xlim = c(min(ggggg[,1], na.rm = T), max(ggggg[,1], na.rm = T)),
           xlab = "cumulative predicted benefit function quantiles",
           ylab = "",
           main = paste("Tail Oriented STEPP", splitnr), xaxt = "n")
      axis(1, at = 1:(nrquantilesxx + (nrquantilesxx - 1)), labels = weights)
      abline(h = 0, lty = 1)
      abline(h = ggggg[(ggggg[,1] == nrquantilesxx),],lty = 3, lwd = 2, col = 2)
      lines(ggggg[,1], hh3$fit, col = 1, lty = 1,lwd = 2)
      lines(ggggg[,1], hh3$fit + 2 * hh3$se.fit, col = 1,lty = 2, lwd = 1.5)
      lines(ggggg[,1], hh3$fit - 2 * hh3$se.fit, col = 1,lty = 2, lwd = 1.5)
    }
    if (var(trteffectxx, na.rm = T) ==0) {
      plot(c(1:10), rep(0,10), 
           col = NULL,
           ylim = c(-1,1),
           xlab = "predicted benefit function quantiles",
           ylab = "",
           main = paste("Tail Oriented STEPP", splitnr))
      text(x = c(1,1), y = c(0,-0.1),labels = c("There is no variation in the risk-differences;", "   probably both prediction models are empty."), adj = 0)
    }
  }
}

# TAIL ORIENTED STEPP PLOT: RD - AVERAGE
avgtailOrientedSTEPPPlotF <- function(nrquantilesxx, Intermedresults1xx, trt2xx,
                                      outcome2xx, splinedfxx){
  par(mfrow = c(1,1))
  par(mar = c(4, 3, 4, 2) + 0.1)
  MeanResults1 <- apply(Intermedresults1xx, c(2,3), mean, na.rm = T)
  MeanResults2 <- apply(MeanResults1, 1, mean, na.rm = T)
  MeanResults1 <- MeanResults1[is.na(MeanResults2) == FALSE,]
  MeanResults2 <- data.frame(trt2 = trt2xx, outcome2 = outcome2xx, MeanResults1 = MeanResults1[,c(3,4)])
  trt2xx       <- MeanResults2[,1]
  outcome2xx   <- MeanResults2[,2]
  trteffectxx  <- MeanResults2[,3] - MeanResults2[,4]
  if (var(trteffectxx, na.rm = T) > 0) {
    cutpoints <- unique(c(-1, quantile(trteffectxx, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
    subpops   <- as.numeric(cut(trteffectxx, cutpoints))
    ggggg <- c()
    for (jj in 1:(nrquantilesxx - 1)) {
      hh10  <- table(trt2xx[(subpops <= jj)], outcome2xx[(subpops <= jj)])
      hh11  <- hh10[,2] /apply(hh10, 1, sum)
      ggggg <- rbind(ggggg, c(jj, hh11[1] - hh11[2]))
      hh10  <- table(trt2xx[(subpops > jj)], outcome2xx[(subpops > jj)])
      hh11  <- hh10[,2] / apply(hh10,1,sum)
      ggggg <- rbind(ggggg, c(jj+10, hh11[1]- hh11[2]))
    }
    jj        <- nrquantilesxx
    hh10      <- table(trt2xx, outcome2xx)
    hh11      <- hh10[,2] / apply(hh10,1,sum)
    ggggg     <- rbind(ggggg,c(jj,hh11[1]-hh11[2]))
    rangorde  <- order(ggggg[,1])
    ggggg     <- ggggg[rangorde,]
    weights   <- c(1:nrquantilesxx, (nrquantilesxx - 1):1)
    hh3       <- predict(lm(ggggg[,2] ~ ns(ggggg[,1], df = splinedfxx), weights = weights), se.fit = TRUE)
    maxy      <- max(ggggg[,2], hh3$fit + 2 * hh3$se.fit, na.rm = TRUE)
    miny      <- min(ggggg[,2], hh3$fit - 2 * hh3$se.fit, na.rm = TRUE)
    par(mar = c(4,3,4,2) + 0.1)
    plot(ggggg[,1], ggggg[,2],
         ylim = c(miny, maxy), 
         xlim = c(min(ggggg[,1], na.rm = T), max(ggggg[,1], na.rm = T)),
         xlab = "cumulative predicted benefit function quantiles",
         ylab = "",
         main = "Average Tail Oriented STEPP", xaxt = "n")
    axis(1, at = 1:(nrquantilesxx + (nrquantilesxx - 1)), labels = weights)
    abline(h = 0, lty = 1)
    abline(h = ggggg[(ggggg[,1] == nrquantilesxx),],lty = 3, lwd = 2, col = 2)
    lines(ggggg[,1], hh3$fit, col = 1, lty = 1,lwd = 2)
    lines(ggggg[,1], hh3$fit + 2 * hh3$se.fit, col = 1,lty = 2, lwd = 1.5)
    lines(ggggg[,1], hh3$fit - 2 * hh3$se.fit, col = 1,lty = 2, lwd = 1.5)
  }
  if (var(trteffectxx, na.rm = T) ==0) {
    plot(c(1:10), rep(0,10), 
         col = NULL,
         ylim = c(-1,1),
         xlab = "predicted benefit function quantiles",
         ylab = "",
         main = "")
    text(x = c(1,1), y = c(0,-0.1),labels = c("There is no variation in the risk-differences;", "   probably both prediction models are empty."), adj = 0)
  }
}

# PLOT FOR CONTROL GROUP BIOMARKER COMBINATION 
biomarkercombinationplotF <- function(splitnr, intermediateResults2axx,
                                      intermediateResults1xx, nrquantilesxx, splinedfxx) {
  selectedGR0 <- intermediateResults2axx[splitnr,]
  selectedGR0 <- selectedGR0[is.na(selectedGR0) == F]
  if (length(selectedGR0) >= 1) {
    trt2xx     <- intermediateResults1xx[splitnr,,1]
    outcome2xx <- intermediateResults1xx[splitnr,,2]
    outcome2xx <- outcome2xx[is.na(trt2xx) == F]
    probsxx    <- intermediateResults1xx[splitnr,,3]
    probsxx    <- probsxx[is.na(trt2xx) == F]
    trt2xx     <- trt2xx[is.na(trt2xx) == F]
    cutpoints  <- unique(c(0, quantile(probsxx, probs = c((1:nrquantilesxx)/nrquantilesxx), na.rm = TRUE)))
    subpops    <- cut(probsxx, cutpoints)
    hh2        <- table(trt2xx, outcome2xx, subpops)[,2,] / apply(table(trt2xx, outcome2xx, subpops), c(1,3), sum)
    diffs      <- hh2[1,] - hh2[2,]
    hh3        <- predict(lm(diffs ~ ns(as.numeric(cutpoints[-1]), df = splinedfxx)), se.fit = TRUE)
    plot(cutpoints[-1], diffs,
         col = NULL,
         main = trtA,
         xlab = paste("Quantiles of", trtA, ":", selectedGR0),
         ylab = paste("Predicted Risk Difference(",trtA, "-", trtB,")"),
         ylim = c(min(c(diffs, hh3$fit - 2 * hh3$se.fit), na.rm = T),
                  max(c(diffs, hh3$fit + 2 * hh3$se.fit), na.rm=T)))
    abline(h=0)
    points(as.numeric(cutpoints[-1]), diffs)
    lines(as.numeric(cutpoints[-1]), hh3$fit, col = 1, lty = 1, lwd = 2)
    lines(as.numeric(cutpoints[-1]), hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
    lines(as.numeric(cutpoints[-1]), hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 1.5)
  }
  if (length(selectedGR0) == 0) {
    plot(c(1:10), rep(0,10),
         col = NULL,
         ylim = c(-1,1),
         xlab = paste("quantiles of", trtA, ":", selectedGR0),
         ylab = paste("Predicted Risk Difference(",trtA, "-", trtB,")"), 
         main = "trtA")
    text(x = 1, y = 0, labels = c("The prediction model for the control group is empty."), adj = 0)
  }
}

# INDIVIDUAL PLOT FOR THE CONTROL GROUP BIOMARKER COMBINATION 
individualbiomarkercombinationplotF <- function(splitnr, intermediateResults2axx, intermediateResults1xx, splinedfxx) {
  selectedGR0 <- intermediateResults2axx[splitnr,]
  selectedGR0 <- selectedGR0[is.na(selectedGR0) == F]
  if (length(selectedGR0) >= 1) {
    trt2xx     <- intermediateResults1xx[splitnr,,1]
    probs0xx   <- intermediateResults1xx[splitnr,,3]
    probs0xx   <- probs0xx[is.na(trt2xx) == F]
    probs1xx   <- intermediateResults1xx[splitnr,,4]
    probs1xx   <- probs1xx[is.na(trt2xx) == F]
    trt2xx     <- trt2xx[is.na(trt2xx) == F]
    trteffectxx <- probs0xx - probs1xx
    rangordexx  <- order(probs0xx)
    trteffectxx <- trteffectxx[rangordexx]
    probs0xx    <- probs0xx[rangordexx]
    probs0xx <- probs0xx[is.na(probs0xx) == F]
    trteffectxx <- probs0xx[is.na(probs0xx) == F]
    hh3         <- predict(lm(trteffectxx ~ ns(probs0xx, df = splinedfxx)), se.fit = TRUE)
    plot(probs0xx, trteffectxx,
         xlab = paste("predicted Risk Biomarker =",  selectedGR0),
         ylab = paste("Predicted Risk Difference(",trtA, "-", trtB,")"),
         main = "Control Group")
    abline(h = 0)
    lines(probs0xx, hh3$fit, col = 1, lty = 1, lwd = 3)
    lines(probs0xx, hh3$fit + 2 * hh3$se.fit, col = 1, lty = 2, lwd = 2.5)
    length(probs0xx)
    length(hh3$fit)
    lines(probs0xx, hh3$fit - 2 * hh3$se.fit, col = 1, lty = 2, lwd = 2.5)
  }
  if (length(selectedGR0) == 0) {
    plot(c(1:10), rep(0,10),
         col = NULL,
         ylim = c(-1,1),
         xlab = "Predicted Risk - control-group-biomarker-combination values",
         ylab = paste("Predicted Risk Difference(",trtA, "-", trtB,")"))
    text(x = 1, y = 0, labels = c("The prediction model for the control group is empty."), adj = 0)
  }
}

# AVERAGE BRIER AND AUC STATISTICS
fitstats <- function(nrsplitsxx, intermediateResults1xx, splinedfxx) {
  brier <- c()
  auc   <- c()
  for (splitnr in 1:nrsplitsxx){
    trt2xx     <- intermediateResults1xx[splitnr,,1]
    outcome2xx <- intermediateResults1xx[splitnr,,2]
    outcome2xx <- outcome2xx[is.na(trt2xx) == F]
    probs0xx   <- intermediateResults1xx[splitnr,,3]
    probs0xx   <- probs0xx[is.na(trt2xx) == F]
    probs1xx   <- intermediateResults1xx[splitnr,,4]
    probs1xx   <- probs1xx[is.na(trt2xx) == F]
    trt2xx     <- trt2xx[is.na(trt2xx) == F]
    brier[splitnr] <- (sum((outcome2xx[trt2xx == 0] - probs0xx[trt2xx == 0])^2, na.rm = T) + 
                         sum((outcome2xx[trt2xx == 1] - probs1xx[trt2xx == 1])^2, na.rm = T))/length(trt2xx)
    probsxx    <- probs0xx * (trt2xx == 0) + probs1xx  * (trt2xx == 1)
    cc         <- roc(outcome2xx ~ probsxx, ci = T, plot = FALSE, auc = T, na.rm = T)
    auc        <- rbind(auc, as.numeric(cc$ci))
  }
  colnames(auc) <- c("lower", "auc", "upper")
  return(list(brier = brier, auc = auc))
}

# PERFORMANCE MEASURES
summarymeasures <- function(nrsplitsxx, intermediateResults1xx) {   
  markerpositivityrate                <- c() 
  averagebenefitoftreatment           <- c()
  averagebenefitofnotreatment         <- c()
  decreaseinoutcomewithmarkerstrategy <- c()
  for (splitnr in 1:nrsplitsxx){
    trt2xx     <- intermediateResults1xx[splitnr,,1]
    outcome2xx <- intermediateResults1xx[splitnr,,2]
    outcome2xx <- outcome2xx[is.na(trt2xx) == F]
    probs0xx   <- intermediateResults1xx[splitnr,,3]
    probs0xx   <- probs0xx[is.na(trt2xx) == F]
    probs1xx   <- intermediateResults1xx[splitnr,,4]
    probs1xx   <- probs1xx[is.na(trt2xx) == F]
    trt2xx     <- trt2xx[is.na(trt2xx) == F]
    trteffect  <- probs0xx - probs1xx
    markerpositivityrate[splitnr] = sum((trteffect <= 0), na.rm = TRUE)/sum((trteffect <= 9999), na.rm = TRUE)
    averagebenefitoftreatment[splitnr] = mean(trteffect[(trteffect <= 0)], na.rm = TRUE)
    averagebenefitofnotreatment[splitnr] = mean(trteffect[(trteffect >= 0)], na.rm = TRUE)
    decreaseinoutcomewithmarkerstrategy[splitnr] = markerpositivityrate[splitnr] * averagebenefitoftreatment[splitnr]
  }
  return(cbind(markerpositivityrate, averagebenefitoftreatment, averagebenefitofnotreatment, decreaseinoutcomewithmarkerstrategy))
}

# PERFORMANCE MEASURES PLOTS
plotSummaryMeasures <- function(summeasxx, nrsplitsxx){
  
  layout(matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = T))
  
  plot(summeasxx[,1], 
       xlab = paste("Nr of Splits:", nrsplits),
       ylab = "marker positivity rate",
       ylim = c(0,1))
  abline(h = mean(summeasxx[,1], na.rm = TRUE), lty = 2, col = 2, lwd = 2)
  
  plot(summeasxx[,2], 
       xlab = paste("Nr of Splits:", nrsplits),
       ylab = "average benefit of treatment",
       ylim = c(min(summeasxx[,2], na.rm=T) - 0.1,max(summeasxx[,2], na.rm = T) + 0.1))
  abline(h = mean(summeasxx[,2], na.rm=TRUE), lty = 2, col = 2, lwd = 2)
  abline(h = 0, lty = 1, col = 1, lwd = 1)
  
  plot(summeasxx[,3],
       xlab = paste("Nr of Splits:", nrsplits),
       ylab = "average benefit of no treatment",
       ylim = c(min(summeasxx[,3], na.rm = T) - 0.1, max(summeasxx[,3], na.rm = T) + 0.1))
  abline(h = mean(summeasxx[,3], na.rm = TRUE), lty = 2, col = 2, lwd = 2)
  abline(h = 0, lty = 1, col = 1, lwd = 1)
  
  plot(summeasxx[,4],
       xlab = paste("Nr of Splits:", nrsplits),
       ylab = "outcome change with a model-based strategy",
       ylim = c(min(summeasxx[,4], na.rm=T) - 0.1, max(summeasxx[,4], na.rm = T) + 0.1))
  abline(h = mean(summeasxx[,4], na.rm = TRUE), lty = 2, col = 2, lwd = 2)
  abline(h = 0, lty = 1, col = 1, lwd = 1)
  
  title("population impact measures estimated in the test sets", outer = T, line =-2)
  
}


