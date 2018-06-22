#######################################################################################
### RUN FUNCTIONS 
#########################################################################################
# checks for sufficient variation of biomarkers
length(which(apply(marker, 2, sd, na.rm = T) <= 0))
check1(trt)
check1(outcome)
table(trt, outcome)

# OVERALL EFFECT COMPARISON
OverallTrtEffect <- OverallTrtEffectF(trtxx = trt, outcomexx = outcome)
OverallTrtEffect

# DESCRIPTIVE PRESENTATION OF MARKERS
describe(marker)
summary(marker)

# UNIVARIATE ANALYSIS
InteractionsBM <- InteractionsF(trtxx = trt, outcomexx = outcome, biomarkersxx = marker)
InteractionsBM
MostSignificantBM <- MostSignificantBMF(resxx = InteractionsBM, biomarkersxx = marker, outcomexx = outcome, trtxx = trt, nbootxx = 100, yaxislimitsxx = yaxislimits)
MostSignificantBM
#############
### PLOTS ###
#############
# UNIVARIATE ANALYSIS
# INTERACTIONS
ORplot <- ORplotF(resxx = InteractionsBM)
# SINGLE MARKER VS OBSERVED RISK DIFFERECE
RDSingleMarkerForAllTypes <- RDSingleMarkerForAllTypesF(selectedbiomarker = "age")
RDSingleMarkerForAllTypes
# SINGLE MARKER VS OBSERVED RISK
RiskIndividualTrtsSingleMarker <-RiskIndTrt1MarkerF(selectedbiomarker = "age")
RiskIndividualTrtsSingleMarker

# MULTIVARIATE ANALYSIS
nrmarkers = length(marker)/2
splitsdata <- splitdata(biomarkersxx = marker, trtxx = trt, outcomexx = outcome, caditxx = cadit, propindevpsetxx = propindevpset)
crossvalidation <- CV(nrsplits = nrsplits)
crossvalidation$IntermediateResults1[1,1:10,]
crossvalidation$IntermediateResults2axx[,1:nrmarkers]
crossvalidation$IntermediateResults2bxx[,1:nrmarkers]
crossvalidation$IntermediateResults2cxx[,1:nrmarkers]
#############
### PLOTS ###
#############
# FREQUENCY OF SELECTED BIOMARKERS AFTER CROSSVALIDATION
chosenMarkersPlotRD  <- chosenMarkersPlotRDF (intermediateresults2a = crossvalidation$IntermediateResults2axx, intermediateresults2b = crossvalidation$IntermediateResults2bxx, bmsnames = names(marker), plotsizexx = plotsize)
chosenMarkersPlotCadit <- chosenMarkersPlotCaditF(intermediateresults2c = crossvalidation$IntermediateResults2cxx, bmsnames = names(marker), plotsizexx = plotsize)
# PREDICTED VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED RISK DIFFERENCE): 
# OPTION TO CHOOSE PLOT YOU WANT TO SEE
subpopplotQ(splitnr = 1, IntermediateResults1 = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf, thresholdxx = threshold)
subpopplotQ(splitnr = 2, IntermediateResults1 = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf, thresholdxx = threshold)
# CADIT VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED RISK DIFFERENCE): 
# OPTION TO CHOOSE PLOT YOU WANT TO SEE
subpopplotQCADIT(splitnr = 1, IntermediateResults1 = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf, thresholdxx = threshold)
subpopplotQCADIT(splitnr = 2, IntermediateResults1 = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf, thresholdxx = threshold)          
# PREDICTED VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED RISK DIFFERENCE) 
# SEE ALL PLOTS IN CROSSVALIDATION
subpopplotQAll(nrsplitsxx = nrsplits, IntermediateResults1 = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf, thresholdxx = threshold)
# CADIT VS OBSERVED RISK DIFFERENCE (SUBPOPS BASED ON PREDICTED RISK DIFFERENCE) 
# SEE ALL PLOTS IN CROSSVALIDATION
subpopplotQAllCADIT(nrsplitsxx = nrsplits, IntermediateResults1 = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf, thresholdxx = threshold)
# SUBPOPULATION ANALYSIS USING TREATMENT EFFECT: AVERAGE OF ALL PLOTS RD
avgsubpopplotRDQ(nrquantilesxx = nrquantiles, Intermedresults1xx = crossvalidation$IntermediateResults1xx, trt2xx = splitsdata$TrtV01, outcome2xx = splitsdata$OutcomeV01, splinedfxx = splinedf)
# SUBPOPULATION ANALYSIS USING TREATMENT EFFECT: AVERAGE OF ALL PLOTS CADIT
avgsubpopplotCADITQ(nrquantilesxx = nrquantiles, Intermedresults1xx = crossvalidation$IntermediateResults1xx, trt2xx = splitsdata$TrtV01, outcome2xx = splitsdata$OutcomeV01, splinedfxx = splinedf)
# TAIL ORIENTED STEPP PLOT: RD - SEE INDIVIDUAL PLOTS
tailOrientedSTEPPPlot <- tailOrientedSTEPPPlotF(splitnr = 1, intermediateResults1xx = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf)
# TAIL ORIENTED STEPP PLOT: RD - SEE ALL PLOTS
tailOrientedSTEPPPlotAll <- tailOrientedSTEPPPlotAllF(nrsplitsxx = nrsplits, intermediateResults1xx = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf)
# TAIL ORIENTED STEPP PLOT: RD - AVERAGE
avgtailOrientedSTEPPPlotF(nrquantilesxx = nrquantiles, Intermedresults1xx = crossvalidation$IntermediateResults1xx, trt2xx = splitsdata$TrtV01, outcome2xx = splitsdata$OutcomeV01, splinedfxx = splinedf)
# PLOT FOR CONTROL GROUP BIOMARKER COMBINATION 
biomarkercombinationplotF(splitnr = 1, intermediateResults2axx = crossvalidation$IntermediateResults2axx, intermediateResults1xx = crossvalidation$IntermediateResults1xx, nrquantilesxx = nrquantiles, splinedfxx = splinedf)
# PREDICTED RISK FOR CONTROL GROUP VS PREDICTED RISK DIFFERENCE - INDIVIDUAL MARKERS
individualbiomarkercombinationplotF(splitnr = 1, intermediateResults2axx = crossvalidation$IntermediateResults2axx, intermediateResults1xx = crossvalidation$IntermediateResults1xx, splinedfxx = splinedf)

# AVERAGE BRIER AND AUC STATISTICS
stats <- fitstats(nrsplitsxx = nrsplits, intermediateResults1xx = crossvalidation$IntermediateResults1xx, splinedfxx = splinedf)
stats 
# Brier 
hist((stats$brier), xlab = "Brier scores in the test sets", main = paste(nrsplits, "data splits"))
c(MeanBrierScore = mean(stats$brier), SDBrierScore = sd(stats$brier))
# AUC
hist(stats$auc[,2], xlab = "C-statistics in the test sets", main = paste(nrsplits, "data splits"))
c(Mean = apply(stats$auc, 2, mean)[2], SD = apply(stats$auc, 2, sd)[2])

# PERFORMANCE MEASURES
summeas <- summarymeasures(nrsplitsxx = nrsplits, intermediateResults1xx = crossvalidation$IntermediateResults1xx)
# PERFORMANCE MEASURES PLOTS
plotSummaryMeasures(summeasxx = summeas, nrsplitsxx = nrsplits)













