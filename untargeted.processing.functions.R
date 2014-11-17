getBestMS2PeakNum <- function(xset, groupID) {
  peakGroup <- as.data.frame(xset@peaks[xset@groupidx[[groupID]],]);
  numClasses <- length(table(xset@phenoData));
  samplesPerClass <- mean(table(xset@phenoData));
  
  peakGroup$sampleClass <- trunc((peakGroup$sample-1)/samplesPerClass+1);
  peakGroup$replicateNum <- apply(peakGroup,1,function(y) {y[11]-samplesPerClass*(y[12]-1)})
  class <- 1:numClasses;
  #Create list of peak intensities divided by sample class, for sample classes with peaks in every replicate
  tL <- lapply(1:numClasses,function(y) {
    if(sum(1:samplesPerClass %in% peakGroup[(peakGroup$sampleClass == y) & (peakGroup$maxo != 0),"replicateNum"]) == samplePerClass) {
      return(peakGroup[peakGroup$sampleClass==y,"maxo"]);
    } else {
      return(NA);
    } 
  } )
  maxClass <- which.max(sapply(tL,min)); #find sample class with highest minimum maxo value
  if(invalid(maxClass)) {
    return(NA); #in case none of the classes have representatives in every group: shouldn't be a problem for our purposes, but makes error detection clear
  } else {
    testRef <- which(peakGroup$maxo == max(tL[[maxClass]]), arr.ind=TRUE);
    testRef <- testRef[1]; #in case there's a tie
    
    peakNumber <- xset@groupidx[[groupID]][testRef];
    return(peakNumber);
  }
}

#This function only works if xset has a constant number of samples per sample class

mergeMaxoandFormat <- function (annDiffReport, xset, tag, sampleSequence = 1:sum(phenoData)) {
  #initializing merged dataframe
  maxo.Data <- as.data.frame(cbind(xset@groups,groupval(xset,value="maxo",method="maxint")));
  annDiffReport$ID <- rownames(annDiffReport);
  calcFrame <- merge(annDiffReport,maxo.Data,by="mzmed");
  rownames(calcFrame) <- calcFrame$ID;
  
  # Setting xset-specific parameters that determine location of column variables
  metlin <- "metlin" %in% names(calcFrame);
  phenoData <- table(xset@phenoData);
  if (var(phenoData) != 0) stop("xset does not have consistent number of samples per class");
  numClasses <- length(phenoData);
  samplesPerClass <- mean(phenoData);
  
  # Removing redundant columns
  lengthOfGarbage <- 6 + numClasses;
  startOfGarbage <- 12 + numClasses + as.integer(metlin) + (sum(phenoData)) + 4 + 1;
  for (i in 1:(lengthOfGarbage)) {calcFrame[,startOfGarbage] <- NULL}
  
  #Renaming columns to improve clarity
  colnames(calcFrame)[7:(12 + numClasses)] <- vapply(colnames(calcFrame)[7:(12 + numClasses)], 
                                                     FUN.VALUE = "", function(y) gsub('.x', '', y, fixed = TRUE));
  sampleNames <- formatC(sampleSequence, width=2, format="d", flag="0");
  tag1 <- tag;
  tag2 <- c("into","maxo")
  startOfInto <- 12 + numClasses + as.integer(metlin) + 1;
  startOfMaxo <- startOfInto + (sum(phenoData)) + 4;
  colnames(calcFrame)[startOfInto:(startOfInto + sum(phenoData) - 1)] <- sapply(sampleNames,function(x) paste(c(tag1,x,tag2[1]),collapse='.'));
  colnames(calcFrame)[startOfMaxo:(startOfMaxo + sum(phenoData) - 1)] <- sapply(sampleNames,function(x) paste(c(tag1,x,tag2[2]),collapse='.'));
  
  return(calcFrame);
}

quadFilterC18 <- function(calcFrame, xset) {
  metlin <- "metlin" %in% names(calcFrame);
  phenoData <- table(xset@phenoData);
  if (var(phenoData) != 0) stop("xset does not have consistent number of samples per class");
  numClasses <- length(phenoData);
  samplesPerClass <- mean(phenoData);
  startOfInto <- 12 + numClasses + as.integer(metlin) + 1;
  startOfMaxo <- startOfInto + (sum(phenoData)) + 4;
  
  # A) Feature is not in void volume
  calcFrameA <- subset(calcFrame, rtmed > 300);
  
  # B) Minimum intensity on at least one day is > 8000	
  maxo.Matrix <- data.matrix(calcFrameA[,startOfMaxo:(startOfMaxo + sum(phenoData) -1)]);
  minIntensityPerClass <- t(apply(maxo.Matrix, 1, function(y) sapply(1:numClasses, function(z) min(y[((z - 1) * samplesPerClass + 1):(z * samplesPerClass)]))));
  colnames(minIntensityPerClass) <- sapply(levels(sampclass(xset)), function(y) paste(c(y,"min"),collapse='.'));
  calcFrameB <- data.frame(calcFrameA, minIntensityPerClass);
  calcFrameB$maxMin <- apply(calcFrameB[,(dim(calcFrameA)[2] + 1):(dim(calcFrameA)[2] + dim(minIntensityPerClass)[2])], 1, max);
  calcFrameB$MS2TargetClass <- apply(calcFrameB[,(dim(calcFrameA)[2] + 1):(dim(calcFrameA)[2] + dim(minIntensityPerClass)[2])], 1, function (y) levels(sampclass(xset))[which.max(y)]);
  calcFrameB <- subset(calcFrameB, maxMin > 8000);
  
  # C) peak is real: at least one class should have peaks in all samples
  calcFrameB$classMaxN <- apply(calcFrameB[,13:(12 + numClasses)], 1, max);
  calcFrameC <- subset(calcFrameB, maxPresence == samplesPerClass);
  
  # D) removing isotopic peaks
  calcFrameC$notIsotope <- sapply(calcFrameC[, "isotopes"], function(x) (x == "" | grepl("[M]", x, fixed = TRUE)));
  calcFrameD <- calcFrameC[calcFrameC$notIsotope==TRUE,];
  
  return(calcFrameD)
}

rmNonMetlin <- function(filterFrame, metlinBatch) {
  metlinBatch <- subset(metlinBatch, msms == "YES"); #removing hits without MS2 data from metlin batch file
  filterFrame$ms2_avail <- trunc(filterFrame$mzmed*10^5) %in% trunc(metlinBatch$inputmass*10^5); #apparently %in% doesn't work well with floating point numbers
  filterFrame$multiCharge <- sapply(filterFrame[, "isotopes"], function(x) (grepl("[M]2+", x, fixed = TRUE) | grepl("[M]3+", x, fixed = TRUE) | grepl("[M]4+", x, fixed = TRUE)));
  filterFrame <- subset(filterFrame, multiCharge | ms2_avail);
  return(filterFrame);
}

readyParallel <- function() {
  library(plyr)
  library(foreach)
  library(doParallel)
  
  nodes <- detectCores() - 1
  cl <- makeCluster(nodes)
  registerDoParallel(cl)
}

loadStdLibs <- function() {
  library(xcms);
  library(multtest);
  library(CAMERA);
  library(gtools);
}

installLibs <- function(update = FALSE) {
  libList <- c("xcms","multtest","CAMERA","gtools","plyr","foreach","doParallel");
  for (i in 1:length(libList)) {
    if(!require(package = libList[i], character.only = TRUE)) {
      install.packages(libList[i]);
    }
  }
  if(update) {
    install.packages(libList);
  }
}

updateNateLibs <- function() {
  remove.packages("metab.tools")
  remove.packages("credential")
  remove.packages("creDBle.tools")
  
  library(devtools)
  
  install_github("nathanielguy/metab.tools", username="jlspaldi", auth_token="649a72de733657a2e99ae8f026bb4470dc8e2435") # PRIVATE
  install_github("pattilab/credential", username="jlspaldi", auth_token="649a72de733657a2e99ae8f026bb4470dc8e2435") # PRIVATE
  install_github("pattilab/creDBle.tools", username="jlspaldi", auth_token="649a72de733657a2e99ae8f026bb4470dc8e2435") # PRIVATE
}

loadNateLibs <- function() {
  library(metab.tools)
  library(credential)
  library(creDBle.tools)
}

runXCMS <- function(mode = 4, path, fillp = TRUE) {
  # Load libraries, set file paths
  loadStdLibs();
  originalwd <- getwd();
  setwd(path);
  files <- as.list(list.files(recursive=TRUE, full.names=TRUE));
  
  # Mode-specific settings
  bandwidth1 <- 35;
  bandwidth2 <- 15;
  if ((mode == 1) | (mode == 2)) {
    bandwidth1 <- 30;
    bandwidth2 <- 10;
  }
  
  #using crude version of parallel computing
  readyParallel();
  
  xset_holder <- llply(.data = files, .parallel = TRUE, 
                       .paropts = list(.packages = "xcms"), .progress = "text", function(y) {
                         return(xcmsSet(y, method = "centWave", ppm = 12, peakwidth = c(10,120), 
                                        snthresh = 6, integrate = 2, mzdiff = .01));
                       })
  
  # Actual xcms Calculations
  xset1 <- xset_holder[[1]];
  for (i in 2:length(files)) {
    xset1 <- c(xset1, xset_holder[[i]]); #pasting everything back together into one xcmsSet
  }
  
  xset1 <- group(xset1, mzwid = .025, max = 100, bw = bandwidth1);
  xset2 <- retcor(xset1, method = "obiwarp", profStep = .5);
  xset2 <- group(xset2, bw = bandwidth2, mzwid = .025, max = 100);
  
  if(!fillp) {
    registerDoSEQ();
    return(xset2);
  } else {
    xset3 <- fillPeaks(xset2);
    registerDoSEQ();
    return(xset3);
  }
}