make an error   # do not source this file.
##################################################################################################################################################
# the within-subjects analyses for Alan's data, without the permutation testing.
# Averaged timepoints rather than the offsets. "precue" is 4-12 seconds (offsets -4:0), "postcue" is 16-24 seconds (offsets 2:6).
# the permutation tests are contained in ClassificationScript_HalfSplit.R (also in git).
# Jo Etzel, 17 July 2013.
##################################################################################################################################################
# Alan's half-splitting within-subjects analyses; not permutation tests. 

library(e1071);  # R interface to libsvm

rm(list=ls()); where.run <- "cluster";
# rm(list=ls()); where.run <- "Jo"

SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");
ROIS <- c("PFC_mask_native", "BG_LR_CaNaPu_native", "Parietal_mask_native");

#OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5); # timepoints to classify, in TR from the first PAIR1 or PAIR2 in each event.
# in the datatables, eventType follows a pattern: ITI - memset - delay - upgreen - ITI - probe - ITI - ITI . Offset 0 is the first upgreen (or whatever).
# ITI is used as a filler in eventType, indicating a pause, and occurs during trials, not just in between trials.
# also, the eventType column is in real time, not adjusted for any lag in the BOLD.

doSVM <- function(train, test, DO_DEFAULT_SCALING) {  # train <- useTest; test <- useTrain;
  test <- subset(test, select=c(-subID, -run, -TR));  # get rid of non-classify or voxel columns
  train <- subset(train, select=c(-subID, -run, -TR));
  if (colnames(test)[2] != "v1" | colnames(train)[2] != "v1") { stop("v1 not found where expected"); }
  
  fit <- svm(eventType~., data=train, type="C-classification", kernel="linear", cost=1, scale=DO_DEFAULT_SCALING);  
  tree <- table(test$eventType, predict(fit, test));
  if (dim(tree)[2]==1 | dim(tree)[1]==1) { wrT <- 0.5; } else { wrT <- sum(diag(tree))/sum(tree); }
  
  return(wrT);
}

if (where.run == "cluster") {
  inpath <- "~/tmp/input/";
  outpath <- "~/tmp/output/"; 
  
  cA <- commandArgs();   
  num  <- as.numeric(cA[5]);  # specifies the pair to run.
  if (num == 1) { PAIR1 <- "upempty"; PAIR2 <- "upgreen"; } 
  if (num == 2) { PAIR1 <- "upempty"; PAIR2 <- "upred"; } 
  if (num == 3) { PAIR1 <- "upgreen"; PAIR2 <- "upred";  }
}
if (where.run == "NIL") {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
  permpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/permInputFiles/";   # location of 6eachTable.txt, 8eachTable.txt, etc.
  PAIR1 <- "upempty"; PAIR2 <- "upgreen";
}
if (where.run == "Jo") {
  inpath <- "d:/temp/Alan/for_halfSplit/";
  outpath <- "d:/temp/Alan/halfSplit_out/";
  PAIR1 <- "upempty"; PAIR2 <- "upgreen";
}
  
SEEDS <- c(51591, 36414, 56347, 38442, 20176, 51348, 89727, 67106, 23543, 52663);  # 10 random seeds, from sample(1:100000)[1:10]
NUMSPLITS <- 10;    # there will be 10, but we're running them one at a time.
if (PAIR1 == PAIR2) { stop("PAIR1 == PAIR2"); }

DO_RUN_COLUMN_MS <- FALSE;   # these flags indicate which sort of scaling to do.
DO_RUN_COLUMN_SC <- FALSE; 
DO_ROW_SCALING <- FALSE;
DO_DEFAULT_SCALING <- TRUE;
sc.lbl <- "_";  # make a label for the output files showing the type of scaling used for this classification
if (DO_RUN_COLUMN_MS == TRUE) { sc.lbl <- paste(sc.lbl, "RunColMS", sep=""); }
if (DO_RUN_COLUMN_SC == TRUE) { sc.lbl <- paste(sc.lbl, "RunColSC", sep=""); }
if (DO_ROW_SCALING == TRUE) { sc.lbl <- paste(sc.lbl, "RowSc", sep=""); }
if (DO_DEFAULT_SCALING == TRUE) { sc.lbl <- paste(sc.lbl, "DefaultSc", sep=""); 
} else { sc.lbl <- paste(sc.lbl, "Only", sep=""); }

for (ROI in ROIS) {
for (OFFSET in c("precue", "postcue")) {   # OFFSET <- "precue";
  rtbl <- data.frame(array(NA, c(NUMSPLITS*length(SUBS), 6)));   # results table
  colnames(rtbl) <- c("splitNum", "timePoint", "subID", "firstTest", "secondTest", "avgProp");
  ctr <- 1;
  for (sid in 1:length(SUBS)) {  # sid <- 1;
    tbl <- read.table(gzfile(paste(inpath, SUBS[sid], "_", ROI, "_meanSub.gz", sep="")), comment.char=""); # read in the data
    TRS <- unique(tbl$TR);
    if (length(TRS) != 210) { stop("too many or few TRs"); }
    FIRSTVOXEL <- which(colnames(tbl) == "v1");  # column number of the first voxel column; all larger-index columns are voxels.
    # change some of the eventTypes so can classify
    ind <- which(levels(tbl$eventType) == "smain"); levels(tbl$eventType)[ind] <- "upempty";
    ind <- which(levels(tbl$eventType) == "smainnp"); levels(tbl$eventType)[ind] <- "upempty";
    ind <- which(levels(tbl$eventType) == "delayempty"); levels(tbl$eventType)[ind] <- "upempty";
    ind <- which(levels(tbl$eventType) == "rmain"); levels(tbl$eventType)[ind] <- "upred";
    ind <- which(levels(tbl$eventType) == "rmainup"); levels(tbl$eventType)[ind] <- "upred";
    ind <- which(levels(tbl$eventType) == "rmainnp"); levels(tbl$eventType)[ind] <- "upred";
    ind <- which(levels(tbl$eventType) == "delayred"); levels(tbl$eventType)[ind] <- "upred";
    ind <- which(levels(tbl$eventType) == "update"); levels(tbl$eventType)[ind] <- "upgreen";
    ind <- which(levels(tbl$eventType) == "updateop"); levels(tbl$eventType)[ind] <- "upgreen";
    ind <- which(levels(tbl$eventType) == "updatenp"); levels(tbl$eventType)[ind] <- "upgreen";
    ind <- which(levels(tbl$eventType) == "delaygreen"); levels(tbl$eventType)[ind] <- "upgreen";
    
    # find the **start** of each trial that's type pair1 or pair2: offset goes forwards and backwards from each first TR labeled 'PAIR'
    inds <- which(tbl$eventType == PAIR1 | tbl$eventType == PAIR2);
    tinds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
    
    # ************ new stuff *************************
    # generate the examples we want to classify.
    if (exists("do.offsets")) { rm(do.offsets); }
    if (OFFSET == "precue") { do.offsets <- c(-4, -3, -2, -1, 0); }  # old-style offsets corresponding to the two time windows: which TR to average.
    if (OFFSET == "postcue") { do.offsets <- c(2, 3, 4, 5, 6); } 
    base.inds <- inds[tinds];   # start finding data to classify; offset from time point 0: start of each trial

    t0tbl <- data.frame(array(NA, c(length(tinds), ncol(tbl))));  # make an empty table to hold the averaged examples
    colnames(t0tbl) <- colnames(tbl);  # add the column names
    t0tbl[,1:(FIRSTVOXEL-1)] <- tbl[base.inds,1:(FIRSTVOXEL-1)];   # and the label columns
    # and now the averaged voxels
    for (i in 1:length(base.inds)) {  # i <- 1;
      these.inds <- base.inds[i] + do.offsets;
      t0tbl[i,FIRSTVOXEL:ncol(tbl)] <- apply(tbl[these.inds,FIRSTVOXEL:ncol(tbl)], 2, mean);  # actually do the averaging.
    }
    # ************ new stuff *************************
    
    RUNS <- sort(unique(t0tbl$run));
    if (length(RUNS) > 12) { stop("too many runs"); }   # try to catch if the input data is wrong
    
    if (DO_RUN_COLUMN_MS == TRUE) {  # do run-column scaling (each run separately) AFTER separating into offsets, but before separating classes
      for (RUN in RUNS) {   # RUN <- RUNS[1];
          thisRun <- subset(t0tbl, run == RUN);
          others <- subset(t0tbl, run != RUN);
          scaled <- thisRun[,FIRSTVOXEL:ncol(thisRun)];   # get voxel columns; cols 1 to 3 are labels
          scaled <- scale(scaled, center=TRUE, scale=FALSE);  # mean-subtract column-wise
          scaled <- cbind(thisRun[,1:(FIRSTVOXEL-1)], scaled);  # put label columns back on
          colnames(scaled) <- colnames(thisRun);
          t0tbl <- rbind(scaled, others);
      }
      rm(others, scaled, thisRun);   # clean up
    }
    
    if (DO_RUN_COLUMN_SC == TRUE) {  # do run-column scaling (each run separately) AFTER separating into offsets, but before separating classes
      for (RUN in RUNS) {   # RUN <- RUNS[1];
          thisRun <- subset(t0tbl, run == RUN);
          others <- subset(t0tbl, run != RUN);
          scaled <- thisRun[,FIRSTVOXEL:ncol(thisRun)];   # get voxel columns; cols 1 to 3 are labels
          scaled <- scale(scaled, center=TRUE, scale=TRUE);  # scale column-wise
          scaled <- cbind(thisRun[,1:(FIRSTVOXEL-1)], scaled);  # put label columns back on
          colnames(scaled) <- colnames(thisRun);
          t0tbl <- rbind(scaled, others);
      }
      rm(others, scaled, thisRun);   # clean up
    }
    
    if (DO_ROW_SCALING == TRUE) {  # do row scaling (all voxels in each summary volume)
      scaled <- t0tbl[,FIRSTVOXEL:ncol(t0tbl)];   # get voxel columns; cols 1 to 3 are labels
      scaled <- scale(t(scaled));  # scale does columns, so transpose
      t0tbl <- cbind(t0tbl[,1:(FIRSTVOXEL-1)], t(scaled));  # put label columns back on, transpose back
      rm(scaled)
    }
    
    numInTrain <- round(length(RUNS)/2);
    trainRuns <- RUNS[1:numInTrain];
    testRuns <- RUNS[(numInTrain+1):length(RUNS)];
    if (length(intersect(trainRuns, testRuns)) > 0) { stop("length(intersect(trainRuns, testRuns)) > 0"); }
    if (length(trainRuns) + length(testRuns) != length(RUNS)) { stop("length(trainRuns + testRuns) != length(RUNS)"); }
    
    testinds <- which(is.element(t0tbl$run, testRuns));  # rows with test-set runs
    traininds <- (1:nrow(t0tbl))[-testinds];  # all the others are training-set runs
    if (length(testinds) == 0) { stop("no testing data"); }  # can't be empty
    if (length(traininds) == 0) { stop("no training data"); }
    if (length(intersect(testinds, traininds)) > 0) { stop("overlap in test and traininds"); }
    
    t0tbl$eventType <- factor(t0tbl$eventType);   # gets rid of empty factor levels (simplifies balancing)
    allTest <- t0tbl[testinds,];  # subset into training and testing sets
    allTrain <- t0tbl[traininds,];
    
    for (splt in 1:NUMSPLITS) {   # splt <- 1;
      set.seed(SEEDS[splt]);  # set random seed so splits are repeat-able
      # balance number of pair1 and pair2 entries in the training and testing sets by omitting rows of bigger class
      # check if classes are balanced in the training data, and balance if not
      oneCount <- summary(allTrain$eventType)[[1]];
      twoCount <- summary(allTrain$eventType)[[2]];
      if (twoCount == oneCount) {   # balanced already.
        useTrain <- allTrain;
      } else {
        if (oneCount > twoCount) {   # too many of PAIR1
          rows <- sample(which(allTrain$eventType==PAIR1));
          rows <- rows[1:(oneCount-twoCount)];
          useTrain <- allTrain[-rows,];   # remove extra rows by index
        } 
        if (twoCount > oneCount) {  # too many of PAIR2
          rows <- sample(which(allTrain$eventType==PAIR2));
          rows <- rows[1:(twoCount-oneCount)];
          useTrain <- allTrain[-rows,];
        }    
      }
      if (summary(useTrain$eventType)[[1]] != summary(useTrain$eventType)[[2]]) { stop("training data not balanced"); }
      if (summary(useTrain$eventType)[[1]] == 0) { stop("no training data after balancing"); }
      
      # check if classes are balanced in the testing data, and balance if not
      oneCount <- summary(allTest$eventType)[[1]];
      twoCount <- summary(allTest$eventType)[[2]];
      if (twoCount == oneCount) {
        useTest <- allTest;
      } else {
        if (oneCount > twoCount) {   # too many of PAIR1
          rows <- sample(which(allTest$eventType==PAIR1));
          rows <- rows[1:(oneCount-twoCount)];
          useTest <- allTest[-rows,];
        } else if (twoCount > oneCount) {  # too many of PAIR2
          rows <- sample(which(allTest$eventType==PAIR2));
          rows <- rows[1:(twoCount-oneCount)];
          useTest <- allTest[-rows,];
        }    
      }
      if (summary(useTest$eventType)[[1]] != summary(useTest$eventType)[[2]]) { stop("testing data not balanced"); }  
      if (summary(useTest$eventType)[[1]] == 0) { stop("no testing data after balancing"); }
      
      # do the classification
      rtbl[ctr,1] <- splt;
      rtbl[ctr,2] <- OFFSET;
      rtbl[ctr,3] <- SUBS[sid];
      rtbl[ctr,4] <- doSVM(useTrain, useTest, DO_DEFAULT_SCALING);  
      rtbl[ctr,5] <- doSVM(useTest, useTrain, DO_DEFAULT_SCALING); 
      ctr <- ctr + 1; 
    }
  }
  
  # save the results file
  rtbl[,6] <- apply(rtbl[,4:5], 1, mean);  # calculate the average for each row
  rownames(rtbl) <- 1:nrow(rtbl);
  write.table(rtbl, paste(outpath, "halfSplit_", ROI, "_", PAIR1, "_", PAIR2, "_", OFFSET, sc.lbl, ".txt", sep=""));    
}
}

#############################################################################################################################################################
# results tables. all offsets and people are in the files; one file per pair and ROI.

rm(list=ls());

inpath <- "d:/temp/Alan/halfSplit_out/";
SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");
PAIRS <- c("upempty_upgreen","upgreen_upred","upempty_upred")
OFFSETS <- c("precue", "postcue");
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native", "Parietal_mask_native");
source("d:/svnFiles/other/R_code/formatNumberOutput.R");  # for doformat

suffix <- "DefaultSc";  # suffix <- "RowScOnly"

all.out <- data.frame(array(NA, c(length(OFFSETS)*length(ROIS), length(PAIRS)+2)))
all.means <- data.frame(array(NA, c(length(OFFSETS)*length(ROIS), length(PAIRS)+2)))
colnames(all.out) <- c("ROI", "timebin", PAIRS)
colnames(all.means) <- colnames(all.out);

for (pr in 1:length(PAIRS)) {
  ctr <- 1;
  for (r in 1:length(ROIS)) {    
    for (OFFSET in OFFSETS) {   # OFFSET <- OFFSETS[1]; pr <- 1; r <- 1;
      # average over the splits, so just one mean per person and offset.
      fname <- paste(inpath, "halfSplit_", ROIS[r], "_", PAIRS[pr], "_", OFFSET, "_", suffix, ".txt", sep="")
      if (file.exists(fname)) {
        tbl <- read.table(fname);
        means <- rep(NA, length(SUBS))          
        for (sn in 1:length(SUBS)) {    # sn <- 1;
          stbl <- subset(tbl, subID==SUBS[sn] & timePoint==OFFSET);
          if (dim(stbl)[1] > 10) { stop("ran more than 10 splits?"); }
          means[sn] <- mean(stbl$avgProp);
        }   
        
        # return the mean and a p-value from a t-test for mean > 0.5; an estimate of significance.
        ttest <- t.test(means, alternative="greater", mu=0.5);
        all.out[ctr,1] <- ROIS[r];    all.means[ctr,1] <- ROIS[r];
        all.out[ctr,2] <- OFFSET;     all.means[ctr,2] <- OFFSET;
        all.out[ctr,pr+2] <- paste(doformat(mean(means),2), " (", doformat(ttest$p.value, 3), ")", sep="");
        all.means[ctr,pr+2] <- mean(means);  # average over the people so just one mean per offset.
        ctr <- ctr + 1;
      }
    }
  }
}
all.out

#############################################################################################################################################################

















































#
