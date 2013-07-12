make an error
# alan's half-splitting within-subjects analyses; not permutation tests.


library(e1071);  # R interface to libsvm

rm(list=ls()); ONNIL <- FALSE; ONCLUSTER <- TRUE; JOCOMPUTER <- FALSE;
# rm(list=ls()); ONNIL <- FALSE; ONCLUSTER <- FALSE; JOCOMPUTER <- TRUE;

SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");
#ROI <- "PFC_mask_native"     
ROI <- "BG_LR_CaNaPu_native"
#ROI <- "Parietal_mask_native"

OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5); # timepoints to classify, in TR from the first PAIR1 or PAIR2 in each event.
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

if (ONCLUSTER == TRUE) {
  inpath <- "~/tmp/input/";
  outpath <- "~/tmp/output/"; 
  
  cA <- commandArgs();   
  num  <- as.numeric(cA[5]);  # specifies the pair to run.
  if (num == 1) { PAIR1 <- "upempty"; PAIR2 <- "upgreen"; } 
  if (num == 2) { PAIR1 <- "upempty"; PAIR2 <- "upred"; } 
  if (num == 3) { PAIR1 <- "upgreen"; PAIR2 <- "upred";  }
}
if (ONNIL == TRUE) {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
  permpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/permInputFiles/";   # location of 6eachTable.txt, 8eachTable.txt, etc.
  SUB <- "sub1003";
  PAIR1 <- "upempty"; PAIR2 <- "upgreen";
}
if (JOCOMPUTER == TRUE) {
  inpath <- "d:/temp/Alan/for_halfSplit/";
  outpath <- "d:/temp/Alan/halfSplit_out/";
  PAIR1 <- "upempty"; PAIR2 <- "upgreen";
}
  
SEEDS <- c(51591, 36414, 56347, 38442, 20176, 51348, 89727, 67106, 23543, 52663);  # 10 random seeds, from sample(1:100000)[1:10]
NUMSPLITS <- 10;    # there will be 10, but we're running them one at a time.
if (PAIR1 == PAIR2) { stop("PAIR1 == PAIR2"); }

DO_RUN_COLUMN_MS <- FALSE;   # these flags indicate which sort of scaling to do.
DO_RUN_COLUMN_SC <- TRUE; 
DO_ROW_SCALING <- FALSE;
DO_DEFAULT_SCALING <- FALSE;
sc.lbl <- "_";  # make a label for the output files showing the type of scaling used for this classification
if (DO_RUN_COLUMN_MS == TRUE) { sc.lbl <- paste(sc.lbl, "RunColMS", sep=""); }
if (DO_RUN_COLUMN_SC == TRUE) { sc.lbl <- paste(sc.lbl, "RunColSC", sep=""); }
if (DO_ROW_SCALING == TRUE) { sc.lbl <- paste(sc.lbl, "RowSc", sep=""); }
if (DO_DEFAULT_SCALING == TRUE) { sc.lbl <- paste(sc.lbl, "DefaultSc", sep=""); 
} else { sc.lbl <- paste(sc.lbl, "Only", sep=""); }


for (OFFSET in rev(OFFSETS)) {   # OFFSET <- OFFSETS[11]
  rtbl <- data.frame(array(NA, c(NUMSPLITS*length(SUBS), 6)));   # results table
  colnames(rtbl) <- c("splitNum", "timePoint", "subID", "firstTest", "secondTest", "avgProp");
  ctr <- 1;
  for (sid in 1:length(SUBS)) {  # sid <- 1;
    tbl <- read.table(gzfile(paste(inpath, SUBS[sid], "_", ROI, ".txt.gz", sep="")), comment.char=""); # read in the data
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
    
    useinds <- inds[tinds] + OFFSET;  # data to classify; offset from time point 0: start of each trial
    if (length(which(useinds < 0)) > 0) { stop("yes, have negative rows"); }
    t0tbl <- tbl[useinds,];
    if (OFFSET != 0) {  # need to get labels for which events these go with; offset labels are other things than what we're classifying.
      uselbls <- tbl$eventType[inds[tinds]];     
      if (length(uselbls) != dim(t0tbl)[1]) { stop("very wrong lengths for uselbls"); } else { t0tbl$eventType <- uselbls; }
    }
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

#############################################################################################################################################################
# results tables. all offsets and people are in the files; one file per pair and ROI.

rm(list=ls());

inpath <- "d:/temp/Alan/halfSplit_out/";
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");
SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");
PAIRS <- c("upempty_upgreen","upgreen_upred","upempty_upred")
OFFSETS <- -5:5

ROI <- ROIS[1];

#> head(tbl)
#  splitNum timePoint   subID firstTest secondTest   avgProp
#1        1        -5 sub1003 0.5000000  0.5384615 0.5192308
#2        2        -5 sub1003 0.4444444  0.6923077 0.5683761
#3        3        -5 sub1003 0.5277778  0.6153846 0.5715812

all.out <- array(NA, c(length(OFFSETS), length(PAIRS)));
all.means <- array(NA, c(length(OFFSETS), length(PAIRS)));
colnames(all.out) <- PAIRS;
colnames(all.means) <- PAIRS;

for (pr in 1:length(PAIRS)) {   # pr <- 1;
  # average over the splits, so just one mean per person and offset.
  mtbl <- data.frame(array(NA, c(length(SUBS)*length(OFFSETS), 3)))
  colnames(mtbl) <- c("subID", "offset", "meanProp");
  ctr <- 1;  # row counter for output table
  for (OFFSET in OFFSETS) {   # SUB <- SUBS[1]; OFFSET <- OFFSETS[1];
    fname <- paste(inpath, "halfSplit_", ROI, "_", PAIRS[pr], "_", OFFSET, "_RowSCOnly.txt", sep="")
    if (file.exists(fname)) {
      tbl <- read.table(fname);
      for (SUB in SUBS) {
        stbl <- subset(tbl, subID==SUB & timePoint==OFFSET);
        if (dim(stbl)[1] > 10) { stop("ran more than 10 splits?"); }
        mtbl[ctr,1] <- SUB;
        mtbl[ctr,2] <- OFFSET;
        mtbl[ctr,3] <- mean(stbl$avgProp);
        ctr <- ctr + 1;
      }
    }
  }
  
  # now average over the people so just one mean per offset.
  # return the mean and a p-value from a t-test for mean > 0.5. this is *just* an estimate of significance.
  out <- rep(NA, length(OFFSETS));
  means <- rep(NA, length(OFFSETS));
  for (i in 1:length(OFFSETS)) {   # i <- 1;
    stbl <- subset(mtbl, offset==OFFSETS[i]);
    if (nrow(stbl) > 0) {
      if (nrow(stbl) > length(SUBS)) { stop("too many rows"); }
      ttest <- t.test(stbl$meanProp, alternative="greater", mu=0.5);
      out[i] <- paste(round(mean(stbl$meanProp),2), " (", round(ttest$p.value, 3), ")", sep="");
      means[i] <- mean(stbl$meanProp)
    }
  }
  all.out[,pr] <- out;
  all.means[,pr] <- means;
}
cbind(OFFSETS, all.out)

boxplot(meanProp~offset, data=mtbl)

plot(x=OFFSETS, y=means, type='l', ylab="classification accuracy", xlab=ROI, main="")

seq(from=2,to=22, by=2)
#############################################################################################################################################################
# plots showing the means

library(ggplot2);

rm(list=ls());


inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/";
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");
SUBS <- paste("sub", c(1003:1009, 1011:1018), sep="");
PAIRS <- c("_upempty_upred","_upempty_upgreen","_upgreen_upred")

ROI <- ROIS[1];
PAIR <- PAIRS[1}

tbl <- read.table(paste(inpath, ROI, PAIR, sep=""));

OFFSETS <- unique(tbl$timePoint);
# average over the splits, so just one mean per person and offset.
means <- rep(NA, length(SUBS)*length(OFFSETS));
lbls <- array(NA, c(length(SUBS)*length(OFFSETS), 2));  # label table
ctr <- 1;  # row counter for output table
for (SUB in SUBS) {
    for (OFFSET in OFFSETS) {   # SUB <- SUBS[1]; OFFSET <- OFFSETS[1];
        stbl <- subset(tbl, subID==SUB & timePoint==OFFSET);
        if (dim(stbl)[1] > 10) { stop("ran more than 10 splits?"); }
        means[ctr] <- mean(stbl$avgProp);
        lbls[ctr,1] <- SUB;
        lbls[ctr,2] <- OFFSET;
        ctr <- ctr + 1;
    }
}
mtbl <- data.frame(lbls, means);
colnames(mtbl) <- c("subID", "offset", "meanProp");
mtbl$offset <- factor(mtbl$offset);

# http://stackoverflow.com/questions/8269016/ggplot2-boxplot-horizontal-bar-at-median
f <- function(x, height) {
 ans <- median(x)
 data.frame(ymin = ans-height/2, ymax = ans+height/2, y = ans)
}

ttl <- "put the plot title here";

ggplot(mtbl, aes(x=offset, y=meanProp)) + geom_hline(aes(yintercept=0.5), colour="darkgrey", size=1) + geom_boxplot() + scale_y_continuous(limits=c(0, 1)) +
 stat_summary(fun.data=f, geom="crossbar", height=0.03, colour=NA, fill="skyblue", width=0.8, alpha=0.5) +
 geom_point(aes(colour=subID), position=position_jitter(w=0.1)) + opts(title=ttl, legend.position="none")
# try rug? could if do panels on offset

###############################################################################################################################################################





#############################################################################################################################################################

















































#
