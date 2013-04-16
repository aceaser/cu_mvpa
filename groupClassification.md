cu_mvpa
=======

For scripts to run mvpa for control data.

make an error
##################################################################################################################################################
# Alan Ceasar classification: each time point separately, all people together.
##################################################################################################################################################
# prepare the input datafiles for group classification: average the examples together so one per person, ROI, and type.
# makes one file for each ROI and stimulus type (upgreenCOR, upemptyINCOR, etc); all people and timepoints in each file.
# the code is slow since it reads in each person's datafiles several times, but only has to be run once.
# change the labeling code here (accurate-only or whatever) to make new files as needed.

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
}
if (where.run == "Jo") {
  inpath <- "c:/maile/svnFiles/plein/consulting/Alan/dataFilesFromServer/";
  outpath <- "d:/temp/"; 
}
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");
vox.counts <- c(1098, 2778);   # number of voxels in each ROI, in the same order as ROIS
SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="");   # SUBS <- c(1003:1009, 1011:1019);   # SUBS <- "sub1003";
  #removed 1003, 1004, 1019)
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify

# tbl[1:5,1:10]
# eventType   subID run TR        v1        v2       v3        v4         v5         v6
#       ITI sub1003   1  1 -20.89934 -0.328182 18.35852 14.280652 -0.4373814 -3.9457958
#    memset sub1003   1  2 -22.48345  1.808049 15.40503 15.047131  0.2154506 -0.7131909
#    delay1 sub1003   1  3 -23.16070  5.945988 12.77588  7.411145  0.5202602 -6.1879833
#    delay1 sub1003   1  4 -23.66143  4.240422 21.53211 14.946362 10.2307094  5.1469776
#    delay1 sub1003   1  5 -24.57488 -4.042537 18.86731  8.944958  5.3793910  3.9625294


# what we want the trial types to be called in the output file names.
doing <- "cor_incor";  # doing <- "cor_incor"; #
if (doing == "cor_incor") { do.names <- c("upgreenCOR","upemptyCOR","upredCOR", "upgreenINCOR","upemptyINCOR","upredINCOR");}
if (doing == "first_version") { do.names <- c("upgreen","upempty","upred"); }

for (r in 1:length(ROIS)) {   #   r <- 1; 
  for (p in 1:length(do.names)) {    #  p <- 4;  r <- 1; 
    out.tbl <- data.frame(array(NA, c(length(SUBS)*length(OFFSETS), vox.counts[r]+2)));   # set up output table (to be written out)
    colnames(out.tbl) <- c("subID", "offset", paste("v", 1:vox.counts[r], sep=""));
    rowctr <- 1;  # counter for rows in out.tbl
    for (SUB in SUBS) {     # SUB <- SUBS[1];
      tbl <- read.table(gzfile(paste(inpath, SUB, "_", ROIS[r], "_meanSub.gz", sep="")), comment.char=""); # read in the data
      if (length(unique(tbl$TR)) != 210) { stop("too many or few TRs"); }
      if (ncol(tbl) != vox.counts[r]+4) { stop("not the expected number of voxels in the input file."); }
      
      # change the eventTypes in the input file to match the proper labels; these need to match do.names
      if (doing == "cor_incor") {
        ind <- which(levels(tbl$eventType) == "smain"); levels(tbl$eventType)[ind] <- "upemptyCOR";
        ind <- which(levels(tbl$eventType) == "smainnp"); levels(tbl$eventType)[ind] <- "upemptyCOR";
        ind <- which(levels(tbl$eventType) == "delayempty"); levels(tbl$eventType)[ind] <- "upemptyCOR";
        ind <- which(levels(tbl$eventType) == "rmain"); levels(tbl$eventType)[ind] <- "upredCOR";
        ind <- which(levels(tbl$eventType) == "rmainup"); levels(tbl$eventType)[ind] <- "upredCOR";
        ind <- which(levels(tbl$eventType) == "rmainnp"); levels(tbl$eventType)[ind] <- "upredCOR";
        ind <- which(levels(tbl$eventType) == "delayred"); levels(tbl$eventType)[ind] <- "upredCOR";
        ind <- which(levels(tbl$eventType) == "update"); levels(tbl$eventType)[ind] <- "upgreenCOR";
        ind <- which(levels(tbl$eventType) == "updateop"); levels(tbl$eventType)[ind] <- "upgreenCOR";
        ind <- which(levels(tbl$eventType) == "updatenp"); levels(tbl$eventType)[ind] <- "upgreenCOR";
        ind <- which(levels(tbl$eventType) == "delaygreen"); levels(tbl$eventType)[ind] <- "upgreenCOR";
        ind <- which(levels(tbl$eventType) == "smainERR"); levels(tbl$eventType)[ind] <- "upemptyINCOR";
        ind <- which(levels(tbl$eventType) == "smainnpERR"); levels(tbl$eventType)[ind] <- "upemptyINCOR";
        ind <- which(levels(tbl$eventType) == "delayemptyERR"); levels(tbl$eventType)[ind] <- "upemptyINCOR";
        ind <- which(levels(tbl$eventType) == "rmainERR"); levels(tbl$eventType)[ind] <- "upredINCOR";
        ind <- which(levels(tbl$eventType) == "rmainupERR"); levels(tbl$eventType)[ind] <- "upredINCOR";
        ind <- which(levels(tbl$eventType) == "rmainnpERR"); levels(tbl$eventType)[ind] <- "upredINCOR";
        ind <- which(levels(tbl$eventType) == "delayredERR"); levels(tbl$eventType)[ind] <- "upredINCOR";
        ind <- which(levels(tbl$eventType) == "updateERR"); levels(tbl$eventType)[ind] <- "upgreenINCOR";
        ind <- which(levels(tbl$eventType) == "updateopERR"); levels(tbl$eventType)[ind] <- "upgreenINCOR";
        ind <- which(levels(tbl$eventType) == "updatenpERR"); levels(tbl$eventType)[ind] <- "upgreenINCOR";
        ind <- which(levels(tbl$eventType) == "delaygreenERR"); levels(tbl$eventType)[ind] <- "upgreenINCOR";
      }
      
      if (doing == "first_version") {
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
      }    
      
      # find the start of each trial that's of the type we want (do.names) 
      inds <- which(tbl$eventType == do.names[p]);
      if (length(inds) < 1) { stop("too few trials of this type"); }
      start.inds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
      if (length(start.inds) < 1) { stop("too few starting trials of this type"); }
      
      for (OFFSET in OFFSETS) {    # OFFSET <- OFFSETS[1]
        useinds <- inds[start.inds] + OFFSET;  # data to classify; offset from time point 0: start of each trial
        if (length(which(useinds < 0)) > 0) { stop("have negative rows"); }
        temp.tbl <- tbl[useinds,5:ncol(tbl)]; # voxel columns of the rows we need out of tbl
        temp.tbl <- apply(temp.tbl, 2, mean);  # calculate the over-trials average for each voxel
        if (length(temp.tbl) != vox.counts[r]) { stop("wrong number of voxels in temp.tbl after averaging"); }
        if (length(which(is.na(temp.tbl))) > 0) { stop("NAs in temp.tbl after averaging"); }
        
        # save the average and labels in the output table
        out.tbl[rowctr,1] <- SUB;
        out.tbl[rowctr,2] <- OFFSET;
        out.tbl[rowctr,3:ncol(out.tbl)] <- temp.tbl;
        
        rowctr <- rowctr + 1;  # move counter
      }
    }
    write.table(out.tbl, gzfile(paste(outpath, ROIS[r], "_", do.names[p], "_avg_allSubs.txt.gz", sep=""))); 
  }
}

##################################################################################################################################################
##################################################################################################################################################
# leave-one-person-out cross-validation

library(e1071);  # R interface to libsvm (only needs calling once per R session)

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") {
  #inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
}
if (where.run == "Jo") {
  inpath <- "d:/temp/";
  outpath <- "d:/temp/"; 
}
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");
#SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="")
SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");   # SUBS <- c(1003:1009, 1011:1019);   # SUBS <- "sub1003";
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify

# the things to classify. will classify the first entry of PAIR1S with the first entry of PAIR2S, etc.
#PAIR1S <- c("upgreenCOR","upemptyCOR","upredCOR");  
#PAIR2S <- c("upgreenINCOR","upemptyINCOR","upredINCOR");  
PAIR1S <- c("upgreen","upred","upred");  
PAIR2S <- c("upempty","upempty","upgreen");  


# flags for type of scaling to do.
DO_ROW_SCALING <- FALSE;
DO_DEFAULT_SCALING <- FALSE;

doSVM <- function(train, test) {  # test <- test.set; train <- train.set;
  test <- subset(test, select=c(-subID, -offset));  # get rid of non-classify or voxel columns
  train <- subset(train, select=c(-subID, -offset));
  if (colnames(test)[2] != "v1" | colnames(train)[2] != "v1") { stop("v1 not found where expected"); }
  
  fit <- svm(expType~., data=train, type="C-classification", kernel="linear", cost=1, scale=DO_DEFAULT_SCALING);  
  tree <- table(test$expType, predict(fit, test));
  if (dim(tree)[2]==1 | dim(tree)[1]==1) { wrT <- 0.5; } else { wrT <- sum(diag(tree))/sum(tree); }
  
  return(wrT);
}

SCALING_LABEL <- "_";  # make a label for the output files showing the type of scaling used for this classification
if (DO_ROW_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "RowSc", sep=""); }
if (DO_DEFAULT_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "DefaultSc", sep=""); }

if (length(PAIR1S) != length(PAIR2S)) { stop("length(PAIR1S) != length(PAIR2S)"); }
for (ROI in ROIS) {     # ROI <- ROIS[1]; 
  for (p in 1:length(PAIR1S)) {   # p <- 1; ROI <- ROIS[1]; 
    # make a dataframe with all the data we need to classify this pair of stimuli
    if (PAIR1S[p] == PAIR2S[p]) { stop("PAIR1 == PAIR2"); }
    tbl <- read.table(gzfile(paste(inpath, ROI, "_", PAIR1S[p], "_avg_allSubs.txt.gz", sep="")));
    tbl <- data.frame(PAIR1S[p], tbl);   # add the trial type label - not in the average file, just its filename
    colnames(tbl)[1] <- "expType";   # fix the new column's title
    
    temp.tbl <- read.table(gzfile(paste(inpath, ROI, "_", PAIR2S[p], "_avg_allSubs.txt.gz", sep="")));
    temp.tbl <- data.frame(PAIR2S[p], temp.tbl);   # add the trial type label - not in the average file, just its filename
    colnames(temp.tbl)[1] <- "expType";   # fix the new column's title
    if (nrow(tbl) != nrow(temp.tbl) | ncol(tbl) != ncol(temp.tbl)) { stop("tbl and temp.tbl different sizes"); }
    
    tbl <- rbind(tbl, temp.tbl);  # put both into one dataframe
    rm(temp.tbl);  # take the temporary table out of memory
    if (length(summary(tbl$expType)) != 2) { stop("not two expType in tbl"); }
    if (summary(tbl$expType)[[1]] != summary(tbl$expType)[[2]]) { stop("number of examples not balanced"); }
    if (length(unique(tbl$subID)) != length(SUBS)) { stop("number of subjects don't match in tbl$subID and SUBS"); }

    if (DO_ROW_SCALING == TRUE) {  # do row scaling (all voxels in each volume)
      scaled <- tbl[,4:nrow(tbl)];   # get voxel columns; cols 1 to 4 are labels
      scaled <- scale(t(scaled));  # scale does columns, so transpose
      tbl <- data.frame(tbl[,1:3], t(scaled));  # put label columns back on, transpose back
      rm(scaled);
    }   
    
    # do the leave-one-subject-out cross-validation
    rtbl <- data.frame(array(NA, c(length(OFFSETS), length(SUBS)+3)));   # results table (will be written out)
    colnames(rtbl) <- c("offset", "pair", paste(SUBS, "out", sep=""), "avg.acc");
    for (o in 1:length(OFFSETS)) {   # o <- 1; 
      rtbl[o,1] <- OFFSETS[o];
      rtbl[o,2] <- paste(PAIR1S[p], PAIR2S[p], sep="_");
      for (sn in 1:length(SUBS)) {     # sn <- 1;
        train.set <- subset(tbl, offset == OFFSETS[o] & subID != SUBS[sn]);  # everyone but the testing subject
        test.set <- subset(tbl, offset == OFFSETS[o] & subID == SUBS[sn]);  # just the testing subject
        # double-check the training and testing data a bit
        if (summary(train.set$expType)[[1]] != summary(train.set$expType)[[2]]) { stop("training data not balanced"); }  
        if (summary(test.set$expType)[[1]] != summary(test.set$expType)[[2]]) { stop("testing data not balanced"); }  
        if (nrow(train.set) != (length(SUBS)-1)*2) { stop("not the right number of rows in the training data"); }
        if (nrow(test.set) != 2) { stop("not two rows in the testing data"); }
        
        rtbl[o,sn+2] <- doSVM(train.set, test.set);   
      }
      rtbl[o,ncol(rtbl)] <- mean(as.vector(rtbl[o,3:(length(SUBS)+2)], mode="numeric"));   # average accuracy
    }
    if (SCALING_LABEL == "_") { SCALING_LABEL <- "_noScaling"; }
    rownames(rtbl) <- 1:dim(rtbl)[1];
    write.table(rtbl, paste(outpath, "l1subOut_", ROI, "_", PAIR1S[p], "_", PAIR2S[p], SCALING_LABEL, ".txt", sep=""));    
  }
}

##################################################################################################################################################
# results table-viewing code not really needed: just the avg.acc column in each file.
# but this is a bit of graphing code: accuracy at each offset

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") { path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; }
if (where.run == "Jo") { path <- "d:/temp/"; }
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
ROI <- "BG_LR_CaNaPu_native";  # ROI <- "PFC_mask_native"; # 

scaling <- "RowSc";   # scaling <- "noScaling";

plot(x=OFFSETS, y=rep(0, length(OFFSETS)), col='white', xlab="time offset", ylab="leave-one-subject-out accuracy", ylim=c(0.25,1), las=1);
lines(x=c(-10, 10), y=c(0.5, 0.5), col='grey', lwd=2);  # line at chance

# can just put the file names to plot in here - lots of lines can be plotted
tbl <- read.table(paste(path, "l1subOut_", ROI, "_upred_upgreen_", scaling, ".txt", sep=""));
lines(x=OFFSETS, y=tbl$avg.acc, col='black');

tbl <- read.table(paste(path, "l1subOut_", ROI, "_upred_upempty_", scaling, ".txt", sep=""));
lines(x=OFFSETS, y=tbl$avg.acc, col='red');

tbl <- read.table(paste(path, "l1subOut_", ROI, "_upgreen_upempty_", scaling, ".txt", sep=""));
lines(x=OFFSETS, y=tbl$avg.acc, col='forestgreen');




pair <- "upgreen_upempty"
#

plot(x=OFFSETS, y=rep(0, length(OFFSETS)), col='white', xlab="time offset", ylab="leave-one-subject-out accuracy", ylim=c(0.25,1), las=1);
lines(x=c(-10, 10), y=c(0.5, 0.5), col='grey', lwd=2);  # line at chance

# can just put the file names to plot in here - lots of lines can be plotted
tbl <- read.table(paste(path, "l1subOut_", ROI, "_", pair, "_RowSc.txt", sep=""));
lines(x=OFFSETS, y=tbl$avg.acc, col='blue');

tbl <- read.table(paste(path, "l1subOut_", ROI, "_", pair, "_RowSc.txt", sep=""));
lines(x=OFFSETS, y=tbl$avg.acc, col='red');

tbl <- read.table(paste(path, "l1subOut_", ROI, "_", pair, "_RowSc.txt", sep=""));
lines(x=OFFSETS, y=tbl$avg.acc, col='black');

##################################################################################################################################################


