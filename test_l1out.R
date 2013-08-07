make an error
##################################################################################################################################################
# get the number of examples going into each average. how balanced?

rm(list=ls());

where.run <- "Jo";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/avg_allSubs/"; 
}
if (where.run == "Jo") {
  inpath <- "d:/temp/Alan/for_halfSplit/";
  outpath <- "d:/temp/Alan/"; 
}
ROIS <- c("PFC_mask_native", "BG_LR_CaNaPu_native", "Parietal_mask_native");
vox.counts <- c(2778, 1098, 2150);   # number of voxels in each ROI, in the same order as ROIS
SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="");   # 1003 and 1004 have some categories with no examples (e.g. no errors)
#OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
#OFFSETS <- c("precue", "postcue");
OFFSETS <- -5;

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

r <- 2;  # smallest; doesn't matter which we do.

out.tbl <- data.frame(array(NA, c(length(SUBS)*length(do.names), 3)));   # set up output table (to be written out)
colnames(out.tbl) <- c("subID", "doing", "row.count");
rowctr <- 1;  # counter for rows in out.tbl
for (p in 1:length(do.names)) {    #  p <- 1;  r <- 1; 
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
    
    # find the start of each trial that's of the type we want (do.names) 
    inds <- which(tbl$eventType == do.names[p]);
    if (length(inds) < 1) { stop("too few trials of this type"); }
    start.inds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
    if (length(start.inds) < 1) { stop("too few starting trials of this type"); }
    
    for (OFFSET in OFFSETS) {    # OFFSET <- OFFSETS[1]
      if (OFFSET == "precue" | OFFSET == "postcue") {
        if (OFFSET == "precue") { do.offsets <- c(-4, -3, -2, -1, 0); }  # old-style offsets corresponding to the two time windows: which TR to average.
        if (OFFSET == "postcue") { do.offsets <- c(2, 3, 4, 5, 6); } 
      } else { do.offsets <- OFFSET; }
      
      # useinds holds the row numbers of *all* the trials, all the offsets.
      if (length(do.offsets) > 1) {
        useinds <- inds[start.inds] + do.offsets[1];  # data to classify; offset from time point 0: start of each trial
        for (i in 2:length(do.offsets)) { useinds <- union(useinds, inds[start.inds] + do.offsets[i]); }  # add on all the other timepoints in the timebin
      } else { useinds <- inds[start.inds] + do.offsets; } 
                   
      # save the average and labels in the output table
      out.tbl[rowctr,1] <- SUB;
      out.tbl[rowctr,2] <- do.names[p];
      out.tbl[rowctr,3] <- length(useinds);
      
      rowctr <- rowctr + 1;  # move counter
    }
  }
}
write.table(out.tbl, paste(outpath, "num_examples2.txt", sep="")); 

tbl <- read.table("d:/temp/Alan/num_examples.txt")
tbl

# a LOT of difference in number of correct and incorrect trials. Thankfully, upempty has the most imbalance, but this class wasn't classified,
# making me think that's not the driving factor. Also, all time points have the same amount of imbalance but don't classify the same. (if imbalance
# was the only reason stuff classified, it should be more uniform - similar across timepoints. perhaps just do a few (non-permutation) comparisons
# with balancing the number of examples going into the averages; if similar to existing analysis, no problem.

##################################################################################################################################################
# prepare the input datafiles for group classification: average the examples together so one per person, ROI, and type.
# BALANCED VERSION
# reads in the file made above to know how many examples are needed

rm(list=ls());

inpath <- "d:/temp/Alan/for_halfSplit/";
outpath <- "d:/temp/Alan/balanced/"; 

ROIS <- c("PFC_mask_native", "BG_LR_CaNaPu_native", "Parietal_mask_native");
vox.counts <- c(2778, 1098, 2150);   # number of voxels in each ROI, in the same order as ROIS
SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="");   # 1003 and 1004 have some categories with no examples (e.g. no errors)
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
#OFFSETS <- c("precue", "postcue");


# what we want the trial types to be called in the output file names.
doing <- "cor_incor";  # doing <- "cor_incor"; #
#if (doing == "cor_incor") { do.names <- c("upgreenCOR","upemptyCOR","upredCOR", "upgreenINCOR","upemptyINCOR","upredINCOR");}
if (doing == "cor_incor") { do.names <- c("upgreenCOR","upgreenINCOR");}
if (doing == "first_version") { do.names <- c("upgreen","upempty","upred"); }

ct.tbl <- read.table("d:/temp/Alan/num_examples.txt")
# head(ct.tbl)
#   subID      doing row.count
# sub1005 upgreenCOR        44
# sub1006 upgreenCOR        46
# sub1007 upgreenCOR        46


for (r in 1:length(ROIS)) {   #   r <- 2; 
for (n in 2:5) {
  for (p in 1:length(do.names)) {    #  p <- 4;  r <- 2; 
    out.tbl <- data.frame(array(NA, c(length(SUBS)*length(OFFSETS), vox.counts[r]+2)));   # set up output table (to be written out)
    colnames(out.tbl) <- c("subID", "offset", paste("v", 1:vox.counts[r], sep=""));
    rowctr <- 1;  # counter for rows in out.tbl
    for (SUB in SUBS) {     # SUB <- SUBS[1];
      tbl <- read.table(gzfile(paste(inpath, SUB, "_", ROIS[r], "_meanSub.gz", sep="")), comment.char=""); # read in the data
      if (length(unique(tbl$TR)) != 210) { stop("too many or few TRs"); }
      if (ncol(tbl) != vox.counts[r]+4) { stop("not the expected number of voxels in the input file."); }
      
      # figure out how many examples we need for this class; minimum of incorrect and correct trial counts
      if (substr(do.names[p], 1,3) == "upg") { ct.inds <- which(ct.tbl$subID == SUB & (ct.tbl$doing == "upgreenCOR" | ct.tbl$doing == "upgreenINCOR")); }
      if (substr(do.names[p], 1,3) == "upr") { ct.inds <- which(ct.tbl$subID == SUB & (ct.tbl$doing == "upredCOR" | ct.tbl$doing == "upredINCOR")); }
      if (substr(do.names[p], 1,3) == "upe") { ct.inds <- which(ct.tbl$subID == SUB & (ct.tbl$doing == "upemptyCOR" | ct.tbl$doing == "upemptyINCOR")); }
      need.ct <- min(ct.tbl$row.count[ct.inds]);
      
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
     
      # find the start of each trial that's of the type we want (do.names) 
      inds <- which(tbl$eventType == do.names[p]);
      if (length(inds) < 1) { stop("too few trials of this type"); }
      start.inds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
      if (length(start.inds) < 1) { stop("too few starting trials of this type"); }
      
      for (OFFSET in OFFSETS) {    # OFFSET <- OFFSETS[1]
        if (OFFSET == "precue" | OFFSET == "postcue") {
          if (OFFSET == "precue") { do.offsets <- c(-4, -3, -2, -1, 0); }  # old-style offsets corresponding to the two time windows: which TR to average.
          if (OFFSET == "postcue") { do.offsets <- c(2, 3, 4, 5, 6); } 
        } else { do.offsets <- OFFSET; }
        
        if (length(do.offsets) > 1) {
          useinds <- inds[start.inds] + do.offsets[1];  # data to classify; offset from time point 0: start of each trial
          for (i in 2:length(do.offsets)) { useinds <- union(useinds, inds[start.inds] + do.offsets[i]); }  # add on all the other timepoints in the timebin
        } else { useinds <- inds[start.inds] + do.offsets; } 
        
        if (length(useinds) > need.ct) { useinds <- sample(useinds)[1:need.ct]; }  # only keep a random subset of the examples if too many
        if (length(useinds) < need.ct) { stop("length(useinds) < need.ct"); }
        temp.tbl <- apply(tbl[useinds,5:ncol(tbl)], 2, mean); # voxel columns of the rows we need out of tbl; calculate the over-trials average for each voxel
        if (length(temp.tbl) != vox.counts[r]) { stop("wrong number of voxels in temp.tbl after averaging"); }
        if (length(which(is.na(temp.tbl))) > 0) { stop("NAs in temp.tbl after averaging"); }
        
        # save the average and labels in the output table
        out.tbl[rowctr,1] <- SUB;
        out.tbl[rowctr,2] <- OFFSET;
        out.tbl[rowctr,3:ncol(out.tbl)] <- temp.tbl;
        
        rowctr <- rowctr + 1;  # move counter
      }
    }
    write.table(out.tbl, gzfile(paste(outpath, ROIS[r], "_", do.names[p], "_avg_allSubs_timebins_bal", n, ".txt.gz", sep=""))); 
  }
}

##################################################################################################################################################
##################################################################################################################################################
# run the classification on balanced data, no permutations.

library(e1071);  # R interface to libsvm (only needs calling once per R session)

rm(list=ls());

where.run <- "Jo";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") {
  #inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/avg_allSubs/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/l1subOut_result/"; 
}
if (where.run == "Jo") {
  inpath <- "d:/temp/Alan/balanced/";
  outpath <- "d:/temp/Alan/balanced_out/"; 
}
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native","Parietal_mask_native");
#ROIS <- c("Parietal_mask_native");
SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="")
#SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");   # SUBS <- c(1003:1009, 1011:1019);   # SUBS <- "sub1003";
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify

# the things to classify. will classify the first entry of PAIR1S with the first entry of PAIR2S, etc.
#PAIR1S <- c("upgreenCOR","upemptyCOR","upredCOR");  
#PAIR2S <- c("upgreenINCOR","upemptyINCOR","upredINCOR");  
PAIR1S <- c("upgreenCOR")
PAIR2S <- c("upgreenINCOR")
#PAIR1S <- c("upgreen","upred","upred");  
#PAIR2S <- c("upempty","upempty","upgreen");  


# flags for type of scaling to do.
DO_DEFAULT_SCALING <- TRUE;

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
if (DO_DEFAULT_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "DefaultSc", sep=""); }

if (length(PAIR1S) != length(PAIR2S)) { stop("length(PAIR1S) != length(PAIR2S)"); }
#for (ROI in ROIS) {     # 
ROI <- ROIS[1]; 
p <- 1;
  #for (p in 1:length(PAIR1S)) {   # p <- 1; ROI <- ROIS[1]; 
for (n in 1:5) {   # n <- 1;
    # make a dataframe with all the data we need to classify this pair of stimuli
    if (PAIR1S[p] == PAIR2S[p]) { stop("PAIR1 == PAIR2"); }
    tbl <- read.table(gzfile(paste(inpath, ROI, "_", PAIR1S[p], "_avg_allSubs_timebins_bal", n, ".txt.gz", sep="")));
    tbl <- data.frame(PAIR1S[p], tbl);   # add the trial type label - not in the average file, just its filename
    colnames(tbl)[1] <- "expType";   # fix the new column's title
    
    temp.tbl <- read.table(gzfile(paste(inpath, ROI, "_", PAIR2S[p], "_avg_allSubs_timebins_bal", n, ".txt.gz", sep="")));
    temp.tbl <- data.frame(PAIR2S[p], temp.tbl);   # add the trial type label - not in the average file, just its filename
    colnames(temp.tbl)[1] <- "expType";   # fix the new column's title
    if (nrow(tbl) != nrow(temp.tbl) | ncol(tbl) != ncol(temp.tbl)) { stop("tbl and temp.tbl different sizes"); }
    
    tbl <- rbind(tbl, temp.tbl);  # put both into one dataframe
    rm(temp.tbl);  # take the temporary table out of memory
    if (length(summary(tbl$expType)) != 2) { stop("not two expType in tbl"); }
    if (summary(tbl$expType)[[1]] != summary(tbl$expType)[[2]]) { stop("number of examples not balanced"); }
    if (length(unique(tbl$subID)) != length(SUBS)) { stop("number of subjects don't match in tbl$subID and SUBS"); }
    
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
    write.table(rtbl, paste(outpath, "l1subOut_", ROI, "_", PAIR1S[p], "_", PAIR2S[p], SCALING_LABEL, "bal", n, ".txt", sep=""));    
  }
}

##################################################################################################################################################
# plot the results of the subsetted averages for comparison with the no-subsetting averages

rm(list=ls());

path <- "d:/temp/Alan/balanced_out/"; 

ROI <- "BG_LR_CaNaPu_native"
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
PAIR1 <- "upgreenCOR"
PAIR2 <- "upgreenINCOR"

tbl <- read.table("d:/temp/Alan/forplotting/upgreenCOR_upgreenINCOR_BG_LR_CaNaPu_native_nullDists.txt"); 
reals <- tbl[1,];  # real accuracy
accs <- array(NA, c(11,5))

windows(10, 5);  # specify the plot window size so everything will look like it should.
layout(matrix(1:6, c(2,3), byrow=TRUE));
par(mar=c(3, 3.75, 2, 1), mgp=c(3, 1, 0)) 
# mar: c(bottom, left, top, right)’ gives the number of lines of margin to be specified on the four sides of the plot. default is ‘c(5, 4, 4, 2) + 0.1’.
# mgp: margin line (in ‘mex’ units) for the axis title, axis labels and axis line. The default is ‘c(3, 1, 0)’.

for (n in 1:5) {   # n <- 1;
  tbl <- read.table(paste(path, "l1subOut_", ROI, "_", PAIR1, "_", PAIR2, "_DefaultScbal", n, ".txt", sep=""))
  accs[,n] <- tbl$avg.acc;
  ttl <- paste("subsetting #", n, sep="");
  plot(x=seq(from=2,to=22,by=2), y=reals, ylim=c(0.2,0.8), type='l', lwd=2, col='darkorange', main=ttl, xlab="", ylab="", 
    cex.axis=1.5, cex.lab=1.5, cex.main=1.5);   # no-subsetting accuracy
  lines(x=c(-1,25), y=rep(0.5,2), col='darkgrey');   # horizontal line at chance
    
  lines(x=seq(from=2,to=22,by=2), y=tbl$avg.acc, lwd=2, col='black')  # accuracy line for this subset
}

ttl <- "average of the 5 replications";
mn <- apply(accs, 1, mean);
plot(x=seq(from=2,to=22,by=2), y=reals, ylim=c(0.2,0.8), type='l', lwd=2, col='darkorange', main=ttl, xlab="", ylab="",
   cex.axis=1.5, cex.lab=1.5, cex.main=1.5);   # no-subsetting accuracy
lines(x=c(-1,25), y=rep(0.5,2), col='darkgrey');   # horizontal line at chance 
lines(x=seq(from=2,to=22,by=2), y=mn, lwd=2, col='cadetblue')


##################################################################################################################################################
















































#
