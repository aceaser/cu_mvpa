make an error
##################################################################################################################################################
# Alan Ceasar classification: each time point separately, all people together.
# COMPLETE PERMUTATION version.
library(e1071);  # R interface to libsvm (only needs calling once per R session)

rm(list=ls());

where.run <- "cluster";   # where.run <- "Jo";  # set the paths according to which computer this is being run on

if (where.run == "cluster") {
   inpath <- "~/AccPerm/";
   outpath <- "~/compPerms/output/"; 
   perm.path <- "~/compPerms/";   # location of l1outPermlbls.txt
    
   cA <- commandArgs();   # get the numbers sent in the job file
   o  <- as.numeric(cA[5]);  # which timepoints to run 
   
   num  <- as.numeric(cA[6]);  # which region to run
   if (num == 1) { ROI <- "BG_LR_CaNaPu_native"; }
   if (num == 2) { ROI <- "PFC_mask_native"; }
   
   num  <- as.numeric(cA[7]);  # which pair to run
   if (num == 1) { PAIR1 <- "upgreenCOR"; PAIR2 <- "upgreenINCOR"; }
   if (num == 2) { PAIR1 <- "upemptyCOR"; PAIR2 <- "upemptyINCOR"; }
   if (num == 3) { PAIR1 <- "upredCOR"; PAIR2 <- "upredINCOR"; }
   
   num  <- as.numeric(cA[8]);   # which permutations to do
   if (num == 1) { do.perms <- 1:1000; }
   if (num == 2) { do.perms <- 1001:2000; }
   if (num == 3) { do.perms <- 2001:3000; }
   if (num == 4) { do.perms <- 3001:4000; }
   if (num == 5) { do.perms <- 4001:5000; }
   if (num == 6) { do.perms <- 5001:6000; }
   if (num == 7) { do.perms <- 6001:7000; }
   if (num == 8) { do.perms <- 7001:8000; }
   if (num == 9) { do.perms <- 8001:8191; }
   perm.key <- letters[num]
}

if (where.run == "Alan") {
  #inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
  perm.path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
  o <- 1; 
  ROI <- "BG_LR_CaNaPu_native";
  PAIR1 <- "upredCOR"; PAIR2 <- "upredINCOR";
  do.perms <- 1:5; perm.key <- "a1"; 
}
if (where.run == "Jo") {
  inpath <- "d:/temp/";
  outpath <- "d:/temp/"; 
  perm.path <- "c:/maile/svnFiles/plein/consulting/Alan/setupPerms/"; 
  o <- 11; 
  ROI <- "BG_LR_CaNaPu_native";
  PAIR1 <- "upredCOR"; PAIR2 <- "upredINCOR";
  do.perms <- 1:5; perm.key <- "a1"; 
}
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="")
#SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");   # SUBS <- c(1003:1009, 1011:1019);   # SUBS <- "sub1003";
perm.lbls <- read.table(paste(perm.path, "l1outPermlbls_complete.txt", sep=""));  # all possible relabelings, plus the true one (first row)
if (max(do.perms) > nrow(perm.lbls)) { stop("do.perms doesn't match perm.lbls"); }

# flags for type of scaling to do.
DO_ROW_SCALING <- FALSE;
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
if (DO_ROW_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "RowSc", sep=""); }
if (DO_DEFAULT_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "DefaultSc", sep=""); }

# make a dataframe with all the data we need to classify this pair of stimuli
tbl <- read.table(gzfile(paste(inpath, ROI, "_", PAIR1, "_avg_allSubs.txt.gz", sep="")));
tbl <- data.frame(PAIR1, tbl);   # add the trial type label - not in the average file, just its filename
colnames(tbl)[1] <- "expType";   # fix the new column's title

temp.tbl <- read.table(gzfile(paste(inpath, ROI, "_", PAIR2, "_avg_allSubs.txt.gz", sep="")));
temp.tbl <- data.frame(PAIR2, temp.tbl);   # add the trial type label - not in the average file, just its filename
colnames(temp.tbl)[1] <- "expType";   # fix the new column's title
if (nrow(tbl) != nrow(temp.tbl) | ncol(tbl) != ncol(temp.tbl)) { stop("tbl and temp.tbl different sizes"); }

tbl <- rbind(tbl, temp.tbl);  # put both into one dataframe
rm(temp.tbl);  # take the temporary table out of memory
if (length(summary(tbl$expType)) != 2) { stop("not two expType in tbl"); }
if (summary(tbl$expType)[[1]] != summary(tbl$expType)[[2]]) { stop("number of examples not balanced"); }
if (length(unique(tbl$subID)) != length(SUBS)) { stop("number of subjects don't match in tbl$subID and SUBS"); }

tbl <- subset(tbl, offset == OFFSETS[o]);  # all the data for one offset: where we need to put the new labels.
# check o.tbl: the rows should already be sorted by expType, and the number of rows should match the number of new labels in perm.lbls
if (nrow(tbl) != ncol(perm.lbls)) { stop("nrow(tbl) != ncol(perm.lbls)"); }
if (length(unique(tbl$expType[1:length(SUBS)])) != 1) { stop("length(unique(tbl$expType[1:length(SUBS)])) != 1"); }

if (DO_ROW_SCALING == TRUE) {  # do row scaling (all voxels in each volume)
  scaled <- tbl[,4:nrow(tbl)];   # get voxel columns; cols 1 to 4 are labels
  scaled <- scale(t(scaled));  # scale does columns, so transpose
  tbl <- data.frame(tbl[,1:3], t(scaled));  # put label columns back on, transpose back
  rm(scaled);
}   

# do the leave-one-subject-out cross-validation
rtbl <- data.frame(array(NA, c(nrow(perm.lbls), length(SUBS)+4)));   # results table (will be written out)
colnames(rtbl) <- c("perm.num", "offset", "pair", paste(SUBS, "out", sep=""), "avg.acc");
ctr <- 1;  # row counter for rtbl

# put on the labels and classify this permutation.
for (perm.num in do.perms) {    # perm.num <- 1;
  rtbl[ctr,2] <- OFFSETS[o];
  rtbl[ctr,3] <- paste(PAIR1, PAIR2, sep="_");
  rtbl[ctr,1] <- perm.num - 1;   # perm.num - 1 so the true labeling gets called 0 in the results file (my convention)
  tbl$expType <- as.factor(t(perm.lbls[perm.num,]));  # now a and b instead of the real labels, but it doesn't matter.
  for (sn in 1:length(SUBS)) {     # sn <- 1;
    train.set <- subset(tbl, subID != SUBS[sn]);  # everyone but the testing subject
    test.set <- subset(tbl, subID == SUBS[sn]);  # just the testing subject
    # double-check the training and testing data a bit
    if (summary(train.set$expType)[[1]] != summary(train.set$expType)[[2]]) { stop("training data not balanced"); }  
    if (summary(test.set$expType)[[1]] != summary(test.set$expType)[[2]]) { stop("testing data not balanced"); }  
    if (nrow(train.set) != (length(SUBS)-1)*2) { stop("not the right number of rows in the training data"); }
    if (nrow(test.set) != 2) { stop("not two rows in the testing data"); }
    
    rtbl[ctr,sn+3] <- doSVM(train.set, test.set);    # actually classify!
  }
  rtbl[ctr,ncol(rtbl)] <- mean(as.vector(rtbl[ctr,4:(length(SUBS)+3)], mode="numeric"));   # average accuracy
  ctr <- ctr + 1;
}
if (SCALING_LABEL == "_") { SCALING_LABEL <- "_noScaling"; }
rownames(rtbl) <- 1:dim(rtbl)[1];
write.table(rtbl, gzfile(paste(outpath, "l1subOut_", OFFSETS[o], "_", ROI, "_", PAIR1, "_", PAIR2, SCALING_LABEL, "_comPerms_", perm.key, ".txt.gz", sep="")));

##################################################################################################################################################
# make the batch files

rm(list=ls());

path <- "d:/temp/batches/";  

bs <- c("pfc", "bg")
cs <- c("g", "e", "r");
for (a in 1:11) {
  for (b in 1:2) {
    for (i in 1:3) {
      for (p in 1:9) {
        fout <- file(paste(path, "job_", a, "_", bs[b], "_", cs[i], "_", letters[p], sep=""), "wt");
        cat("#!/bin/bash\n", file=fout)
        cat(paste("#PBS -N job_", a, "_", bs[b], "_", cs[i], "_", letters[p], "\n", sep=""), file=fout)
        cat("\n", file=fout);
        cat("#PBS -l nodes=1:ppn=1\n", file=fout);
        cat(paste("#PBS -l walltime=15:00:00\n", sep=""), file=fout);
        cat("\n", file=fout);
        cat(paste("cd $HOME/compPerms\n", sep=""), file=fout);
        cat("\n", file=fout);
        cat(paste("/export/R-2.14.2/bin/R --no-restore --no-save --no-readline ", a, " ", b, " ", i, " ", p, " < completePerms.R\n", sep=""), file=fout);
        close(fout); unlink(fout);
      }
    }
  }
}


# make the shell script to start all the jobs
b <- 2;
fout <- file(paste(path, "startAll_", bs[b], "jobs.sh", sep=""), "wt");
cat("# submit all the jobs listed (they should be in the directory cd to here)\n", file=fout)
cat("# by typing ./startAllJobs.sh at the command prompt in putty.\n", file=fout)
cat("# if it complains about Permission denied, type\n", file=fout)
cat("# chmod 755 startAllJobs.sh at the command prompt then try running it again.\n", file=fout)
cat("\n", file=fout)

cat(paste("cd jobFiles\n", sep=""), file=fout);
cat("\n", file=fout)
for (a in 1:11) {
  for (i in 1:3) {
    for (p in 1:9) {
      cat(paste("set jid=`qsub job_", a, "_", bs[b], "_", cs[i], "_", letters[p], "`\n", sep=""), file=fout); 
      cat("sleep 5\n", file=fout);   # wait five seconds before submitting the next job
      cat(paste('echo "submitted job_', a, "_", bs[b], "_", cs[i], "_", letters[p], '"\n', sep=""), file=fout)
    }
  }
}
cat("\n", file=fout)
cat('echo "all done submitting the jobs"\n', file=fout)
close(fout); unlink(fout);


##################################################################################################################################################

















































#
