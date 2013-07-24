make an error   # do not source this file
##################################################################################################################################################################
# Read and convert the event/onset log files. These files are in .fidl format. Offset 0 occurs at the onset of the update cue. No time shift
# was done to account for lag. For half tr's the script rounds up for .5n and randomly rounds up or down for .5 even.
##################################################################################################################################################################

rm(list=ls())

inpath <- "C:/maile/svnFiles/plein/consulting/Alan/";
outpath <- inpath;

SUBS <- c(1003:1009,1011:1019)
#SUBS <- 1003;   # just for Jo; put SUB into a loop

for (SUB in SUBS) {    
  tbl <- read.delim(paste(inpath, SUB, "errTrltype.fidl", sep="")); 
  #           X2 delayempty memset delay1 smain smainnp rmain rmainup rmainnp update updateop updatenp delaygreen delayred probe delayemptyERR memsetERR delay1ERR smainERR
  # 1      3.206          1      3     NA    NA      NA    NA      NA      NA     NA       NA       NA         NA       NA    NA            NA        NA        NA       NA
  # 2      6.311          2      7     NA    NA      NA    NA      NA      NA     NA       NA       NA         NA       NA    NA            NA        NA        NA       NA
  types <- colnames(tbl);
  types <- types[2:length(types)]
  tbl <- tbl[,1:3];
  colnames(tbl) <- c("onsetTime", "eventType", "duration");
  tbl$eventType <- types[tbl$eventType + 1]   # change number codes to the text descriptions
  
  # change the onsetTime to TR (volume number). 210 TR/run; 12 runs; 2 second TR.
  timeInTR <- tbl$onsetTime/2
  tbl <- data.frame(timeInTR, tbl);
  durationInTR <- tbl$duration/2
  tbl <- data.frame(durationInTR, tbl);
  diffInSec <- c(0,diff(tbl$onsetTime));
  tbl <- data.frame(diffInSec, tbl);
  
  # 11.537 or so in diffInSec corresponds to the start of a new run
  runNum <- rep(1, dim(tbl)[1]);
  runTimeTR <- rep(NA, dim(tbl)[1]);
  ctr <- 0;
  for (i in 1:dim(tbl)[1]) { 
    if (tbl$diffInSec[i] > 11 & tbl$diffInSec[i] < 13) { 
      ctr <- ctr + 1; 
    }
    runNum[i] <- runNum[i] + ctr;
    runTimeTR[i] <- tbl$onsetTime[i]/2 - 210*ctr;
  }
  tbl <- data.frame(runTimeTR, tbl);
  cts <- summary(as.factor(runNum))
  for (i in 1:length(cts)) { if (cts[[i]] != 60) { stop("not 60 rows in a run"); } }
  tbl <- data.frame(runNum, tbl);
  #cbind(round(tbl$runTimeTR), tbl)
  
  # try to make the repeated labels
  alllbls <- rep("ITI", 210*12);
  for (i in 1:dim(tbl)[1]) { 
    thisLabel <- tbl$eventType[i];
    thisDuration <- floor(tbl$durationInTR[i]);
    thisStart <- round(tbl$timeInTR[i]);
    if (thisDuration == 1) { alllbls[thisStart] <- thisLabel; }
    if (thisDuration == 3) {
      alllbls[thisStart] <- thisLabel;
      alllbls[thisStart + 1] <- thisLabel;
      alllbls[thisStart + 2] <- thisLabel;
    }
    if (thisDuration != 1 & thisDuration != 3) { stop("unexpected thisDuration"); }
  }
  alllbls <- data.frame(alllbls, rep(1:210,12), c(rep(1,210), rep(2,210), rep(3,210), rep(4,210), rep(5,210), rep(6,210), rep(7,210), rep(8,210), rep(9,210), rep(10,210), rep(11,210), rep(12,210)))
  colnames(alllbls) <- c("eventType", "TRnumber", "runNumber");   
  timeInMsec <- (alllbls$TRnumber * 2000) - 3000;  # 10 June 2013, Jo note: I don't see why we - 3000 here, but I don't think this is column is not used.
  alllbls <- data.frame(alllbls, timeInMsec)
  
  write.table(alllbls, paste(outpath, SUB, "_rewrittenLog.txt", sep=""));
}
####################################################################################################################################################################
# prepare the functional images for MVPA
# convert all of the images from nifti to arrays (just numbers). this converts all of the files in specific directories.
####################################################################################################################################################################

library(oro.nifti);

rm(list=ls());  # clear out R's memory

inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";   # nii files
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";   # file output to here

SUBS <- c(1003:1009,1011:1019); # subject ID codes, matching the way the files are labeled
RUNS <- 1:12;  # run names

NUMVOLUMES <- 210;  # number of volumes in each run
XDIM <- 64;   # must be same order as dim(img); should be same ordering as in mricron
YDIM <- 64;
ZDIM <- 35;
NUMVOXELS <- XDIM * YDIM * ZDIM;

doUnwrap <- function(intbl) {    # unwrap the 3D matrix intbl into a long 1d array
    newVol <- rep(NA, NUMVOXELS);
    ctr <- 1;            
    for (k in 1:ZDIM) {
        for (j in 1:YDIM) { # k <- 1; j <- 1;   
            add <- intbl[,j,k];   # get the numbers
            add0 <- gsub("NaN", 0, add);   # check for NaN, replace with zero if found
            add0 <- as.numeric(add0);    # make sure all numbers
            newVol[ctr:(ctr + XDIM - 1)] <- add0;
            ctr <- ctr + XDIM;  # move the row counter
        }
    }
    
    return(newVol);
}

for (SUB in SUBS) {    
      #if (SUB==1003) {RUNS <- c(1:6,8:12)}
      #if (SUB==1004) {RUNS <- c(1:7,10)}
      #if (SUB==1019) {RUNS <- 1:6}

    for (RUN in RUNS) {   # RUN <- RUNS[1]; SUB <- SUBS[1];
        #img <- readNIfTI(paste(inpath, SUB, "/", SUB, "_brun", RUN, "_STD.nii", sep=""));   # by default this will adjust orientation as specified in the headers
        img <- readNIfTI(paste(inpath, SUB, "_brun", RUN, "_STD.nii", sep=""));   # by default this will adjust orientation as specified in the headers
        if (XDIM != dim(img)[1]) { stop("x dims do not match"); }
        if (YDIM != dim(img)[2]) { stop("y dims do not match"); }
        if (ZDIM != dim(img)[3]) { stop("z dims do not match"); }
        if (NUMVOLUMES != dim(img)[4]) { stop("number of volumes do not match"); }

        # turn into a big flat array
        bigtbl <- array(0.0, c(NUMVOXELS, NUMVOLUMES));  # a very big array to hold all the volumes for this person and run
        for (i in 1:NUMVOLUMES) {   # i <- 1;
            bigtbl[,i] <- doUnwrap(img[,,,i]);
        }
        write.table(bigtbl, gzfile(paste(outpath, "sub", SUB, "_run", RUN, ".gz", sep="")));  # this saves a text file, gzipped to make it smaller.
        # write.table works but can be very slow. the code below is often much faster, but saves the files as R-only binaries, so they aren't readable in anything else.
        #save(bigtbl, file=paste(outpath, "sub", SUB, "_run", r, ".Robj", sep=""), ascii=FALSE);
    }
}

####################################################################################################################################################################
# prepare the MASKS images for MVPA
# convert all of the images from nifti to arrays (just numbers). this converts all of the files in specific directories.
####################################################################################################################################################################

library(oro.nifti);

rm(list=ls());  # clear out R's memory

inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";   # nii files
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";   # file output to here
ROIS <- c("BG_LR_CaNaPu_native","PFC_mask_native","Parietal_mask_native"); # subject ID codes, matching the way the files are labeled
#ROI <- "Parietal_mask_native"; # subject ID codes, matching the way the files are labeled

NUMVOLUMES <- 1;  # 183 volumes in each run
XDIM <- 64;   # must be same order as dim(img); should be same ordering as in mricron
YDIM <- 64;
ZDIM <- 35;
NUMVOXELS <- XDIM * YDIM * ZDIM;

doUnwrap <- function(intbl) {    # unwrap the 3D matrix intbl into a long 1d array
    newVol <- rep(NA, NUMVOXELS);
    ctr <- 1;            
    for (k in 1:ZDIM) {
        for (j in 1:YDIM) { # k <- 1; j <- 1;   
            add <- intbl[,j,k];   # get the numbers
            add0 <- gsub("NaN", 0, add);   # check for NaN, replace with zero if found
            add0 <- as.numeric(add0);    # make sure all numbers
            newVol[ctr:(ctr + XDIM - 1)] <- add0;
            ctr <- ctr + XDIM;  # move the row counter
        }
    }
    
    return(newVol);
}


img <- readNIfTI(paste(inpath, ROI, ".nii", sep=""));   # by default this will adjust orientation as specified in the headers
        if (XDIM != dim(img)[1]) { stop("x dims do not match"); }
        if (YDIM != dim(img)[2]) { stop("y dims do not match"); }
        if (ZDIM != dim(img)[3]) { stop("z dims do not match"); }


# turn into a big flat array
bigtbl <- doUnwrap(img);
write.table(bigtbl, gzfile(paste(outpath, ROI, ".gz", sep="")));  # this saves a text file, gzipped to make it smaller.
        # write.table works but can be very slow. the code below is often much faster, but saves the files as R-only binaries, so they aren't readable in anything else.
        #save(bigtbl, file=paste(outpath, "sub", SUB, "_run", r, ".Robj", sep=""), ascii=FALSE);

####################################################################################################################################################################
# separate the functional images into ROIs and add label columns. normalized data, with one set of ROI masks for everyone.
# do not add event type labels yet, but just label by TR.
# combine files so just one output per ROI (not one per run), adding subID and run columns
##################################################################################################################################################

rm(list=ls());

ROIpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";   # path to the ROI mask nii files
inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step1/";   # gz-text files from previous step
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/";   # new files output to here
#ROIpath <- "c:/maile/svnFiles/plein/consulting/Alan/AlanCode/";   # path to the ROI mask nii files
#inpath <- "c:/maile/svnFiles/plein/consulting/Alan/AlanCode/";   # gz-text files from previous step
#outpath <- "c:/maile/svnFiles/plein/consulting/Alan/AlanCode/";   # new files output to here

SUBS <- c(1003:1009,1011:1019); # subject ID codes, matching the way the files are labeled
RUNS <- 1:12;  # run names
#ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");   # ROI names (of the text files)
ROIS <- c("Parietal_mask_native");

NUMVOLUMES <- 210;  # number of volumes in each run
TRs <- 1:NUMVOLUMES;  # to use for the TR column in the output file
XDIM <- 64;   # must be same order as dim(img); should be same ordering as in mricron
YDIM <- 64;
ZDIM <- 35;
NUMVOXELS <- XDIM * YDIM * ZDIM;

#for (ROI in ROIS) {    # To do just one roi select that and use this instead of loop ROI <- ROIS[1];
    region <- read.table(paste(ROIpath, ROI, ".gz", sep=""));
    if (dim(region)[1] != NUMVOXELS) { stop("number of voxels in roi file not NUMVOXELS"); }  # safety code
    inds <- which(region[,1] > 0, arr.ind=TRUE);  # important line: finds all the locations with a value > 0 and stores the indicies
    # FYI: NUMROIVOXELS <- length(inds);  # how many voxels will be in the ROI files.
    if (exists('outtbl') == TRUE) { rm(outtbl); }  # we will make this file, but need to make a new one for each ROI. this gets rid of the variable.
    for (SUB in SUBS) { 
      RUNS <- 1:12;  # run names
      if (SUB==1003) {RUNS <- c(1:6,8:12)}
      if (SUB==1004) {RUNS <- c(1:7,10)}
      if (SUB==1019) {RUNS <- 1:6}
        for (RUN in RUNS) {    # RUN <- RUNS[1]; SUB <- SUBS[1];
            img <- read.table(gzfile(paste(inpath, "sub", SUB, "_run", RUN, ".gz", sep="")), comment.char="", header=TRUE);  # since saved with write.table; can be slow
            if (dim(img)[1] != NUMVOXELS) { stop("number of voxels in img file not NUMVOXELS"); }
            img <- img[inds,];   # just keep the voxels for this ROI; this line gets rid of the other rows.
            
            # add the label columns and row/column names
            if (length(TRs) != dim(img)[2]) { stop("not expected number of volumes"); }
            img <- data.frame(paste("sub", SUB, sep=""), RUN, TRs, t(img)); # now have subID, run, voxel1, voxel2, ... voxeln
            colnames(img) <- c("subID", "run", "TR", paste("v", 1:(dim(img)[2] - 3), sep=""));  # put the column names on
            if (exists('outtbl') == FALSE) { outtbl <- img; } else { outtbl <- rbind(outtbl, img); }   # puts all the peoples' data together. not the fastest way.
            rm(img);
        }
        rownames(outtbl) <- 1:(dim(outtbl)[1]);   # number the rows to be neat
        write.table(outtbl, gzfile(paste(outpath, "allData_", ROI, ".gz", sep=""))); 
    }
 }

####################################################################################################################################################################
# remove invariant voxels, all subjects together since normalized. aka "cleaning" the files.
##################################################################################################################################################

rm(list=ls());

ROIpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";
inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/";   # gz-text files from previous step
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/";   # new files output to here
#inpath <- "c:/maile/svnFiles/plein/consulting/Alan/AlanCode/";   # gz-text files from previous step
#outpath <- "c:/maile/svnFiles/plein/consulting/Alan/AlanCode/";   # new files output to here

SUBS <- c(1003:1009,1011:1019); # subject ID codes, matching the way the files are labeled
RUNS <- 1:12;  # run names
#ROIS <- c("Parietal_mask_native")
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");   # ROI names (of the text files)

#for (ROI in ROIS) {    # To do just one roi select that and use this instead of loop ROI <- ROIS[1];
   intbl <- read.table(gzfile(paste(inpath, "allData_", ROI, ".gz", sep="")), comment.char="", header=TRUE);
   
    FIRSTVOXEL <- which(colnames(intbl) == "v1");  # index at which "clean" voxels should start; column number of the v1 column in the input files.
    NUMVOXELS <- dim(intbl)[2] - FIRSTVOXEL + 1;  # the number of voxels in the ROI

    NUMVOLUMES <- max(summary(intbl$subID));  # number of volumes in each run; should be 210
    inds <- array(NA, c(length(SUBS), NUMVOLUMES));
    for (s in 1:length(SUBS)) {     # here we store the rows where each subjects' data is located         
  tmp <- which(intbl$subID == paste("sub", SUBS[s], sep=""));
  if (length(tmp) < 1) { stop(paste("no rows for sub", SUBS[s], "in intbl")); }
        inds[s,1:length(tmp)] <- tmp; 
    }
    
    # loop through all the people to see if they are zero-variance for the voxel
    vec <- rep(0, NUMVOXELS);   # saving this vector is useful for making images from the classification results later
    voxinds <- FIRSTVOXEL:(dim(intbl)[2]);
    for (i in 1:length(voxinds)) {    # i <- FIRSTVOXEL;
        include <- "yes";
        for (s in 1:length(SUBS)) { 
            sinds <- inds[s,];  # get this subject's indicies
            sinds <- subset(sinds, is.na(sinds)==FALSE);  # some subs have a missing run
            if (var(intbl[sinds,voxinds[i]])==0) { include <- "no"; break; }  
        }
        if (include=="yes") { vec[i] <- 1; }  # change this entry to 1 if include is still "yes"
    }
    if (length(vec) != NUMVOXELS) { stop("length(vec) != NUMVOXELS"); }
    vec <- c(rep(2, FIRSTVOXEL-1), vec);   # want to keep the label columns so add them on the front
    write.table(vec, paste(outpath, "vec_", ROI, ".txt", sep=""));
    inds <- which(vec == 0, arr.ind=TRUE);  # output voxels if the ROI is > zero
    
    # now output the cleaned files
    if (length(inds) > 0) {     # take out if some invariant voxels
        intbl <- intbl[,-inds]; 
        tds <- 1:(dim(intbl)[2] - FIRSTVOXEL + 1);
        tds <- paste("cv", tds, sep="");   # make names for the voxel columns
        colnames(intbl)[FIRSTVOXEL:dim(intbl)[2]] <- tds;
    } 
    write.table(intbl, paste(outpath, ROI, "_groupCleaned2.txt", sep="")); 
}

##################################################################################################################################################
# add on the eventType column from the log files and save one file per person and ROI. This WAS mean-subtraction, but now just does the labeling.
# these are the files used in the half-split (within each person individually) classification.
##################################################################################################################################################

rm(list=ls()); on.computer <- "Alan";
# rm(list=ls()); on.computer <- "Jo";

if (on.computer == "Alan") {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/";    # where the voxel-output files are from the 'cleaning'
  logpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/logs/";   # where the log files are (R-written version)
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/"; 
}
if (on.computer == "Jo") {
  inpath <- "d:/temp/Alan/";    # where the voxel-output files are from the 'cleaning'
  logpath <- "d:/temp/Alan/logs/";   # where the log files are (R-written version)
  outpath <- "d:/temp/Alan/for_halfSplit/"; 
}
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native","Parietal_mask_native");
SUBS <- c(1003:1009, 1011:1019);

FIRSTVOXEL <- 4;   # column number of first voxel column (v1)
NUMVOLUMES <- 210;  # 210 volumes for everyone in each run (no missing volumes, just whole runs missing)

for (ROI in ROIS) {     # ROI <- ROIS[1]; 
  # get all of the data for everyone into memory
  tbl <- read.table(gzfile(paste(inpath, ROI, "_groupCleaned2.txt.gz", sep="")), comment.char=""); # read in the data file
  if (colnames(tbl)[FIRSTVOXEL] != "v1") { stop("v1 column not at FIRSTVOXEL"); }
  for (SUB in SUBS) {     # SUB <- SUBS[1];
    stbl <- subset(tbl, tbl$subID == paste("sub", SUB, sep=""));

    # add the log file columns
    ltbl <- read.table(paste(logpath, SUB, "_rewrittenLog.txt", sep=""));
    if (nrow(ltbl) != nrow(stbl)) { stop("row counts for ltbl and stbl don't match"); }

    inds_match <- which(ltbl$runNumber == stbl$run & ltbl$TRnumber == stbl$TR)
    if (length(inds_match) != nrow(stbl)) {  # not all match, so fix the row order to match ltbl
      stbl <- stbl[order(stbl$run, stbl$TR),]; # check help for order
      inds_match2 <- which(ltbl$runNumber == stbl$run & ltbl$TRnumber == stbl$TR)
      if (length(inds_match2) != nrow(stbl)) { stop("inds_match2 don't!"); }
    }
    stbl <- cbind(ltbl$eventType, stbl);
    colnames(stbl)[1] <- "eventType";  # fix the column name; don't want it it be ltbl$eventType
    stbl$subID <- paste("sub", SUB, sep="");  # put back to the correct string; turned into a level in earlier steps
    write.table(stbl, gzfile(paste(outpath, "sub", SUB, "_", ROI, ".txt.gz", sep="")));  # one file per subject and ROI
    rm(stbl, ltbl, inds_match);
  }
}

##################################################################################################################################################
# prepare the input datafiles for group classification: average the examples together so one per person, ROI, and type.
# makes one file for each ROI and stimulus type (upgreenCOR, upemptyINCOR, etc); all people and timepoints in each file.
# the code is slow since it reads in each person's datafiles several times, but only has to be run once.
# change the labeling code here (accurate-only or whatever) to make new files as needed.
##################################################################################################################################################

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/avg_allSubs/"; 
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
# Half-Split classification
# note: cannot handle subjects that don't have all runs
#
# Prepare Perms
# First, figure out how many of each type are in each run for each person.
# this is a ROI-based analysis, classifying three types of stimuli (upempty, upred, upgreen), pairwise.
# always default scaling. within-subjects classification, each timepoint separately. There are missings and imbalance.
# partition on the runs, but not each run separately (not all stimuli in all runs); rather groups of runs at once.
# only need to permute labels once, not for each pair separately.

# should permute the upempty, upred, upgreen trial-type labels only; keep the time structure intact. probably best to permute labels within
# the partitions (several runs) if possible. Might need a different scheme for every person ...
##################################################################################################################################################
# figure out the counts: what permutations do we need?

rm(list=ls());

#inpath <- "c:/maile/svnFiles/plein/consulting/Alan/dataFilesFromServer/"; 
inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/"
#outpath <- "c:/maile/svnFiles/plein/consulting/Alan/setupPerms/";

SUBS <- c(1003:1009, 1011:1019);
#SUBS <- c(1003, 1005);
ROI <- "BG_LR_CaNaPu_native";
PAIRS <- c("upempty","upred","upgreen"); 
#> tbl[1:5,1:7]
#  eventType   subID run TR         v1          v2         v3
#1       ITI sub1008   1  1  -8.387116 -10.2185471 -18.805070
#2    memset sub1008   1  2 -11.601106  -2.8148606  -8.486344
#3    delay1 sub1008   1  3  -8.180573  -0.7729905  -7.959977

rtbl <- array(NA, c(length(SUBS)*length(PAIRS), 4));   # results table
rctr <- 1;
for (sn in 1:length(SUBS)) {     # sn <- 1;
  tbl <- read.table(gzfile(paste(inpath, "sub", SUBS[sn], "_", ROI, "_meanSub.gz", sep="")), comment.char=""); # read in the data    
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
  
  RUNS <- sort(unique(tbl$run));
  for (p in 1:length(PAIRS)) {   # p <- 1;
    rtbl[rctr,1] <- paste("sub", SUBS[sn]);
    rtbl[rctr,2] <- PAIRS[p];
    # find the start of each trial that's type pair1 or pair2
    inds <- which(tbl$eventType == PAIRS[p]);
    tinds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
    useinds <- inds[tinds];  # data to classify; offset from time point 0: start of each trial. this is just OFFSET == 0
    if (length(which(useinds < 0)) > 0) { stop("yes, have negative rows"); }
    t0tbl <- tbl[useinds,];
    
    numInTrain <- round(length(RUNS)/2);
    trainRuns <- RUNS[1:numInTrain];
    testRuns <- RUNS[(numInTrain+1):length(RUNS)];
    if (length(intersect(trainRuns, testRuns)) > 0) { stop("length(intersect(trainRuns, testRuns)) > 0"); }
    if (length(trainRuns) + length(testRuns) != length(RUNS)) { stop("length(trainRuns + testRuns) != length(RUNS)"); }
    
    testinds <- which(is.element(t0tbl$run, testRuns));  # rows with test-set runs
    traininds <- (1:dim(t0tbl)[1])[-testinds];  # all the others are training-set runs
    if (length(testinds) == 0) { stop("no testing data"); }  # can't be empty
    if (length(traininds) == 0) { stop("no training data"); }
    if (length(intersect(testinds, traininds)) > 0) { stop("overlap in test and traininds"); }
    allTest <- t0tbl[testinds,];  # subset into training and testing sets
    allTrain <- t0tbl[traininds,];
    allTest$eventType <- factor(allTest$eventType);   # gets rid of empty factor levels (simplifies balancing)
    allTrain$eventType <- factor(allTrain$eventType);
    rtbl[rctr,3] <- dim(allTrain)[1];
    rtbl[rctr,4] <- dim(allTest)[1];
    
    rctr <- rctr + 1;
  }
}
colnames(rtbl) <- c("subID", "item", "train1", "test1");
rtbl
write.table(rtbl, "c:/maile/svnFiles/plein/consulting/Alan/setupPerms/countTable.txt");

tbl <- read.table("c:/maile/svnFiles/plein/consulting/Alan/setupPerms/countTable.txt");
#> tbl
#      subID    item train1 test1
#1  sub 1003 upempty     13    18
#2  sub 1003   upred     31    18
#3  sub 1003 upgreen     19    22
#4  sub 1004 upempty     12     2

# these people are included: paste("sub", c(1003, 1005:1009, 1011:1019), sep="");

# probably easiest just to read the values out of tbl.
#           1003    1004     1005       1006        1007         1008     1009         1011     1012          1013      1014     1015
cts <- c(13,19,18, 11,14,9, 18,16,23, 13,20,18,25, 14,21,15,24, 9,14,19, 12,21,17,24, 12,19,11, 13,22,12,17, 14,21,24, 10,20,14, 13,20,12,19,
         #    1016     1017      1018     1019
         13,17,14, 11,16,17, 14,22,23, 8,11,6,9);
sort(unique(cts))   # 6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25

# these are the number of examples of each class in the training set, so smallest is 6 "a" and 6 "b".

# smallest number of each row.
#           1003    1004     1005       1006        1007         1008     1009         1011     1012          1013      1014     1015
cts <- c(13,18,19, 2,15,11, 9,23,18,  13,25,20,   14,24,21,    9,17,14, 12,24,21,   11,19,19,  12,17,22,  14,24,21,   10,20,20, 12,20,19,
         #    1016     1017      1018     1019
         13,14,17, 11,17,16, 14,22,14,   6,11,9);
sort(unique(cts))   #  2  6  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25


##################################################################################################################################################
# enumerate all of the permutations for the case of 6 examples of each class

library(gtools)

# datafile for classification. half the runs into training, half into testing.  13 of each for this person. Ignore run structure during permutations.
#> useTrain[,1:10]
#     eventType   subID run  TR         v1         v2          v3          v4           v5          v6
#160    upempty sub1003   1 160   4.872750 -17.588802 -21.8477757 -11.1878543  -8.25525251  -4.1752880
#178    upempty sub1003   1 178  19.613595  -2.468319 -24.3555882 -17.8549686  -6.42877546  -4.5097607

# will need to sort by eventType before the permutations.

rm(list=ls());

# 6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
NUMS <- c(6, 8, 9, 10);  # , 18, 19, 20, 21, 22, 23, 24, 25
for (NUM in NUMS) {     #    NUM <- 6;   # how many in each class
  ans <- 0;  # figure out how many relabelings are possible. this is by how many of the true-labeled data of each type have their label changed.
  for (i in 0:NUM) { ans <- ans + choose(NUM,i); }
  ans
  
  newlbls <- array(NA, c(ans, NUM*2));
  ctr <- 1;
  # row 1: 0 relabelings (1)
  newlbls[ctr,] <- c(rep("a", NUM), rep("b", NUM));
  ctr <- ctr + 1;
  
  for (j in 1:NUM) {    # 1 relabeling
    inds <- combinations(NUM,j)
    for (i in 1:dim(inds)[1]) {   # i <- 1;
      a <- rep("a", NUM);
      a[inds[i,]] <- "b";
      
      b <- rep("b", NUM);
      b[inds[i,]] <- "a";
      newlbls[ctr,] <- c(a, b);
      ctr <- ctr + 1;
    }
  }
  
  # confirm no repeated rows
  for (i in 1:(dim(newlbls)[1] -1)) {
    for (j in (i+1):dim(newlbls)[1]) {   # i <- 1; j <- 2;
      if (length(which(newlbls[i,] == newlbls[j,])) == (NUM*2)) { stop("found a match"); break; }
    }
  }
  write.table(newlbls, paste("c:/maile/svnFiles/plein/consulting/Alan/setupPerms/", NUM, "eachTable.txt", sep=""));
}

##################################################################################################################################################
# make a thousand permutations at random if there are too many to enumerate them all.

rm(list=ls());

hasMatch <- function(thislbl, otherlbls) {    # confirm no repeated rows
  ans <- FALSE;   
  for (j in 1:dim(otherlbls)[1]) {
    if (length(which(thislbl == otherlbls[j,])) == (NUM*2)) { ans <- TRUE; break; }
  }
  
  return(ans)
}

NUMS <- 11:25;
for (NUM in NUMS) {     #    NUM <- 15;   # how many in each class
  lbls <- c(rep("a", NUM), rep("b", NUM))
  newlbls <- array(NA, c(1001, NUM*2));
  newlbls[1,] <- lbls;    # row 1: 0 relabelings
  for (i in 2:1001) {    # i <- 2;
    albl <- sample(lbls);
    while (hasMatch(albl, newlbls[1:i,]) == TRUE) { albl <- sample(lbls); }    # could run forever ...
    newlbls[i,] <- albl;
  }
  write.table(newlbls, paste("c:/maile/svnFiles/plein/consulting/Alan/setupPerms/", NUM, "eachTable.txt", sep=""));
}

##################################################################################################################################################
# Half-Split Class and Perm
# C L U S T E R 
# Need "permpath" files (from previous step), job files (doBigPerms/jobFiles/), inputfiles from /step2/MeanSub/
# Run by submitting ./startAll_pair1.sh (or "pair2.sh", etc). Runs pair for every subject for individual roi. Must manually change roi
# in script and re-run. Subject 1004 can't be run.
##################################################################################################################################################

library(e1071);  # R interface to libsvm

rm(list=ls()); ONNIL <- FALSE; ONCLUSTER <- TRUE; JOCOMPUTER <- FALSE;
# rm(list=ls()); ONNIL <- FALSE; ONCLUSTER <- FALSE; JOCOMPUTER <- TRUE;

SUBS <- paste("sub", c(1003, 1005:1009, 1011:1019), sep="");
#ROI <- "PFC_mask_native"     
#ROI <- "BG_LR_CaNaPu_native"
ROI <- "Parietal_mask_native"

OFFSETS <- c("precue", "postcue"); # timebins to classify
# in the datatables, eventType follows a pattern: ITI - memset - delay - upgreen - ITI - probe - ITI - ITI . Offset 0 is the first upgreen (or whatever).
# ITI is used as a filler in eventType, indicating a pause, and occurs during trials, not just in between trials.
# also, the eventType column is in real time, not adjusted for any lag in the BOLD.

if (ONCLUSTER == TRUE) {
  inpath <- "/scratch/aceaser/input/";
  outpath <- "/scratch/aceaser/output/"; 
  permpath <- "/scratch/aceaser/permInputFiles/";   # location of 6eachTable.txt, 8eachTable.txt, etc.
  
  cA <- commandArgs();   
  num  <- as.numeric(cA[5]);  # which timepoints to run
  OFFSET <- OFFSETS[num];    # 2 total
  
  num  <- as.numeric(cA[6]);   # second number sent specifies the pair to run.
  if (num == 1) { PAIR1 <- "upempty"; PAIR2 <- "upgreen"; } 
  if (num == 2) { PAIR1 <- "upempty"; PAIR2 <- "upred"; } 
  if (num == 3) { PAIR1 <- "upgreen"; PAIR2 <- "upred";  }
  
  splt  <- as.numeric(cA[7]);   # and third the split, 1:10.
}
if (ONNIL == TRUE) {
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
  permpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/permInputFiles/";   # location of 6eachTable.txt, 8eachTable.txt, etc.
  SUB <- "sub1003";
  PAIR1 <- "upempty"; PAIR2 <- "upgreen";
}
if (JOCOMPUTER == TRUE) {
  inpath <- "d:/temp/"; 
  outpath <- "d:/temp/";
  permpath <- "c:/maile/svnFiles/plein/consulting/Alan/setupPerms/";
  SUB <- "sub1003";
  PAIR2 <- "upred"; PAIR1 <- "upgreen";
}


SEEDS <- c(51591, 36414, 56347, 38442, 20176, 51348, 89727, 67106, 23543, 52663);  # 10 random seeds, from sample(1:100000)[1:10]
NUMSPLITS <- 10;  
MAXPERMS <- 1024;  # the most permutations for any person and pair. there will usually be fewer than this, but this sets up the results table.
DO_DEFAULT_SCALING <- TRUE;    # flags for type of scaling to do.

doSVM <- function(train, test, DO_DEFAULT_SCALING) {  # test <- useTest; train <- useTrain;
  test <- subset(test, select=c(-subID, -run, -TR));  # get rid of non-classify or voxel columns
  train <- subset(train, select=c(-subID, -run, -TR));
  if (colnames(test)[2] != "v1" | colnames(train)[2] != "v1") { stop("v1 not found where expected"); }
  
  fit <- svm(eventType~., data=train, type="C-classification", kernel="linear", cost=1, scale=DO_DEFAULT_SCALING);  
  tree <- table(test$eventType, predict(fit, test));
  if (dim(tree)[2]==1 | dim(tree)[1]==1) { wrT <- 0.5; } else { wrT <- sum(diag(tree))/sum(tree); }
  
  return(wrT);
}

if (PAIR1 == PAIR2) { stop("PAIR1 == PAIR2"); }
for (SUB in SUBS) {  #SUB <- "sub1005";
  fname <- paste(inpath, SUB, "_", ROI, "_meanSub.gz", sep="");
  if (!file.exists(fname)) { stop(paste("missing file", fname)); }
  tbl <- read.table(gzfile(fname), comment.char=""); # read in the data
  RUNS <- sort(unique(tbl$run));
  if (length(RUNS) > 12) { stop("too many runs"); }   # try to catch if the input data is wrong
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
  
  # old code: use if running individual timepoints (one TR)
 # useinds <- inds[tinds] + OFFSET;  # data to classify; offset from time point 0: start of each trial
 # if (length(which(useinds < 0)) > 0) { stop("yes, have negative rows"); }
 # t0tbl <- tbl[useinds,];
 # if (OFFSET != 0) {  # need to get labels for which events these go with; offset labels are other things than what we're classifying.
 #   uselbls <- tbl$eventType[inds[tinds]];     
 #   if (length(uselbls) != dim(t0tbl)[1]) { stop("very wrong lengths for uselbls"); } else { t0tbl$eventType <- uselbls; }
 # }
  
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

  rtbl <- array(NA, c(NUMSPLITS*MAXPERMS, 7));   # results table
  subID <- rep(NA, NUMSPLITS*MAXPERMS);   # subject ID column of results table
  rowctr <- 1;
  
  numInTrain <- round(length(RUNS)/2);
  trainRuns <- RUNS[1:numInTrain];
  testRuns <- RUNS[(numInTrain+1):length(RUNS)];
  if (length(intersect(trainRuns, testRuns)) > 0) { stop("length(intersect(trainRuns, testRuns)) > 0"); }
  if (length(trainRuns) + length(testRuns) != length(RUNS)) { stop("length(trainRuns + testRuns) != length(RUNS)"); }
  
  testinds <- which(is.element(t0tbl$run, testRuns));  # rows with test-set runs
  traininds <- (1:dim(t0tbl)[1])[-testinds];  # all the others are training-set runs
  if (length(testinds) == 0) { stop("no testing data"); }  # can't be empty
  if (length(traininds) == 0) { stop("no training data"); }
  if (length(intersect(testinds, traininds)) > 0) { stop("overlap in test and traininds"); }
  
  allTest <- t0tbl[testinds,];  # subset into training and testing sets
  allTrain <- t0tbl[traininds,];
  allTest$eventType <- factor(allTest$eventType);   # gets rid of empty factor levels (simplifies balancing)
  allTrain$eventType <- factor(allTrain$eventType);
  
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
  
  
  # now the permutation testing. I'm permuting the labels on BOTH the training and testing sets. If there are the same number of permuted
  # labels for both than I just use them all. If there are too few of one then I bootstrap the smaller one to match the bigger.
  # read in the precomputed label files
  permlblsTrain <- read.table(paste(permpath, dim(useTrain)[1]/2, "eachTable.txt", sep=""));  # labels go by number of examples
  if (dim(permlblsTrain)[2] != dim(useTrain)[1]) { stop("dim(permlblsTrain)[2] != dim(useTrain)[1]"); }
  permlblsTest <- read.table(paste(permpath, dim(useTest)[1]/2, "eachTable.txt", sep=""));  # labels go by number of examples
  if (dim(permlblsTest)[2] != dim(useTest)[1]) { stop("dim(permlblsTest)[2] != dim(useTrain)[1]"); }
  
  if (dim(permlblsTrain)[1] < dim(permlblsTest)[1]) {  # too few training labels
    needrows <- array(NA, c((dim(permlblsTest)[1] - dim(permlblsTrain)[1]), dim(permlblsTrain)[2]))
    for (k in 1:dim(needrows)[1]) {   # k <- 1;
      needrows[k,] <- t(permlblsTrain[sample(2:dim(permlblsTrain)[1])[1],]);
    }
    permlblsTrain <- rbind(permlblsTrain, needrows);
    if (dim(permlblsTrain)[1] != dim(permlblsTest)[1]) { stop("bootstrapping permlblsTrain went wrong"); }
  }
  if (dim(permlblsTest)[1] < dim(permlblsTrain)[1]) {  # too few testing labels
    needrows <- array(NA, c((dim(permlblsTrain)[1] - dim(permlblsTest)[1]), dim(permlblsTest)[2]))
    for (k in 1:dim(needrows)[1]) {   # k <- 1;
      needrows[k,] <- t(permlblsTest[sample(2:dim(permlblsTest)[1])[1],]);
    }
    permlblsTest <- rbind(permlblsTest, needrows);
    if (dim(permlblsTrain)[1] != dim(permlblsTest)[1]) { stop("bootstrapping permlblsTest went wrong"); }
  }
  if (dim(permlblsTrain)[1] != dim(permlblsTest)[1]) { stop("dim(permlblsTrain)[1] != dim(permlblsTest)[1]"); }
  if (dim(permlblsTrain)[1] > MAXPERMS | dim(permlblsTrain)[1] < 64) { stop("dim(permlblsTrain)[1] > MAXPERMS | dim(permlblsTrain)[1] < 64"); }
  
  
  # sort by eventType so the permutation labels can be put on the rows.
  useTrain <- useTrain[order(useTrain$eventType),];
  useTest <- useTest[order(useTest$eventType),];
  
  # finally, do the permutations
  for (pnum in 1:dim(permlblsTrain)[1]) {   # pnum <- 1;
    rtbl[rowctr,1] <- splt;
    rtbl[rowctr,2] <- OFFSET;
    subID[rowctr] <- SUB;
    useTrain$eventType <- factor(t(permlblsTrain[pnum,]));  # put on the train labels for this perm
    useTest$eventType <- factor(t(permlblsTest[pnum,]));    # and the test ones.
    rtbl[rowctr,3] <- doSVM(useTrain, useTest, DO_DEFAULT_SCALING);  
    rtbl[rowctr,4] <- doSVM(useTest, useTrain, DO_DEFAULT_SCALING);  
    
    rtbl[rowctr,5] <- (dim(useTrain)[1])/2;
    rtbl[rowctr,6] <- (dim(useTest)[1])/2;
    rtbl[rowctr,7] <- pnum - 1;   # first row is the true labeling, which I call permutation zero from habit.
    rowctr <- rowctr + 1;
  }
  
  # save the results file
  rtbl <- rtbl[1:(rowctr-1),];   # take off extra rows
  subID <- subID[1:(rowctr-1)];
  avgProp <- apply(rtbl[,3:4], 1, mean, na.rm=TRUE);  # calculate the average for each row
  colnames(rtbl) <- c("splitNum", "timePoint", "firstTest", "secondTest", "1stTrainSize", "2ndTrainSize", "permNum");
  rtbl <- data.frame(subID, rtbl, avgProp);
  rownames(rtbl) <- 1:dim(rtbl)[1];
  write.table(rtbl, gzfile(paste(outpath, "splt", splt, "_", SUB, "_", OFFSET, "_", ROI, "_", PAIR1, "_", PAIR2, "_defSc_halfSplt_perms.txt.gz", sep="")));    
}

##################################################################################################################################################
# C L U S T E R
# Combine split files. Looks for all files from previous step and creates subject files per roi and pair combo.
##################################################################################################################################################

rm(list=ls());

inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/outputSplits/";
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/output/"; 
#inpath <- "/scratch/aceaser/output/";
#outpath <- "/scratch/aceaser/output/"; 
#SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");
SUBS <- paste("sub", c(1005), sep="");
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);
PAIRS <- c("upempty_upgreen", "upempty_upred", "upgreen_upred");

# ROI <- "PFC_mask_native"     
# ROI <- "BG_LR_CaNaPu_native"
 ROI <- "Parietal_mask_native"

for (SUB in SUBS) {  #SUB <- "sub1019";
  
  for (PAIR in PAIRS) {
    for (OFFSET in OFFSETS) {
      if (exists('bigtbl')) { rm(bigtbl); }
      for (splt in 1:10) {
        intbl <- read.table(gzfile(paste(inpath, "splt", splt, "_", SUB, "_", OFFSET, "_", ROI, "_", PAIR, "_defSc_halfSplt_perms.txt.gz", sep="")));    
        if (exists('bigtbl')) { bigtbl <- rbind(bigtbl, intbl); } else { bigtbl <- intbl; }
      }
      write.table(bigtbl, gzfile(paste(outpath, SUB, "_", OFFSET, "_", ROI, "_", PAIR, "_defSc_halfSplt_perms.txt.gz", sep="")));    
    }
  }
}

##################################################################################################################################################
# Combine subject output files; make results table output file giving group accuracy per time point, 
# significance from permutation, and indiv sub accuracy.
##################################################################################################################################################

rm(list=ls()); 

path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/output/";
set.seed(234984);  # so can reproduce the permutation p-values, hopefully they don't jitter much with different seeds.

OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5); 
SUBS <- paste("sub", c(1003, 1005:1009, 1011:1019), sep="");

PAIR1 <- "upempty"; PAIR2 <- "upgreen";
#PAIR1 <- "upempty"; PAIR2 <- "upred";
#PAIR1 <- "upgreen"; PAIR2 <- "upred";

#ROI <- "BG_LR_CaNaPu_native"    # "PFC_mask_native"
#ROI <- "PFC_mask_native"
ROI <- "Parietal_mask_native"

#> head(intbl)
#    subID splitNum timePoint firstTest secondTest X1stTrainSize X2ndTrainSize permNum   avgProp
#1 sub1003        1        -5 0.5277778  0.4615385            13            18       0 0.4946581
#2 sub1003        1        -5 0.3055556  0.4230769            13            18       1 0.3643162
#3 sub1003        1        -5 0.5277778  0.5384615            13            18       2 0.5331197

ptbl_all <- array(NA, c(length(OFFSETS),2));  # mean & p-value averaged across included subjects
ptbl_eachSub <- array(NA, c(length(OFFSETS), length(SUBS)));
for (o in 1:length(OFFSETS)) {   # o <- 2;
  temptbl <- array(NA, c(1001, length(SUBS)));  # put all the subjects' permutation accuracies into here
  for (s in 1:length(SUBS)) { # s <- 1;
    fname <- paste(path, SUBS[s], "_", OFFSETS[o], "_", ROI, "_", PAIR1, "_", PAIR2, "_defSc_halfSplt_perms.txt.gz", sep="")
    if (file.exists(fname)) {
      intbl <- read.table(gzfile(fname));
      # intbl has accuracies for each of the splits. average so one average accuracy per permNum.
      # permNum == 0 is the real-labeled accuracy
      PERMS <- unique(intbl$permNum);
      PERMS <- PERMS[which(PERMS != 0)];
      if (max(PERMS) > 1000) { PERMS <- sample(PERMS)[1:1000]; } # a few people have more than 1000 perms; subset in this case.
      means <- rep(NA, length(PERMS)+1)
      # real labeling first
      stbl <- subset(intbl, intbl$permNum == 0);
      if (dim(stbl)[1] != 10) { stop("(dim(stbl)[1] != 10)"); }
      means[1] <- mean(stbl$avgProp);
      for (p in 1:length(PERMS)) {   # p <- 1;
        stbl <- subset(intbl, intbl$permNum == PERMS[p])
        if (dim(stbl)[1] != 10) { stop("(dim(stbl)[1] != 10)"); }
        means[p+1] <- mean(stbl$avgProp);
      }
      temptbl[1:(length(PERMS)+1),s] <- means;
      
      # calculate the p-value for this person
      ptbl_eachSub[o,s] <- length(which(means[2:length(means)] > means[1]))/(length(PERMS)-1); 
    } else { print(paste("missing file:", fname)); }
  }
  
  # bootstrap the people with missings so that they have 1000 as well. this was recommended in Good (resampling methods pdf book), pg 105.
  for (s in 1:length(SUBS)) { # s <- 16;
    numMissings <- length(which(is.na(temptbl[,s])))
    if (numMissings > 0 & numMissings != dim(temptbl)[1]) {  # have a missing, so bootstrap. don't try if no data at all.
      means <- temptbl[2:1001,s];
      means <- means[which(is.na(means)==FALSE)];
      startRow <- min(which(is.na(temptbl[,s])));  # where to start putting in the bootstrapped rows
      for (i in startRow:1001) {
        if (!is.na(temptbl[i,s])) { stop("(!is.na(temptbl[i,s]))"); }
        temptbl[i,s] <- sample(means)[1]; 
      }
    }
  }
  
  # calculate the p-value averaging across everyone
  rm(means);
  means <- apply(temptbl, 1, mean, na.rm=TRUE);  # na.rm to take out all-NA columns
  if (length(means) != 1001 | length(which(is.na(means))) > 0) { stop("wrong means"); }
  ptbl_all[o,1] <- means[1]; 
  ptbl_all[o,2] <- length(which(means[2:length(means)] > means[1]))/(length(means)); 
}
colnames(ptbl_eachSub) <- SUBS;
colnames(ptbl_all) <- c("allSubsMean", "allSubsP");
cbind(OFFSETS, round(ptbl_all,3), round(ptbl_eachSub,3));  # first column is the 'all' p-value

apply(ptbl_eachSub, 1, mean)

out.tbl <- cbind(OFFSETS, ptbl_all, ptbl_eachSub, apply(ptbl_eachSub, 1, mean))
# either write or set an outpath, write out what the file name should be called.txt
#write.table(out.tbl, some.path)
     # write.table(out.tbl, gzfile(paste(outpath, ROI, "_", PAIR1, "_", PAIR2, "_defSc_halfSplt_perms.txt.gz", sep="")));    

######################################################################################################################################################################
