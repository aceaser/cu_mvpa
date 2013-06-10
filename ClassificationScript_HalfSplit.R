make an error   # do not source this file
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
  # intbl2 <-   then intbl <- rbind(intbl, intbl2)  # then rm(intbl2); to clean up
   
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
# perform mean-subtraction. also add on the eventType column from the log files and save one file per person and ROI.
##################################################################################################################################################

rm(list=ls());

ROIpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/";
inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/";    # where the voxel-output files are from the 'cleaning'
logpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/logs/";   # where the log files are (R-written version)
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/"; 
#ROIS <- "Parietal_mask_native";
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native","Parietal_mask_native");
# SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="")
SUBS <- c(1003:1009, 1011:1019);

FIRSTVOXEL <- 4;   # column number of first voxel column (v1)
NUMVOLUMES <- 210;  # 210 volumes for everyone in each run (no missing volumes, just whole runs missing)

#for (ROI in ROIS) {     # ROI <- ROIS[1]; 
    # get all of the data for everyone into memory
    tbl <- read.table(paste(inpath, ROI, "_groupCleaned2.txt", sep=""), comment.char=""); # read in the data file
    #tbl2 <- read.table(the rest of the data)
    #tbl <- rbind(tbl, tbl2);  # this will be very slow, but should be ok since it only has to be done once per ROI
    #rm(tbl2);  # get the extra table out of memory
    if (colnames(tbl)[FIRSTVOXEL] != "v1") { stop("v1 column not at FIRSTVOXEL"); }
    for (SUB in SUBS) {     # SUB <- 1007;
        stbl <- subset(tbl, tbl$subID == paste("sub", SUB, sep=""));
        RUNS <- unique(stbl$run);  # set RUNS here since it varies for different people
        outtbl <- data.frame(array(NA, dim(stbl)));  # the output file will be the same size as the input (but different voxel values)
        for (r in 1:length(RUNS)) {      # r <- 1;
            intbl <- subset(stbl, stbl$run == RUNS[r]);  # just the data for one person and one run
            if (dim(intbl)[1] != NUMVOLUMES) { stop(paste("not expected number of volumes for", ROI, SUB, RUN)); }
            
            vtbl <- intbl[,FIRSTVOXEL:dim(intbl)[2]];  # just the voxel columns (no label columns)
            means <- apply(vtbl, 2, mean);  # calculate the mean, column-wise (an average for each voxel)
            means <- matrix(data=rep(means, NUMVOLUMES), nrow=NUMVOLUMES, ncol=length(means), byrow=TRUE);  # repeat means to make a matrix
           
            vtbl <- vtbl - means;   # want mean-subtraction VOXEL-wise, not ROW-wise; transposing ensures the subtraction is what we want
            if (length((((r-1)*210)+1):(r*210)) != dim(vtbl)[1] | length((((r-1)*210)+1):(r*210)) != dim(intbl)[1]) { stop("row counts don't match"); }
            outtbl[(((r-1)*210)+1):(r*210),1:(FIRSTVOXEL-1)] <- intbl[,1:(FIRSTVOXEL-1)];
            outtbl[(((r-1)*210)+1):(r*210),FIRSTVOXEL:dim(intbl)[2]] <- vtbl;
        }
        colnames(outtbl) <- colnames(stbl);
        rownames(outtbl) <- 1:dim(outtbl)[1];
        
        # add the log file columns
        ltbl <- read.table(paste(logpath, SUB, "_rewrittenLog.txt", sep=""));
        if (dim(ltbl)[1] != dim(outtbl)[1]) { stop("row counts for ltbl and outtbl don't match"); }

        inds_match <- which(ltbl$runNumber == outtbl$run & ltbl$TRnumber == outtbl$TR)
	if (length(inds_match) != dim(outtbl)[1]) {  # not all match, so fix the row order to match ltbl
		outtbl <- outtbl[order(outtbl$run, outtbl$TR),]; # check help for order
        	inds_match2 <- which(ltbl$runNumber == outtbl$run & ltbl$TRnumber == outtbl$TR)
		if (length(inds_match2) != dim(outtbl)[1]) { stop("inds_match2 don't!"); }
	}
        outtbl <- data.frame(ltbl$eventType, outtbl);
        colnames(outtbl)[1] <- "eventType";  # fix the column name; don't want it it be ltbl$eventType
        outtbl$subID <- paste("sub", SUB, sep="");  # put back to the correct string; turned into a level in earlier steps
        write.table(outtbl, gzfile(paste(outpath, "sub", SUB, "_", ROI, "_meanSub.gz", sep="")));  # one file per subject and ROI
    }
}

# spot-check the mean-subtraction
vtbl[1:5,1:5]
intbl[1:5, FIRSTVOXEL:(FIRSTVOXEL+5)]
means[1,1:5]
intbl[1,FIRSTVOXEL] - means[1,1]
intbl[2,FIRSTVOXEL] - means[1,1]

intbl[1, (FIRSTVOXEL+100-1)]
vtbl[1,100]
means[1,100]
intbl[1,(FIRSTVOXEL+100-1)] - means[1,100]


intbl[10, (FIRSTVOXEL+100-1)]
vtbl[10,100]
means[1,100]
intbl[10,(FIRSTVOXEL+100-1)] - means[1,100]

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
# Cannot handle subjects that don't have all runs
# Outputs all subjects into one file per comparison per roi
##################################################################################################################################################

library(e1071);  # R interface to libsvm

rm(list=ls());

inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native","Parietal_mask_native");
#ROIS <- c("Parietal_mask_native");
#ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native");
SUBS <- paste("sub", c(1005:1009,1011:1018), sep="");   
#SUBS <- c(1003:1009, 1011:1019);
#SUBS <- c(1005:1009, 1011:1019);
# SUBS <- "sub1019";
SEEDS <- c(51591, 36414, 56347, 38442, 20176, 51348, 89727, 67106, 23543, 52663);  # 10 random seeds, from sample(1:100000)[1:10]
NUMSPLITS <- 10;

#FIRSTVOXEL <- 4;   # column number of first voxel column (v1)
#NUMVOLUMES <- 210;  # 210 volumes for everyone in each run (no missing volumes, just whole runs missing)

PAIR1S <- c("upempty","upempty","upgreen");  # the things to classify. will classify the first entry of PAIR1S with the first entry of PAIR2S, etc.
PAIR2S <- c("upgreen","upred","upred");  # each pair1 needs to be alphabetically before each corresponding pair2 for the balancing code to work.
#PAIR1S <- "upempty"
#PAIR2S <- "upgreen"

# flags for type of scaling to do.
DO_ROW_SCALING <- FALSE;
DO_DEFAULT_SCALING <- TRUE;
DO_COLUMN_SCALING <- FALSE;

doSVM <- function(train, test, DO_DEFAULT_SCALING) {  # test <- useTest; train <- useTrain;
  test <- subset(test, select=c(-subID, -run, -TR));  # get rid of non-classify or voxel columns
  train <- subset(train, select=c(-subID, -run, -TR));
  if (colnames(test)[2] != "v1" | colnames(train)[2] != "v1") { stop("v1 not found where expected"); }
  
  fit <- svm(eventType~., data=train, type="C-classification", kernel="linear", cost=1, scale=DO_DEFAULT_SCALING);  
  tree <- table(test$eventType, predict(fit, test));
  if (dim(tree)[2]==1 | dim(tree)[1]==1) { wrT <- 0.5; } else { wrT <- sum(diag(tree))/sum(tree); }
  
  return(wrT);
}


#> tbl[1:5,1:7]
#  eventType   subID run TR         v1          v2         v3
#1       ITI sub1008   1  1  -8.387116 -10.2185471 -18.805070
#2    memset sub1008   1  2 -11.601106  -2.8148606  -8.486344
#3    delay1 sub1008   1  3  -8.180573  -0.7729905  -7.959977


SCALING_LABEL <- "_";  # make a label for the output files showing the type of scaling used for this classification
if (DO_COLUMN_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "ColSc", sep=""); }
if (DO_ROW_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "RowSc", sep=""); }
if (DO_DEFAULT_SCALING == TRUE) { SCALING_LABEL <- paste(SCALING_LABEL, "DefaultSc", sep=""); }

STEPBY <- 4;  # one less than how many runs to include in each cross-validation fold.
MAXTESTS <- 4;  # maximum number of cross-validation folds that will be done for a person; 12/3. can be larger, but not smaller, than actual number of folds.
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify

if (length(PAIR1S) != length(PAIR2S)) { stop("length(PAIR1S) != length(PAIR2S)"); }
for (p in 1:length(PAIR1S)) {   # p <- 1;
  PAIR1 <- PAIR1S[p];  # get what we're classifying this time. 
  PAIR2 <- PAIR2S[p];
  if (PAIR1 == PAIR2) { stop("PAIR1 == PAIR2"); }
  for (ROI in ROIS) {     # ROI <- ROIS[1]; 
    rtbl <- array(NA, c(length(SUBS) * length(OFFSETS) * NUMSPLITS, (MAXTESTS+2)));   # results table
    subID <- rep(NA, length(SUBS) * length(OFFSETS) * NUMSPLITS);   # subject ID column of results table
    rowctr <- 1;
    for (SUB in SUBS) {     # SUB <- "sub1003";
      # NOTE: to set different cross-validation schemes for different people, set STEPBY here (if (SUB == WHOEVER) { STEPBY <- 4; })
      # but MAXTESTS needs to be set to the most cross-validation folds for everyone (or rtbl will be too small).
      tbl <- read.table(gzfile(paste(inpath, SUB, "_", ROI, "_meanSub.gz", sep="")), comment.char=""); # read in the data
      RUNS <- sort(unique(tbl$run));
      if (length(RUNS) > 12) { stop("too many runs"); }   # try to catch if the input data is wrong
      TRS <- unique(tbl$TR);
      if (length(TRS) != 210) { stop("too many or few TRs"); }
      FIRSTVOXEL <- which(colnames(tbl) == "v1");  # column number of the first voxel column; bigger columns are voxels.
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
      
      if (DO_ROW_SCALING == TRUE) {  # do row scaling (all voxels in each volume)
        scaled <- tbl[,FIRSTVOXEL:(dim(tbl)[2])];   # get voxel columns; cols 1 to LABELCOLS are labels
        scaled <- scale(t(scaled));  # scale does columns, so transpose
        tbl <- cbind(tbl[,1:(FIRSTVOXEL-1)], t(scaled));  # put label columns back on, transpose back
        rm(scaled);
      }    
      
      # find the start of each trial that's type pair1 or pair2
      inds <- which(tbl$eventType == PAIR1 | tbl$eventType == PAIR2);
      tinds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
      
      for (OFFSET in OFFSETS) {    # OFFSET <- OFFSETS[1]
        colctr <- 3;
        useinds <- inds[tinds] + OFFSET;  # data to classify; offset from time point 0: start of each trial
        #useinds <- useinds[which(useinds > 0)];  # don't want any negative rows
        if (length(which(useinds < 0)) > 0) { stop("yes, have negative rows"); }
        t0tbl <- tbl[useinds,];
        if (OFFSET != 0) {  # need to get labels for which events these go with; offset labels are other things than what we're classifying.
          uselbls <- tbl$eventType[inds[tinds]];     
          if (length(uselbls) != dim(t0tbl)[1]) { stop("very wrong lengths for uselbls"); } else { t0tbl$eventType <- uselbls; }
        }
        
        for (r in seq(from=1, to=length(RUNS), by=STEPBY)) {   # r <- 1;
          testRuns <- RUNS[r];  # runs to go into test set for this cross-validation fold
          for (stepby in 1:(STEPBY-1)) {  # fiddly since not all people have all runs: need to make sure doesn't try to include a run that isn't there.
            if ((r+stepby) < length(RUNS)) { testRuns <- c(testRuns, RUNS[r+stepby]); } 
          }
          #if ((r+1) < length(RUNS)) { testRuns <- c(testRuns, RUNS[r+1]); }  # this is fiddly since not all people have all runs.
          #if ((r+2) < length(RUNS)) { testRuns <- c(testRuns, RUNS[r+2]); }
          testinds <- which(is.element(t0tbl$run, testRuns));  # rows with test-set runs
          traininds <- (1:dim(t0tbl)[1])[-testinds];  # all the others are training-set runs
          if (length(testinds) == 0) { stop("no testing data"); }  # can't be empty
          if (length(traininds) == 0) { stop("no training data"); }
          if (length(intersect(testinds, traininds)) > 0) { stop("overlap in test and traininds"); }
          allTest <- t0tbl[testinds,];  # subset into training and testing sets
          allTrain <- t0tbl[traininds,];
          allTest$eventType <- factor(allTest$eventType);   # gets rid of empty factor levels (simplifies balancing)
          allTrain$eventType <- factor(allTrain$eventType);
          
          for (splt in 1:NUMSPLITS) {    # splt <- 1;
            set.seed(SEEDS[splt]);  # set random seed so splits are repeat-able
            doSplits <- TRUE; # marker to recognize when the examples are balanced, so don't split when unnecessary.
            
            # balance number of pair1 and pair2 entries in the training and testing sets by omitting rows of bigger class
            # check if classes are balanced in the training data, and balance if not
            oneCount <- summary(allTrain$eventType)[[1]];
            twoCount <- summary(allTrain$eventType)[[2]];
            if (twoCount == oneCount) {   # balanced already.
              useTrain <- allTrain;
              doSplits <- FALSE;
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
              # setting splt to 10 should ensure that this person is only run through once, since they're training and testing sets are balanced
              # there's no need to do the splits.
              if (doSplits == FALSE) { splt <- 10; }  
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
            
            rtbl[rowctr,1] <- splt;
            rtbl[rowctr,2] <- OFFSET;
            subID[rowctr] <- SUB;
            rtbl[rowctr,colctr] <- doSVM(useTrain, useTest, DO_DEFAULT_SCALING);  
            rm(useTest, useTrain);  # take these out of memory
            rowctr <- rowctr + 1;
          }
          rowctr <- rowctr - NUMSPLITS;
          colctr <- colctr + 1;
        }
        rowctr <- rowctr + NUMSPLITS;
      }
    }
    avgProp <- apply(rtbl[,3:(MAXTESTS+2)], 1, mean, na.rm=TRUE);  # calculate the average for each row
    colnames(rtbl) <- c("splitNum", "timePoint", paste("runSet", 1:MAXTESTS, "out", sep=""));
    rtbl <- data.frame(subID, rtbl, avgProp);
    rownames(rtbl) <- 1:dim(rtbl)[1];
    if (SCALING_LABEL == "_") { SCALING_LABEL <- "_noScaling"; }
    write.table(rtbl, paste(outpath, SUB, "_", ROI, "_", PAIR1, "_", PAIR2, SCALING_LABEL, ".txt", sep=""));    
  }
}

##################################################################################################################################################
# Prepare Perms
# First, figure out how many of each type are in each run for each person.
# this is a ROI-based analysis, classifying three types of stimuli (upempty, upred, upgreen), pairwise.
# always default scaling. within-subjects classification, each timepoint separately. There are missings and imbalance.
# partition on the runs, but not each run separately (not all stimuli in all runs); rather groups of runs at once.
# only need to permute labels once, not for each pair separately.

# should permute the upempty, upred, upgreen trial-type labels only; keep the time structure intact. probably best to permute labels within
# the partitions (several runs) if possible. Might need a different scheme for every person ...
# NEED TO CHANGE PATHS FOR ALL STEPS TO RUN (currently set up to run from Jo's directories)
##################################################################################################################################################

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

rtbl <- array(NA, c(length(SUBS)*length(PAIRS), 10));   # results table
rctr <- 1;
for (sn in 1:length(SUBS)) {     # sn <- 1;
  tbl <- read.table(gzfile(paste(inpath, "sub", SUBS[sn], "_", ROI, "_meanSub.gz", sep="")), comment.char=""); # read in the data    
  # STEPBY <- 4;
  #if (SUBS[sn] == 1003 | SUBS[sn] == 1009 | SUBS[sn] == 1019) { STEPBY <- 3; } 
  # if (SUBS[sn] == 1019) { STEPBY <- 3; }
  
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
    
    #for (r in seq(from=1, to=length(RUNS), by=STEPBY)) {   # r <- 1;
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
colnames(rtbl) <- c("subID", "item", "train1", "test1", "train2", "test2", "train3", "test3", "train4", "test4");
rtbl
rtbl <- rtbl[,1:4];
write.table(rtbl, "c:/maile/svnFiles/plein/consulting/Alan/setupPerms/countTable.txt");

##################################################################################################################################################
# figure out the counts: what permutations do we need?

rm(list=ls());

tbl <- read.table("c:/maile/svnFiles/plein/consulting/Alan/setupPerms/countTable.txt");
#> tbl
#      subID    item train1 test1
#1  sub 1003 upempty     13    18
#2  sub 1003   upred     31    18
#3  sub 1003 upgreen     19    22
#4  sub 1004 upempty     12     2

ftable(tbl$item, tbl$train1)
#         8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 27 28 29 30 31                                                                        
#upempty  1 2  1  1  3  5  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
#upgreen  0 0  0  1  0  0  3  0  1  1  1  2  3  3  1  0  0  0  0  0  0  0
#upred    0 0  0  1  0  0  0  1  0  1  0  0  0  2  1  1  1  1  2  1  1  3

# probably easiest just to read the values out of tbl.
#           1003    1004     1005       1006        1007         1008     1009         1011     1012          1013      1014     1015
cts <- c(13,19,18, 11,14,9, 18,16,23, 13,20,18,25, 14,21,15,24, 9,14,19, 12,21,17,24, 12,19,11, 13,22,12,17, 14,21,24, 10,20,14, 13,20,12,19,
         #    1016     1017      1018     1019
         13,17,14, 11,16,17, 14,22,23, 8,11,6,9);
cts <- unique(cts);
sort(cts)   # 6  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25

# these are the number of examples of each class in the training set, so smallest is 6 "a" and 6 "b".

##################################################################################################################################################
# enumerate all of the permutations for the case of 6 examples of each class

library(gtools)

# datafile for classification. half the runs into training, half into testing.  13 of each for this person. Ignore run structure during permutations.
#> useTrain[,1:10]
#     eventType   subID run  TR         v1         v2          v3          v4           v5          v6
#160    upempty sub1003   1 160   4.872750 -17.588802 -21.8477757 -11.1878543  -8.25525251  -4.1752880
#178    upempty sub1003   1 178  19.613595  -2.468319 -24.3555882 -17.8549686  -6.42877546  -4.5097607
#196    upempty sub1003   1 196  18.503121  -8.132381 -19.7862523  -8.1713138  -5.24597517  -2.6767529
#229    upgreen sub1003   2  19  14.557575  50.837628  17.8342445  11.0546977   7.35696644  16.6106152
#282    upgreen sub1003   2  72 -18.238080 -14.006488   0.8774574  -0.1830953  -0.38198620  -8.9532276
#317    upempty sub1003   2 107  10.783649  11.108746  -8.0478356   0.6866557   3.78232044   0.8348583
#335    upgreen sub1003   2 125   8.351642  17.503033   2.3585365   2.5915629  -3.77834850   9.0567822
#353    upempty sub1003   2 143   6.681842   8.543316  -6.8422081  -2.9763082   7.32596058  15.9610180
#388    upempty sub1003   2 178  10.084308   7.738995 -11.3744957   5.6793925  -3.60024792   1.0660595
#406    upgreen sub1003   2 196 -12.293377 -10.609149   2.5709388   3.4344584  -0.70101696   3.9908642
#422    upempty sub1003   3   2  34.621426  60.294717   3.3604033   0.6563276   1.68389980  19.8985180
#527    upgreen sub1003   3 107 -14.027500 -27.213950   9.3878691  12.8508466  -7.93359288 -12.0057789
#545    upempty sub1003   3 125 -14.289951 -22.920737   1.7990630   6.0984663  -2.22344884 -13.7485767
#598    upgreen sub1003   3 178 -12.016879 -13.781577  17.7986357  14.0569013   4.41546718  -6.7927662
#649    upempty sub1003   4  19  11.392789 -19.649489 -32.6713640 -12.1860334 -13.08576718 -13.6444731
#667    upgreen sub1003   4  37  -1.230014 -26.352126 -15.2751848  -2.7804548  -7.28053037  -6.6732207
#684    upempty sub1003   4  54  -7.499668 -32.327834  -2.3118059 -10.2841169 -10.67713681  -8.7518340
#702    upgreen sub1003   4  72  -6.497959 -23.920973  -5.8902972  -3.5492536  -7.97498837  -8.2425566
#720    upempty sub1003   4  90 -11.374423 -29.782912 -14.0034563 -13.6246320  -3.96729794  -2.4583159
#1018   upgreen sub1003   5 178  10.652058  36.371021  12.4823338   5.0550177  16.78656645  10.1619963
#1035   upgreen sub1003   5 195  -1.001385  38.072803   5.5977513  11.1405890   7.65418120  12.8994231
#1069   upgreen sub1003   6  19   5.884363  30.548124  -8.9742170  -7.4789475   4.00917039   4.9442970
#1122   upgreen sub1003   6  72   9.731409  16.643949  -6.2416120  -4.9791306  -2.28306594  -5.5044945
#1210   upgreen sub1003   6 160 -11.884558 -26.946261  -6.7518049   2.6082717   0.02369675  -2.1811303
#1228   upempty sub1003   6 178 -11.757239 -29.026217   0.4969134   3.8113357  -3.48435989  -3.2608422
#1245   upempty sub1003   6 195 -14.832800 -46.463961  -2.2309309  -0.1621140  -6.66783157  -2.5804222

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
# confirm randomization repeatable with the same seeds

rm(list=ls());

set.seed(234234);
lbls <- c(rep("a", 4), rep("b", 4))
sample(lbls)    # [1] "a" "b" "a" "b" "a" "b" "a" "b"
sample(lbls)    # [1] "a" "b" "b" "b" "a" "a" "a" "b"
# same thing comes up even if r closed and restarted.


# but not with a different seed:
set.seed(281236);
lbls <- c(rep("a", 4), rep("b", 4))
sample(lbls)    # [1] "a" "b" "b" "a" "a" "a" "b" "b"
sample(lbls)    # [1] "b" "a" "b" "b" "a" "b" "a" "a"

##################################################################################################################################################
# Half-Split Class and Perm
# C L U S T E R 
# Need "permpath" files (from previous step), job files (doBigPerms/jobFiles/), inputfiles from /step2/MeanSub/
# Run by submitting ./startAll_pair1.sh (or "pair2.sh", etc). Runs pair for every subject for individual roi. Must manually change roi
# in script and re-run. Subject 2004 can't be run.
##################################################################################################################################################

library(e1071);  # R interface to libsvm

rm(list=ls()); ONNIL <- FALSE; ONCLUSTER <- TRUE; JOCOMPUTER <- FALSE;
# rm(list=ls()); ONNIL <- FALSE; ONCLUSTER <- FALSE; JOCOMPUTER <- TRUE;

SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");
#ROI <- "PFC_mask_native"     
ROI <- "BG_LR_CaNaPu_native"
#ROI <- "Parietal_mask_native"

OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);

for (SUB in SUBS) {  #SUB <- "sub1005";
  
  if (ONCLUSTER == TRUE) {
    inpath <- "/scratch/aceaser/input/";
    outpath <- "/scratch/aceaser/output/"; 
    permpath <- "/scratch/aceaser/permInputFiles/";   # location of 6eachTable.txt, 8eachTable.txt, etc.
    
    cA <- commandArgs();   
    num  <- as.numeric(cA[5]);  # which timepoints to run
    OFFSET <- OFFSETS[num];    # 11 total
    
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
    inpath <- "c:/maile/svnFiles/plein/consulting/Alan/dataFilesFromServer/"; 
    outpath <- "d:/temp/";
    permpath <- "c:/maile/svnFiles/plein/consulting/Alan/setupPerms/";
    SUB <- "sub1003";
    PAIR1 <- "upred"; PAIR2 <- "upgreen";
  }
  
  
  SEEDS <- c(51591, 36414, 56347, 38442, 20176, 51348, 89727, 67106, 23543, 52663);  # 10 random seeds, from sample(1:100000)[1:10]
  NUMSPLITS <- 1;    # there will be 10, but we're running them one at a time.
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
  tbl <- read.table(gzfile(paste(inpath, SUB, "_", ROI, "_meanSub.gz", sep="")), comment.char=""); # read in the data
  RUNS <- sort(unique(tbl$run));
  if (length(RUNS) > 12) { stop("too many runs"); }   # try to catch if the input data is wrong
  TRS <- unique(tbl$TR);
  if (length(TRS) != 210) { stop("too many or few TRs"); }
  FIRSTVOXEL <- which(colnames(tbl) == "v1");  # column number of the first voxel column; bigger columns are voxels.
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
  
  # find the start of each trial that's type pair1 or pair2
  inds <- which(tbl$eventType == PAIR1 | tbl$eventType == PAIR2);
  tinds <- c(1, (1 + which(diff(inds) > 2)));  # these are the rows that start a trial of the type we want
  
  rtbl <- array(NA, c(NUMSPLITS*MAXPERMS, 7));   # results table
  subID <- rep(NA, NUMSPLITS*MAXPERMS);   # subject ID column of results table
  rowctr <- 1;
  
  useinds <- inds[tinds] + OFFSET;  # data to classify; offset from time point 0: start of each trial
  if (length(which(useinds < 0)) > 0) { stop("yes, have negative rows"); }
  t0tbl <- tbl[useinds,];
  if (OFFSET != 0) {  # need to get labels for which events these go with; offset labels are other things than what we're classifying.
    uselbls <- tbl$eventType[inds[tinds]];     
    if (length(uselbls) != dim(t0tbl)[1]) { stop("very wrong lengths for uselbls"); } else { t0tbl$eventType <- uselbls; }
  }
  
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
outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
set.seed(234984);  # so can reproduce the permutation p-values, hopefully they don't jitter much with different seeds.

OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5); 
SUBS <- paste("sub", c(1003, 1005:1009, 1011:1019), sep="");

# PAIR1 <- "upempty"; PAIR2 <- "upgreen";
# PAIR1 <- "upempty"; PAIR2 <- "upred";
 PAIR1 <- "upgreen"; PAIR2 <- "upred";

ROI <- "BG_LR_CaNaPu_native"    # "PFC_mask_native"
#ROI <- "PFC_mask_native"
#ROI <- "Parietal_mask_native"

#> head(intbl)
#    subID splitNum timePoint firstTest secondTest X1stTrainSize X2ndTrainSize permNum   avgProp
#1 sub1003        1        -5 0.5277778  0.4615385            13            18       0 0.4946581
#2 sub1003        1        -5 0.3055556  0.4230769            13            18       1 0.3643162
#3 sub1003        1        -5 0.5277778  0.5384615            13            18       2 0.5331197

ptbl_all <- array(NA, c(length(OFFSETS),2));  # mean & p-value averaged across included subjects
ptbl_eachSub <- array(NA, c(length(OFFSETS), length(SUBS)));
for (o in 1:length(OFFSETS)) {   # o <- 1;
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
  ptbl_all[o,2] <- length(which(means[2:length(means)] > means[1]))/(length(PERMS)-1); 
}
colnames(ptbl_eachSub) <- SUBS;
colnames(ptbl_all) <- c("allSubsMean", "allSubsP");
cbind(OFFSETS, round(ptbl_all,3), round(ptbl_eachSub,3));  # first column is the 'all' p-value

apply(ptbl_eachSub, 1, mean)

out.tbl <- cbind(OFFSETS, ptbl_all, ptbl_eachSub, apply(ptbl_eachSub, 1, mean))
# either write or set an outpath, write out what the file name should be called.txt
#write.table(out.tbl, some.path)
      write.table(out.tbl, gzfile(paste(outpath, ROI, "_", PAIR1, "_", PAIR2, "_defSc_halfSplt_perms.txt.gz", sep="")));    

######################################################################################################################################################################
