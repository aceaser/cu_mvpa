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
<<<<<<< HEAD
	# intbl2 <-   then intbl <- rbind(intbl, intbl2)  # then rm(intbl2); to clean up
=======
  # intbl2 <-   then intbl <- rbind(intbl, intbl2)  # then rm(intbl2); to clean up
>>>>>>> cb7aff6f45a9705a73849e254f8915362362c947
   
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
# leave-one-person-out cross-validation
##################################################################################################################################################

library(e1071);  # R interface to libsvm (only needs calling once per R session)

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") {
  #inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/";
  inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/avg_allSubs/";
  outpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/l1subOut_result/"; 
}
if (where.run == "Jo") {
  inpath <- "d:/temp/";
  outpath <- "d:/temp/"; 
}
ROIS <- c("BG_LR_CaNaPu_native", "PFC_mask_native","Parietal_mask_native");
#ROIS <- c("Parietal_mask_native");
#SUBS <- paste("sub", c(1005:1009, 1011:1018), sep="")
SUBS <- paste("sub", c(1003:1009, 1011:1019), sep="");   # SUBS <- c(1003:1009, 1011:1019);   # SUBS <- "sub1003";
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify

# the things to classify. will classify the first entry of PAIR1S with the first entry of PAIR2S, etc.
PAIR1S <- c("upgreenCOR","upemptyCOR","upredCOR");  
PAIR2S <- c("upgreenINCOR","upemptyINCOR","upredINCOR");  
#PAIR1S <- c("updateopCOR")
#PAIR2S <- c("updateopINCOR")
#PAIR1S <- c("upgreen","upred","upred");  
#PAIR2S <- c("upempty","upempty","upgreen");  


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
##### C L U S T E R #####
# Alan Ceasar classification: each time point separately, all people together.
# COMPLETE PERMUTATION (8000) version.
# Need files "startAll_bgjobs.sh" for BG roi and "startAll_pfcjobs.sh" for PFC roi, as well
# as "l1outPermlbls_complete.txt" (randomized perm labels).
# Files from previous step need to be on the cluster, as well as the above mentioned files.  
# Execute these jobs on the cluster using "./startAll_pfcjobs.sh". Calls all job files
# created that exist in the "jobs" folder. Need a script to create these files for new 
# study.
##################################################################################################################################################

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
   if (num == 1) { ROI <- "PFC_mask_native"; }
   if (num == 2) { ROI <- "BG_LR_CaNaPu_native"; }
   
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
# combine the multiple output files made by the complete permutation test (complete_l1outPerms.R).
#gzfile(paste(outpath, "l1subOut_", OFFSETS[o], "_", ROI, "_", PAIR1, "_", PAIR2, SCALING_LABEL, "_comPerms_", perm.key, ".txt.gz", sep="")));
# this won't overwrite an existing output file; delete it manually.
##################################################################################################################################################

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") { 
  in.path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/AllPerm/perm_multi_output/"; 
  out.path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/AllPerm/";
}
if (where.run == "Jo") { 
  in.path <- "d:/gitFiles/cu_mvpa/"; 
  out.path <- "d:/temp/";
}

for (do.offset in c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) {
 # for (ROI in c("BG_LR_CaNaPu_native", "PFC_mask_native","Parietal_mask_native")) {
  for (ROI in c("Parietal_mask_native")) { 
   for (pair in c("upgreenCOR_upgreenINCOR", "upredCOR_upredINCOR", "upemptyCOR_upemptyINCOR")) { 
      # do.offset <- -5; ROI <- "BG_LR_CaNaPu_native"; pair <- "upredCOR_upredINCOR";
      # check to see if all the permutation files for this ROI, pair, and offset are done before trying to combine them
      cts <- 0;  # use this as a counter to see how many permutation files are done; need 9.
      for (perm.key in letters[1:9]) {
        fname <- paste(in.path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_DefaultSc_comPerms_", perm.key, ".txt.gz", sep="")
        if (file.exists(fname)) { cts <- cts + 1; }
      }
      
      out.fname <- paste(out.path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_DefaultSc_comPerms.txt.gz", sep="")
      if (cts == 9) {   # got them all, so make the combined version if the output file isn't already there.
        out.tbl <- data.frame(array(NA, c(8191,17)));  # hard-coded the dimensions (and all over in this bit of code)
        for (i in 1:9) {    # i <- 1;
          tbl <- read.table(gzfile(paste(in.path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_DefaultSc_comPerms_", letters[i], ".txt.gz", sep="")));
          if (i == 1) { do.perms <- 1:1000; last.perm <- 1000; }   # copied this from the code in complete_l1outPerms.R; specifies which job is which perms
          if (i == 2) { do.perms <- 1001:2000; last.perm <- 1000; }
          if (i == 3) { do.perms <- 2001:3000; }
          if (i == 4) { do.perms <- 3001:4000; }
          if (i == 5) { do.perms <- 4001:5000; }
          if (i == 6) { do.perms <- 5001:6000; }
          if (i == 7) { do.perms <- 6001:7000; }
          if (i == 8) { do.perms <- 7001:8000; }
          if (i == 9) { do.perms <- 8001:8191; last.perm <- 191; }

	  if (length(do.perms) != last.perm) { stop("wrong lengths"); }
          out.tbl[do.perms,] <- tbl[1:last.perm,];
          # move the rows we need over to the all-together file. code is stupid because the output files
# were made always starting at the first row, but have 8191 rows.
        }
	
	# confirm no NA in out.tbl
	if (length(which(is.na(out.tbl))) > 0) { stop("NAs???"); }
	if (max(diff(out.tbl[,1])) != 1) { stop("missed one??"); }
	if (min(diff(out.tbl[,1])) != 1) { stop("missed one??"); }
	colnames(out.tbl) <- colnames(tbl);
        write.table(out.tbl,gzfile(out.fname));
      }
    }
  }
}

##################################################################################################################################################
# graphing and p-values
##################################################################################################################################################

rm(list=ls());

where.run <- "Alan";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") { path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/AllPerm/"; }
if (where.run == "Jo") { path <- "d:/gitFiles/cu_mvpa/"; }
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
#ROI <- "BG_LR_CaNaPu_native";  
#ROI <- "PFC_mask_native"; 
ROI <- "Parietal_mask_native";

scaling <- "DefaultSc";   # scaling <- "noScaling";
#pair <- "upgreenCOR_upgreenINCOR"
#pair <- "upredCOR_upredINCOR"
pair <- "upemptyCOR_upemptyINCOR"
yttl <- "frequency, in 1000 permutations";
xttl <- "accuracy of the permuted-label datasets"
    
layout(matrix(1:12, c(3,4), byrow=TRUE));
for (do.offset in OFFSETS) {  # do.offset <- 5;
  fname <- paste(path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_", scaling, "_comPerms.txt.gz", sep="")
  if (file.exists(fname)) {
    in.tbl <- read.table(fname)
  
    if (var(in.tbl$avg.acc) == 0) {   # recalculate the avg.acc column - it was messed up in the first file
      if (do.offset != -5 | scaling != "DefaultSc") { stop("no variance in avg.acc?"); } 
      for (i in 1:nrow(in.tbl)) { in.tbl$avg.acc[i] <- mean(as.vector(in.tbl[i,4:16], mode="numeric")); }
    }
  
    p.value <- (1+length(which(in.tbl$avg.acc[2:nrow(in.tbl)] > in.tbl$avg.acc[1])))/nrow(in.tbl)
    
    ttl <- paste("offset= ", do.offset, " p=", round(p.value, 3), "\n", pair, sep="");
    hist(in.tbl$avg.acc[2:nrow(in.tbl)], xlim=c(0,1), breaks=seq(from=0, to=1, by=0.05), xlab=xttl, main=ttl, ylab=yttl, col='cornsilk')
    lines(x=c(0.5,0.5), y=c(-10,9000), col='grey', lwd=2);  # chance
    lines(x=rep(in.tbl$avg.acc[1], 2), y=c(-10,9000), col='salmon', lwd=2);  # true-labeled accuracy
  }
}
#write.table(in.tbl, gzfile(paste(outpath, ROI, "_", pair, "_defSc_l1Out_perms.txt.gz", sep="")));    
