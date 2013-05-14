make an error
##################################################################################################################################################
# 7 May 2013. code for figuring out the leave-one-person-out cross-validation permutation testing results
##################################################################################################################################################
# graphing and p-values

rm(list=ls());

where.run <- "Jo";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") { path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; }
if (where.run == "Jo") { path <- "d:/gitFiles/cu_mvpa/"; }
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5);  # offset in TR from the trial starts that we will classify
ROI <- "BG_LR_CaNaPu_native";  # ROI <- "PFC_mask_native"; # 

scaling <- "DefaultSc";   # scaling <- "noScaling";
pair <- "upgreenCOR_upgreenINCOR"
yttl <- "frequency, in 1000 permutations";
xttl <- "accuracy of the permuted-label datasets"
    
layout(matrix(1:12, c(3,4), byrow=TRUE));
for (do.offset in OFFSETS) {  # do.offset <- -5;
  fname <- paste(path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_", scaling, "_perms.txt", sep="")
  if (file.exists(fname)) {
    in.tbl <- read.table(fname)
    if (nrow(in.tbl) != 1001 | ncol(in.tbl) != 17) { stop("not expected dims"); }
  
    if (var(in.tbl$avg.acc) == 0) {   # recalculate the avg.acc column - it was messed up in the first file
      if (do.offset != -5 | scaling != "DefaultSc") { stop("no variance in avg.acc?"); } 
      for (i in 1:nrow(in.tbl)) { in.tbl$avg.acc[i] <- mean(as.vector(in.tbl[i,4:16], mode="numeric")); }
    }
  
    p.value <- (1+length(which(in.tbl$avg.acc[2:1001] > in.tbl$avg.acc[1])))/1001
    
    ttl <- paste("offset= ", do.offset, " p=", round(p.value, 3), "\n", pair, sep="");
    hist(in.tbl$avg.acc[2:1001], xlim=c(0,1), breaks=seq(from=0, to=1, by=0.05), xlab=xttl, main=ttl, ylab=yttl, col='cornsilk')
    lines(x=c(0.5,0.5), y=c(-10,1000), col='grey', lwd=2);  # chance
    lines(x=rep(in.tbl$avg.acc[1], 2), y=c(-10,1000), col='salmon', lwd=2);  # true-labeled accuracy
  }
}

##################################################################################################################################################
# combine the multiple output files made by the complete permutation test (complete_l1outPerms.R).
#gzfile(paste(outpath, "l1subOut_", OFFSETS[o], "_", ROI, "_", PAIR1, "_", PAIR2, SCALING_LABEL, "_comPerms_", perm.key, ".txt.gz", sep="")));
# this won't overwrite an existing output file; delete it manually.

rm(list=ls());

where.run <- "Jo";   # where.run <- "Alan";  # set the paths according to which computer this is being run on

if (where.run == "Alan") { 
  in.path <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/classify/"; 
  out.path <- whatever
}
if (where.run == "Jo") { 
  in.path <- "d:/gitFiles/cu_mvpa/"; 
  out.path <- "d:/temp/";
}

for (do.offset in c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)) {
  for (ROI in c("BG_LR_CaNaPu_native", "PFC_mask_native")) {
    for (pair in c("upgreenCOR_upgreenINCOR", "upredCOR_upredINCOR", "upemptyCOR_upemptyINCOR")) { 
      # do.offset <- -5; ROI <- "BG_LR_CaNaPu_native"; pair <- "upredCOR_upredINCOR";
      # check to see if all the permutation files for this ROI, pair, and offset are done before trying to combine them
      cts <- 0;  # use this as a counter to see how many permutation files are done; need 9.
      for (perm.key in letters[1:9]) {
        fname <- paste(in.path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_DefaultSc_comPerms_", perm.key, ".txt.gz", sep="")
        if (file.exists(gzfile(fname))) { cts <- cts + 1; }
      }
      
      out.fname <- paste(out.path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_DefaultSc_comPerms.txt.gz", sep="")
      if (cts == 9 & !file.exists(out.fname)) {   # got them all, so make the combined version if the output file isn't already there.
        out.tbl <- data.frame(array(NA, c(8191,17)));  # hard-coded the dimensions (and all over in this bit of code)
        for (i in 1:9) {
          tbl <- read.table(gzfile(paste(in.path, "l1subOut_", do.offset, "_", ROI, "_", pair, "_DefaultSc_comPerms_", letters[i], ".txt.gz", sep="")));
          if (i == 1) { do.perms <- 1:1000; }   # copied this from the code in complete_l1outPerms.R; specifies which job is which perms
          if (i == 2) { do.perms <- 1001:2000; }
          if (i == 3) { do.perms <- 2001:3000; }
          if (i == 4) { do.perms <- 3001:4000; }
          if (i == 5) { do.perms <- 4001:5000; }
          if (i == 6) { do.perms <- 5001:6000; }
          if (i == 7) { do.perms <- 6001:7000; }
          if (i == 8) { do.perms <- 7001:8000; }
          if (i == 9) { do.perms <- 8001:8191; }
          out.tbl[do.perms,] <- tbl[do.perms,];   # move the rows we need over to the all-together file
        }
        write.table(gzfile(out.fname));
      }
    }
  }
}
   
##################################################################################################################################################

















































#
