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

















































#
