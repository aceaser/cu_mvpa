make an error
#############################################################################################################################################
# 'prettier' pictures of Alan's MVPA results
# Jo Etzel, 11 July 2013
#############################################################################################################################################
# get the null distribution for each pair and offset; within-subjects

rm(list=ls()); 

in.path <- "d:/temp/Alan/output_fromServer/";
out.path <- "d:/temp/Alan/forplotting/"; 

seeds <- c(5307, 6716, 3170, 7532, 2169, 9275, 1931, 9598, 3288, 6123, 3329);    # random seeds for reproducibility. from sample(1:10000)[1:11]
OFFSETS <- c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5); 
SUBS <- paste("sub", c(1003, 1005:1009, 1011:1019), sep="");

# PAIR <- "upempty_upgreen";
# PAIR <- "upempty_upred";
#PAIR <- "upgreen_upred";

#ROI <- "BG_LR_CaNaPu_native"    # "PFC_mask_native"
#ROI <- "PFC_mask_native"
ROI <- "Parietal_mask_native"

#> head(intbl)
#    subID splitNum timePoint firstTest secondTest X1stTrainSize X2ndTrainSize permNum   avgProp
#1 sub1003        1        -5 0.5277778  0.4615385            13            18       0 0.4946581
#2 sub1003        1        -5 0.3055556  0.4230769            13            18       1 0.3643162
#3 sub1003        1        -5 0.5277778  0.5384615            13            18       2 0.5331197

for (PAIR in c("upempty_upgreen", "upempty_upred", "upgreen_upred")) {

ptbl_all <- array(NA, c(length(OFFSETS),2));  # mean & p-value averaged across included subjects
colnames(ptbl_all) <- c("allSubsMean", "allSubsP");
ptbl_eachSub <- array(NA, c(length(OFFSETS), length(SUBS)));
colnames(ptbl_eachSub) <- SUBS;
null.tbl <- array(NA, c(1001, length(OFFSETS)));
colnames(null.tbl) <- paste("offset", OFFSETS, sep="");

for (o in 1:length(OFFSETS)) {   # o <- 1;
  set.seed(seeds[o]);  # so can reproduce the permutation p-values, hopefully they don't jitter much with different seeds.
  temptbl <- array(NA, c(1001, length(SUBS)));  # put all the subjects' permutation accuracies into here
  for (s in 1:length(SUBS)) { # s <- 1;
    fname <- paste(in.path, SUBS[s], "_", OFFSETS[o], "_", ROI, "_", PAIR, "_defSc_halfSplt_perms.txt.gz", sep="")
    if (file.exists(fname)) {
      intbl <- read.table(gzfile(fname));
      # intbl has accuracies for each of the splits. average so one average accuracy per permNum.
      PERMS <- unique(intbl$permNum);
      PERMS <- PERMS[which(PERMS != 0)];     # permNum == 0 is the real-labeled accuracy
      if (max(PERMS) > 1000) { PERMS <- sample(PERMS)[1:1000]; } # a few people have more than 1000 perms; subset in this case.
      means <- rep(NA, length(PERMS)+1);  # to hold the across-splits accuracies (one per PERM)
      # real labeling first
      stbl <- subset(intbl, intbl$permNum == 0);
      if (nrow(stbl) != 10) { stop("(dim(stbl)[1] != 10)"); }
      means[1] <- mean(stbl$avgProp);
      for (p in 1:length(PERMS)) {   # p <- 1;
        stbl <- subset(intbl, intbl$permNum == PERMS[p])
        if (nrow(stbl) != 10) { stop("(dim(stbl)[1] != 10)"); }   # always 10 splits
        means[p+1] <- mean(stbl$avgProp);
      }
      temptbl[1:(length(PERMS)+1),s] <- means;
      
      # calculate the p-value for this person
      ptbl_eachSub[o,s] <- length(which(means[2:length(means)] > means[1]))/length(means); 
    } else { print(paste("missing file:", fname)); }
  }
  
  # bootstrap the people with missings so that they have 1000 as well. this was recommended in Good (resampling methods pdf book), pg 105.
  for (s in 1:length(SUBS)) { # s <- 6;
    numMissings <- length(which(is.na(temptbl[,s])))
    if (numMissings > 0) {  # have a missing, so bootstrap.
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
  means <- apply(temptbl, 1, mean);  # across-subjects means for each of the 1001 rows
  if (length(means) != 1001 | length(which(is.na(means))) > 0) { stop("wrong means"); }
  ptbl_all[o,1] <- means[1]; 
  ptbl_all[o,2] <- length(which(means[2:length(means)] > means[1]))/length(means); 
  
  null.tbl[,o] <- means;  # save the null distribution for plotting
}
cbind(OFFSETS, round(ptbl_all,3), round(ptbl_eachSub,3));  # first column is the 'all' p-value

write.table(null.tbl, paste(out.path, PAIR, "_", ROI, "nullDists.txt", sep=""));
}


##################################################################################################################################################
# within-subjects graphs: shading for significance
# nice color charts at http://gucky.uni-muenster.de/cgi-bin/rgbtab-en

rm(list=ls());

path <- "d:/temp/Alan/forplotting/";    # these files are on the server: /data/nil-external/ccp/ALAN_CU/FORMVPA/classify/permResults_plotting/
pair.ids <- c("upempty_upgreen", "upgreen_upred", "upempty_upred");
pair.ttls <- c("Upempty vs. Upgreen", "Upgreen vs. Upred", "Upred vs. Upempty");  # plot titles; same order as pair.ids
roi.ids <- c("BG_LR_CaNaPu_native", "PFC_mask_native", "Parietal_mask_native");
roi.ttls <- c("Striatum", "PFC", "Parietal");  # plot titles; same order as roi.ids
timepoints <- seq(from=2,to=22,by=2)

clr.center <- "steelblue";
clr.mid <- "lightsteelblue";
clr.edge <- "lightblue";

windows(10,8);  # specify the plot window size so everything will look like it should.
layout(matrix(1:9,c(3,3),byrow=TRUE));
par(mar=c(2, 3, 1, 1.5), mgp=c(1.5, 0.5, 0)) 
# mar: c(bottom, left, top, right)’ gives the number of lines of margin to be specified on the four sides of the plot. default is ‘c(5, 4, 4, 2) + 0.1’.
# mgp: margin line (in ‘mex’ units) for the axis title, axis labels and axis line. The default is ‘c(3, 1, 0)’.

for (r in 1:length(roi.ids)) {
  for (p in 1:length(pair.ids)) {   # p <- 1; r <- 1;
    tbl <- read.table(paste(path, pair.ids[p], "_", roi.ids[r], "nullDists.txt", sep=""));
    # if (r == 1) { ttl <- pair.ttls[p]; } else { ttl <- ""; }    # only put main title in the first row
    # if (p == 1) { yttl <- roi.ttls[r]; } else { yttl <- ""; }   # and y-axis title in first column
    # if (r == 3) { xttl <- "TIME (sec)"; } else { xttl <- ""; }   # and x-axis title in last row
    ttl <- ""; yttl <- ""; xttl <- "";
    
    # start with a blank plot
    plot(x=timepoints, y=rep(0,11), col='white', ylab=yttl, xlab=xttl, main=ttl, xlim=c(2,22), ylim=c(0.43,0.72), las=1, cex.axis=1.3, cex.main=1.5, cex.lab=1.3)
    for (yval in c(0.4, 0.45, 0.55, 0.6, 0.65, 0.7, 0.75)) { lines(x=c(-1,25), y=rep(yval,2), lty='dotted', col='darkgrey'); }
    
    # add shading to show the boundaries of the null distribution (and so the p-values)
    for (tp in 1:10) {    # tp <- 1;
      polygon(c(timepoints[tp],timepoints[tp],timepoints[tp+1],timepoints[tp+1]), 
        c(max(tbl[2:1001,tp]),min(tbl[2:1001,tp]),min(tbl[2:1001,(tp+1)]),max(tbl[2:1001,(tp+1)])), col=clr.edge, lty='blank')

      bounds1 <- quantile(tbl[2:1001,tp], probs=c(0.01, 0.99));
      bounds2 <- quantile(tbl[2:1001,(tp+1)], probs=c(0.01, 0.99));
      polygon(c(timepoints[tp],timepoints[tp],timepoints[tp+1],timepoints[tp+1]), 
        c(bounds1[[1]],bounds1[[2]], bounds2[[2]],bounds2[[1]]), col=clr.mid, lty='blank')
        
      bounds1 <- quantile(tbl[2:1001,tp], probs=c(0.05, 0.95));
      bounds2 <- quantile(tbl[2:1001,(tp+1)], probs=c(0.05, 0.95));
      polygon(c(timepoints[tp],timepoints[tp],timepoints[tp+1],timepoints[tp+1]), 
        c(bounds1[[1]],bounds1[[2]], bounds2[[2]],bounds2[[1]]), col=clr.center, lty='blank')
    }
    
    # horizontal and vertical lines
    lines(x=c(-1,25), y=rep(0.5,2), col='darkgrey');   # horizontal line at chance
    lines(x=rep(10,2), y=c(0,1), col='darkgrey', lwd=2);   # vertical lines for timepoint boundaries
    lines(x=rep(13,2), y=c(0,1), col='darkgrey', lwd=2);
    lines(x=rep(20,2), y=c(0,1), col='darkgrey', lwd=2);

    # put on significance stars.
    for (tp in 1:11) {    # tp <- 11;
      star.str <- " ";
      p.val <- (1 + length(which(tbl[2:1001,tp] > tbl[1,tp])))/1001;
      if (p.val < 0.05) { star.str <- "*"; }
      if (p.val < 0.01) { star.str <- "**"; }
      text(x=timepoints[tp], y=(tbl[1,tp] + 0.02), label=star.str, cex=2);
    }
 
    lines(x=seq(from=2,to=22,by=2), y=tbl[1,], lwd=3, col='orangered1');   # classification accuracy line
  }
}

#############################################################################################################################################
##################################################################################################################################################
# within-subjects graphs: lines for significance
# nice color charts at http://gucky.uni-muenster.de/cgi-bin/rgbtab-en

rm(list=ls());

path <- "d:/temp/Alan/forplotting/"; 
pair.ids <- c("upempty_upgreen", "upgreen_upred", "upempty_upred");
pair.ttls <- c("Upempty vs. Upgreen", "Upgreen vs. Upred", "Upred vs. Upempty");  # plot titles; same order as pair.ids
roi.ids <- c("BG_LR_CaNaPu_native", "PFC_mask_native", "Parietal_mask_native");
roi.ttls <- c("Striatum", "PFC", "Parietal");  # plot titles; same order as roi.ids
timepoints <- seq(from=2,to=22,by=2)

y.lim <- c(0.4,0.65)

windows(10,8);  # specify the plot window size so everything will look like it should.
layout(matrix(1:9,c(3,3),byrow=TRUE));
par(mar=c(2.5, 3.5, 2, 1.5), mgp=c(2, 0.5, 0)) 
# 'mgp The margin line (in ‘mex’ units) for the axis title, axis labels and axis line. The default is ‘c(3, 1, 0)’.
# mar: c(bottom, left, top, right)’ gives the number of lines of margin to be specified on the four sides of the plot. default is ‘c(5, 4, 4, 2) + 0.1’.

for (r in 1:length(roi.ids)) {
  for (p in 1:length(pair.ids)) {   # p <- 1; r <- 1;
    tbl <- read.table(paste(path, pair.ids[p], "_", roi.ids[r], "nullDists.txt", sep=""));
    # head(tbl)
    #    offset.5  offset.4  offset.3  offset.2  offset.1   offset0   offset1   offset2   offset3   offset4   offset5
    #   0.5181228 0.4965174 0.5019507 0.4951193 0.4931519 0.4925066 0.5278038 0.5088556 0.5410935 0.5771908 0.5664269
    #   0.4772542 0.4761773 0.4924951 0.4963407 0.4883575 0.5073966 0.5094105 0.4874892 0.4921644 0.4807036 0.4946070
    
    if (r == 1) { ttl <- pair.ttls[p]; } else { ttl <- ""; }    # only put main title in the first row
    if (p == 1) { yttl <- paste("Classification Accuracy:", roi.ids[r]); } else { yttl <- ""; }   # and y-axis title in first column
    if (r == 3) { xttl <- "TIME (sec)"; } else { xttl <- ""; }   # and x-axis title in last row
    
    plot(x=timepoints, y=rep(0,11), col='white', ylab=yttl, xlab=xttl, main=ttl, ylim=y.lim, cex.axis=1.1, cex.main=1.5, cex.lab=1.5)
   # lines(x=c(-1,25), y=rep(0.45,2), col='grey', lty='dashed');
    lines(x=c(-1,25), y=rep(0.5,2), col='grey', lwd=2);   # horizontal line at chance, lty='dashed'
   # lines(x=c(-1,25), y=rep(0.55,2), col='grey', lty='dashed');
   # lines(x=c(-1,25), y=rep(0.6,2), col='grey', lty='dashed');
   # lines(x=c(-1,25), y=rep(0.65,2), col='grey', lty='dashed');
    
    lines(x=rep(10,2), y=c(0,1), col='darkgrey', lwd=2);   # vertical lines for timepoint boundaries
    lines(x=rep(13,2), y=c(0,1), col='darkgrey', lwd=2);
    lines(x=rep(20,2), y=c(0,1), col='darkgrey', lwd=2);

    # lines showing the boundaries of the null distribution (and so the p-values)
    for (tp in 1:10) {    # tp <- 1;
      lines(x=c(timepoints[tp], timepoints[tp+1]), y=c(max(tbl[2:1001,tp]), max(tbl[2:1001,(tp+1)])), col='steelblue')
      lines(x=c(timepoints[tp], timepoints[tp+1]), y=c(min(tbl[2:1001,tp]), min(tbl[2:1001,(tp+1)])), col='steelblue')

      bounds1 <- quantile(tbl[2:1001,tp], probs=c(0.01, 0.99));
      bounds2 <- quantile(tbl[2:1001,(tp+1)], probs=c(0.01, 0.99));
      lines(x=c(timepoints[tp], timepoints[tp+1]), y=c(bounds1[[1]], bounds2[[1]]), col='lightsteelblue')
      lines(x=c(timepoints[tp], timepoints[tp+1]), y=c(bounds1[[2]], bounds2[[2]]), col='lightsteelblue')
      
      bounds1 <- quantile(tbl[2:1001,tp], probs=c(0.05, 0.95));
      bounds2 <- quantile(tbl[2:1001,(tp+1)], probs=c(0.05, 0.95));
      lines(x=c(timepoints[tp], timepoints[tp+1]), y=c(bounds1[[1]], bounds2[[1]]), col='lightblue')
      lines(x=c(timepoints[tp], timepoints[tp+1]), y=c(bounds1[[2]], bounds2[[2]]), col='lightblue')      
    }
    
    lines(x=seq(from=2,to=22,by=2), y=tbl[1,], lwd=3, col='orangered1');   # classification accuracy line
  }
}

#############################################################################################################################################





















































#
