make an error
##################################################################################################################################################
# this is a ROI-based analysis, classifying three types of stimuli (upempty, upred, upgreen), pairwise.
# always default scaling. within-subjects classification, each timepoint separately. There are missings and imbalance.
# partition on the runs, but not each run separately (not all stimuli in all runs); rather groups of runs at once.
# only need to permute labels once, not for each pair separately.

# should permute the upempty, upred, upgreen trial-type labels only; keep the time structure intact. probably best to permute labels within
# the partitions (several runs) if possible. Might need a different scheme for every person ...
##################################################################################################################################################
# first, figure out how many of each type are in each run for each person.

rm(list=ls());

inpath <- "c:/maile/svnFiles/plein/consulting/Alan/dataFilesFromServer/"; 
#inpath <- "/data/nil-external/ccp/ALAN_CU/FORMVPA/step2/MeanSub/"
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
##################################################################################################################################################
# permutation labels for the leave-one-subject-out analyses. 13 subjects, two examples per person (one of each class). Flip the labels pairwise
# sub1 flipped, sub 1 & 3, etc. to keep balance. True labeling is sub1-class1 sub2-class1 ... sub1-class2 sub2class2 ...

# how many possible?   choose(13,1) + choose(13,2) + choose(13,3) + choose(13,4) + choose(13,5) + choose(13,6) + choose(13,7) + choose(13,8) + choose(13,9) +
# choose(13,10) + choose(13,11) + choose(13,12)   # 8190

rm(list=ls());

set.seed(3421);
num.subs <- 13;
num.perms <- 1000;

hasMatch <- function(thislbl, done.lbls) {    # confirm no repeated rows
  ans <- FALSE;   
  
  if (done.lbls > 0) {
    for (j in 1:done.lbls) {   # j <- 2;
      if (length(which(thislbl == out.lbls[j,])) == (num.subs*2)) { ans <- TRUE; break; }
    }
  }
    
  return(ans)
}

makeNew <- function() {
  num.to.flip <- sample(1:(num.subs))[1];  # can't flip all 13: true labeling
  a.s <- rep("a", num.subs)
  b.s <- rep("b", num.subs)
  for (i in 1:num.to.flip) {   # i <- 1;
    flip <- sample(1:(num.subs))[1]; 
    a.s[flip] <- "b";
    b.s[flip] <- "a";
  }
  
  return(c(a.s, b.s));
}

out.lbls <- array(NA, c((num.perms+1), num.subs*2))
out.lbls[1,] <- c(rep("a", num.subs), rep("b", num.subs));  # true labeling first
for (p in 1:num.perms) {  # p <- 1;
  this.lbl <- makeNew();
  while(hasMatch(this.lbl, (p-1)) == TRUE) { this.lbl <- makeNew(); }   # add if new; will loop forever.
  out.lbls[(p+1),] <- this.lbl;
}

write.table(out.lbls, paste("c:/maile/svnFiles/plein/consulting/Alan/setupPerms/l1outPermlbls.txt", sep=""));

##################################################################################################################################################
##################################################################################################################################################
# permutation labels for the leave-one-subject-out analyses. 13 subjects, two examples per person (one of each class). Flip the labels pairwise
# sub1 flipped, sub 1 & 3, etc. to keep balance. True labeling is sub1-class1 sub2-class1 ... sub1-class2 sub2class2 ...
# COMPLETE set of permutations

# how many possible?   choose(13,1) + choose(13,2) + choose(13,3) + choose(13,4) + choose(13,5) + choose(13,6) + choose(13,7) + choose(13,8) + choose(13,9) +
# choose(13,10) + choose(13,11) + choose(13,12)   # 8190

library(gtools)
rm(list=ls());

num.subs <- 13;
all.perms <- array(NA, c(8191, num.subs*2));  # 13 subjects * 2 classes
all.perms[1,] <- c(rep("a", num.subs), rep("b", num.subs));  # true labeling first
ctr <- 2;

# flip 1
tmp <- combinations(num.subs, 1, 1:num.subs);  # which people to flip
for (i in 1:nrow(tmp)) {   # i <- 1;
  flip <- tmp[i]
  a.s <- rep("a", num.subs)
  a.s[flip] <- "b";
  b.s <- rep("b", num.subs)
  b.s[flip] <- "a";
  all.perms[ctr,] <- c(a.s, b.s)
  ctr <- ctr + 1;
}

# flip 2 to 11
for (nf in 2:11) {  # nf <- 11;
  tmp <- combinations(num.subs, nf, 1:num.subs);  # which people to flip
  for (i in 1:nrow(tmp)) {   # i <- 1;
    flip <- tmp[i,]
    a.s <- rep("a", num.subs)
    a.s[flip] <- "b";
    b.s <- rep("b", num.subs)
    b.s[flip] <- "a";
    all.perms[ctr,] <- c(a.s, b.s)
    ctr <- ctr + 1;
  }
}

# flip 12
tmp <- combinations(num.subs, 12, 1:num.subs);  # which people to flip
for (i in 1:nrow(tmp)) {   # i <- 1;
  flip <- tmp[i,]
  a.s <- rep("a", num.subs)
  a.s[flip] <- "b";
  b.s <- rep("b", num.subs)
  b.s[flip] <- "a";
  all.perms[ctr,] <- c(a.s, b.s)
  ctr <- ctr + 1;
}

write.table(all.perms, paste("c:/maile/svnFiles/plein/consulting/Alan/setupPerms/l1outPermlbls_complete.txt", sep=""));

# double-check for repeats (shouldn't be any!)
for (i in 1:(nrow(all.perms)-1)) { 
  for (j in (i+1):nrow(all.perms)) {   # i <- 1; j <- 2;
    if (length(which(all.perms[i,] == all.perms[j,])) == 26) { stop("found a match!!!"); } 
  }
}

##################################################################################################################################################






















































#
