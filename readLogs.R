##################################################################################################################################################################
# read the log files

rm(list=ls())

inpath <- "C:/maile/svnFiles/plein/consulting/Alan/";
outpath <- inpath;

SUB <- 1003;   # just for Jo; put SUB into a loop

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

##################################################################################################################################################################
