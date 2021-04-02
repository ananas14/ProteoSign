options(warn=1)
# DEPLOYMENT VERSION
# ======================================
# WARNING: Make sure to install the following packages as administrator/root (for them to be available to all users)
# ======================================
DEBUG <- FALSE;
# # latest limma version might not be compatible!! limma 3.18.13 suits best for PS
# if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
# BiocManager::install("limma")
# if(!require("statmod")){ install.packages("statmod") }
# if(!require("ggplot2")){ install.packages("ggplot2", repos="http://cran.fhcrc.org") }
# if(!require("stringr")){ install.packages("stringr", repos="http://cran.fhcrc.org") }
# if(!require("reshape")){ install.packages("reshape", repos="http://cran.fhcrc.org") }
# if(!require("plyr")){ install.packages("plyr", repos="http://cran.fhcrc.org") }
# if(!require("gtools")){ install.packages("gtools", repos="http://cran.fhcrc.org") }
# if(!require("gtools")){ install.packages("labeling", repos="http://cran.fhcrc.org") }
# if(!require("data.table")){ install.packages("data.table", repos="http://cran.fhcrc.org") }
# if(!require("data.table")){ install.packages("data.table", repos="http://cran.fhcrc.org") }
# if(!require("outliers")){ install.packages("outliers", repos="http://cran.fhcrc.org") }
# if(!require("pryr")){ install.packages("pryr", repos="http://cran.fhcrc.org") }
# if(!require("gprofiler2")){  install.packages("gprofiler2") }


library(limma)
library(statmod)
library(stringr)
library(reshape)
library(plyr)
library(ggplot2)
library(labeling)
library(gtools)
library(data.table)
library(outliers)
library(pryr)
library(gprofiler2)
#library(devtools)
#library(lineprof)


# DEBUGGING log flag/level (0 translates to no debugging log at all)
debuglog <- 10
# DEBUGGING log indentation level (psoitive int, global) for levellog function
loglvl <- 0
lastSysTime <- NA
firstSysTime <- NA
#  Print msg with specific indentation (loglvl*3)
#  - change [int]: change previous indentation level by <change> amount
#  - reset [bool]: reset indentation level to 0 (no indentation)
#  - after [bool]: change indentation after printing msg
levellog <- function (msg, change=0, reset=F,after=F, supression=debuglog){
  currSysTime<-Sys.time()
  if(is.na(firstSysTime)){
    firstSysTime <<- currSysTime
  }
  if(is.na(lastSysTime)){
    lastSysTime <<- currSysTime
    td <- 0
  }else{
    td <- difftime(currSysTime, lastSysTime, unit='secs')
    if(td < 0.1){
      td <- 0
    }
    lastSysTime <<- currSysTime
  }
  if(td > 0){
    tdstr <- sprintf('+%.1fs, ', td)
  }else{
    tdstr <- ''
  }
  tdt <- difftime(currSysTime, firstSysTime, unit='secs')
  if(tdt < 0.1){
    tdt <- 0
  }
  if(tdt > 0){
    tdstr<-paste0(tdstr, sprintf('%.1fs elapsed, ', tdt))
  }
  
  prev_loglvl <- loglvl
  if(reset){
    loglvl <<- 0
  }
  else if(change < 0){
    if(-change > loglvl){
      loglvl <<- 0
    }else{
      loglvl <<- loglvl + change
    }
  }else if (change > 0){
    loglvl <<- loglvl + change
  }
  if(nchar(msg)>0 && ifelse(after , prev_loglvl<supression , loglvl<supression)){
    if(after){
      cat(paste0(paste0(rep(" ",times=prev_loglvl*3),collapse="")," +-- ",msg," [ ", tdstr, round(mem_used()/(1000*1000), digits=1)," MB used ]\n"))
    }else{
      cat(paste0(paste0(rep(" ",times=loglvl*3),collapse="")," +-- ",msg," [ ", tdstr, round(mem_used()/(1000*1000), digits=1)," MB used ]\n"))
    }
  }
}

# FROM: http://musicroamer.com/blog/2011/01/16/r-tips-and-tricks-modified-pairs-plot/
# Scatterplot matrix functionality

panel.cor.scale <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  if(is.na(r))
  {
    txt="NA"
    text(0.5, 0.5, txt, cex = cex * 0.25)
  }
  else
  {
    text(0.5, 0.5, txt, cex = cex * abs(r))
  }
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  
  # This function computes the Pearsons's r for each paair of replicates for a specific pair of conditions
  # the x and y vectors are two columns of the allratios data frame. All possible pairs will be parsed through panel.cor
  # allratios is a data frame containing u columns
  # (i.e. as many replicates we have) each for one replicate for a specific condition pair with the log ratios of the two conditions (e.g. L vs H) for each
  # protein
  # see panel.lmline for an example of allratios and a debugging note
  
  # First get a backup of the usr parameter and set it to a specific value
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  # Compute the Pearson's R and print it on the plot, the greater the value the larger it will be printed
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex )
}

panel.hist <- function(x, ...){
  
  # This function draws the histograms in the diagonal of the panels, the panel.hist.breaks stores all values
  # that will break the log ratios into intervals for the histogram
  # The x vector is a column from allratios data frame. allratios is a data frame containing u columns
  # (i.e. as many replicates we have) each for one replicate for a specific condition pair with the log ratios of the two conditions (e.g. L vs H) for each
  # protein
  # see panel.lmline for an example of allratios and a debugging note
  
  # First get a backup of the usr parameter and set it to a specific value
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  
  # generate the data and rectangles of the histogram
  h <- hist(x, breaks=panel.hist.breaks,plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  #If all values are 0 create a simple rectangle in the middle:
  non_zero_values <- x != 0
  if(any(non_zero_values))
  {
    rect(breaks[-nB], 0, breaks[-1], y, col=ratios.hist.colour, ...)
  }
  else
  {
    rect(-0.25, 0, 0.25, max(y), col=ratios.hist.colour, ...)
  }
}


panel.lmline = function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = reps.scatter.lmline.colour, ...){
  
  # FROM: http://www-personal.umich.edu/~ladamic/presentations/Rtutorial/Rtutorial.R
  # The following function creates the scatterplots located in the bottom panels of the Reproducibility plot
  # the x argument is a column of allratios and the y another one. allratios is a data frame containing u columns
  # (i.e. as many replicates we have) each for one replicate for a specific condition pair with the log ratios of the two conditions (e.g. L vs H) for each
  # protein
  
  # all possible pairs of all_ratios columns will be parsed to panel.lmline
  
  # an example of all_ratios
  # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
  # |                                                                                    | log2.L.H 1 | log2.L.H 2 | log2.L.H 3 | log2.L.H 4 | log2.L.H 5 | log2.L.H 6 |
  # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
  # | ABV71576.1;ACZ79967.1;D3YRY4;P16783   [Triplex capsid protein UL46 ...]            | -2.56223   | -0.75275   | 0.30299    | 0.352557   | -0.71676   | -1.25695   |
  # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
  # | ABV71585.1;D3YRZ3;P08546;ACZ79976.1 [DNA   polymerase ...]                         | 0.585262   | -0.61003   | 0.047927   | -0.20727   | -0.43676   | -1.21811   |
  # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
  # | ABV71588.1;D3YRZ6; [Single-stranded   DNA-binding protein;Major DNA-bind ...]      | -1.18307   | -1.52482   | -0.48778   | -5.17598   | 1.473903   | 2.505526   |
  # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
  # | ABV71602.1;ACZ79991.1;ABV71603.1;ACZ79992.1;D3YS04;D3YS05;P16753;P16753-2   [ ...] | -0.8204    | -0.09631   | 1.174011   | 1.432664   | 1.765752   | 1.724972   |
  # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
  
  # debugging note:
  # notice that even if you manually set x to an allratios column and y to another the following code will fail to exevute manually because of two reasons
  # first, the panel.lmline is designed to be executed from pairs.panels thus pairs.panels running the command pairs handles creating a plot area beforehand
  # to avoid this execute plot.new() first. Also, the line abline(lm(y[ok] ~ x[ok]), col = col.smooth, ...) can not be executed manually since it inherits
  # the ellipsis arguments by the parent function. Simply erase them when testing using R console
  
  # The following lines plot a scatterplot of x vs y and computes a simple line of best fit and draws it
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  unequal_values <- x != y
  if (any(ok) && any(unequal_values))
  {
    lm_slope = coef(lm(y[ok] ~ x[ok]))[2]
    if (!is.na(lm_slope))
    {
      abline(lm(y[ok] ~ x[ok]), col = col.smooth, ...)
    }
    else
    {
      levellog("Warning!: panel.lmline: found abline with NA slope, the regression line will not be drawn")
    }
  }
}


# Called by do_results_plot (with smooth=TRUE,scale=TRUE,lm=TRUE), the x argument is all_ratios, a data frame containing u columns
# (i.e. as many replicates we have) each for one replicate for a specific condition pair with the log ratios of the two conditions (e.g. L vs H) for each
# protein
pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE,lm=FALSE){
  if (smooth){
    if (scale) {
      if(lm){
        # This line is by default executed:
        # the pairs function will produce a matrix of plots we will call "reproducibility plot" as a whole
        
        # The reproducibility plot is quite complicated, it is broken into u * u panels where the diagonal from top-left to bottom-right
        # show histograms (one per replicate), the panels above them show Pearson's R values - one for each replicate pair (i.e. u(u-1)/2 panels)
        # and the ones below show scatter plots one for each replicate pair as well
        # the function that produces the histograms is panel.hist, the one that produces the Pearson's Rs is panel.cor.scale
        # and the one that produces the scatterplots is panel.lmline
        
        pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.lmline)
      }else{
        pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth)
      }
    }else{
      if(lm){
        pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.lmline)
      }else{
        pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
      }
    }
  }else{
    if(scale){
      pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale)
    }else{
      pairs(x,diag.panel=panel.hist,upper.panel=panel.cor)
    }
  }
}



# to be used with apply on a limma-results data.frame, called by do_result_plots
# calculate mean, sd and N of ratios between available channels
# conds_cols_idxs the column indexes of quant channels intensities
# x a limma-results data.frame row (respective to one protein)
calcRowStats<-function(x,conds_cols_idxs,ratio_combs){
  ratios<-c()
  m<-c()
  std<-c()
  N<-c()
  avg.I<-c()
  # TODO: what should we do here in case there is replication mismatch?
  
  # for the following comments let u be the amount of different replicates and k the number of combinations of conditions choose 2 ({conditions}C2)
  
  for(i in 1:nrow(ratio_combs)){
    tmp1<-as.numeric(x[conds_cols_idxs[ratio_combs[i,2],]]) #columns of "lighter" label - vector of u log transformed intensities (one per replicate)
    tmp2<-as.numeric(x[conds_cols_idxs[ratio_combs[i,1],]]) #columns of "heavier" label - vector of u log transformed intensities (one per replicate)
    ratio_i<-(tmp1-tmp2) #ratios of log-2 transformed data (vector of u ratios - one per replicate)
    
    # Example of ratio_i:
    #
    #              [,1]       [,2]      [,3]       [,4]      [,5]        [,6]
    # ratio_i -1.102544 -0.3613926 -1.365139 -0.6276527 0.4658263 -0.09864526
    
    ratios<-rbind(ratios, ratio_i) # ratios will have one line of data per condition pair
    
    m<-rbind(m, mean(ratio_i,na.rm=T)) # average of all ratios calculated - one line per condition pair - this data frame has by definition one column
    std<-rbind(std, sd(ratio_i,na.rm=T)) #standard deviation of all ratios calculated - one line per condition pair - this data frame has by definition one column
    N<-rbind(N, length(which(!is.na(ratio_i) & ratio_i != 0))) # number of all ratios calculated - one line per condition pair - this data frame has by definition one column
    avg.I<-rbind(avg.I, mean(c(tmp1,tmp2),na.rm=T))  # average log 2 intensity of the protein between all conditions and replicates of the experiment
  }
  # the columns of the returned vector are: averages of ratios calculated (one column per condition combination)/standard deviation (one column per condition combination)
  #/number of ratios calculated (one column per condition combination) /average of log2 quant data between all replicates and conditions (one column per condition combination)
  # The rest of the columns are ratios calculated for the protein (one per replicate - condition pair)
  # Thus the returned vector should have k(4 + u) elements
  
  # The columns returned are very important and will be appended to the results data frame. The headers that will be given to each one of the columns
  # will be: (given we have 2 conditions and they are called L and H)
  # log2.avg.L.H         for the average of the u log ratios calculated (it gives us a measure of how much the protein is overexpressed in L compared to H)
  # log2.sd.L.H          for the standard deviation of the u log ratios calculated (a measure of deviation of the relative expression of the protein in L over H between the replicates)
  # log2.N.L.H           for the amount of the u log ratios calculated
  # log2.avg.I.L.H       for the average of all intensities measured in both conditions in all replicates (a measure of how overally abundant the protein is)
  # log2.L.H 1 .. u      u columns with the log ratios of L over H across the u different replicates
  
  return(c(m[,1],std[,1],N[,1],avg.I[,1],as.vector(t(ratios))))
}

# Produces Reproducibility, Volcano plot, MA plot and Scatterplot (matrix). Called by do_limma_analysis subroutine.
do_results_plots<-function(norm.median.intensities,time.point,exportFormat="pdf",outputFigsPrefix=""){
  levellog("",change=1)
  
  ratio_combs<-combinations(nConditions,2,1:nConditions)
  levellog("Preparing data ...")
  
  # The following line loads the results table that was created by the last function in MSdiffexp.R (do limma analysis)
  # This table must contain all information regarding the linear models applied per protein These are: the Amean that is the mean intensity of the protein
  # between all replicates, the coefficient that is the coefficient of the respective linear model, the t value of a t test that test against
  # the protein not being differentialy expressed between 2 conditions,the p value from this t statistic, an adjusted p value for multiple testing (that
  # is most of the times done by Bonferroni correction), an F statistic also testing against the protein not being differentially expressed and
  # the p value of that. The comment below describes such a typical truncated table
  
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # |                                                                                                           | A           | Coef         | t            | P.value     | P.value.adj | F           | F.p.value   |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1   [55 kDa immediate-early protein 1;Regulatory protei ...] | 24.10433767 | -0.735770791 | -0.914260979 | 0.380577872 | 0.636484278 | 0.835873138 | 0.380577872 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | AAR31362.1;ABV71635.1;P19893;Q6SWJ2;ACZ80024.1   [45 kDa immediate-early protein 2;Protein UL122;Reg ...] | 21.3205914  | -0.870745297 | -0.751712727 | 0.468321054 | 0.702943609 | 0.565072024 | 0.468321054 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | ABV71513.1;ACZ80070.1;D3YS51;P09722   [Early nuclear protein HWLF1;Tegument protein US22 ...]             | 25.90074158 | 0.446692684  | 0.602978241  | 0.555664536 | 0.750648511 | 0.363582759 | 0.555664536 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | ABV71543.1;ACZ79936.1;E2RU82;P16755   [Uncharacterized protein UL13 ...]                                  | 22.1478784  | 1.302623384  | 0.756100333  | 0.469392312 | 0.702943609 | 0.571687714 | 0.469392312 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | ABV71555.1;P16761;ACZ79948.1;D3YRW5   [Phosphoprotein 85;Phosphoprotein UL25;Tegument pro ...]            | 25.58216923 | -0.190863261 | -0.273146993 | 0.788520393 | 0.87708609  | 0.07460928  | 0.788520393 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | ABV71559.1;P16764 [Uncharacterized   protein UL29 ...]                                                    | 24.69084533 | -0.514924453 | -0.636687892 | 0.534066825 | 0.746459413 | 0.405371472 | 0.534066825 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  # | ABV71561.1;ACZ79953.1;D3YRX0;P16848   [Protein UL31;Uncharacterized protein UL31 ...]                     | 18.81611641 | -0.084411176 | -0.050037077 | 0.961212836 | 0.979920248 | 0.002503709 | 0.961212836 |
  # +-----------------------------------------------------------------------------------------------------------+-------------+--------------+--------------+-------------+-------------+-------------+-------------+
  
  # The table above is a typical table emerging from write fit in do limma analysis. In this display we have moved all headers one cell to the right for presentation purposes, in fact the first line has one header less than
  # the columns of the data showing tht the first column contains the row names, however some versions of limma have altered the write.fit function that produces the table above to add the leading \t character and
  # create the table as expected and shown above. Read the comment below the next table to see how we deal with this heterogeneity
  
  
  # In contrast when there are more than 2 conditions and all conditions are tested against all others, the columns are somewhat different providing the statistics for all comparisons
  # For example, the following table is typical for a 3plex experiment
  
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # |                                                                                                 | A           | Coef.L-H     | Coef.M-H     | Coef.M-L     | t.L-H        | t.M-H        | t.M-L        | p.value.L-H | p.value.M-H | p.value.M-L | p.value.adj.L-H | p.value.adj.M-H | p.value.adj.M-L | F           | F.p.value   |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A0AV56;B7ZLP6;Q15424;Q15424-2;B7Z2F6;F5GZU3;K7EII0;K7ES42   [Scaffold attachment factor B1 ...] | 23.60686297 | 0.270861295  | 0.053446983  | -0.217414312 | 0.516894882  | 0.101994904  | -0.414899978 | 0.611499633 | 0.919884241 | 0.683100478 | 0.999791281     | 0.999895105     | 0.99962642      | 0.149908424 | 0.861846471 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A0AVT1;A0AVT1-2;A0AVT1-3;A0AVT1-4;H0Y8S8   [Ubiquitin-like modifier-activating enzyme 6 ...]    | 22.24456697 | 0.558433305  | -0.291823418 | -0.850256723 | 0.774796695  | -0.404889567 | -1.179686262 | 0.448482308 | 0.690312526 | 0.253426552 | 0.999791281     | 0.999895105     | 0.99962642      | 0.718635052 | 0.500833966 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A0FGR8-2;H7BXI1;A0FGR8-6;A0FGR8;A0FGR8-4;A0FGR8-5;F2Z3K9   [Extended synaptotagmin-2 ...]       | 20.25866894 | 0.795171585  | 0.405121807  | -0.390049777 | 0.873132877  | 0.44484131   | -0.428291567 | 0.396293776 | 0.66275852  | 0.67449857  | 0.999791281     | 0.999895105     | 0.99962642      | 0.381226159 | 0.689427632 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A0MZ66-3;A0MZ66;A0MZ66-6;A0MZ66-5;A0MZ66-4;B7Z7Z9;A0MZ66-2;A0MZ66-7   [Shootin-1 ...]           | 22.10586558 | 0.503531284  | 0.494434774  | -0.009096509 | 0.471405061  | 0.462888925  | -0.008516135 | 0.653827745 | 0.659583917 | 0.993478348 | 0.999791281     | 0.999895105     | 0.99962642      | 0.145520471 | 0.867501947 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A1KXE4;A1KXE4-2 [Protein FAM168B ...]                                                           | 18.40201357 | -2.381544854 | 0.016774271  | 2.398319125  | -3.424623696 | 0.024121135  | 3.448744831  | 0.007490873 | 0.981278544 | 0.007209552 | 0.221396911     | 0.999895105     | 0.181096036     | 7.874156732 | 0.010417083 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A1L0T0;E9PJS0;E9PL44;E9PNL1 [Acetolactate   synthase-like protein ...]                          | 20.76969805 | 0.381252074  | -0.032625841 | -0.413877916 | 0.420595561  | -0.03599268  | -0.456588241 | 0.68144408  | 0.971876787 | 0.65607695  | 0.999791281     | 0.999895105     | 0.99962642      | 0.12888964  | 0.880265087 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  # | A1L188 [Uncharacterized protein C17orf89   ...]                                                 | 17.93232776 | 0.025989213  | -0.21841115  | -0.244400363 | 0.015707934  | -0.13200815  | -0.147716084 | 0.987971279 | 0.899247482 | 0.887354377 | 0.999791281     | 0.999895105     | 0.99962642      | 0.013164311 | 0.986950072 |
  # +-------------------------------------------------------------------------------------------------+-------------+--------------+--------------+--------------+--------------+--------------+--------------+-------------+-------------+-------------+-----------------+-----------------+-----------------+-------------+-------------+
  
  # Notice that there is a column for the {coefficient, t stat, p value, p value adjusted} for all the possible combinations of conditions. The F statistic (and its p value) being always an overall statistic for each protein is of course displayed in one column alone
  # To deal with this heterogeinity of the result tables depending on the amount of conditions per experiment we alter slightly the column names in case we have solely 2 conditions
  # (see if (nrow(ratio_combs) == 1) (...) below)
  
  # Now since other limma versions add the leading \t in the headers and others do not we have to make sure that when the table is loaded from read.table
  # the table will load the first column as row names and not as a column with undefined name. To predict the behaviour of read.table in both circumstances we read
  # in the documentation of read.table, argument row.names:
  
  # row.names:  a vector of row names. This can be a vector giving the actual row names, or a single number giving the column of the table which contains the row names, or character string giving the name of the table column containing the row names.
  #             If there is a header and the first row contains one fewer field than the number of columns, the first column in the input is used for the row names. Otherwise if row.names is missing, the rows are numbered.
  #             Using row.names = NULL forces row numbering. Missing or NULL row.names generate row names that are considered to be ‘automatic’ (and not preserved by as.matrix)
  
  # Thus, in case the first line contains one fewer header than the columns, the file will be loaded succesfully and the first column will be considered row names correctly
  # If - though - the first line has a blank first header and the headers are as much as the columns, the first column will be considered as an unnamed column
  # that **will be named "X" while loading** and an index will be applied as row names. To deal with this we will search for a possible blank first header
  # (that will be a \t leading character in the file) and if this is the case we will load the file with row.names argument = "X". Otherwise we will simply load the file
  
  # Search for the leading tab character
  
  FileHandle <- file(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""),"r")
  first_line <- readLines(FileHandle,n=1)
  hasLeadingTab <- (substr(first_line, 0, 1) == "\t")
  close(FileHandle)
  if (hasLeadingTab)
  {
    results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "", row.names = "X")
  } else {
    results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "")
  }
  
  
  
  
  if (nrow(ratio_combs) == 1) {
    # limma does not indicate the conditions compared in case only two conditions are compared in the column names so if ratio_combs = 1 add the conditions manually to the colnames of results
    colnames(results)[grep("p\\.value\\.adj",colnames(results), ignore.case = T)]<-paste("p.value.adj.",conditions.labels[2],".",conditions.labels[1],sep="")
  }
  
  
  tmp<-as.data.frame(t(norm.median.intensities))
  
  #rownames(tmp)<-colnames(norm.median.intensities)
  if (.GlobalEnv[["replicate_mismatch"]] == F) {
    # Here each column is renamed to "<condition> index" for example if the columns where Light.b1t1, Light.b1t2, Heavy.b1t1... they will be renamed to Light 1, Light 2, Heavy 1...
    nsamples<-length(colnames(tmp))/nConditions
    colnames(tmp)<-apply(data.frame(cbind(rep(conditions.labels,each=nsamples),rep(1:nsamples))),1,function(x) paste(x['X1'],x['X2']))
  } else {
    # in case not all conditions have the same amount of replicates a special way to achieve the desirable renaming should be performed:
    # first get all colnames of tmp and then make the replacement as needed
    allcols = colnames(tmp)
    allcols = str_extract(allcols, "\\..*")
    allcols = unique(allcols)
    colsreplacements = cbind(allcols, replacements = paste0(" ", 1:length(allcols)))
    for (i in 1:length(allcols)) {
      colnames(tmp) <- gsub(colsreplacements[i, "allcols"], colsreplacements[i, "replacements"], colnames(tmp), fixed = T)
    }
  }
  
  
  results<-cbind(results,tmp)
  if (.GlobalEnv[["replicate_mismatch"]] == F) {
    results$N <- apply(results[, colnames(tmp)], 1, function(x)(nsamples * nConditions) - length(which(is.na(x))))
  } else {
    results$N <- apply(results[, colnames(tmp)], 1, function(x)(length(colnames(tmp))) - length(which(is.na(x))))
  }
  if (!is.null(results$X))
  {
    results <- subset(results, select=-c(X))
  }
  levellog("Filtering data based on P-value(s) ...")
  signTruth<-rep(FALSE,nrow(results))
  for(i in 1:nrow(ratio_combs)){
    col_desc_<-paste("p.value.adj.",paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep=""),sep="")
    na_indexes<-which(is.na(results[,col_desc_]))
    if(length(na_indexes)>0){
      results[na_indexes,col_desc_]<-1
      signTruth<-(signTruth | results[,col_desc_]<pThreshold)
      results[na_indexes,col_desc_]<-NA
    }else{
      signTruth<-(signTruth | results[,col_desc_]<pThreshold)
    }
  }
  
  ndiffexp<-nrow(results[signTruth,])
  # conds_cols_idxs contain the indices on the columns for each condition. That is the first row will contain all results indices that contain quantification values deriving from this condition (but from different experiment replicates), the second will contain the same information for the second condition etc.
  conds_cols_idxs <- c()
  if (.GlobalEnv[["replicate_mismatch"]] == F) {
    # when there is not replicate mismatch creating conds_cols_idxs is easy:
    for (lbl_i in conditions.labels) {
      conds_cols_idxs <- rbind(conds_cols_idxs, grep(paste("^", lbl_i, sep = ""), colnames(results)))
    }
  } else {
    # when there is replicate mismatch not all conditions match all replicates, that means that some replicate - condition combinations should have an NA value in conds_cols_idxs
    maxcols <- 0
    for (lbl_i in conditions.labels) {
      if (length(grep(paste("^", lbl_i, sep = ""), colnames(results))) > maxcols) {
        maxcols = length(grep(paste("^", lbl_i, " ", sep = ""), colnames(results)))
      }
    }
    nsamples <- maxcols
    for (lbl_i in conditions.labels) {
      row_to_add <- grep(paste("^", lbl_i, sep = ""), colnames(results))
      row_to_add_2 <- rep(NA, maxcols)
      for (row_index in row_to_add) {
        row_to_add_2[as.integer(str_extract(colnames(results)[row_index], "\\d*$"))] <- row_index
      }
      conds_cols_idxs <- rbind(conds_cols_idxs, row_to_add_2)
    }
  }
  
  
  levellog("Calculating data frame-row statistics ...")
  # in case there is only one replicate for one condition the statistics for the proteins can not be calculated, in this case warn the user that plots concerning this conditiojn can not be calculated since R uses the std formula with "n-1" as denominator, so all standard deviations can not be calculated
  for (i in 1:nrow(conds_cols_idxs)) {
    if (sum(!is.na(conds_cols_idxs[i, ])) == 1) {
      levellog(paste0("Warn User: Condition ", conditions.labels[i], " was found replicated just once, plots displaying comparisons between this condition and others might not be available!"))
    }
  }
  
  # The following line parses all rows of results one by one to calcrowstats that calculates all statistics for this row.
  # For example the first vextor that will be sen to calcRowStats might be:
  
  #          A            Coef               t         P.value p.value.adj.L.H               F       F.p.value             H 1             H 2             H 3 
  # 24.1043377      -0.7357708      -0.9142610       0.3805779       0.6364843       0.8358731       0.3805779      24.9853920      24.2196678      24.6194861 
  # H 4             H 5             H 6             L 1             L 2             L 3             L 4             L 5             L 6               N 
  # 24.0643463              NA              NA      24.3807884      24.2406723      23.4744687      22.8498796              NA              NA       8.0000000
  
  # that correspond to a specific protein, the conds_cols_idxs are the indices of the columns that show the intensities of conditions in all replicates
  # in the example above conds_cols_idxs is:
  
  #       	[,1] [,2] [,3] [,4] [,5] [,6]
  # 	[1,]    8    9   10   11   12   13
  # 	[2,]   14   15   16   17   18   19
  
  # that means that the columns 8-13 correspond to the heavy intensities and 14-19 are the light intensities in 6 different replicates each
  
  # see calcRowStats for more information
  
  
  d<-data.frame(t(apply(results,1, function(x) calcRowStats(x,conds_cols_idxs,ratio_combs))))
  
  # d has all calculated statistics for each protein. One row corresponds to a protein and a column to a statistisc. See the comments on calcRowStats
  # to see which statistics are returned.
  
  # let u be the amount of different replicates and k the number of combinations of conditions choose 2 ({conditions}C2)
  # d now has k(4 + u) columns that are untitled but the informative titles of them can be easily produced with the following lines that store
  # the appropriate titles to colnames_d_
  
  levellog("Performing final formatting operations ...")
  colnames_d_<-c()
  # d contains all calculated statistics for ratios calculated by calcrowstats. the columns will be renamed to a more informative name. For a specific combination of conditions the columns contain the following info: mean/standard deviation/number of ratios calculated/average of log2 quant data/the rest of the columns are ratios calculated (one per replicate)
  # the first columns are standard:
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep="")
    colnames_d_<-rbind(colnames_d_, c(
      paste("log2.avg.",ratio_i_str,sep=""),
      paste("log2.sd.",ratio_i_str,sep=""),
      paste("log2.N.",ratio_i_str,sep=""),
      paste("log2.avg.I.",ratio_i_str,sep="")))
  }
  colnames_d_ <- as.vector(colnames_d_)
  
  
  # the following columns are named depending on the experimental structure
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep="") # for example Heavy.Light
    colnames_d_<-c(colnames_d_, paste(paste("log2.",ratio_i_str,sep=""),1:nsamples))
  }
  
  # as an example in a 2plex experiment of 6 differente replicates we have:
  # colnames_d_:"log2.avg.L.H"   "log2.sd.L.H"    "log2.N.L.H"     "log2.avg.I.L.H" "log2.L.H 1"     "log2.L.H 2"     "log2.L.H 3"     "log2.L.H 4"     "log2.L.H 5"    "log2.L.H 6"    
  
  # the columns of d will be renamed to colnames_d_
  
  colnames(d)<-colnames_d_
  
  # in case of replicate mismatch there will be some columns that contain NA ratios since some replicates do not contain all the conditions
  
  # Due to the mysterious bug, the following was added ...
  results$ID <- rownames(results)
  
  # all this information of statistics contained in d will be added to results creating a large data frame that has the proteins as rownames, the limma
  # analysis results (A, coefficient, pvallues and adjusted p values F statistic and its p value), the intensities of the protein per condition-replicate
  # a column called ID that has the name of the protein as well, and after that all the statistics that are contained in d as described above
  
  results<-cbind(results,d)
  
  
  results$ID <- factor(results$ID, levels=unique(as.character(results$ID)))
  results$nID<-1:nrow(results)
  
  levellog("Plotting time ...")
  theme_set(theme_bw())
  # customized colorblind-friendly palette from http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2)/
  cbPalette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")
  
  # Save the results and other parameters to a RData file, this binary copy of data will be used by PlotGenerator, an other Rscript given to the user
  # to help him rerun a part of the analysis (specifically the procedures below to create the plots) customized as needed
  
  save(results, pThreshold, quantitated_items_lbl, nConditions, calcRowStats, time.point, outputFigsPrefix, conditions.labels, IsobaricLabel, PDdata, log.intensities, norm.intensities, fit2.coefficients, file = "Plot_Generator.RData")
  
  # Now the results data frame has all the necessary information to create all target plots, notice that results will be expanded
  # even more while generating the plots, go to the comments after plot generation to see the new columns
  
  # An example of the first rows of results at the moment is:
  # 
  #   +-----------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+-----------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+-----------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+-----+
  #   |                                                                                                           | A        | Coef     | t        | P.value  | p.value.adj.L.H | F        | F.p.value | H 1      | H 2      | H 3      | H 4      | H 5      | H 6      | L 1      | L 2      | L 3      | L 4      | L 5      | L 6      | N  | ID                                                                                                        | log2.avg.L.H | log2.sd.L.H | log2.N.L.H | log2.avg.I.L.H | log2.L.H 1 | log2.L.H 2 | log2.L.H 3 | log2.L.H 4 | log2.L.H 5 | log2.L.H 6 | nID |
  #   +-----------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+-----------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+-----------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+-----+
  #   | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1   [55 kDa immediate-early protein 1;Regulatory protei ...] | 24.10434 | -0.73577 | -0.91426 | 0.380578 | 0.636484        | 0.835873 | 0.380578  | 24.98539 | 24.21967 | 24.61949 | 24.06435 | NA       | NA       | 24.38079 | 24.24067 | 23.47447 | 22.84988 | NA       | NA       | 8  | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1 [55 kDa immediate-early   protein 1;Regulatory protei ...] | -0.73577     | 0.573453    | 4          | 24.10434       | -0.6046    | 0.021004   | -1.14502   | -1.21447   | NA         | NA         | 1   |
  #   +-----------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+-----------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+-----------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+-----+
  #   | AAR31362.1;ABV71635.1;P19893;Q6SWJ2;ACZ80024.1   [45 kDa immediate-early protein 2;Protein UL122;Reg ...] | 21.32059 | -0.87075 | -0.75171 | 0.468321 | 0.702944        | 0.565072 | 0.468321  | 23.30772 | 23.70334 | 20.22017 | 19.79263 | NA       | NA       | 22.38749 | 22.90319 | 21.32653 | 16.92367 | NA       | NA       | 8  | AAR31362.1;ABV71635.1;P19893;Q6SWJ2;ACZ80024.1 [45 kDa immediate-early   protein 2;Protein UL122;Reg ...] | -0.87075     | 1.623705    | 4          | 21.32059       | -0.92022   | -0.80015   | 1.106357   | -2.86897   | NA         | NA         | 2   |
  #   +-----------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+-----------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+-----------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+-----+
  #   | ABV71513.1;ACZ80070.1;D3YS51;P09722   [Early nuclear protein HWLF1;Tegument protein US22 ...]             | 25.90074 | 0.446693 | 0.602978 | 0.555665 | 0.750649        | 0.363583 | 0.555665  | 25.85704 | 26.30705 | 26.76279 | 26.47704 | 25.52645 | 23.13399 | 26.81915 | 25.76735 | 26.76279 | 27.31044 | 26.68177 | 23.40302 | 12 | ABV71513.1;ACZ80070.1;D3YS51;P09722 [Early nuclear protein HWLF1;Tegument   protein US22 ...]             | 0.446693     | 0.651371    | 5          | 25.90074       | 0.962109   | -0.5397    | 0          | 0.833401   | 1.15532    | 0.269025   | 3   |
  #   +-----------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+-----------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+-----------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+-----+
  
  
  # for the comments in plot generation let u be the amount of different replicates and k the number of combinations of conditions choose 2 ({conditions}C2)
  
  
  # Plot generation:
  for(i in 1:nrow(ratio_combs)){
    levellog(paste("Generating plots for combination #",i," ..."),change=1,after=T)
    ratio_i_str <- paste(conditions.labels[ratio_combs[i, 2]], ".", conditions.labels[ratio_combs[i, 1]], sep = "") # for example H.L
    
    # 1 - Volcano plot (-log10 P-value vs log ratio)
    result <- tryCatch({
      levellog("Making volcano plot ...")
      figsuffix<-paste("_",ratio_i_str,"-volcano","_",sep="")
      if(exportFormat == "pdf"){
        pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      
      ratio_i_<-paste("log2.",ratio_i_str,sep="")
      ratio_i_sd_col<-paste("log2.sd.",ratio_i_str,sep="")
      
      # Before starting the volcano plot generation we will compute all values needed for our plots. Forst goal is to find per replicate the (log2 intensity) +- (sd of log intensities) for the current condition pair
      # tmp2 is a data frame that gets for the current condition pair and for each protein the value (log2 intensity) + (sd of log intensities) for all the replicates
      
      tmp2<-results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_, " ")),colnames(results))]]+results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_sd_col, "$")),colnames(results))]]
      
      # an example of tmp2 is:
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # |                                                                                    | log2.L.H 1 | log2.L.H 2 | log2.L.H 3 | log2.L.H 4 | log2.L.H 5 | log2.L.H 6 |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71576.1;ACZ79967.1;D3YRY4;P16783   [Triplex capsid protein UL46 ...]            | -1.47929   | 0.330185   | 1.385926   | 1.435492   | 0.366172   | -0.17401   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71585.1;D3YRZ3;P08546;ACZ79976.1 [DNA   polymerase ...]                         | 1.196955   | 0.001665   | 0.65962    | 0.404419   | 0.17493    | -0.60642   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71588.1;D3YRZ6;P17147; [Single-stranded DNA-binding protein;Major DNA-bind ...] | 1.49977    | 1.158022   | 2.19506    | -2.49314   | 4.156746   | 5.188369   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71602.1;ACZ79991.1;ABV71603.1;ACZ79992.1;D3YS04;D3YS05;P16753;P16753-2   [ ...] | 0.250421   | 0.974516   | 2.244834   | 2.503487   | 2.836576   | 2.795796   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # 
      # notice that the column names are not informative, in fact log2.L.H 1 could better be named (log2.L.H 1 + log2.sd.L.H) etc.
      
      # tmp1 stores the values (log2 intensity) - (sd of log intensities)
      tmp1<-results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_, " ")),colnames(results))]]-results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_sd_col, "$")),colnames(results))]]
      
      
      # an example is:
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # |                                                                                                | log2.L.H 1 | log2.L.H 2 | log2.L.H 3 | log2.L.H 4 | log2.L.H 5 | log2.L.H 6 |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71576.1;ACZ79967.1;D3YRY4;P16783   [Triplex capsid protein UL46 ...]                        | -3.64516   | -1.83569   | -0.77995   | -0.73038   | -1.7997    | -2.33988   |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71585.1;D3YRZ3;P08546;ACZ79976.1 [DNA   polymerase ...]                                     | -0.02643   | -1.22172   | -0.56377   | -0.81897   | -1.04846   | -1.82981   |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71588.1;D3YRZ6; [Single-stranded   DNA-binding protein;Major DNA-bind ...]                  | -3.86591   | -4.20766   | -3.17062   | -7.85882   | -1.20894   | -0.17732   |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71602.1;ACZ79991.1;ABV71603.1;ACZ79992.1;D3YS04;D3YS05;P16753;P16753-2   [ ...]             | -1.89123   | -1.16713   | 0.103187   | 0.36184    | 0.694928   | 0.654148   |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71605.1;ACZ79994.1;D3YS07;P06725   [Tegument protein pp65;65 kDa matrix phosphoprotein ...] | -0.27818   | -0.70204   | -1.46206   | -1.0021    | -0.75168   | -0.4241    |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71607.1;D3YS09;P16728;ACZ79996.1   [Probable capsid protein VP23;Protein UL85 ...]          | -3.802     | -4.14723   | -0.04019   | -0.14354   | -0.52758   | -0.53625   |
      # +------------------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # 
      
      # it can be easily noted that the corresponding cells in tmp1 are 2sd lower than the ones in tmp2 as expected (in the case of Triplex capsid protein UL46 e.g. 2sd = 2* 1.083 = 2.166 that is the difference of -3.645 and -1.479)
      
      # to create our histograms we need to find the maximum absolute value stored in both tmp1 and tmp2
      # one easy to understand way would be to bind these data frames compute their abs values and get the maximum one:
      
      # > mytmp3 = cbind(tmp1, tmp2)
      # > mytmp3 = abs(mytmp3)
      # > max(mytmp3, na.rm=T)
      
      # the following line does the same thing but is more effective and returns the ceiling of the respective result
      
      ratiolim<-ceiling(max(max(range(tmp1,na.rm=T),range(tmp2,na.rm=T)),abs(min(range(tmp1,na.rm=T),range(tmp2,na.rm=T)))))
      
      # This will help break all values in intervals of size 1 from -ratiolim to +ratiolim
      
      if(ratiolim == 0)
      {
        ratiolim <- 5
      }
      panel.hist.breaks<<-(-ratiolim:ratiolim)
      
      
      ratio_i_p.value.adj<-paste("p.value.adj.",paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep=""),sep="")
      ratio_i_avg_col<-paste("log2.avg.",ratio_i_str,sep="")
      mlog10_ratio_i_p.value.adj<-paste("mlog10_",ratio_i_p.value.adj,sep="")
      diffexp_ratio_i<-paste("diffexp_",ratio_i_str,sep="")
      
      # -log(p_val) [where p_val is the adjusted p value for the differential expression of a protein in the current condition pair] is the y axis
      # in a typical volcano plot. The following line computes this and adds it to the results data frame as new column (named e.g. mlog10_p.value.adj.L.H)
      results[,mlog10_ratio_i_p.value.adj]<-(-log10(results[,ratio_i_p.value.adj]))
      
      # the following lines will add another (boolean) column to results indicating if the respective protein is differentially expressed or not
      # it simply marks as differentially expressed each one that has an adjusted p value < pThreshold and NA if the pValue was not calculated for some reason
      na_indexes<-which(is.na(results[,ratio_i_p.value.adj]))
      if(length(na_indexes)>0){
        results[na_indexes,ratio_i_p.value.adj]<-1
        results[,diffexp_ratio_i]<-results[,ratio_i_p.value.adj]<pThreshold
        results[na_indexes,ratio_i_p.value.adj]<-NA
      }else{
        results[,diffexp_ratio_i]<-results[,ratio_i_p.value.adj]<pThreshold
      }
      
      # Prepare the labels of the axes for the volcano plot
      
      if(!IsobaricLabel)
      {
        myxlab <- paste("average log2 ",sub("\\.","/",ratio_i_str),sep="")
      }else{
        if(!PDdata)
        {
          myxlab <- paste("average log2 ", ratio_i_str, sep="")
          myxlab <- gsub("Reporter\\.intensity\\.", "Reporter ", myxlab)
        }else{
          myxlab <- paste("average log2 ",ratio_i_str ,sep="")
          myxlab <- gsub("X([[:digit:]])", "\\1", myxlab)
        }
      }
      myxlab <- gsub("\\.", "/", myxlab)
      
      
      # Volcano plot is nothing more than a plot of the proteins ratios (the average in all replicates) vs -log10 of the p-value for this ratio:
      p<-ggplot(data=results, aes_string(x=ratio_i_avg_col, y=mlog10_ratio_i_p.value.adj, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=1.75) +
        theme(legend.position = "none", axis.title.y=element_text(vjust=0.2), axis.title.x=element_text(vjust=0), plot.title = element_text(vjust=1.5, lineheight=.8, face="bold")) +
        xlim(c(-ratiolim, ratiolim)) + ylim(c(0, 6)) + scale_colour_manual(values=cbPalette) +
        xlab(myxlab) + ylab("-log10 P-value") + ggtitle("P-value vs Fold change") +
        geom_hline(aes(yintercept=-log10(pThreshold)), colour="#990000", linetype="dashed") +
        geom_text(size=2.5, hjust=1, vjust=-0.5,aes(x=-4.2, y=-log10(pThreshold)), label=paste0("P-value=", pThreshold),colour="#990000")
      
      
      print(p)
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,time.point,".emf",sep=""),type="emf")
      }
      dev.off()
      png(paste("../",outputFigsPrefix,figsuffix,time.point,".png",sep=""), width = 1500, height = 1050)
      
      p<-ggplot(data=results, aes_string(x=ratio_i_avg_col, y=mlog10_ratio_i_p.value.adj, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=4.25) +
        theme(legend.position = "none", axis.title.y=element_text(size = 22.5, vjust=0.2), axis.title.x=element_text(size = 22.5, vjust=0), plot.title = element_text(size = 30, vjust=1.5, lineheight=.8, face="bold"), axis.text.x = element_text(size = 22.5), axis.text.y = element_text(size = 22.5)) +
        xlim(c(-ratiolim, ratiolim)) + ylim(c(0, 6)) + scale_colour_manual(values=cbPalette) +
        xlab(myxlab) + ylab("-log10 P-value") + ggtitle("P-value vs Fold change") +
        geom_hline(aes(yintercept=-log10(pThreshold)), colour="#990000", linetype="dashed") +
        geom_text(size=5, hjust=1, vjust=-0.5,aes(x=-4.2, y=-log10(pThreshold)), label=paste0("P-value=", pThreshold),colour="#990000") 
      print(p)
      dev.off()#script for volcano plot END
    }, error = function(err){
      levellog(paste0("Warn User: ", ratio_i_str, " volcano plot failed"))
    })
    
    # 2 - value-ordered - log ratio
    result <- tryCatch({
      levellog("Making value-ordered plot ...")
      figsuffix<-paste("_",ratio_i_str,"-value-ordered-log-ratio","_",sep="")
      
      if(exportFormat == "pdf"){
        pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      
      # The value ordered plot plots the average log ratio of the proteins for the specific condition pair(a measure of how much the protein is overexpressed in L compared to H)
      # for all proteins sorted from the one with the smallest average log ratio to the one with the highest (from the most underexpressed protein in e.g. L over H to the most overexpressed one)
      # We also mark the +- sd of the log ratios and we mark the differentially expressed proteins with orange. It is quite common in these cases
      # to find differentially expressed proteins in areas of significan under or overexpression i.e. to the left or to the right of the graph
      
      
      results<-results[with(results, order(results[,c(ratio_i_avg_col)])),] # sort the data from the most under to the most overexpressed protein
      results$nID<-1:nrow(results) # add an ID to help the plotting procedure
      
      # Add 2 new columns to results to help the plotting procedure with the formulae {log ratio + sd} and {log ratio - sd}
      
      ratio_i_avg_col_ymax<-paste(ratio_i_avg_col,".ymax",sep="")
      ratio_i_avg_col_ymin<-paste(ratio_i_avg_col,".ymin",sep="")
      results[,ratio_i_avg_col_ymax]<-results[,ratio_i_avg_col]+results[,ratio_i_sd_col]
      results[,ratio_i_avg_col_ymin]<-results[,ratio_i_avg_col]-results[,ratio_i_sd_col]
      
      # Prepare the labels of the plot and generate it
      
      if(!IsobaricLabel)
      {
        myylab <- paste("average log2 ",sub("\\.","/",ratio_i_str),sep="")
      }else{
        if(!PDdata)
        {
          myylab <- paste("average log2 ", ratio_i_str, sep="")
          myylab <- gsub("Reporter\\.intensity\\.", "Reporter ", myylab)
        }else{
          myylab <- paste("average log2 ", ratio_i_str, sep="")
          myylab <- gsub("X([[:digit:]])", "\\1", myylab)
        }
      }
      myylab <- gsub("\\.", "/", myylab)
      p<-ggplot(data=results, aes_string(x="nID", y=ratio_i_avg_col, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=1.5) +
        geom_errorbar(aes_string(ymin=ratio_i_avg_col_ymin, ymax=ratio_i_avg_col_ymax), width=1.5) +
        theme(legend.position = "none", axis.title.y=element_text(vjust=0.2), axis.title.x=element_text(vjust=0), plot.title = element_text(vjust=1.5, lineheight=.8, face="bold")) +
        ylim(c(-ratiolim, ratiolim)) + scale_colour_manual(values=cbPalette) +
        xlab(paste(quantitated_items_lbl,"ID")) + ylab(myylab) + ggtitle("Value-ordered fold change")
      print(p)
      ggsave(paste(outputFigsPrefix,figsuffix,time.point,".png",sep=""), plot = p, device = "png", path = "..")
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,time.point,".emf",sep=""),type="emf")
      }
      dev.off()    
    }, error = function(err){
      levellog(paste0("Warn User: ", ratio_i_str, " value-ordered plot failed"))
    })
    
    # 3 - MA plot
    result <- tryCatch({
      levellog("Making MA plot ...")
      
      # The MA plot show the average log intensity of a protein between all the intensity mesaures that happened for a condition pair (i.e. across both
      # conditions and all u replicates at once) - that is a measure of how generally abundant the protein is against the (average) log ratio of the protein
      # for our current condition pair (how much the protein is over or under expressed). We also mark the differentially expressed proteins with orange.
      # As described in the value ordered plot the DE proteins will tend to be found in areas of the plot where significantly over or under expressed proteins lie
      # i.e. in this plot in the highest and lowest values on the y axis
      
      figsuffix<-paste("_",ratio_i_str,"-MA","_",sep="")
      ratio_i_avgI_col<-paste("log2.avg.I.",ratio_i_str,sep="")
      
      if(exportFormat == "pdf"){
        pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      
      if(!IsobaricLabel)
      {
        myylab <- paste("A (average log2 ",sub("\\.","/",ratio_i_str),")",sep="")
      }else{
        if(!PDdata)
        {
          myylab <- paste("A (average log2 ", ratio_i_str, ")", sep="")
          myylab <- gsub("Reporter\\.intensity\\.", "Reporter ", myylab)
        }else{
          myylab <- paste("A (average log2 ",ratio_i_str,")",sep="")
          myylab <- gsub("X([[:digit:]])", "\\1", myylab)
        }
      }
      myylab <- gsub("\\.", "/", myylab)
      
      p<-ggplot(data=results, aes_string(x=ratio_i_avgI_col, y=ratio_i_avg_col, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=1.75) +
        theme(legend.position = "none", axis.title.y=element_text(vjust=0.2), axis.title.x=element_text(vjust=0), plot.title = element_text(vjust=1.5, lineheight=.8, face="bold")) +
        ylim(c(-ratiolim, ratiolim)) + scale_colour_manual(values=cbPalette) +
        xlab("M (average log2 Intensity)") + ylab(myylab) + ggtitle("MA plot")
      print(p)
      ggsave(paste(outputFigsPrefix,figsuffix,time.point,".png",sep=""), plot = p, device = "png", path = "..")
      
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,time.point,".emf",sep=""),type="emf")
      }
      dev.off()
    }, error = function(err){
      levellog(paste0("Warn User: ", ratio_i_str, " MA plot failed"))
    })
    
    # 4 - Reproducibility plots & histograms
    result <- tryCatch({
      levellog("Making reproducibility plot ...")
      figsuffix<-paste("_",ratio_i_str,"-reproducibility","_",sep="")
      
      # The reproducibility plot is quite complicated, it is broken into u * u panels where the diagonal from top-left to bottom-right
      # show histograms (one per replicate), the panels above them show Pearson's R values - one for each replicate pair (i.e. u(u-1)/2 panels)
      # and the ones below show scatter plots one for each replicate pair as well
      
      # the function that produces the histograms is panel.hist, the one that produces the Pearson's Rs is panel.cor.scale
      # and the one that produces the scatterplots is panel.lmline, they are called by the pairs.panels function below that sends all_ratios as argument to them
      # see the respective functions for more information
      
      # First lets isolate the log ratios of all proteins for the current condition pair across all replicates
      allratios <- results[, colnames(results)[grep(paste0(ratio_i_, " "), colnames(results))]]
      
      # allratios now has u columns and an example of it is shown below:
      
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # |                                                                                    | log2.L.H 1 | log2.L.H 2 | log2.L.H 3 | log2.L.H 4 | log2.L.H 5 | log2.L.H 6 |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71576.1;ACZ79967.1;D3YRY4;P16783   [Triplex capsid protein UL46 ...]            | -2.56223   | -0.75275   | 0.30299    | 0.352557   | -0.71676   | -1.25695   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71585.1;D3YRZ3;P08546;ACZ79976.1 [DNA   polymerase ...]                         | 0.585262   | -0.61003   | 0.047927   | -0.20727   | -0.43676   | -1.21811   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71588.1;D3YRZ6; [Single-stranded   DNA-binding protein;Major DNA-bind ...]      | -1.18307   | -1.52482   | -0.48778   | -5.17598   | 1.473903   | 2.505526   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      # | ABV71602.1;ACZ79991.1;ABV71603.1;ACZ79992.1;D3YS04;D3YS05;P16753;P16753-2   [ ...] | -0.8204    | -0.09631   | 1.174011   | 1.432664   | 1.765752   | 1.724972   |
      # +------------------------------------------------------------------------------------+------------+------------+------------+------------+------------+------------+
      
      
      # drop all the columns that have all values equal to NA
      allratios <- allratios[, colSums(is.na(allratios)) < nrow(allratios)]
      
      # Prepare the labels of the plot
      
      if(!IsobaricLabel)
      {
        colnames(allratios)<-sub(ratio_i_,paste("log2(",sub("\\.","/",ratio_i_str),") ",sep=""),colnames(allratios))
      }else{
        if(!PDdata){
          colnames(allratios)<-sub(ratio_i_,paste("log2(",ratio_i_str,") ",sep=""),colnames(allratios))
          colnames(allratios) <- gsub("Reporter\\.intensity\\.", "Reporter ", colnames(allratios))
        }else{
          colnames(allratios)<-sub(ratio_i_,paste("log2(",ratio_i_str,") ",sep=""),colnames(allratios))
          colnames(allratios) <- gsub("X([[:digit:]])", "\\1", colnames(allratios))
        }
      }
      
      
      colnames(allratios) <- gsub("\\.", "/", colnames(allratios))
      
      # Call pairs.panels to generate the plot. See the comments in pairs.panel for more information
      
      if(exportFormat == "pdf"){
        pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      pairs.panels(allratios,scale=T,lm=T)
      dev.off()
      png(paste("../",outputFigsPrefix,figsuffix,time.point,".png",sep=""), width = 1500, height = 1050)
      pairs.panels(allratios,scale=T,lm=T)
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,time.point,".emf",sep=""),type="emf")
      }
      dev.off()
    }, error = function(err){
      levellog(paste0("Warn User: ", ratio_i_str, " reproducibility plot failed"))
    })
    levellog("",change=-1)
  }
  
  # While generating the plots, results was expanded with new columns giving more information per protein. these columns are (for 2 conditions L and H):
  # mlog10_p.value.adj.L.H, diffexp_L.H,                                                     log2.avg.L.H.ymax, log2.avg.L.H.ymin that correspond to:
  # the -log(p_val),        a boolean (if the protein was differentially expressed or not),  log ratio + sd,    log ratio - sd
  # See the comments above, just before generating the volcano plot for more information
  
  levellog("Saving plots results to files ...")
  
  # The last thing to do in this function is to rename some column names of the results data frame to make them more informative
  # An example of column names in case we have two conditions i.e. L and H is:
  #
  # "A"                      "Coef"                   "t"                      "P.value"                "p.value.adj.L.H"        "F"                     
  # "F.p.value"              "H 1"                    "H 2"                    "H 3"                    "H 4"                    "H 5"                   
  # "H 6"                    "L 1"                    "L 2"                    "L 3"                    "L 4"                    "L 5"                   
  # "L 6"                    "N"                      "ID"                     "log2.avg.L.H"           "log2.sd.L.H"            "log2.N.L.H"            
  # "log2.avg.I.L.H"         "log2.L.H 1"             "log2.L.H 2"             "log2.L.H 3"             "log2.L.H 4"             "log2.L.H 5"            
  # "log2.L.H 6"             "nID"                    "mlog10_p.value.adj.L.H" "diffexp_L.H"            "log2.avg.L.H.ymax"      "log2.avg.L.H.ymin"   
  # 
  # The content of each one of the columns is described in comments above
  
  colnames(results)<-gsub("[\\.]"," ",colnames(results))
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste("(",conditions.labels[ratio_combs[i,2]],") (",conditions.labels[ratio_combs[i,1]],")",sep="")
    colnames(results)<-gsub(ratio_i_str,"\\1/\\2",colnames(results))
  }
  colnames(results)<-gsub("p value adj","P-value adjusted",colnames(results))
  colnames(results)<-gsub("p value","P-value",colnames(results))
  colnames(results)<-gsub("mlog10","-log10",colnames(results))
  colnames(results)<-gsub("log2 avg","avg log2",colnames(results))
  colnames(results)<-gsub("log2 sd","sd log2",colnames(results))
  colnames(results)<-gsub("log2 N","N log2",colnames(results))
  colnames(results)<-gsub("ymax","+sd",colnames(results))
  colnames(results)<-gsub("ymin","-sd",colnames(results))
  
  # Erase unnecessary columns
  
  results<-results[,which(!grepl("^-log10",colnames(results)))]
  results<-results[,which(!grepl("sd$",colnames(results)))]
  results<-results[,which(!grepl("^nID$",colnames(results)))]
  results<-results[,which(!grepl("^diffexp",colnames(results)))]
  
  # Log the amount of quantified proteins
  
  quant_species<-"proteins"
  if(!ProteinQuantitation){
    quant_species<-"peptides"
  }  
  
  levellog(paste("do_results_plots: Quantified ",quant_species,": ",nrow(results)," (",time.point,")",sep=""))
  
  conditions.labels <- gsub("\\.", " ", conditions.labels)
  
  
  # Log amount of differentially expressed proteins
  
  for(i in 1:nrow(ratio_combs)){
    col_desc_<-paste("P-value adjusted ",paste(conditions.labels[ratio_combs[i,2]],"/",conditions.labels[ratio_combs[i,1]],sep=""),sep="")
    ndiffexp_tmp<-length(which(results[,col_desc_]<pThreshold))
    levellog(paste("do_results_plots: Differentially expressed for ",conditions.labels[ratio_combs[i,2]]," vs ",conditions.labels[ratio_combs[i,1]]," : ",ndiffexp_tmp,sep=""))
  }
  
  if(nrow(ratio_combs) > 1){
    levellog(paste("do_results_plots: Differentially expressed in at least one combination of conditions: ",ndiffexp,sep=""))
  }
  
  # Create the diffexp data frame, a data frame that will contain information only for the differentially expressed proteins
  # and will be printed to the user
  
  # Isolate the pvalue adjusted, average log ratio, average of intensities in all proteins and standard deviation of log ratios 
  diffexp<-results[,c(grep("^P-value adjusted",colnames(results)),grep("avg log2 [^I]+",colnames(results)),grep("sd log2 ",colnames(results)),grep("N log2 ",colnames(results)),grep("avg log2 I ",colnames(results)))]
  
  # Add a column names "Proteins" (or "Peptides") with the proteins' names
  diffexp[,quantitated_items_lbl]<-rownames(diffexp)
  
  # Isolate the average log ration and the p value adjusted
  
  diffexp<-diffexp[,c(grep(quantitated_items_lbl,colnames(diffexp)),grep("avg log2 [^I]+",colnames(diffexp)),grep("log2 sd ",colnames(diffexp)),grep("P-value ",colnames(diffexp)), grep("log2 N ",colnames(diffexp)),grep("log2 avg I ",colnames(diffexp)))]
  
  # Get a copy of the protein groups data frame that contains the intensity of each protein in all replicates, the unique sequences and the ratio counts
  # i.e. amount of peptides quantified for a protein per replicate of these proteins as well as other information
  # and rename the "Protein.ID" column to Protein to help merge with the data of diffexp above
  
  tmp_protein_groups<-protein_groups
  colnames(tmp_protein_groups)[grep(paste(quantitated_items_lbl,".IDs",sep=""),colnames(tmp_protein_groups))]<-quantitated_items_lbl
  
  # Merge diffexp and protein groups and keep only the ratio counts for all the replicates as information from the protein groups data frame
  
  diffexp<-merge(diffexp,tmp_protein_groups[,c(quantitated_items_lbl,sort(colnames(tmp_protein_groups)[grep("Ratio\\.counts",colnames(tmp_protein_groups))]))],by=quantitated_items_lbl,all.x=T)
  
  # Add a new column with the sum of all ratio counts (from al replicates) per protein and rename it to Ratio.counts.total
  diffexp$newcol<-rowSums(diffexp[,colnames(diffexp)[grep("Ratio\\.counts",colnames(diffexp))]],na.rm=T)
  colnames(diffexp)[length(colnames(diffexp))]<-"Ratio.counts.total"  
  
  
  # Create a boolean vector showing which lines in diffexp are respective to differentialy expressed proteins
  signTruth<-rep(FALSE,nrow(diffexp))
  col_desc_<-grep("P-value ",colnames(diffexp))
  for(cond_i_col in col_desc_){
    na_indexes<-which(is.na(diffexp[,cond_i_col]))
    if(length(na_indexes)>0){
      diffexp[na_indexes,cond_i_col]<-1
      signTruth<-(signTruth | diffexp[,cond_i_col]<pThreshold)
      diffexp[na_indexes,cond_i_col]<-NA
    }else{
      signTruth<-(signTruth | diffexp[,cond_i_col]<pThreshold)
    }    
  }
  
  colnames(diffexp)<-gsub("\\."," ",colnames(diffexp))
  
  if(!ProteinQuantitation){
    diffexp<-merge(diffexp,protein_groups[c("Peptide.IDs","Protein.IDs")],by.x=c("Peptide"),by.y=c("Peptide.IDs"),all.x=T)
    colnames(diffexp)[grep("^Protein\\.IDs$",colnames(diffexp))]<-"Protein"
  }
  
  
  dec <- "."
  # Restore original rep description for output (in case bio or tech reps were not successive i.e. we have 3 bio reps that are b1, b3, b4)
  newcolumns <- names(diffexp[signTruth,])
  oldcolumns <- newcolumns
  for(my_column in newcolumns){
    for(my_repdesc in .GlobalEnv[["rep_structure"]]$rep_desc){
      if (grepl(my_repdesc, my_column)){
        temp_name <- .GlobalEnv[["original_rep_structure"]]$rep_desc[match(my_repdesc, .GlobalEnv[["rep_structure"]]$rep_desc)]
        newcolumns[match(my_column, newcolumns)] <- sub(my_repdesc, temp_name, my_column)
      }
    }
  }
  
  
  colnames(diffexp) <- newcolumns
  
  # Write the DE proteins from diffexp to a table
  write.table(diffexp[signTruth,],dec=dec,file=paste(outputFigsPrefix,"_diffexp_",time.point,".txt",sep=""),sep="\t",row.names=F,quote=F)
  
  colnames(diffexp) <- oldcolumns
  
  # Merge diffexp and results to create a large data frame containing in fact all the information created from ProteoSign up to now,
  # i.e. all information stored in results and the ratio counts per protein
  diffexp<-merge(diffexp,results[,-grep("^(avg log2|P-value adjusted)",colnames(results))],by.x=c(quantitated_items_lbl),by.y=c("ID"),all.x=T)
  
  # Save this information to a file
  write.table(diffexp,dec=dec,file=paste(outputFigsPrefix,"_results_",time.point,".txt",sep=""),sep="\t",row.names=F,quote=F)
  
  levellog("",change=-1)
  return(results)
}

# Performs the differential expression analysis through limma, after quantile normalization.
do_limma_analysis<-function(working_pgroups,time.point,exp_design_fname,exportFormat="pdf",outputFigsPrefix=""){
  levellog("",change=1)
  levellog("Preparing limma input data frame ...")
  # Read the sample key
  # Assigns sample names (from the data file) to groups
  # Sample order must be the same as the main data file, but excludes technical
  # replicates as we will aggregate into one value per sample.
  sample.key <- read.delim(exp_design_fname, header=TRUE,row.names=1,colClasses="character")
  sample.key$Category <- factor(sample.key$Category, levels=unique(sample.key$Category))
  
  # Extract protein/peptide quantitation columns only from quantitation input file
  # The quantitation data is in columns 10 to 90.
  
  #prot.intensities <- quantitation[,10:90]
  prot.intensities <- working_pgroups
  
  # Extract the protein names/peptide sequences (imported from the data file) into a
  # separate list for future reference
  
  prot.names <- rownames(prot.intensities)
  
  # Take log2 of intensities
  
  log.intensities  <- log2(prot.intensities)
  
  setwd(limma_output)
  
  #rename the column names to their original ones to display to the graph
  names(log.intensities) <- rownames(sample.key)
  levellog("Saving limma input frame ...")
  #also restore original repdesc to output to the limma-input file
  newcolumns <- names(working_pgroups)
  oldcolumns <- newcolumns
  for(my_column in newcolumns){
    for(my_repdesc in .GlobalEnv[["rep_structure"]]$rep_desc){
      if (grepl(my_repdesc, my_column)){
        temp_name <- .GlobalEnv[["original_rep_structure"]]$rep_desc[match(my_repdesc, .GlobalEnv[["rep_structure"]]$rep_desc)]
        newcolumns[match(my_column, newcolumns)] <- sub(my_repdesc, temp_name, my_column)
      }
    }
  }
  colnames(working_pgroups) <- newcolumns
  write.table(working_pgroups,file=paste(outputFigsPrefix,"_limma-input_",quantitated_items_lbl,"Groups.txt",sep=""),sep="\t",row.names = T, col.names=NA)
  colnames(working_pgroups) <- oldcolumns
  if(exportFormat == "pdf"){
    pdf(file=paste(outputFigsPrefix,"_limma-graphs_",time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
  }
  
  #Remove infinite records
  
  for(i in 1:ncol(log.intensities)){
    mi <- NA
    mi <- which(log.intensities[,i] == -Inf)
    if(length(mi) != 0){
      log.intensities <- log.intensities[-mi,]
    }
    mi <- which(log.intensities[,i] == Inf)
    if(length(mi) != 0){
      log.intensities <- log.intensities[-mi,]
    }
  }
  
  if (IsobaricLabel){
    if(!PDdata){
      colnames(log.intensities) <- sub('Reporter\\.intensity\\.', "Reporter.", colnames(log.intensities))
    }else{
      colnames(log.intensities) <- sub('X([[:digit:]])', "\\1", colnames(log.intensities))
    }
  }
  
  # Box plot before normalisation
  boxplot(log.intensities)
  title(main="Intensities Before Normalisation")
  log.intensities <<- log.intensities
  
  # Perform quantile normalisation
  levellog("Performing quantile normalisation ...")
  norm.intensities <- normalizeBetweenArrays(data.matrix(log.intensities), method="quantile");
  
  
  # Box plot after normalisation
  boxplot(norm.intensities)
  title(main="Intensities After Normalisation")
  norm.intensities <<- norm.intensities
  
  norm.median.intensities<-as.data.frame(t(as.matrix(norm.intensities)))
  
  # Assign row names to our aggregated intensities from the sample key
  
  row.names(norm.median.intensities) <- row.names(sample.key)
  #Assign a blocking variable that groups the data based on technical replication
  blocking_var<-c()
  if(.GlobalEnv[["n_techreps"]] > 1){
    #In Label-Free projects each raw file is derived from a different experiment. Thus, data from different conditions
    #must be in different groups
    if(.GlobalEnv[["LabelFree"]]){
      max_bioreps = max(rep_structure$biorep)
      for(j in 1:nConditions)
      {
        for(i in unique(rep_structure$biorep)){
          blocking_var<-c(blocking_var, rep((j-1)*max_bioreps + i,length(unique(rep_structure[rep_structure$biorep == i,]$techrep))))
        }
      }
    }
    else{
      #All other supported experiment methods are multi-plexed, so all conditions derived from the same
      #technical replicate must be in the same group
      for(i in unique(rep_structure$biorep)){
        blocking_var<-c(blocking_var, rep(i,length(unique(rep_structure[rep_structure$biorep == i,]$techrep))))
      }
      blocking_var<-rep(blocking_var, nConditions)
    }
  }else{
    #When there is no technical replication, each sample.key record must be in its own group
    blocking_var<-1:nrow(sample.key)
  }
  
  # Setup design matrix
  # This specifies the design of the experiment for limma, replicating
  # the info in the sample key, but representing it in a matrix format
  levellog("Constructing the design matrix ...")
  design <- model.matrix(~0 + factor(sample.key$Category))
  colnames(design) <- levels(sample.key$Category)
  write.table(design,file=paste(outputFigsPrefix,"_limma-design-matrix_",quantitated_items_lbl,"Groups.txt",sep=""),sep="\t",row.names = T, col.names=NA)
  write.table(blocking_var,file=paste(outputFigsPrefix,"_limma-blocking-variable_",quantitated_items_lbl,"Groups.txt",sep=""),sep="\t",row.names = T, col.names=NA)
  fit<-""
  
  if(.GlobalEnv[["LabelFree"]]){
    fitMethod <- 'robust'
  }else{
    fitMethod <- 'ls'
  }
  
  levellog("Fitting the model ...")
  if(.GlobalEnv[["n_bioreps"]] > 1 & .GlobalEnv[["n_techreps"]] > 1){
    # technical replication specification
    corfit <- duplicateCorrelation(t(norm.median.intensities), design=design, block = blocking_var, trim = duplicateCorrelation_trim)
    # Fit the limma model to the data
    # Pass the protein names/peptide sequences to limma as the genes option
    suppressWarnings(fit <- lmFit(t(norm.median.intensities), design, genes=prot.names, block = blocking_var, cor = corfit$consensus, method=fitMethod))
  }else{
    suppressWarnings(fit <- lmFit(t(norm.median.intensities), design, genes=prot.names, method=fitMethod))
  }
  
  # Setup contrast matrix
  # The contrast matrix specifies what comparisons we want to make between groups.
  # We pass the design matrix as the levels option to the makeContrasts function.
  
  levellog("Constructing the contrast matrix ...")
  ratio_combs<-combinations(nConditions,2,1:nConditions)
  contrasts<-c()
  for(i in 1:nrow(ratio_combs)){
    contrasts<-c(contrasts,paste(conditions.labels[ratio_combs[i,2]],"-",conditions.labels[ratio_combs[i,1]],sep=""))
  }
  contrasts <- makeContrasts(contrasts=contrasts, levels=design)
  write.table(contrasts,file=paste(outputFigsPrefix,"_limma-contrasts-matrix_",quantitated_items_lbl,"Groups.txt",sep=""),sep="\t",row.names = T, col.names=NA)  
  # Apply contrast matrix and do empirical bayes analysis to get p-values etc.
  
  levellog("Performing hypothesis testing ...")
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  
  # Plot a Histogram of co-efficients (log2 ratio)
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],"/",conditions.labels[ratio_combs[i,1]],sep="")
    hist(fit2$coefficients[,i],main=paste("Log2 Fold Change ",ratio_i_str,sep=""), xlab="Log2 Fold Change", breaks=50 )
  }   
  fit2.coefficients <<- fit2$coefficients
  if(exportFormat == "emf"){
    savePlot(filename=paste(outputFigsPrefix,"_limma-graphs_",time.point,"_hist.emf",sep=""),type="emf")
  }
  
  dev.off()
  
  
  #save the graphs to pngs
  png(paste("../",outputFigsPrefix,"_limma-graphs_",time.point,"_intensities-before-normalization.png",sep=""), width = 1500, height = 1050)
  boxplot(log.intensities)
  title(main="Intensities Before Normalisation")
  
  dev.off()
  
  png(paste("../",outputFigsPrefix,"_limma-graphs_",time.point,"_intensities-after-normalization.png",sep=""), width = 1500, height = 1050)
  boxplot(norm.intensities)
  title(main="Intensities After Normalisation")
  
  dev.off()
  
  tmp_conditions.labels <- conditions.labels
  if (IsobaricLabel)
  {
    if(!PDdata){
      conditions.labels <- sub("Reporter\\.intensity\\.", "Reporter.", conditions.labels)
    }else{
      conditions.labels <- sub('X([[:digit:]])', "\\1", conditions.labels)
    }
  }
  
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],"/",conditions.labels[ratio_combs[i,1]],sep="")
    ratio_i_str_with_dash <-paste(conditions.labels[ratio_combs[i,2]],"-",conditions.labels[ratio_combs[i,1]],sep="")
    png(paste("../",outputFigsPrefix,"_limma-graphs_",time.point,"_log2-fold-change-histogram", ratio_i_str_with_dash ,".png",sep=""), width = 1500, height = 1050)
    hist(fit2$coefficients[,i],main=paste("Log2 Fold Change ",ratio_i_str,sep=""), xlab="Log2 Fold Change", breaks=50 )
    dev.off()
  } 
  
  if (IsobaricLabel)
  {
    conditions.labels <- tmp_conditions.labels
  }
  
  # Output analysis details to file
  # adjust="BH" means adjust the calculated p-values for multiple testing using
  # the Benjamini Hochberg method (FDR)
  levellog("Saving analysis results to file ...")
  write.fit(fit2, file=paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), adjust="BH")
  
  # Notice that write.fit output will give us a table with average log2 intensity (A), fold changes for each comparison between conditions, p-values and adjusted p-values, F and F p values for each protein (the protein names are not displayed in this file)
  # Warning! latest releases of limma have changed the way the columns are displayed! this might not work in latest limma version! 3.18.13 suits best for PS
  if (packageVersion("limma") > "3.30.0" & nConditions == 2)
  {
    
    # Opening the results file might not be trivial, see the comment in do results plot for more infomation
    # and why it will be opened in such a complicated manner:
    
    FileHandle <- file(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""),"r")
    first_line <- readLines(FileHandle,n=1)
    hasLeadingTab <- (substr(first_line, 0, 1) == "\t")
    close(FileHandle)
    if (hasLeadingTab)
    {
      results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "", row.names = "X")
    } else {
      results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "")
    }
    
    colnames(results)[1] <- "A"
    colnames(results)[2] <- "Coef"
    colnames(results)[3] <- "t"
    colnames(results)[4] <- "p.value"
    colnames(results)[5] <- "p.value.adj"
    write.table(results, file=paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), sep = "\t", quote = FALSE)
  }
  return(norm.median.intensities)
}

read.pgroups_v3<-function(fname,evidence_fname,time.point,keepEvidenceIDs=F){
  levellog("",change=1)
  levellog("Reading data file ...");
  evidence<-read.table(evidence_fname, header = T, sep = "\t",quote="",stringsAsFactors=F,comment.char = "")
  
  if(PDdata) {
    if ('Protein.Group.Accessions' %in% colnames(evidence)) {
      pgroups_colname<-'Protein.Group.Accessions'
    }
    else if ('Protein.Accessions' %in% colnames(evidence)) {
      #pgroups_colname<-'Protein.Accessions'
      ###Ismini edit: Pattern must starts with Protein.Accessions because id PD 2.4 there is also another column as Master.Protein.Accessions
      pgroups_colname<-'^Protein.Accessions'
    } else {
      levellog("Error User: The dataset does not contain the columns 'Protein Group Accessions' or 'Protein Accessions'")
    }
  } else {
    pgroups_colname<-'^Proteins$'
  }
  colnames(evidence)[grepl(pgroups_colname,colnames(evidence))]<-'Protein.IDs'
  if(!PDdata){
    ## For MaxQuant correct protein groups in the evidence file using the protein groups file.
    pgroups<-read.table(fname, header = T, sep = "\t",quote="",stringsAsFactors=F,comment.char = "")
    # If there isn't a Protein.Names or Protein.names column (depends on MQ version), create one from the Fasta Headers column
    col_Protein.names <- length(grep('Protein.Names',colnames(pgroups))) > 0
    col_Protein.namesLC <- length(grep('Protein.names',colnames(pgroups))) > 0
    if(! col_Protein.names & ! col_Protein.namesLC)
    {
      pgroups$Protein.Names <- str_match(pgroups$Fasta.headers, '>[:alnum:]+[^[:alnum:]]+([^;>]+)')[,2]
    }
    if (col_Protein.namesLC)
    {
      colnames(pgroups)[grepl('Protein.names',colnames(pgroups))] <- 'Protein.Names'
    }
    # Construct a table, mapping the correct protein groups IDs (and the corresponding proteins names) to the evidence IDs
    #First check if there are any blank lines and remove them:
    mi<-which(pgroups$Evidence.IDs == "")
    if (length(mi)>0)
    {
      pgroups <- pgroups[-mi,]
    }
    if (!'Protein.IDs' %in% colnames(pgroups) & 'Peptide.IDs' %in% colnames(pgroups))
    {
      colnames(pgroups)[colnames(pgroups) == 'Peptide.IDs'] <- 'Protein.IDs'
    }
    tmp.table.1<-data.table(do.call(rbind, apply(pgroups[,c('Protein.IDs','Protein.Names','Evidence.IDs')], 1, function(x){return(cbind(x['Protein.IDs'], x['Protein.Names'], unlist(strsplit(x['Evidence.IDs'], ';'))))})))
    setnames(tmp.table.1, colnames(tmp.table.1), c('Protein.IDs', 'Protein.Names', 'id'))
    class(tmp.table.1$id)<-'integer'
    setkey(tmp.table.1, id)
    # Get the evidence records
    tmp.table.2<-data.table(evidence)
    setkey(tmp.table.2, id)
    # Remove the incorrect Protein.IDs column
    tmp.table.2[, c('Protein.IDs', 'Protein.Names') := NULL]
    # Inner join the mapping table with the evidence table and return the data frame that we ought to have in the first place
    evidence<-data.frame(tmp.table.1[tmp.table.2])
    if(length(grep("Contaminant", colnames(evidence))) > 0){
      evidence<-evidence[evidence$Contaminant == '', ]
    }
    if(length(grep("Reverse", colnames(evidence))) > 0){
      evidence<-evidence[evidence$Reverse == '', ]
    }
  }
  
  #In the case of Isobaric labeling we should reformat the table before proceeding, afterwards we will treat the data as
  #if they were label-free data
  if(IsobaricLabel)
  {
    if(!PDdata)
    {
      evidence$Intensity <- NULL
      varcolnames <- grep("^Reporter.intensity.[[:digit:]]", colnames(evidence), value = TRUE)
      evidence <- reshape(evidence, varying = varcolnames, v.names = "Intensity", timevar = "Labeling.State", times = varcolnames, direction = "long", new.row.names=sequence(prod(length(varcolnames), nrow(evidence))))
      conditions.labels<<-sub("^X", "Reporter.intensity.", conditions.labels)
      if (AllowLabelRename == T)
      {
        Rename_Array$old_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1", Rename_Array$old_label)
        Rename_Array$new_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1", Rename_Array$new_label)
      }
      if (AllowLS == T)
      {
        Ls_array$first_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  Ls_array$first_label)
        Ls_array$second_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  Ls_array$second_label)
      }
      LabelFree<-T;
      filterL_lbl <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  filterL_lbl)
      if(RMisused){
        RMtagsdata$name <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  RMtagsdata$name)
      }
    }
    else{
      evidence$Intensity <- NULL
      if (any(grepl("^Abundance..", colnames(evidence)))){
        colnames(evidence) <- sub("^Abundance..", "X", colnames(evidence))
      }
      varcolnames <- grep("^X[[:digit:]]+[[:alpha:]]?$", colnames(evidence), value = TRUE)
      evidence <- reshape(evidence, varying = varcolnames, v.names = "Intensity", timevar = "Modifications", times = varcolnames, direction = "long", new.row.names=sequence(prod(length(varcolnames), nrow(evidence))))
      LabelFree<-T;
      if (AllowLabelRename == T)
      {
        Rename_Array$old_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Rename_Array$old_label)
        Rename_Array$new_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Rename_Array$new_label)
      }
      if (AllowLS == T)
      {
        Ls_array$first_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Ls_array$first_label)
        Ls_array$second_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Ls_array$second_label)
      }
      filterL_lbl <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", filterL_lbl)
      if(RMisused){
        RMtagsdata$name <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", RMtagsdata$name)
      }
    }
  }
  # #Erase all rows in rename_array that try to rename a label to an already existing
  # mi<-which(Rename_Array$new_label %in% conditions.labels)
  # Rename_Array <- Rename_Array[-mi,]
  # Rename_Array <<- Rename_Array
  
  levellog(paste0("read.pgroups_v3: Identified proteins: ",length(unique(evidence$Protein.IDs))," (",time.point,")"))
  
  n1<-nrow(evidence)
  evidence<-evidence[nchar(evidence$Protein.IDs) > 0,]
  levellog(paste0("read.pgroups_v3: Discarded PSM records due to unassigned protein group: ",(n1-nrow(evidence))))
  ## Make Protein.IDs human-readable
  if(PDdata){
    pgroups_colname<-'Protein.Descriptions'
    if (!'Protein.Descriptions' %in% colnames(evidence)) {
      evidence[, c('Protein.Descriptions')] <- ""
    }
  }else{ pgroups_colname<-'Protein.Names' }
  tmp.table<-data.table(cbind(evidence[, c('Protein.IDs', pgroups_colname)], i=1:nrow(evidence)))
  setkey(tmp.table, Protein.IDs)
  # Generate data.table with unique Protein.IDs
  tmp.table2<-tmp.table[,.(n=.N),by=Protein.IDs]
  setkey(tmp.table2, Protein.IDs)
  # Make a new protein description column in other data.table
  tmp.table[, pdesc := paste0(Protein.IDs, ' [',strtrim(get(pgroups_colname), 50), ' ...]')]
  # set the Protein.IDs in the original data frame
  evidence$Protein.IDs<-tmp.table2[tmp.table][order(i),pdesc]
  ## Assign defined labels (conditions), one for each PSM record
  if(PDdata){ rawfile_col<-'Spectrum.File' }else{
    if(length(grep("Raw.File", colnames(evidence))) > 0)
    {
      rawfile_col<-'Raw.File'
    }
    else
    {
      rawfile_col<-'Raw.file'
    }
  }
  
  
  if(LabelFree){
    cond_spec_col<-rawfile_col
    if(IsobaricLabel){
      if(PDdata){ cond_spec_col<-'Modifications' }else{ cond_spec_col<-'Labeling.State' }
    }
  }else{
    if(PDdata){ cond_spec_col<-'Modifications' }else{ cond_spec_col<-'Labeling.State' }
  }
  
  if(RMisused){
    levellog("read.pgroups_v3: Transforming data for Replication Multiplexing ...")
    #when RM is chosen, in this line evidence has two main columns, one called Labeling.State or Modifications
    #that tells us what tag it was tagged and one called 
    #Spectrum File or Raw File that says from which raw file did the psm come from
    #in case of RM breps treps and conditions may come from either raw files or tags but in this data format creating
    #two new columns describing the structure correctly is not difficult
    #first create the column for conditions
    RMrawfilesdata <- RMrawfilesdata[!RMrawfilesdata$used == 'false',]
    RMtagsdata <- RMtagsdata[!RMtagsdata$used == 'false',]
    if(RMconditionsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'cond')], by.x = rawfile_col, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'cond')], by.x = cond_spec_col, by.y = 'name')
    }
    #the conditions.labels array is already set from the front end
    #Now we will initialize the new_raw_file column that will contain pseudo-raw files describing our bioreps, techreps and fracs
    colnames(RMrawfilesdata)[3:5] <- c('new_brep', 'new_trep', 'new_frac')
    colnames(RMtagsdata)[3:5] <- c('new_brep', 'new_trep', 'new_frac')
    if(RMbrepsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_brep')], by.x = rawfile_col, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'new_brep')], by.x = cond_spec_col, by.y = 'name')
    }
    #do the same thing for treps
    if(RMtrepsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_trep')], by.x = rawfile_col, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'new_trep')], by.x = cond_spec_col, by.y = 'name')
    }
    if(RMbrepsinrawfiles | RMtrepsinrawfiles) {
      #do the same thing for fracs
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_frac')], by.x = rawfile_col, by.y = 'name')
    }
    if(!RMbrepsinrawfiles & !RMtrepsinrawfiles) {
      evidence$new_raw_file <- paste0('b', evidence$new_brep, 't', evidence$new_trep)
    } else {
      evidence$new_raw_file <- paste0('b', evidence$new_brep, 't', evidence$new_trep, 'f', evidence$new_frac)
    }
    #Now lets refresh the rep_structure array the pseudo raw files we created are descriptive and contain the breps treps and fracs
    pseudo_raw_files <- unique(evidence$new_raw_file)
    new_rep_structure <- rep_structure
    new_rep_structure <- new_rep_structure[0,]
    for (i in 1:length(pseudo_raw_files)) {
      levels(new_rep_structure$raw_file) <- c(levels(new_rep_structure$raw_file), pseudo_raw_files[i])
      new_rep_structure[i,] <- c(as.character(pseudo_raw_files[i]), NA, NA, NA, NA)
    }
    colnames(new_rep_structure) <- c('raw_file','biorep','techrep','fraction', 'rep_desc')
    if(RMbrepsinrawfiles | RMtrepsinrawfiles) {
      for(i in 1:nrow(new_rep_structure))
      {
        new_rep_structure[i, c('raw_file', 'biorep','techrep','fraction')] <- str_match_all(new_rep_structure$raw_file[i], "b(.*?)t(.*?)f(.*)")[[1]][1,]
      }
    } else{
      for(i in 1:nrow(new_rep_structure))
      {
        new_rep_structure[i, c('raw_file', 'biorep','techrep')] <- str_match_all(new_rep_structure$raw_file[i], "b(.*?)t(.*)")[[1]][1,]
      }
      new_rep_structure[, 'fraction'] <- "1"
    }
    
    if(length(unique(new_rep_structure$techrep)) > 1){
      if(length(unique(new_rep_structure$fraction)) > 1){
        # we have techreps and fractions
        new_rep_structure$rep_desc<-paste(paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep=''),'f',new_rep_structure$fraction,sep='')
      }else{
        #we have bioreps and techreps
        new_rep_structure$rep_desc<-paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep='')
      }
    }else{
      if(length(unique(new_rep_structure$fraction)) > 1){
        # we have fractions but not techreps
        new_rep_structure$rep_desc<-paste(paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep=''),'f',new_rep_structure$fraction,sep='')
      }else{
        # we just have bioreps
        new_rep_structure$rep_desc<-paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep='')
      }
    }
    .GlobalEnv[["rep_structure"]]<-new_rep_structure
    #now erase the old columns for raw files and conditions and replace with the new ones
    evidence[, c(rawfile_col, cond_spec_col)] <- list(NULL)
    colnames(evidence)[colnames(evidence) == "new_raw_file"] <- rawfile_col
    colnames(evidence)[colnames(evidence) == "cond"] <- cond_spec_col
    if(length(unique(rep_structure$biorep)) == 1){
      levellog("Error User: Cannot accept dataset with just one biological replicate. Aborting ...")
      return(F)
    }
  }
  
  levellog("read.pgroups_v3: Assigning labels ...")
  levellog("",change=1)
  evidence$label_<-NA
  background_species_lbl<-NA
  
  for(i in 1:length(conditions.labels)){
    if(PDdata){
      if(LabelFree){
        if(!IsobaricLabel)
        {
          mi<-which(grepl(conditions.labels[i], LFQ_conds[, "condition"]))
          mi2<-which(grepl(paste(LFQ_conds[mi,]$raw_file, collapse="|"), evidence[, cond_spec_col]))
          evidence[mi2,]$label_<-conditions.labels[i]
        }
        else
        {
          evidence$label_<-evidence$Modifications
        }
      }else{
        evidence$label_<-evidence$Quan.Channel
      }
    }else{
      if(LabelFree){
        if(!IsobaricLabel){
          mi<-which(grepl(conditions.labels[i], LFQ_conds[, "condition"]))
          mi2<-which(grepl(paste(LFQ_conds[mi,]$raw_file, collapse="|"), evidence[, cond_spec_col]))
          evidence[mi2,]$label_<-conditions.labels[i]
        }
        else{
          evidence$label_<-evidence$Labeling.State
        }
        
      }else{
        # MQ nomenclature for labels: 0 the first label, 1 the second etc ...
        # Since the user might opted for excluding some labels
        # find the index of the included labels and parse them to the label_ column
        mi0<-which(All_MQ_Labels == conditions.labels[i])
        mi<-which(grepl((mi0[1]-1), evidence[, cond_spec_col]))
        evidence[mi,]$label_<-conditions.labels[i]
      }
    }
    levellog(paste0("read.pgroups_v3: Assigned label '", conditions.labels[i],"'."))
  }
  #Rename any labels if necessary
  if (AllowLabelRename == T)
  {
    if (length(Rename_Array$old_label) != 0)
    {
      for(i in 1:length(Rename_Array$old_label))
      {
        if(Rename_Array$old_label[i] != Rename_Array$new_label[i])
        {
          #The case where old label = new label is common since if the user did not ask for a label rename (merge) for a label
          #this label is sent to be renamed in R to itself
          #in any other case rename the labels that are merged to the same label so that they become indistinguishable
          #and refresh the conditions labels by erasing the old label and adding the new if necessary
          mi<-which(evidence$label_ == Rename_Array$old_label[i])
          if (LabelFree == FALSE & IsobaricLabel == FALSE)
          {
            #in case of precursor ion data add you need to create a new column with the new cond name containing
            #the intensities of the old labels to the line where the label is found:
            prefix<-NA
            if (PDdata)
            {
              prefix<-""
            }
            else
            {
              prefix<-"Intensity."
            }
            evidence[mi, paste0(prefix, Rename_Array$new_label[i])] <- evidence[mi,  paste0(prefix, Rename_Array$old_label[i])]
          }
          evidence$label_[mi] <- Rename_Array$new_label[i]
          mi<-which(conditions.labels == Rename_Array$old_label[i])
          if (length(mi) > 0)
          {
            conditions.labels <- conditions.labels[-mi]
          }
          mi<-which(conditions.labels == Rename_Array$new_label[i])
          if (length(mi) == 0)
          {
            conditions.labels <- c(conditions.labels, Rename_Array$new_label[i])
          }
        }
      } 
      conditions.labels <<- conditions.labels
      nConditions<<-length(conditions.labels)
    }
    #No condition can be named "N" in ProteoSign so take care of this scarce situation:
    for(i in 1:length(conditions.labels)){
      if(conditions.labels[i] == "N")
      {
        conditions.labels[i] <- "condN"
        conditions.labels <<- conditions.labels
        evidence$label_ <- sub("^N$", "condN", evidence$label_)
      }
    }
  }
  
  #Rename again the labels in case of label swapping
  if (AllowLS == T)
  {
    for(i in 1:length(Ls_array$first_label))
    {
      mi1<-which(evidence[, rawfile_col] == Ls_array$selected_raw_file[i] & evidence$label_ == Ls_array$first_label[i])
      mi2<-which(evidence[, rawfile_col] == Ls_array$selected_raw_file[i] & evidence$label_ == Ls_array$second_label[i])
      evidence$label_[mi1] <- as.character(Ls_array$second_label[i])
      evidence$label_[mi2] <- as.character(Ls_array$first_label[i])
      if (LabelFree == FALSE & IsobaricLabel == FALSE)
      {
        prefix<-NA
        if (PDdata)
        {
          prefix<-""
        }
        else
        {
          prefix<-"Intensity."
        }
        #in case of precursor ion we should also swap the values between the columns "Intensity" of the first and second label
        if (length(mi1) > 0)
        {
          evidence[mi1, "temp_LS_Intensities"] <- evidence[mi1, paste0(prefix, Ls_array$second_label[i])]
          evidence[mi1, paste0(prefix, Ls_array$second_label[i])] <- evidence[mi1, paste0(prefix, Ls_array$first_label[i])]
          evidence[mi1, paste0(prefix, Ls_array$first_label[i])] <- evidence[mi1, "temp_LS_Intensities"]
        }
        if (length(mi2) > 0)
        {
          evidence[mi2, "temp_LS_Intensities"] <- evidence[mi2, paste0(prefix, Ls_array$second_label[i])]
          evidence[mi2, paste0(prefix, Ls_array$second_label[i])] <- evidence[mi2, paste0(prefix, Ls_array$first_label[i])]
          evidence[mi2, paste0(prefix, Ls_array$first_label[i])] <- evidence[mi2, "temp_LS_Intensities"]
        }
      }
    }
  }
  levellog("",change=-1)
  mi<-which(is.na(evidence$label_))
  if(is.na(background_species_lbl)){
    if(length(mi) > 0){
      evidence<-evidence[-mi,]
      levellog(paste("read.pgroups_v3: Discarded PSM records due to unassigned label: ",length(mi),sep=""))
    }
  }else{
    evidence[mi,]$label_<-background_species_lbl
  }
  # Now add the experimental structure information
  evidence<-merge(evidence, .GlobalEnv[["rep_structure"]], by.x=c(rawfile_col), by.y=c('raw_file'))
  new_cond_labels <- NULL
  for (cond_i in conditions.labels)
  {
    if (!(cond_i %in% evidence$label_))
    {
      levellog(paste0("Warn User: ", cond_i, " label was not found in the selected raw files and was not used in comparisons!"))
      if(filterL_lbl == cond_i)
      {
        filterL<-F
        levellog("Warning!: the filter label was not found in the selected raw files. Filtering will not take place!")
      }
    }
    else
    {
      new_cond_labels <- c(new_cond_labels, cond_i)
    }
  }
  if(length(new_cond_labels)>1)
  {
    conditions.labels <- new_cond_labels
    conditions.labels <<- conditions.labels
    nConditions<<-length(conditions.labels)
  }else{
    levellog(paste0("Error User: Not enough labels left, aborting..."))
    return(F)
  }
  ## If we have fractionation, remake the rep_desc column and don't take into account the fraction number
  if(length(unique(.GlobalEnv[["rep_structure"]]$fraction)) > 1){
    evidence$rep_desc <- paste0('b',evidence$biorep,'t',evidence$techrep)
  }
  
  ## Generate Venn data for the identified proteins and output to a file
  levellog("read.pgroups_v3: Generating ID Venn data ...")
  tmp.table<-data.table(evidence[, c('Protein.IDs', 'biorep', 'techrep', 'fraction')])
  setkey(tmp.table,  Protein.IDs, biorep, techrep, fraction)
  setwd(limma_output)
  write.table(tmp.table[, .(n=.N), by=.(Protein.IDs,rep=biorep)][,.(Protein.IDs,rep)],file=paste0(outputFigsPrefix,"_id_venn3-data_",time.point,".txt"),sep="\t",row.names=F)
  setwd("..")    
  
  # Bring Labeled or Label-free data to the following common format (table headers):
  # rep_desc Protein.IDs UniqueSequences.Intensity.condition_1 ... UniqueSequences.Intensity.condition_N Intensity.condition_1 ... Intensity.condition_N
  
  levellog("read.pgroups_v3: Standarizing data format ...")
  if(!PDdata){
    colnames(evidence)[grepl('Peptide.ID',colnames(evidence))]<-'Unique.Sequence.ID'
    if (!IsobaricLabel){
      # colnames(evidence)[grepl('Intensity\\..+',colnames(evidence))]<-conditions.labels
      colnames(evidence) <- sub('Intensity\\.(.+)', "\\1", colnames(evidence))
    }
    # else{
    # evidence[,conditions.labels]<-NA
    # for (my_cond in conditions.labels){
    #   mi<-which(grepl(my_cond, evidence$Labeling.State))
    #   evidence[mi, my_cond] <- evidence[mi, "Intensity"]
    # }
    # }
    
  }
  
  if (!'Unique.Sequence.ID' %in% colnames(evidence)){
    if ('Annotated.Sequence' %in% colnames(evidence)){
      colnames(evidence)[colnames(evidence) == 'Annotated.Sequence'] <- 'Unique.Sequence.ID'
      evidence$Unique.Sequence.ID <- sub(".*?\\.(.*?)\\..*", "\\1", evidence$Unique.Sequence.ID)
    }
  }
  
  #Here we have a column containing a unique sequence ID or the unique sequence of a peptide i.e. all PSMs that correspond
  #to the same peptide will have the same identifier in this column in evidence data table
  
  #Below we will buid the evidence.dt data table so that it contains the following information per PROTEIN:
  if(LabelFree){
    if(PDdata){
      # Precursor Area is unfortunately buggy (sometimes 0/NA), so we are left with Intensity to work with
      #intensityCol <- 'Precursor.Area'
      #intensityCol <- 'Intensity'
      #Ismini edit: in PD 2.4 there is no Intensity column, instead there is the "PRECURSOR.ABUNDANCE" column
      #intensityCol <- 'Intensity'
      if(any(grepl('Intensity',colnames(evidence)))){
        intensityCol <- 'Intensity'
      }else if(any(grepl('Precursor.Abundance',colnames(evidence)))){
        intensityCol <- 'Precursor.Abundance'
      }else{
        levellog("Error User: The dataset does not contain the columns 'Intensity' or 'Precursor Abundance'")
      }
      
      ###Ismini edit: Must handle PD's ver.2.4 colnames. There is no "Unique Sequence ID", instead there is a "PSMs Peptide ID" column name
      if(any(grepl('Peptide.ID',colnames(evidence)))){
        colnames(evidence)[grepl('Peptide.ID',colnames(evidence))]<-'Unique.Sequence.ID'
      }
    }else{
      intensityCol <- 'Intensity'
    }
    # Retrieve the following information for all PSMs from evidence: Protein ID, Unique sequence ID (so one column for protein and one
    # for peptide information), its Intensity, the condition it derives from (label_) and the repetition it derived from (rep_desc)
    evidence.dt<-data.table(evidence[, c('Protein.IDs', 'Unique.Sequence.ID', intensityCol,'label_', 'rep_desc')])
    setkey(evidence.dt, rep_desc, Protein.IDs, Unique.Sequence.ID, label_)
    # In case a peptide has been identified in a single MS run more than one time this means that during LC elution it has been detected multiple times
    # get the maximum of these intensities as the most representative quantification value for these peptides
    # Get maximum PSM intensity per peptide/protein/[(rep_desc/label) = raw_file]
    suppressWarnings(evidence.dt<-evidence.dt[, .(maxI=as.double(max(get(intensityCol), na.rm = T))), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)][maxI != -Inf])
  }else{
    if(PDdata){
      #get only the PSMs that PD suggested as eligible for quantification
      if (!'Quan.Usage' %in% colnames(evidence) & 'Peptide.Quan.Usage' %in% colnames(evidence)){
        colnames(evidence)[colnames(evidence) == 'Peptide.Quan.Usage'] <- 'Quan.Usage'
      }
      ###Ismini edit: Must handle PD's ver.2.4 colnames. There is no "Unique Sequence ID", instead there is a "PSMs Peptide ID" column name
      if(any(grepl('Peptide.ID',colnames(evidence)))){
        colnames(evidence)[grepl('Peptide.ID',colnames(evidence))]<-'Unique.Sequence.ID'
      }
      # Retrieve the following information for all PSMs from evidence: Protein ID, Unique sequence ID (so one column for protein and one
      # for peptide information), its Intensity for each condition seperately, the condition it derives from (label_) and the repetition it derived from (rep_desc)
      # For PD data also retrieve the quantification info column describing the eligibility for quantification
      evidence.dt<-data.table(evidence[, c('Quan.Usage','Protein.IDs', 'Unique.Sequence.ID', conditions.labels,'rep_desc', 'label_')])
    }else{
      # Retrieve the following information for all PSMs from evidence: Protein ID, Unique sequence ID (so one column for protein and one
      # for peptide information), its Intensity for each condition seperately, the condition it derives from (label_) and the repetition it derived from (rep_desc)
      evidence.dt<-data.table(evidence[, c('Protein.IDs', 'Unique.Sequence.ID', conditions.labels,'rep_desc', 'label_')])
    }
    setkey(evidence.dt, rep_desc, Protein.IDs, Unique.Sequence.ID)    
  }
  
  # Here evidence.dt contains PSMs as rows. Each row contains info that matches the respective PSM to a specific peptide, to each most probable parent protein
  # and to each intensity.
  
  ## Calculate identified peptide counts per protein for each condition/label and replicate in the following three steps
  # 1. For each condition (per sequnce, protein and replicate), set a corresponding column to TRUE if there are > 0 evidence.dt (PSMs) records, FALSE otherwise
  evidence.dt.seqCounts<-dcast.data.table(evidence.dt[, .(n=.N > 0), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)], rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=FALSE)
  # 2. Add a column flagging the common, between conditions/labels, sequences.
  # In case of more than two conditions/labels, the flag designates that there are at least two conditions/labels where the peptide is common
  evidence.dt.seqCounts[, 'common' := rowSums(.SD) > 1,.SDcols=conditions.labels]    
  # 3. Collapse the records for each protein (per replicate) and count the TRUEs.
  # evidence.dt[, .(n.Unique.Sequence.IDs=.N), by=.(rep_desc, Protein.IDs)]
  evidence.dt.seqCounts<-evidence.dt.seqCounts[,c(n.Unique.Sequence.IDs=.N,lapply(.SD, function(x){return(length(which(x)))})), by=.(rep_desc,Protein.IDs),.SDcols=c(conditions.labels, 'common')]
  # 4. Calculate the percentage columns
  evidence.dt.seqCounts[, paste0(conditions.labels,'p') := lapply(.SD, function(x){return((x/sum(.SD))*100)}), by=.(rep_desc,Protein.IDs),.SDcols=c(conditions.labels)]
  ## Rename the peptide counts columns
  setnames(evidence.dt.seqCounts,colnames(evidence.dt.seqCounts)[which(colnames(evidence.dt.seqCounts) %in% conditions.labels)],paste('UniqueSequences',conditions.labels,sep='.'))    
  ## Calculate the protein intensity = (sum of unique peptide intensities) for each condition/label and replicate in the following two steps
  if(LabelFree){
    # 1. Cast the data so that we have columns for each label and intensity separately
    evidence.dt<-dcast.data.table(evidence.dt, rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=0)    
  }else{
    if(PDdata){
      # 1. Take the (Quan.Usage == 'Used') records and for each peptide keep only the PSM record with the highest intensity
      evidence.dt<-evidence.dt[Quan.Usage == 'Used' | Quan.Usage == 'Use', lapply(.SD, max), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID), .SDcols=conditions.labels]    
    }else{
      # 2. Take the records with Intensity != NA across labels/conditions and for each peptide keep only the PSM record with the highest intensity
      evidence.dt[, sumI := rowSums(.SD, na.rm = T), .SDcols=conditions.labels]
      evidence.dt<-evidence.dt[sumI > 0, lapply(.SD, max), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID), .SDcols=conditions.labels]    
      evidence.dt[, sumI := NULL]
    }
  }
  
  # Get a vector of unique peptides intensities
  tmp.I<-sort(unique(evidence.dt[,get(conditions.labels)]))
  # If the minimum intensity is zero
  if(tmp.I[1] == 0){
    # Replace 0's with minimum intensity (PD can do this automatically for us)
    minI<-tmp.I[2]
    evidence.dt[, (conditions.labels) := lapply(.SD, function(x){ t<-which(x == 0); if(length(t) > 0){x[t] <- minI}; return(x) }), .SDcols=conditions.labels]
  }else{
    minI<-tmp.I[1]
  }
  
  ## If enabled, do filter out peptides where all 'channels' except filterL_lbl channel have noise-level intensity (peptide-level filtering)
  if(filterL && filterL_lvl){
    evidence.dt[, minIcount := rowSums(.SD == minI), .SDcols=conditions.labels[! conditions.labels %in% filterL_lbl]]
    n1<-nrow(evidence.dt)
    evidence.dt<-evidence.dt[minIcount < (nConditions - 1)]
    n2<-nrow(evidence.dt)
    if(n2 < n1){
      levellog(paste0("read.pgroups_v3: Filtered out ", (n1-n2)," peptides having noise-level intensity in all channels except the '", filterL_lbl,"' channel ..."));
    }
    evidence.dt[, minIcount := NULL]
  }
  
  # 2. Calculate the protein intensity (= sum of unique peptide intensities) for each condition/label and replicate
  # Also count the number of quantifiable peptides (those which do not have intensity NA)
  if(LabelFree){
    # Top three in abundance
    #evidence.dt<-evidence.dt[, lapply(.SD, function(x){x<-x[!is.na(x)]; x<-sort(x, decreasing<-T); if(length(x)<3){return(sum(x))}else{return(sum(x[1:3]))}}), by=.(rep_desc, Protein.IDs), .SDcols=conditions.labels]
    evidence.dt<-evidence.dt[, c(n=.N, nas=length(which(is.na(.SD))) ,lapply(.SD, function(x){x<-x[!is.na(x)]; x<-sort(x, decreasing<-T); if(length(x)<3){return(sum(x))}else{return(sum(x[1:3]))}})), by=.(rep_desc, Protein.IDs), .SDcols=conditions.labels]
  }else{
    # All peptides
    evidence.dt<-evidence.dt[, c(n=.N, nas=length(which(is.na(.SD))) ,lapply(.SD, sum, na.rm = T)), by=.(rep_desc, Protein.IDs), .SDcols=conditions.labels] 
  }
  ## Rename the intensity columns
  setnames(evidence.dt,colnames(evidence.dt)[which(colnames(evidence.dt) %in% conditions.labels)],paste('Intensity',conditions.labels,sep='.'))
  ## Merge with the evidence.dt.seqCounts table
  evidence.dt<-merge(evidence.dt, evidence.dt.seqCounts)
  
  # Add the experimental structure information to evidence.dt based on rep_desc (raw file at this point has no information and is dropped)
  tmp.rep_struct<-.GlobalEnv[["rep_structure"]][! duplicated(.GlobalEnv[["rep_structure"]][,c('biorep','techrep')]), !grepl('raw_file', colnames(.GlobalEnv[["rep_structure"]])) & !grepl('fraction', colnames(.GlobalEnv[["rep_structure"]]) )]
  tmp.rep_struct$rep_desc<-paste0('b',tmp.rep_struct$biorep,'t',tmp.rep_struct$techrep)
  evidence.dt<-merge(evidence.dt ,data.table(tmp.rep_struct), by='rep_desc')
  
  ## If enabled, do filter out proteins based on percentage labeling for the desired label (protein-level filtering)
  if(filterL && !filterL_lvl){
    n1<-length(unique(evidence.dt[get(paste0(filterL_lbl,"p")) == 100.0]$Protein.IDs))
    evidence.dt<-evidence.dt[get(paste0(filterL_lbl,"p")) < 100.0]
    levellog(paste0("read.pgroups_v3: Filtered out ", n1," proteins which where identified solely by '", filterL_lbl, "' peptides ..."));
  }
  
  ## Get protein IDs that were quantified with a total of at least 'nRequiredLeastBioreps' different peptides accross at least 'nRequiredLeastBioreps' biological replicates.
  # E.g. 1: with 3 biological replicates, a protein that was quantified by a single peptide in 'nRequiredLeastBioreps' out of the 3 replicates will be discarded if 'nRequiredLeastBioreps' > 1 (retained otherwise).
  # E.g. 2: with 3 biological replicates, a protein that was quantified by a single peptide in 1 out of the 3 replicates will be discarded if 'nRequiredLeastBioreps' > 1 (retained otherwise).
  # E.g. 3: with 3 biological replicates, a protein that was quantified by two peptides in at one of replicates will be discarded if 'nRequiredLeastBioreps' > 1 (retained otherwise).
  # E.g. 4: with 3 biological replicates, a protein that was quantified by two peptides (in total) in 2 out of the 3 replicates will be discarded if 'nRequiredLeastBioreps' > 2 (retained otherwise).
  
  Protein.IDs.quant <- evidence.dt[, .(c1 = sum(N.x-nas)) , by=.(Protein.IDs, biorep)][, .(nQuantPeps = sum(c1), geqXnRequiredLeastBioreps = .N >= .GlobalEnv[["nRequiredLeastBioreps"]]), by=.(Protein.IDs)][nQuantPeps >= .GlobalEnv[["nRequiredLeastPeps"]] & geqXnRequiredLeastBioreps == T]$Protein.IDs
  levellog(paste0("read.pgroups_v3: Filtered out ", (length(unique(evidence.dt$Protein.IDs)) - length(Protein.IDs.quant))," proteins which were not identified in at least ",nRequiredLeastBioreps," biological replicate(s) with at least a total of ",nRequiredLeastPeps," peptide(s)"));
  evidence.dt[,nQuantPeps := N.x-nas]
  evidence.dt<-evidence.dt[Protein.IDs %in% Protein.IDs.quant]
  
  ## Experimental filter based on outlier removal (grubbs method) based on the first condition specified.
  ## NOTE: It is applied when there are no technical replicates in Label-free data, where variability is expected to be very high.
  # If a protein intensity in condition i and biological replicate j is found to be an outlier based on the distribution
  # of intensities from all biological replicates, then
  # the biological replicate j is removed for that particular protein for all conditions.
  #if(LabelFree && .GlobalEnv[["n_techreps"]] < 2){
  #  evidence.dt.bad <- suppressWarnings(evidence.dt[, lapply(.SD, function(x){p.val = grubbs.test(x)$p.value; if(!is.na(p.val) && p.val < 0.05){outlier.true <- T}else{outlier.true <- F}; if(outlier.true){return(.I[outlier(x, logical=T)][1] )}else{return(as.integer(0))} }),by=.(Protein.IDs),.SDcols=paste0('Intensity.',conditions.labels[1])][,get(paste0('Intensity.',conditions.labels[1]))])
  #  evidence.dt.bad <- evidence.dt.bad[evidence.dt.bad > 0]
  #  evidence.dt<-evidence.dt[! evidence.dt.bad]
  #  levellog(paste0("read.pgroups_v3: Filtered out ", length(evidence.dt.bad)," protein intensities based on outlier detection on condition '",conditions.labels[1],"'."));
  #  Protein.IDs.quant <- evidence.dt[, .(c1 = sum(n-nas)) , by=.(Protein.IDs, biorep)][, .(nQuantPeps = sum(c1), geqXnRequiredLeastBioreps = .N >= .GlobalEnv[["nRequiredLeastBioreps"]]), by=.(Protein.IDs)][nQuantPeps >= .GlobalEnv[["nRequiredLeastBioreps"]] & geqXnRequiredLeastBioreps == T]$Protein.IDs
  #  levellog(paste0("read.pgroups_v3: Filtered out another ", (length(unique(evidence.dt$Protein.IDs)) - length(Protein.IDs.quant))," proteins which were not identified in at least ",nRequiredLeastBioreps," biological replicate(s) with at least a total of ",nRequiredLeastBioreps," peptide(s)"));
  #  evidence.dt[,nQuantPeps := n-nas]
  #  evidence.dt<-evidence.dt[Protein.IDs %in% Protein.IDs.quant]
  #}
  
  
  
  ## Generate Venn data for the identified proteins and output to a file
  levellog("read.pgroups_v3: Generating quant Venn data ...")
  setwd(limma_output)  
  write.table(evidence.dt[, .(Protein.IDs, rep=biorep)],file=paste0(outputFigsPrefix,"_quant_venn3-data-",.GlobalEnv[["nRequiredLeastBioreps"]],"reps_",time.point,".txt"),sep="\t",row.names=F)
  setwd("..")
  
  ## Cast the table to the following format
  # Protein.IDs Intensity.[<rep_desc_X>.<label/condition_Y> ...] [<rep_desc_X>.Ratio.counts ...] [<rep_desc_X>.uniqueSequences ...] time.point [<label/condition_Y> ...] [<label/condition_Y>p ...]
  
  ## Step 1: For each 'rep_desc', add to a growing dataframe the evidence.dt data, renaming the columns accordingly
  # Also, calculate the missing columns required by the target format and drop the unnecessary columns
  setkey(evidence.dt, Protein.IDs)
  pgroups<-data.frame(Protein.IDs = unique(evidence.dt$Protein.IDs))
  setkey(evidence.dt, rep_desc)
  for(rep_desc_i in unique(evidence.dt$rep_desc)){
    rep_desc_i_pgroups<-data.frame(evidence.dt[rep_desc == rep_desc_i,])
    allcols<-colnames(rep_desc_i_pgroups)
    # Rename Intensity cols
    colsl<-grepl('^Intensity' ,allcols)
    colnames(rep_desc_i_pgroups)[colsl]<-gsub("^Intensity(.+)$",paste("Intensity\\1",rep_desc_i,sep='.'), allcols[colsl])
    # Rename UniqueSequences cols
    colsl<-grepl('^UniqueSequences' ,allcols)
    colnames(rep_desc_i_pgroups)[colsl]<-gsub("^UniqueSequences(.+)$",paste(rep_desc_i,"uniqueSequences\\1",sep='.'), allcols[colsl])
    # Add new column <rep_desc_X>.uniqueSequences
    rep_desc_i_pgroups[, paste(rep_desc_i,'uniqueSequences',sep='.')]<-rowSums(rep_desc_i_pgroups[, colnames(rep_desc_i_pgroups)[colsl]])
    # Rename 'p' (percentage) cols
    colsl<-allcols %in% paste0(conditions.labels,'p')
    colnames(rep_desc_i_pgroups)[colsl]<-gsub("^(.+)$",paste("\\1",rep_desc_i,sep='.'), allcols[colsl])
    # Rename the 'nQuantPeps' column to <rep_desc_i>.Ratio.counts
    colsl<-allcols %in% c('nQuantPeps')
    colnames(rep_desc_i_pgroups)[colsl]<-paste(rep_desc_i,'Ratio.counts',sep='.')
    # merge with the growing data frame
    cc<-intersect(names(pgroups), names(rep_desc_i_pgroups))
    pgroups<-merge(pgroups, rep_desc_i_pgroups[, ! colnames(rep_desc_i_pgroups) %in% c('biorep', 'techrep', 'fraction', 'rep_desc', cc[! grepl('Protein.IDs', cc)] )], all.x = T)
  }
  # Step 2: Calculate the columns [<label/condition_Y> ...] containing the number of unique sequences found per condition in all replicates
  allcols<-colnames(pgroups)
  for(cond_i in conditions.labels){
    colsl<-grepl(paste('uniqueSequences', cond_i,sep='\\.') ,allcols)
    pgroups[, cond_i]<-rowSums(pgroups[, allcols[colsl]])
  }
  # Step 3: Calculate the columns [<label/condition_Y>p ...] containing the percentage of unique sequences that were found in a specific condition in all replicates
  allcols<-colnames(pgroups)
  for(cond_i in conditions.labels){
    colsl<-allcols %in% conditions.labels & ! allcols %in% cond_i
    pgroups[, paste0(cond_i,'p')]<-(pgroups[, cond_i]/rowSums(pgroups[, c(cond_i, allcols[colsl])]))*100
  }
  # Step 4: Add time-point column
  pgroups$time.point <- time.point
  # Step 5: in case there is replicate mismatch between conditions i.e. there is at least one replicate of the experiment that contains quantification values not from all conditions, there are some columns in pgroups at the moment that have the format <replicate_description>.uniqueSequences.<condition> e.g. b1t1.uniqueSequences.WildType that contain 0 values solely. These columns should be emmited and the mismatch must be taken into account while sending the data to limma
  replicate_mismatch <<- F
  # We have replicate mismatch in case not all conditions correspond to the same experimental replicates
  for (rep_desc_i in unique(evidence.dt$rep_desc)) {
    for (cond_i in conditions.labels) {
      # For each of the columns of interest check if it contains solely 0 values and if so delete the respective intensity column
      #if (all(pgroups[, paste0(rep_desc_i, ".uniqueSequences.", cond_i)] == 0)) {
      if (all(pgroups[, paste0(rep_desc_i, ".uniqueSequences.", cond_i)] == 0, na.rm = TRUE) || all(is.nan(pgroups[, paste0(rep_desc_i, ".uniqueSequences.", cond_i)]))) {
        
        allcols <- colnames(pgroups)
        pgroups <- pgroups[, - which(grepl(paste0("Intensity", ".", cond_i, ".", rep_desc_i), allcols))]
        replicate_mismatch <<- T
      }
    }
  }
  # Step 6: Remove unnecessary columns (uniqueSequences per rep_desc and percentage unique peptides per rep_desc)
  allcols<-colnames(pgroups)
  pgroups<-pgroups[,-which(grepl('uniqueSequences\\.', allcols) | grepl('p\\.b[0-9]+t[0-9]+$', allcols) | grepl('^common$', allcols))]
  ##
  levellog(paste0("read.pgroups_v3: Quantifiable proteins: ", nrow(pgroups)," (",time.point,")"))
  levellog("",change=-1)
  ## 
  return(pgroups)  
}

#Prepare protein intensity table for differential expression analysis (the format limma requires)
prepare_working_pgroups<-function(working_pgroups){
  rownames(working_pgroups)<-working_pgroups[,paste(quantitated_items_lbl,".IDs",sep="")]
  inten_cols<-c()
  for(cond_i in conditions.labels){
    inten_cols<-c(inten_cols,sort(colnames(working_pgroups)[grep(paste("Intensity.",cond_i,".b",sep=""),colnames(working_pgroups))]))
  }  
  working_pgroups<-working_pgroups[,inten_cols]
  colnames(working_pgroups)<-sub("Intensity\\.","",inten_cols)
  return(working_pgroups)
}

#GLOBAL variables

duplicateCorrelation_trim<-0.15 # use 0.22 for "bad" datasets (too many missing values)

addLabel<-function(lblname, lbl.Modifications){
  levellog("", change=1)
  #If label name is a number some routines won't work, it has to be converted to some acceptable variable name
  lblname<-make.names(lblname)
  labeltxt <- "label";
  if(!LabelFree){
    lbl.Modifications<-gsub("\\(","\\\\(",lbl.Modifications)
    lbl.Modifications<-gsub("\\)","\\\\)",lbl.Modifications)
    unmod_idx<-which(lbl.Modifications == "")
    if(length(unmod_idx) > 0){
      rest_idx<-which(lbl.Modifications[-unmod_idx] != "")
      if(length(rest_idx) > 0){
        lbl.Modifications<-c(lbl.Modifications[unmod_idx],paste(lbl.Modifications[-unmod_idx],"\\)",sep=""))
      }
    }else{
      lbl.Modifications<-paste(lbl.Modifications,"\\)",sep="")
    }
  }else{
    labeltxt <- "condition";
  }
  lblname_i<-which(grepl(paste("^",lblname,"$",sep=""),conditions.labels))
  if(length(lblname_i) != 0){
    levellog(paste("addLabel: Error adding ",labeltxt," '",lblname,"': An existing ",labeltxt," with name '",lblname,"' (specification: ",paste(unlist(conditions.labels.Modifications[lblname_i]),collapse=", "),") already exists. Please try a different name.",sep=""), change=-1)
    return(FALSE)
  }
  
  conditions.labels<<-c(conditions.labels, lblname)
  j<-length(conditions.labels.Modifications)+1
  conditions.labels.Modifications[[j]]<<-lbl.Modifications
  nConditions<<-length(conditions.labels)
  levellog("", change=-1)
}
removeLabel<-function(lblname){
  levellog("", change=1)
  lblname_i<-which(grepl(paste("^",lblname,"$",sep=""),conditions.labels))
  if(length(lblname_i) != 0){
    conditions.labels<<-conditions.labels[-lblname_i]
    conditions.labels.Modifications<<-conditions.labels.Modifications[-lblname_i]
  }
  nConditions<<-length(conditions.labels)
  levellog("", change=-1)
}
clearLabels<-function(){
  levellog("", change=1)
  conditions.labels<<-c()
  conditions.labels.Modifications<<-list()
  nConditions<<-length(conditions.labels)
  levellog("", change=-1)
}

addMod<-function(modname, mod.Modifications){
  levellog("", change=1)
  modname<-make.names(modname)
  mod.Modifications<-gsub("\\(","\\\\(",mod.Modifications)
  mod.Modifications<-gsub("\\)","\\\\)",mod.Modifications)
  unmod_idx<-which(mod.Modifications == "")
  if(length(unmod_idx) > 0){
    rest_idx<-which(mod.Modifications[-unmod_idx] != "")
    if(length(rest_idx) > 0){
      mod.Modifications<-c(mod.Modifications[unmod_idx],paste(mod.Modifications[-unmod_idx],"\\)",sep=""))
    }
  }else{
    mod.Modifications<-paste(mod.Modifications,"\\)",sep="")
  }
  modname_i<-which(grepl(paste("^",modname,"$",sep=""),conditions.Mods))
  if(length(modname_i) != 0){
    levellog(paste("addMod: Error adding modification '",modname,"': An existing modification with name '",modname,"' (specification: ",paste(unlist(conditions.Mods.Modifications[modname_i]),collapse=", "),") already exists. Please try a different name.",sep=""), change=-1)
    return(FALSE)
  }
  i<-1
  for(mod_i in conditions.Mods.Modifications){
    for(mod_i in mod.Modifications){
      mod_i_matches<-which(grepl(paste("^",mod_i,"$",sep=""),mod_i))
      if(length(mod_i_matches) != 0){
        levellog(paste("addMod: Error adding modification '",modname,"': Existing modification with name '",conditions.Mods[i],"' has an identical specification (",mod_i,").",sep=""),change=-1)
        return(FALSE)
      }
    }
    i<-i+1
  }
  conditions.Mods<<-c(conditions.Mods, modname)
  j<-length(conditions.Mods.Modifications)+1
  conditions.Mods.Modifications[[j]]<<-mod.Modifications
  nMods<<-length(conditions.Mods)
  levellog("", change=-1)
}
removeMod<-function(modname){
  levellog("", change=1)
  modname_i<-which(grepl(paste("^",modname,"$",sep=""),conditions.Mods))
  if(length(modname_i) != 0){
    conditions.Mods<<-conditions.Mods[-modname_i]
    conditions.Mods.Modifications<<-conditions.Mods.Modifications[-modname_i]
  }
  nMods<<-length(conditions.Mods)
  levellog("", change=-1)
}
clearMods<-function(){
  levellog("", change=1)
  conditions.Mods<<-c()
  conditions.Mods.Modifications<<-list()
  nMods<<-length(conditions.Mods)
  levellog("", change=-1)
}

# the following functions are for GO analysis


# Pattern of uniprot IDs
my_grep <- function(x){
  #uniprot IDs pattern
  grep('^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',x, value=TRUE)
}

# Get the UNIPROT IDs from the protein names
get_uniprot_ids <- function(results, cond1, cond2){
  
  #get p values associated to DE proteins
  col_desc_<-paste("P-value adjusted ",paste(cond2,"/",cond1,sep=""),sep="")
  
  # Find DE proteins in results data frame
  ind_diffexp_tmp<-which(results[,col_desc_]<pThreshold)
  DE_prot <- rownames(results)[ind_diffexp_tmp]
  
  if (length(DE_prot) == 0)
  {
    uniprot_ids <- c()
    return(uniprot_ids)
  }
  
  uniprot_ids <- c()
  
  # Get the uniprot ID for each one of the DE proteins
  
  for(i in 1:length(DE_prot)){
    a = DE_prot[i]
    b = strsplit(strsplit(a,"\\[")[[1]][1], ";")
    
    # uniprot_ids_iter has the uniprot ID of the specific DE protein
    
    uniprot_ids_iter <-lapply(b, my_grep)[[1]]
    uniprot_ids <- c(uniprot_ids, uniprot_ids_iter)
    
    # MQ store their UniProt IDs differently
    c = str_extract(a, "SWISS-PROT:([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
    if (!is.na(c))
    {
      d = strsplit(c,":")[[1]][2]
      uniprot_ids <- c(uniprot_ids, d)
    }
    
  }
  uniprot_ids <- trimws(uniprot_ids)
  return(uniprot_ids)
}

# Run the enrichment analysis
run_enrichment_analysis <- function(UniProtList, myGOorganism, cond1, cond2)
{
  # Enrichment analysis utilizing gprofiler2 R package
  enrich <- gost(query= UniProtList, organism = myGOorganism, domain_scope = "annotated", significant = T, evcodes = TRUE, sources = c("GO", "KEGG", "REAC", "HPA", "WP"))
  enrich.matrix <- as.matrix(enrich$result[enrich$result$p_value < 0.05,c( "source", "term_name", "term_id", "p_value", "term_size", "query_size", "intersection_size", "intersection")])
  colnames(enrich.matrix) = c("Data Source", "Function", "Term ID", "p-Value", "Term size", "Query size", "Intersection Size", "Intersection")
  
  write.table(enrich.matrix, paste(outputFigsPrefix,"_enrichment_results_" , cond2, ".", cond1 , ".txt",sep=""), row.names=FALSE, sep = "\t", dec = ".", quote = F)
}


# ======= INITIALIZATION =======
unlabeled_peptide_regex<-"^$"
clearLabels()
clearMods()
paramssetfromGUI<-F
working_directory<-getwd()
limma_output<-"msdiffexp_out"
LabelFree<-F
AllowLabelRename<-F
AllowLS<-F
if (!DEBUG)
{
  source("MSdiffexp_definitions.R")
}

perform_analysis<-function(){
  levellog("",change=1)
  setwd(working_directory)
  rep_structure<-read.table(experimental_structure_file,col.names=c('raw_file','biorep','techrep','fraction'), sep="\t")
  rep_structure<-rep_structure[order(rep_structure[,2],rep_structure[,3],rep_structure[,4]),]
  LFQ_conds<-c()
  if(LabelFree)
  {
    #if labelfree load the lfq conditions structure
    LFQ_conds<-read.table(LFQ_conditions_file, col.names=c('raw_file', 'condition'), stringsAsFactors = F)
  }
  if (AllowLabelRename == T)
  {
    Rename_Array <<- read.table(Rename_Array_file, col.names=c('old_label', 'new_label'), stringsAsFactors = F)
  }
  if (AllowLS == T)
  {
    Ls_array <<- read.table(LS_Array_file, col.names=c('selected_raw_file', 'first_label', 'second_label'), stringsAsFactors = F)
  }
  if (RMisused == T)
  {
    RMrawfilesdata <<- read.table(RMrawfilesdata_file, col.names=c('id', 'name', 'brep', 'trep', 'frac', 'cond', 'used', 'selected'), stringsAsFactors = F)
  }
  if (RMisused == T)
  {
    RMtagsdata <<- read.table(RMtagsdata_file, col.names=c('id', 'name', 'brep', 'trep', 'frac', 'cond', 'used', 'selected'), stringsAsFactors = F)
  }
  #Because a condition can not be named "N" in ProteoSign, rename it to condN
  mi <- which(LFQ_conds$condition == "N")
  if(length(mi)>0)
  {
    levels(LFQ_conds$condition) <- c(levels(LFQ_conds$condition), "condN")
    LFQ_conds$condition[which(LFQ_conds$condition == "N")] <- "condN"
    LFQ_conds$condition <- factor(LFQ_conds$condition)
  }
  #take care of the same problem in conditions.labels as well
  #a condition can not be named "N" so just for this case rename it to condN
  for(i in 1:length(conditions.labels)){
    if(conditions.labels[i] == "N")
    {
      conditions.labels[i] <- "condN"
      conditions.labels <<- conditions.labels
    }
  }
  if(filterL_lbl == "N")
  {
    filterL_lbl <<- "condN"
  }
  #we will keep a copy of the original rep_structure to display in the graphs
  original_rep_structure <- rep_structure
  #we are not sure if the biorep and techrep numbers the user typed are sequential, the following code converts them to sequential numbers
  unique_reps <- unique(rep_structure$biorep)
  counter <- 1
  for(rep_i in unique_reps){
    mi <- which(rep_structure$biorep == rep_i)
    rep_structure$biorep[mi] <- counter
    counter <- counter + 1
    unique_techreps <- unique(rep_structure$techrep[mi])
    counter2 <- 1
    for(techrep_i in unique_techreps){
      mi2 <- which(rep_structure$biorep == counter - 1 & rep_structure$techrep == techrep_i)
      rep_structure$techrep[mi2] <- counter2
      counter2 <- counter2 + 1
    }
  }
  original_rep_structure$rep_desc<-paste(paste(paste('b',original_rep_structure$biorep,sep=''),'t',original_rep_structure$techrep,sep=''))
  if (!RMisused)
  {
    if(length(unique(rep_structure$biorep)) == 1){
      levellog("Error User: Cannot accept dataset with just one biological replicate. Aborting ...")
      return(F)
      # single_brep_file <<- T
      # nRequiredLeastBioreps <<- 1
      # } else {
      # single_brep_file <<- F
    }
  }
  if(length(unique(rep_structure$techrep)) > 1){
    if(length(unique(rep_structure$fraction)) > 1){
      # we have techreps and fractions
      rep_structure$rep_desc<-paste(paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep=''),'f',rep_structure$fraction,sep='')
      original_rep_structure$rep_desc<-paste(paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep=''),'f',rep_structure$fraction,sep='')
    }else{
      #we have bioreps and techreps
      rep_structure$rep_desc<-paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep='')
      original_rep_structure$rep_desc<-paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep='')
    }
  }else{
    if(length(unique(rep_structure$fraction)) > 1){
      # we have fractions but not techreps
      rep_structure$rep_desc<-paste(paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep=''),'f',rep_structure$fraction,sep='')
      original_rep_structure$rep_desc<-paste(paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep=''),'f',rep_structure$fraction,sep='')
    }else{
      # we just have bioreps
      rep_structure$rep_desc<-paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep='')
      original_rep_structure$rep_desc<-paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep='')
      # it should be like below, but for backward compatibility with other parts of the code, we keep the convention that in the rep. description, we will always have the terms 'b' (i.e. bio-rep) and 't', even if we don't have tech-reps ...
      # rep_structure$rep_desc<-paste('b',rep_structure$biorep,sep='')
    }
  }
  
  .GlobalEnv[["rep_structure"]]<-rep_structure
  .GlobalEnv[["LFQ_conds"]]<-LFQ_conds
  .GlobalEnv[["original_rep_structure"]]<-original_rep_structure
  .GlobalEnv[["n_bioreps"]]<-max(rep_structure$biorep)
  .GlobalEnv[["n_techreps"]]<-min(ddply(rep_structure[,c("biorep","techrep")],c("biorep"),function(x){return(max(x$techrep))})$V1)
  
  if(ProteinQuantitation){
    quantitated_items_lbl<<-"Protein"
  }else{
    quantitated_items_lbl<<-"Peptide"
  }
  if(file.exists(limma_output)){
    unlink(limma_output, recursive=T, force=T)
  }
  dir.create(limma_output)
  if(grepl("\"",readLines(evidence_fname, n=1))){
    levellog("Removing double quotes from input data file #1 ...")
    tmpdata<-gsub("\"", "", readLines(evidence_fname))
    evidence_fname_cleaned<-file(evidence_fname, open="w")
    writeLines(tmpdata, con=evidence_fname_cleaned)
    close(evidence_fname_cleaned)
  }
  levellog("Reading input data ...")
  if(PDdata){
    protein_groups<<-read.pgroups_v3(pgroups_fname,evidence_fname,time.point,keepEvidenceIDs=T)
  }else{
    if(grepl("\"",readLines(pgroups_fname, n=1))){
      levellog("Removing double quotes from input data file #2 ...")
      tmpdata<-gsub("\"", "", readLines(pgroups_fname))
      pgroups_fname_cleaned<-file(pgroups_fname, open="w")
      writeLines(tmpdata, con=pgroups_fname_cleaned)
      close(pgroups_fname_cleaned)
    }
    protein_groups<<-read.pgroups_v3(pgroups_fname,evidence_fname,time.point,keepEvidenceIDs=T)
  }
  #Restore the original rep descriptions to add to the graph
  if (!RMisused)
  {
    newcolumns <- names(protein_groups)
    oldcolumns = newcolumns
    for(my_column in newcolumns){
      for(my_repdesc in .GlobalEnv[["rep_structure"]]$rep_desc){
        if (grepl(my_repdesc, my_column)){
          temp_name <- .GlobalEnv[["original_rep_structure"]]$rep_desc[match(my_repdesc, .GlobalEnv[["rep_structure"]]$rep_desc)]
          newcolumns[match(my_column, newcolumns)] <- sub(my_repdesc, temp_name, my_column)
        }
      }
    }
    colnames(protein_groups) <- newcolumns
  }
  setwd(limma_output)
  temp_pgroups <- protein_groups
  write.table(temp_pgroups[, -which(names(temp_pgroups) %in% c("N.x","N.y"))],file=paste(outputFigsPrefix,"_proteinGroupsDF.txt",sep=""),row.names=F,sep="\t")
  setwd("..")
  if (!RMisused)
  {
    colnames(protein_groups) <- oldcolumns
  }
  #Create the expdesign table:
  expdesign<-c()
  #Rename conditions labels from condN back to N
  for(i in 1:length(conditions.labels)){
    if(conditions.labels[i] == "condN")
    {
      conditions.labels[i] <- "N"
      conditions.labels <<- conditions.labels
      colnames(protein_groups) <- sub("condN", "N", colnames(protein_groups))
    }
  }
  if(filterL_lbl == "condN")
  {
    filterL_lbl <<- "N"
  }
  #Rename condN back to N in lfq_conds and protein_groups
  mi <- which(LFQ_conds$condition == "condN")
  if(length(mi)>0)
  {
    levels(LFQ_conds$condition) <- c(levels(LFQ_conds$condition), "N")
    LFQ_conds$condition[which(LFQ_conds$condition == "condN")] <- "N"
    LFQ_conds$condition <- factor(LFQ_conds$condition)
    LFQ_conds <<- LFQ_conds
    colnames(protein_groups) <- sub("condN", "N", colnames(protein_groups))
  }
  for(cond_i in conditions.labels){
    expdesign<-rbind(expdesign,cbind(paste(sub("Intensity\\.","",sort(colnames(protein_groups)[grep(paste("Intensity.",cond_i,".b",sep=""),colnames(protein_groups))]))),cond_i))  
  }
  
  colnames(expdesign)<-c("Sample","Category")
  if(!RMisused){
    #the following lines also deal with restoring the original breps and treps numbers
    #temp vector has only the information of the replication (e.g. b1t1 or b1t1f1 if we have fractionation)
    temp_vector <- sub("(.*)\\.","", expdesign[,1])
    temp_vector <- original_rep_structure$rep_desc[match(temp_vector, sub("f.*", "", rep_structure$rep_desc))]
    #Make sure that expdesign (column Sample) contains data in te right format by merging expdesign and tmp_vector:
    tmp_counter <- 0
    for (expdesign_i in expdesign[,1]){
      expdesign[tmp_counter + 1,1] <- sub("(.*)\\..*",paste0("\\1.", temp_vector[tmp_counter + 1]), expdesign_i)
      tmp_counter <- tmp_counter + 1
    }
  }
  #Remove the fractionation information: (if any)
  expdesign[,1] <- sub("(.*\\..*)f.*", "\\1", expdesign[,1], perl = TRUE)
  write.table(expdesign,file="curr_exp_design.txt",row.names=F,quote=F,sep = "\t")
  exp_design_fname<<-"curr_exp_design.txt"
  
  levellog("Performing the analysis ...")
  
  norm.median.intensities <- do_limma_analysis(prepare_working_pgroups(protein_groups),time.point,exp_design_fname,exportFormat="pdf",outputFigsPrefix=outputFigsPrefix)
  
  levellog("Generating analysis plots ...")
  
  
  results<-do_results_plots(norm.median.intensities, time.point, exportFormat="pdf", outputFigsPrefix=outputFigsPrefix)
  
  # The results data frame contains all statistics of the homonymous data frame in the beggining of do results plots function (see the 
  # respective comment) with some columns added that show the intensity of each protein in each one of the replicates and conditions
  # the following table is an example of that:
  # 
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # |                                                                                                                                     | A        | Coef     | t        | P-value  | P-value adjusted L/H | F        | F P-value | H 1      | H 2      | H 3      | H 4      | H 5      | H 6      | L 1      | L 2      | L 3      | L 4      | L 5      | L 6      | N  | ID                                                                                                                                                                                                                     | avg log2 L/H | sd log2 L/H | N log2 L/H | avg log2 I L/H | log2 L/H 1 | log2 L/H 2 | log2 L/H 3 | log2 L/H 4 | log2 L/H 5 | log2 L/H 6 |
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # | H-INV:HIT000007749;Q58EX7   [Pleckstrin homology domain-containing family G mem ...]                                                | 15.44625 | -8.77978 | -7.90621 | 2.81E-05 | 0.005749399          | 62.50818 | 2.81E-05  | 19.17024 | 20.01426 | 20.32391 | NA       | NA       | NA       | 9.764948 | 9.764948 | 13.63917 | NA       | NA       | NA       | 6  | H-INV:HIT000007749;Q58EX7;REFSEQ:NP_001123203;REV__H-INV:HIT000000137;REV__O43304;REV__H-INV:HIT000029273;REV__Q8IWL3;REV__TREMBL:C9J7W5;REV__TREMBL:Q6ZNS1   [Pleckstrin homology domain-containing family G mem ...] | -8.77978     | 1.862788    | 3          | 15.44625       | -9.40529   | -10.2493   | -6.68474   | NA         | NA         | NA         |
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # | H-INV:HIT000038585;P22087;ENSEMBL:ENSP00000339522;ENSP00000397121;A6NHQ2   [34 kDa nucleolar scleroderma antigen;rRNA 2-O-meth ...] | 19.05326 | -8.65061 | -6.69742 | 7.70E-06 | 0.004729014          | 44.85541 | 7.70E-06  | 25.5831  | 24.76394 | 24.37479 | 24.48264 | 20.66157 | 20.40535 | 19.3296  | 22.0828  | 13.63917 | 13.78626 | 9.764948 | 9.764948 | 12 | H-INV:HIT000038585;P22087;ENSEMBL:ENSP00000339522;ENSP00000397121;A6NHQ2   [34 kDa nucleolar scleroderma antigen;rRNA 2-O-meth ...]                                                                                    | -8.65061     | 3.432693    | 6          | 19.05326       | -6.2535    | -2.68114   | -10.7356   | -10.6964   | -10.8966   | -10.6404   |
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # | H-INV:HIT000013642;Q96ME7;   [Zinc finger protein 512;cDNA FLJ52441, highly simi ...]                                               | 23.20329 | -7.21358 | -4.44701 | 0.001709 | 0.045619674          | 19.77591 | 0.001709  | 29.14918 | 29.3402  | NA       | 21.94086 | NA       | NA       | 19.2523  | 19.31348 | NA       | 20.22372 | NA       | NA       | 6  | H-INV:HIT000013642;Q96ME7;ENSEMBL:ENSP00000369040;ENSP00000407038;G3XAG1;B4DSM5;Q86XK6;TREMBL:B4E0X7;B4E0X7   [Zinc finger protein 512;cDNA FLJ52441, highly simi ...]                                                 | -7.21358     | 4.760497    | 3          | 23.20329       | -9.89689   | -10.0267   | NA         | -1.71715   | NA         | NA         |
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # | H-INV:HIT000029498;Q13547;   [Histone deacetylase 1;cDNA FLJ51764, highly simila ...]                                               | 21.75552 | -6.15447 | -3.12955 | 0.00981  | 0.115834791          | 9.794053 | 0.00981   | 24.38079 | 24.31524 | 25.41623 | 25.21878 | NA       | NA       | 24.3675  | 22.92021 | 13.63917 | 13.78626 | NA       | NA       | 8  | H-INV:HIT000029498;Q13547;F5GXM1;TREMBL:B4DSK9;B7Z3S4;ENSEMBL:ENSP00000362642;ENSP00000407859;Q5TEE2   [Histone deacetylase 1;cDNA FLJ51764, highly simila ...]                                                        | -6.15447     | 6.320277    | 4          | 21.75552       | -0.01329   | -1.39503   | -11.7771   | -11.4325   | NA         | NA         |
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # | SWISS-PROT:P55265-4;   [136 kDa double-stranded RNA-binding protein;Double ...]                                                     | 24.43804 | -5.56522 | -5.25566 | 0.000102 | 0.015618143          | 27.62197 | 0.000102  | 30.45977 | 30.4491  | 27.25717 | 25.40813 | 25.04    | 24.70973 | 23.28114 | 20.61311 | 22.99836 | 25.04826 | 19.29469 | 18.69703 | 12 | SWISS-PROT:P55265-4;E7ENU4;REFSEQ:NP_001102;P55265;H-INV:HIT000055014;ENSEMBL:ENSP00000405078;H0YCK3;ENSEMBL:ENSP00000401617   [136 kDa double-stranded RNA-binding protein;Double ...]                                | -5.56522     | 3.158404    | 6          | 24.43804       | -7.17863   | -9.83598   | -4.25881   | -0.35987   | -5.74532   | -6.0127    |
  # +-------------------------------------------------------------------------------------------------------------------------------------+----------+----------+----------+----------+----------------------+----------+-----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------+-------------+------------+----------------+------------+------------+------------+------------+------------+------------+
  # 
  # There are more columns in this data frame that may come in handy in case more functionality should be added to ProteoSign. N is the amount of replicates where the 
  # protein is detected, ID the proteins name as values in a data frame column (they contain the same information as the respective row names) avg log2 cond1/cond2
  # that shows the average log ratio of the two conditions across all replicates, sd log2 cond1/cond2 that show the standard deviation of the log ratios across the replicates,
  # N log2 cond1/cond2 shows how many log ratios we managed to compute, avg log2 I L/H
  
  
  if (GOenrichment)
  {
    # The last step is to perform the GO analysis using the R package gprofiler
    levellog("Perform enrichment analysis.")
    
    # First get the uniprot IDs of the Differentially expressed proteins for each combination of conditions
    ratio_combs<-combinations(nConditions,2,1:nConditions)
    for(j in 1:nrow(ratio_combs))
    {
      result<-tryCatch({
        # In this line conditions.labels[ratio_combs[i, 1]] and conditions.labels[ratio_combs[i, 2]] have the names of the conditions to compare
        uniprot_ids <- c()
        uniprot_ids <- get_uniprot_ids(results, conditions.labels[ratio_combs[j, 1]], conditions.labels[ratio_combs[j, 2]])
        if(length(uniprot_ids) == 0)
        {
          levellog(paste0("Warn User: GO enrichment analysis for conditions: ", conditions.labels[ratio_combs[j, 1]], " and ", conditions.labels[ratio_combs[j, 2]], " failed (", GOorganism, " was selected as target organism) because no UNIprot IDs were found for the differentially expressed proteins (if any)"))
          next;
        }
        # Since uniprot_ids contain the IDs of the DE expressed proteins lets run a GO analysis for them
        run_enrichment_analysis(uniprot_ids,GOorganism, conditions.labels[ratio_combs[j, 1]], conditions.labels[ratio_combs[j, 2]])
      }, error = function(err){
        levellog(paste0("Warn User: GO enrichment analysis for conditions: ", conditions.labels[ratio_combs[j, 1]], " and ", conditions.labels[ratio_combs[j, 2]], " failed (", GOorganism, " was selected as target organism)"))
      })
    }
  }
  
  # The enrichment analysis produces a table file like the following one:
  
  # +-------------+-------------------------------------------------------+------------+----------+-----------+------------+-------------------+----------------------------------------------------------------+
  # | Data Source | Function                                              | Term ID    | p-Value  | Term size | Query size | Intersection Size | Intersection                                                   |
  # +-------------+-------------------------------------------------------+------------+----------+-----------+------------+-------------------+----------------------------------------------------------------+
  # | CORUM       | Nop56p-associated pre-rRNA complex                    | CORUM:3055 | 2.40E-02 | 104       | 10         | 4                 | P22087,Q01081,P15880,P30050                                    |
  # +-------------+-------------------------------------------------------+------------+----------+-----------+------------+-------------------+----------------------------------------------------------------+
  # | GO:BP       | positive regulation of mRNA splicing, via spliceosome | GO:0048026 | 5.03E-06 | 26        | 14         | 4                 | Q13573,P62995,P38159,Q14011                                    |
  # +-------------+-------------------------------------------------------+------------+----------+-----------+------------+-------------------+----------------------------------------------------------------+
  # | GO:BP       | RNA processing                                        | GO:0006396 | 8.24E-06 | 968       | 14         | 9                 | P22087,Q01081,Q13573,P62995,P15880,P38159,Q52LJ0,P14866,Q14011 |
  # +-------------+-------------------------------------------------------+------------+----------+-----------+------------+-------------------+----------------------------------------------------------------+
  # | GO:BP       | positive regulation of mRNA processing                | GO:0050685 | 1.97E-05 | 36        | 14         | 4                 | Q13573,P62995,P38159,Q14011                                    |
  # +-------------+-------------------------------------------------------+------------+----------+-----------+------------+-------------------+----------------------------------------------------------------+
  
  # Which describes the Function detected per line, the term ID etc and all information relative to a GO enrichment analysis
  
  # Return to the main working directory and log the completion of the procedure
  setwd("..")
  levellog("",change=-1)
  levellog("Data analysis finished.")
  levellog("",change=-1)
  return(T)
}

#================ PRODUCTION ===============

if(DEBUG ){  
  if (DEBUG)
  {
    # the following line is for use of RStudio debugging
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    unlabeled_peptide_regex<-"^$"
    clearLabels()
    clearMods()
    paramssetfromGUI<-F
    working_directory<-getwd()
    limma_output<-"msdiffexp_out"
    LabelFree<-F
    AllowLabelRename<-F
    AllowLS<-F
    source("MSdiffexp_definitions.R")
    perform_analysis()
  }
  
}else{
  perform_analysis()
}

