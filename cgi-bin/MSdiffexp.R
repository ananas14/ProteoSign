#
#
# PROTEOSIGN - MSdiffexp.R
# Main Back-end R script for ProteoSign v2: a faster and evolved user-friendly online tool for statistical analyses of differential proteomics
#
# Reference:
# Theodorakis E., Antonakis A. N., Baltsavia I., Pavlopoulos G. A., Samiotaki M., Amoutzias G. D., Theodosiou T., Acuto O., Efstathiou G., Iliopoulos I.
# https://doi.org/10.1093/nar/gkab329
#
#

options(warn=1)
DEBUG <- FALSE

# DEPLOYMENT VERSION
#
# ======================================
# WARNING: Make sure to install the following packages as administrator/root (for them to be available to all users)
# ======================================

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
# if(!require("VennDiagram")){  install.packages("VennDiagram") }

cat("\nStarting ProteoSign...\n")
cat("=== Loading dependencies... ===\n\n")


library(limma, quietly = T, warn.conflicts = F, verbose = F)
library(statmod, quietly = T, warn.conflicts = F, verbose = F)
library(stringr, quietly = T, warn.conflicts = F, verbose = F)
library(reshape, quietly = T, warn.conflicts = F, verbose = F)
library(plyr, quietly = T, warn.conflicts = F, verbose = F)
library(ggplot2, quietly = T, warn.conflicts = F, verbose = F)
library(labeling, quietly = T, warn.conflicts = F, verbose = F)
library(gtools, quietly = T, warn.conflicts = F, verbose = F)
library(data.table, quietly = T, warn.conflicts = F, verbose = F)
library(outliers, quietly = T, warn.conflicts = F, verbose = F)
library(pryr, quietly = T, warn.conflicts = F, verbose = F)
library(gprofiler2, quietly = T, warn.conflicts = F, verbose = F)
library(VennDiagram, quietly = T, warn.conflicts = F, verbose = F)

cat("\n=== Dependencies successfully loaded... ===\n\n")

# levellog variables:
debuglog <- 10
loglvl <- 0
lastSysTime <- NA
firstSysTime <- NA

levellog <- function (msg, change=0, reset=F,after=F, supression=debuglog){
  
  #  Print msg with specific indentation (loglvl*3)
  #  - change [int]: change previous indentation level by <change> amount
  #  - reset [bool]: reset indentation level to 0 (no indentation)
  #  - after [bool]: change indentation after printing msg
  
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
  
  # This function computes the Pearsons's r for each pair of replicates for a specific pair of conditions
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

panel.lmline <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = reps.scatter.lmline.colour, ...){
  
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

pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE,lm=FALSE){
  
  # FROM: http://musicroamer.com/blog/2011/01/16/r-tips-and-tricks-modified-pairs-plot/
  
  # Called by do_results_plot (with smooth=TRUE,scale=TRUE,lm=TRUE), the x argument is all_ratios, a data frame containing u columns
  # (i.e. as many replicates we have) each for one replicate for a specific condition pair with the log ratios of the two conditions (e.g. L vs H) for each
  # protein
  
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

calcRowStats <- function(x,conds_cols_idxs,ratio_combs){
  
  # to be used with apply on a limma-results data.frame, called by do_result_plots
  # calculate mean, sd and N of ratios between available channels
  # conds_cols_idxs the column indexes of quant channels intensities
  # x a limma-results data.frame row (respective to one protein)
  
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

do_results_plots <- function(norm.median.intensities,exp_desc,exportFormat="pdf",outputFigsPrefix=""){
  
  # Produces Reproducibility, Volcano plot, MA plot and Scatterplot (matrix). Called by do_limma_analysis subroutine.
  
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
  
  FileHandle <- file(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""),"r")
  first_line <- readLines(FileHandle,n=1)
  hasLeadingTab <- (substr(first_line, 0, 1) == "\t")
  close(FileHandle)
  if (hasLeadingTab)
  {
    results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "", row.names = "X")
  } else {
    results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "")
  }
  
  
  
  
  if (nrow(ratio_combs) == 1)
  {
    # limma does not indicate the conditions compared in case only two conditions are compared in the column names so if ratio_combs = 1 add the conditions manually to the colnames of results
    colnames(results)[grep("p\\.value\\.adj",colnames(results), ignore.case = T)]<-paste("p.value.adj.",conditions[2],".",conditions[1],sep="")
  } else {
    # In any case make the specific column name case insensitive
    colnames(results) <- sub("p\\.value\\.adj(.*)$", "p.value.adj\\1", colnames(results), ignore.case = T)
  }
  
  
  tmp<-as.data.frame(t(norm.median.intensities))
  
  #rownames(tmp)<-colnames(norm.median.intensities)
  if (.GlobalEnv[["replicate_mismatch"]] == F) {
    # Here each column is renamed to "<condition> index" for example if the columns where Light.b1t1, Light.b1t2, Heavy.b1t1... they will be renamed to Light 1, Light 2, Heavy 1...
    nsamples<-length(colnames(tmp))/nConditions
    colnames(tmp)<-apply(data.frame(cbind(rep(conditions,each=nsamples),rep(1:nsamples))),1,function(x) paste(x['X1'],x['X2']))
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
  for(i in 1:nrow(ratio_combs))
  {
    col_desc_ <- paste("p.value.adj.",paste(conditions[ratio_combs[i,2]],".",conditions[ratio_combs[i,1]],sep=""),sep="")
    # Pick the right clumnin a case insensitive mnner to catch all MQ and PD versions
    col_picker <- grepl(col_desc_, x = colnames(results), ignore.case = T)
    na_indexes<-which(is.na(results[, col_picker]))
    if(length(na_indexes) > 0)
    {
      results[na_indexes,col_picker]<-1
      signTruth<-(signTruth | results[,col_picker]<pThreshold)
      results[na_indexes,col_picker]<-NA
    }else{
      signTruth<-(signTruth | results[,col_picker]<pThreshold)
    }
  }
  
  ndiffexp<-nrow(results[signTruth,])
  # conds_cols_idxs contain the indices on the columns for each condition. That is the first row will contain all results indices that contain quantification values deriving from this condition (but from different experiment replicates), the second will contain the same information for the second condition etc.
  conds_cols_idxs <- c()
  if (.GlobalEnv[["replicate_mismatch"]] == F) {
    # when there is not replicate mismatch creating conds_cols_idxs is easy:
    for (lbl_i in conditions) {
      conds_cols_idxs <- rbind(conds_cols_idxs, grep(paste("^", lbl_i, sep = ""), colnames(results)))
    }
  } else {
    # when there is replicate mismatch not all conditions match all replicates, that means that some replicate - condition combinations should have an NA value in conds_cols_idxs
    maxcols <- 0
    for (lbl_i in conditions) {
      if (length(grep(paste("^", lbl_i, sep = ""), colnames(results))) > maxcols) {
        maxcols = length(grep(paste("^", lbl_i, " ", sep = ""), colnames(results)))
      }
    }
    nsamples <- maxcols
    for (lbl_i in conditions) {
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
      levellog(paste0("Warn User: Condition ", conditions[i], " was found replicated just once, plots displaying comparisons between this condition and others might not be available!"))
    }
  }
  
  # The following line parses all rows of results one by one to calcrowstats that calculates all statistics for this row.
  # For example the first vector that will be sent to calcRowStats might be:
  
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
    ratio_i_str<-paste(conditions[ratio_combs[i,2]],".",conditions[ratio_combs[i,1]],sep="")
    colnames_d_<-rbind(colnames_d_, c(
      paste("log2.avg.",ratio_i_str,sep=""),
      paste("log2.sd.",ratio_i_str,sep=""),
      paste("log2.N.",ratio_i_str,sep=""),
      paste("log2.avg.I.",ratio_i_str,sep="")))
  }
  colnames_d_ <- as.vector(colnames_d_)
  
  
  # the following columns are named depending on the experimental structure
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions[ratio_combs[i,2]],".",conditions[ratio_combs[i,1]],sep="") # for example Heavy.Light
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
  
  save(results, pThreshold, quantitated_items_lbl, nConditions, calcRowStats, exp_desc, outputFigsPrefix, conditions, IsobaricLabel, PDdata, log.intensities, norm.intensities, fit2.coefficients, expdesign, file = "Plot_Generator.RData")
  
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
  for(i in 1:nrow(ratio_combs))
  {
    levellog(paste("Generating plots for combination #",i," ..."),change=1,after=T)
    ratio_i_str <- paste(conditions[ratio_combs[i, 2]], ".", conditions[ratio_combs[i, 1]], sep = "") # for example H.L
    
    # 1 - Volcano plot (-log10 P-value vs log ratio)
    result <- tryCatch({
      levellog("Making volcano plot ...")
      figsuffix<-paste("_",ratio_i_str,"-volcano","_",sep="")
      if(exportFormat == "pdf"){
        pdf(file=paste(outputFigsPrefix,figsuffix,exp_desc,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      
      ratio_i_<-paste("log2.",ratio_i_str,sep="")
      ratio_i_sd_col<-paste("log2.sd.",ratio_i_str,sep="")
      
      # Before starting the volcano plot generation we will compute all values needed for our plots. First goal is to find per replicate the (log2 intensity) +- (sd of log intensities) for the current condition pair
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
      
      ratio_i_p.value.adj<-paste("p.value.adj.",paste(conditions[ratio_combs[i,2]],".",conditions[ratio_combs[i,1]],sep=""),sep="")
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
        savePlot(filename=paste(outputFigsPrefix,figsuffix,exp_desc,".emf",sep=""),type="emf")
      }
      dev.off()
      png(paste("../",outputFigsPrefix,figsuffix,exp_desc,".png",sep=""), width = 1500, height = 1050)
      
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
        pdf(file=paste(outputFigsPrefix,figsuffix,exp_desc,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
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
      ggsave(paste(outputFigsPrefix,figsuffix,exp_desc,".png",sep=""), plot = p, device = "png", path = "..")
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,exp_desc,".emf",sep=""),type="emf")
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
        pdf(file=paste(outputFigsPrefix,figsuffix,exp_desc,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
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
      ggsave(paste(outputFigsPrefix,figsuffix,exp_desc,".png",sep=""), plot = p, device = "png", path = "..")
      
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,exp_desc,".emf",sep=""),type="emf")
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
        pdf(file=paste(outputFigsPrefix,figsuffix,exp_desc,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      pairs.panels(allratios,scale=T,lm=T)
      dev.off()
      png(paste("../",outputFigsPrefix,figsuffix,exp_desc,".png",sep=""), width = 1500, height = 1050)
      pairs.panels(allratios,scale=T,lm=T)
      if(exportFormat == "emf"){
        savePlot(filename=paste(outputFigsPrefix,figsuffix,exp_desc,".emf",sep=""),type="emf")
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
    ratio_i_str<-paste("(",conditions[ratio_combs[i,2]],") (",conditions[ratio_combs[i,1]],")",sep="")
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
  
  levellog(paste("do_results_plots: Quantified ",quant_species,": ",nrow(results)," (",exp_desc,")",sep=""))
  
  conditions <- gsub("\\.", " ", conditions)
  
  
  # Log amount of differentially expressed proteins
  
  for(i in 1:nrow(ratio_combs)){
    col_desc_<-paste("P-value adjusted ",paste(conditions[ratio_combs[i,2]],"/",conditions[ratio_combs[i,1]],sep=""),sep="")
    col_picker <- grepl(col_desc_, x = colnames(results), ignore.case = T)
    ndiffexp_tmp<-length(which(results[,col_picker]<pThreshold))
    levellog(paste("do_results_plots: Differentially expressed for ",conditions[ratio_combs[i,2]]," vs ",conditions[ratio_combs[i,1]]," : ",ndiffexp_tmp,sep=""))
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
  write.table(diffexp[signTruth,],dec=dec,file=paste(outputFigsPrefix,"_diffexp_",exp_desc,".txt",sep=""),sep="\t",row.names=F,quote=F)
  
  colnames(diffexp) <- oldcolumns
  
  # Merge diffexp and results to create a large data frame containing in fact all the information created from ProteoSign up to now,
  # i.e. all information stored in results and the ratio counts per protein
  diffexp<-merge(diffexp,results[,-grep("^(avg log2|P-value adjusted)",colnames(results))],by.x=c(quantitated_items_lbl),by.y=c("ID"),all.x=T)
  
  # Save this information to a file
  write.table(diffexp,dec=dec,file=paste(outputFigsPrefix,"_results_",exp_desc,".txt",sep=""),sep="\t",row.names=F,quote=F)
  
  levellog("",change=-1)
  return(results)
}

do_limma_analysis <- function(working_pgroups,exp_desc,exp_design_fname,exportFormat="pdf",outputFigsPrefix=""){
  
  # Performs the differential expression analysis through limma, after quantile normalization.
  
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
  
  setwd(limma_out_dir)
  
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
    pdf(file=paste(outputFigsPrefix,"_limma-graphs_",exp_desc,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
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
    contrasts<-c(contrasts,paste(conditions[ratio_combs[i,2]],"-",conditions[ratio_combs[i,1]],sep=""))
  }
  contrasts <- makeContrasts(contrasts=contrasts, levels=design)
  write.table(contrasts,file=paste(outputFigsPrefix,"_limma-contrasts-matrix_",quantitated_items_lbl,"Groups.txt",sep=""),sep="\t",row.names = T, col.names=NA)  
  # Apply contrast matrix and do empirical bayes analysis to get p-values etc.
  
  levellog("Performing hypothesis testing ...")
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  
  # Plot a Histogram of co-efficients (log2 ratio)
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions[ratio_combs[i,2]],"/",conditions[ratio_combs[i,1]],sep="")
    hist(fit2$coefficients[,i],main=paste("Log2 Fold Change ",ratio_i_str,sep=""), xlab="Log2 Fold Change", breaks=50 )
  }   
  fit2.coefficients <<- fit2$coefficients
  if(exportFormat == "emf"){
    savePlot(filename=paste(outputFigsPrefix,"_limma-graphs_",exp_desc,"_hist.emf",sep=""),type="emf")
  }
  
  dev.off()
  
  
  #save the graphs to pngs
  png(paste("../",outputFigsPrefix,"_limma-graphs_",exp_desc,"_intensities-before-normalization.png",sep=""), width = 1500, height = 1050)
  boxplot(log.intensities)
  title(main="Intensities Before Normalisation")
  
  dev.off()
  
  png(paste("../",outputFigsPrefix,"_limma-graphs_",exp_desc,"_intensities-after-normalization.png",sep=""), width = 1500, height = 1050)
  boxplot(norm.intensities)
  title(main="Intensities After Normalisation")
  
  dev.off()
  
  tmp_conditions <- conditions
  if (IsobaricLabel)
  {
    if(!PDdata){
      conditions <- sub("Reporter\\.intensity\\.", "Reporter.", conditions)
    }else{
      conditions <- sub('X([[:digit:]])', "\\1", conditions)
    }
  }
  
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions[ratio_combs[i,2]],"/",conditions[ratio_combs[i,1]],sep="")
    ratio_i_str_with_dash <-paste(conditions[ratio_combs[i,2]],"-",conditions[ratio_combs[i,1]],sep="")
    png(paste("../",outputFigsPrefix,"_limma-graphs_",exp_desc,"_log2-fold-change-histogram", ratio_i_str_with_dash ,".png",sep=""), width = 1500, height = 1050)
    hist(fit2$coefficients[,i],main=paste("Log2 Fold Change ",ratio_i_str,sep=""), xlab="Log2 Fold Change", breaks=50 )
    dev.off()
  } 
  
  if (IsobaricLabel)
  {
    conditions <- tmp_conditions
  }
  
  # Output analysis details to file
  # adjust="BH" means adjust the calculated p-values for multiple testing using
  # the Benjamini Hochberg method (FDR)
  levellog("Saving analysis results to file ...")
  write.fit(fit2, file=paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""), adjust="BH")
  
  # Notice that write.fit output will give us a table with average log2 intensity (A), fold changes for each comparison between conditions, p-values and adjusted p-values, F and F p values for each protein (the protein names are not displayed in this file)
  # Warning! latest releases of limma have changed the way the columns are displayed! this might not work in latest limma version! 3.18.13 suits best for PS
  if (packageVersion("limma") > "3.30.0" & nConditions == 2)
  {
    
    # Opening the results file might not be trivial, see the comment in do results plot for more infomation
    # and why it will be opened in such a complicated manner:
    
    FileHandle <- file(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""),"r")
    first_line <- readLines(FileHandle,n=1)
    hasLeadingTab <- (substr(first_line, 0, 1) == "\t")
    close(FileHandle)
    if (hasLeadingTab)
    {
      results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "", row.names = "X")
    } else {
      results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "")
    }
    
    colnames(results)[1] <- "A"
    colnames(results)[2] <- "Coef"
    colnames(results)[3] <- "t"
    colnames(results)[4] <- "p.value"
    colnames(results)[5] <- "p.value.adj"
    write.table(results, file=paste(outputFigsPrefix,"_condition-i_vs_condition-j_",exp_desc,".txt",sep=""), sep = "\t", quote = FALSE)
  }
  return(norm.median.intensities)
}

read.pgroups <- function(fname,evidence_fname,exp_desc){
  
  # read pgroups will read the main quantitation file of the experiment that is the evidence file for MQ or the MultiConsensus file for PD
  levellog("",change=1)
  levellog("Reading data file ...")
  
  # The following line reads the file and stores it to evidence - this is one of the most time consuming lines since these files may be quite big
  # and copying all of them to the memory takes a while
  
  
  evidence<-read.table(evidence_fname, header = T, sep = "\t",quote="",stringsAsFactors=F,comment.char = "")
  
  # Keeping all this information to memory is a waste of memory space so we will drop all unnecessary columns unless we are in DEBUG mode
  # Note that if we decide to use other columns as well in the future we should update nec_evidence_columns to keep them in the process
  
  if (!DEBUG)
  {
    # TODO: Perform drop of unnecessary columns
    levellog("Dropping unnecessary columns ...")
    # nec_evidence_columns stores regular expressions of headers of columns to keep from evidence (case insensitive) - it should contain for both MQ and PD
    nec_evidence_columns <- c('Protein.Group.Accessions', '^Protein.Accessions', '^Proteins$')
  }
  
  # We will start searching for the necessary columns and store their headers to variables. E.g. pgroups_colname will store the regexes that match all headers of columns that have to do
  # with Proteins detected. The information that is retrieved by different preprocessing programs may differ significantly. Lets take an example working on an evidence file of MQ:
  #
  #
  # All lines correspond to a peptide sequence (with or without modifications and in a certain charge state) detected in a single MS run.
  # Each peptide sequence is matched against proteins that it is possible to be derived from. These proteins are stored in our example in the Proteins column and are named
  # after IDs of well established databases. An example is H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446 that shows us 3 ids one from H-invitational database, SwissProt and Ensemble database.
  # This is not the best example since these accessions lead to the same protein. If there will be many accessions pointing to the same protein or not depends on the parameters that wre used in MQ
  # and specifically on the FASTA file of valid proteins that was selected. Most of the times the acessions we will get will be for a single database e.g. UniProt
  #
  # So all these ids refer to the same protein but this is not always the case. For example a line below in our example has this record: H-INV:HIT000091288;P30153;TREMBL:B4DQY1;TREMBL:B4E1Q0;ENSEMBL:ENSP00000391905;C9J9C1
  # Taking the two swissprot accessions (i.e. P30153 and C9J9C1) we can see that they correspond to:
  #
  # Serine/threonine-protein phosphatase 2A 65 kDa regulatory subunit A alpha isoform and
  # cDNA FLJ59603, highly similar to Serine/threonine-protein phosphatase 2A 65 kDa regulatory subunit A alpha isoform respectively that are two very similar protein subunits.
  #
  # Indeed this can be seen in their uniprot database entry. Taking a look at the sequence that was detected in this evidence line we can see that it is
  # AAADGDDSI YPIAVIIDEI R that matches the sequence AAADGDDSL YPIAVLIDEL R in both proteins. It is notable that some leukines are detected as isoleukines but this is a common issue
  # in proteomics so most of the times the preprocessing programs are set so that isoleukines and leukines are treated as equivalent. This is also an example where two proteins can be matched to
  # an evidence line but even in this case these are very similar and it may not be a problem to assign e.g. P30153 and not C9J9C1 to the evidence line
  # Sometimes though 2 quite different proteins might both be matched against an evidence line since a detected sequence can be part of two different proteins in the target organism's
  # proteome. For example the sequence AAAASAAEAGIATTGTEDSDDAIIK that is detected in our example is matched against
  #
  # P55036 - Putative PIP5K1A and PSMD4-like protein and
  # A2A3N6 - 26S proteasome non-ATPase regulatory subunit 4
  #
  # searching the respective uniprot records we see that P55036 is a putative protein - possible product of a gene PIPSL which "appeared in hominoids by L1-mediated retrotransposition in a hominoid ancestor of a readthrough, intergenically spliced transcript between the PIP5K1A and PSMD4 genes."
  # A2A3N6 on the other hand derives from the PSMD4 gene.
  #
  # These kinds of similar proteins are commonly found in the same evidence line and it is the preprocessing program's job to define
  # which one is the most probable or biologically useful Protein to match the detected peptide. This is commonly found in another column that is called "Leading Proteins" "Master Protein" or something similar
  # PD usually gives us one single protein to overcome this issue but MQ shows us the whole protein group. One way to solve this is to always take the first protein that is shown
  # in the Proteins column - this is the one that is detected by the largest number of peptides and usually is the one of interest. This creates a problem if many databases are involved
  # as in the case above since the accession of the first database and only will be fetched or when an unreviewed protein is detected to have more peptide hits than a reviewed one
  # but all these are well known issues that are fixed beforehand using the correct FASTA file when loading MQ.
  #
  # We will take the proteins column or MQ that stores all proteins of the protein group
  
  if(PDdata)
  {
    # The information for proteins matched against evidence lines in PD is stored in the column Protein Group Accessions or Protein Accessions
    if ('Protein.Group.Accessions' %in% colnames(evidence))
    {
      pgroups_colname<-'Protein.Group.Accessions'
    } else if ('Protein.Accessions' %in% colnames(evidence)) {
      # Pattern must start with Protein.Accessions because id PD 2.4 there is also another column as Master.Protein.Accessions
      pgroups_colname<-'^Protein.Accessions'
    } else {
      levellog("Error User: The dataset does not contain the columns 'Protein Group Accessions' or 'Protein Accessions'")
    }
  } else {
    # As for MQ the information for proteins matched against evidence lines is stored in the column "Proteins"
    pgroups_colname<-'^Proteins$'
  }
  
  # Let's rename the aforementioned column to Protein.IDs to make it independable of the preprocessing program
  
  colnames(evidence)[grepl(pgroups_colname,colnames(evidence))]<-'Protein.IDs'
  
  # The next step is to match each protein group to a descriptive name. Sometimes the evidence file already has a column called protein.names that has these descriptions
  # If this is not the case we should create the column using the second file - the protein groups file that has all the info for the protein groups
  
  # Check if protein names exists in evidence and if so rename it to Protein.Names to face case sensitivity ambiguity
  
  col_Protein.names <- length(grep('Protein.Names',colnames(evidence))) > 0
  col_Protein.namesLC <- length(grep('Protein.names',colnames(evidence))) > 0
  
  if (col_Protein.namesLC) colnames(evidence)[grepl('Protein.names',colnames(evidence))] <- 'Protein.Names'
  
  if(!PDdata & !col_Protein.names & !col_Protein.namesLC){
    
    # In MQ, if there are no descriptive names in the evidence file, the protein_groups file is needed to add descriptive names to each protein. These names may be already contained in
    # the Protein Names column of the file but if this is not the case we should create a descriptive column from the protein groups file
    
    # Read the file:
    pgroups<-read.table(fname, header = T, sep = "\t",quote="",stringsAsFactors=F,comment.char = "")
    
    # If both Protein.Names and Protein.names column do not exist on the pgroups file (depends on MQ version), create one from the Fasta Headers column
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
    
    # At this point we have a column called Protein.names that contains the descriptions of all proteins in the pgroups data table
    
    # Construct a table, mapping the correct protein groups IDs (and the corresponding proteins names) to the evidence IDs
    # First check if there are any blank lines and remove them:
    
    mi<-which(pgroups$Evidence.IDs == "")
    if (length(mi)>0)
    {
      pgroups <- pgroups[-mi,]
    }
    
    
    if (!'Protein.IDs' %in% colnames(pgroups) & 'Peptide.IDs' %in% colnames(pgroups))
    {
      colnames(pgroups)[colnames(pgroups) == 'Peptide.IDs'] <- 'Protein.IDs'
    }
    
    # the following lines unfold pgroups evidence column by creating a long table (tmp.table.1) that replicates protein id and protein name for each evidence id
    # matched to a protein. For example for the record in pgroups that has the following protein id - protein name and evidence ids triad:
    # 
    # Protein.IDs Evidence.IDs
    #    7960;9207  77961;90856
    # 
    # Protein.Names
    # GC-box-binding zinc finger protein 1;Transcription factor ZBP-99;Zinc finger DNA-binding prot [...]
    #
    # will be changed to:
    #
    # Protein.IDs Protein.names                                                                                       id
    # 7960;9207	  GC-box-binding zinc finger protein 1;Transcription factor ZBP-99;Zinc finger DNA-binding pro [...]	77961
    # 7960;9207	  GC-box-binding zinc finger protein 1;Transcription factor ZBP-99;Zinc finger DNA-binding pro [...]	90856
    #
    # This way each line corresponds to a single evidence id
    
    tmp.table.1<-data.table(do.call(rbind, apply(pgroups[,c('Protein.IDs','Protein.Names','Evidence.IDs')], 1, function(x){return(cbind(x['Protein.IDs'], x['Protein.Names'], unlist(strsplit(x['Evidence.IDs'], ';'))))})))
    setnames(tmp.table.1, colnames(tmp.table.1), c('Protein.IDs', 'Protein.Names', 'id'))
    
    # Now that we have unfolded pgroups we can inner join it with evidence with common key the evidence ids to add the descriptions to the evidence data table
    class(tmp.table.1$id)<-'integer'
    setkey(tmp.table.1, id)
    # Get the evidence records
    tmp.table.2<-data.table(evidence)
    setkey(tmp.table.2, id)
    
    # Remove the incorrect Protein.IDs column
    tmp.table.2[, c('Protein.IDs', 'Protein.Names') := NULL]
    # Inner join the mapping table with the evidence table and return the data frame that we ought to have in the first place
    evidence<-data.frame(tmp.table.1[tmp.table.2])
    
    

  }
  
  # Get rid of all contaminant and reverse peptides if any (valid only for MQ)
  if (!PDdata)
  {
    if(length(grep("Contaminant", colnames(evidence))) > 0){
      evidence<-evidence[evidence$Contaminant == '', ]
    }
    if(length(grep("Reverse", colnames(evidence))) > 0){
      evidence<-evidence[evidence$Reverse == '', ]
    }
  }

  
  # TODO: test on Iso Labeling:
  
  # In the case of Isobaric labeling we should reformat the table before proceeding, afterwards we will treat the data as
  # if they were label-free data
  if(IsobaricLabel)
  {
    if(!PDdata)
    {
      evidence$Intensity <- NULL
      varcolnames <- grep("^Reporter.intensity.[[:digit:]]", colnames(evidence), value = TRUE)
      evidence <- reshape(evidence, varying = varcolnames, v.names = "Intensity", timevar = "Labeling.State", times = varcolnames, direction = "long", new.row.names=sequence(prod(length(varcolnames), nrow(evidence))))
      conditions<<-sub("^X", "Reporter.intensity.", conditions)
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
    } else {
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
  
  # log the amount of identified proteins
  
  levellog(paste0("read.pgroups: Identified proteins: ",length(unique(evidence$Protein.IDs))," (",exp_desc,")"))
  
  # n1 will hold the number of evidence we have (for logging purposes only):
  n1<-nrow(evidence)
  
  # Drop all evidence that have blank ProteinID
  evidence<-evidence[nchar(evidence$Protein.IDs) > 0,]
  
  levellog(paste0("read.pgroups: Discarded PSM records due to unassigned protein group: ",(n1-nrow(evidence))))
  
  # TODO test for PD data
  ## Make Protein.IDs human-readable
  if(PDdata){
    pgroups_colname<-'Protein.Descriptions'
    if (!'Protein.Descriptions' %in% colnames(evidence)) {
      evidence[, c('Protein.Descriptions')] <- ""
    }
  }else{
    pgroups_colname<-'Protein.Names'
  }
  
  # Let's get only protein information from our evidence table to a data table to make things faster:
  # Setting the i as a primary key will afterwards help us return the edited data back to evidence our main data frame
  tmp.table<-data.table(cbind(evidence[, c('Protein.IDs', pgroups_colname)], i=1:nrow(evidence)))
  setkey(tmp.table, Protein.IDs)
  
  # tmp.table:
  # +-----------------------------------------+----------------------------------+-------+
  # | Protein.IDs                             | Protein.Names                    | i     |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59601 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59602 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59603 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59604 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59605 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59606 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59607 |
  # +-----------------------------------------+----------------------------------+-------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59608 |
  # +-----------------------------------------+----------------------------------+-------+
  # 
  # Notice that setting key to Protein.IDs sorts the respective column
  
  # Make a new protein description column in the tmp.table data table containing both the protein accession IDs and their descriptions. This are separated by
  # a "|:" symbol
  tmp.table[, pdesc := paste0(Protein.IDs, "|: ", get(pgroups_colname))]
  
  # tmp.table: 
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | Protein.IDs                             | Protein.Names                    | i     | pdesc                                                                     |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1                   |                                  | 86714 | AAR31361.1;ABV71636.1:                                                    |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1                   |                                  | 86715 | AAR31361.1;ABV71636.1:                                                    |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59601 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1: 55 kDa immediate-early protein 1 |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59602 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1: 55 kDa immediate-early protein 1 |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59603 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1: 55 kDa immediate-early protein 1 |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59604 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1: 55 kDa immediate-early protein 1 |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59605 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1: 55 kDa immediate-early protein 1 |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1 | 55 kDa immediate-early protein 1 | 59606 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1: 55 kDa immediate-early protein 1 |
  # +-----------------------------------------+----------------------------------+-------+---------------------------------------------------------------------------+
  # 
  # Remember that tmp.table has as many lines as evidence
  
  # The following line orders all data using the i index and returns the detailed descriptions safely to the Protein.ID column of evidence data frame
  
  evidence$Protein.IDs<-tmp.table[order(i),pdesc]
  
  # Save the name of the rawfile column to a variable to make it independent to the preprocessing program
  
  if(PDdata){
    rawfile_col<-'Spectrum.File' 
  } else {
    if(length(grep("Raw.File", colnames(evidence))) > 0)
    {
      rawfile_col<-'Raw.File'
    } else {
      rawfile_col<-'Raw.file'
    }
  }
  
  # There is always a column that defines the condition of the specific evidence line. Save that to cond_spec_col to make PS independent of the preprocessing program and the
  # experiment type
  
  # TODO test for labelfree and iso labeling
  if(LabelFree){
    cond_spec_col<-rawfile_col
    if(IsobaricLabel){
      if(PDdata){
        cond_spec_col<-'Modifications'
      } else {
        cond_spec_col<-'Labeling.State'
      }
    }
  }else{
    if(PDdata){
      cond_spec_col<-'Modifications'
    } else {
      cond_spec_col<-'Labeling.State'
    }
  }
  
  # TODO test RM
  if(RMisused){
    levellog("read.pgroups: Transforming data for Replication Multiplexing ...")
    
    # when RM is chosen, in this line evidence has two main columns, the column that matches an evidence line to a condition
    # (cond_spec_col) and one that matches the line to  raw file (rawfile_col)
    #
    # In case of RM breps treps and conditions may come from either raw files or tags but in this data format creating
    # two new columns describing the structure correctly is not difficult.
    #
    # First create the column for conditions
    
    RMrawfilesdata <- RMrawfilesdata[!RMrawfilesdata$used == 'false',]
    RMtagsdata <- RMtagsdata[!RMtagsdata$used == 'false',]
    if(RMconditionsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'cond')], by.x = rawfile_col, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'cond')], by.x = cond_spec_col, by.y = 'name')
    }
    
    # The conditions array is already set from the front end
    # Now we will initialize the new_raw_file column that will contain pseudo-raw files describing our bioreps, techreps and fracs
    
    colnames(RMrawfilesdata)[3:5] <- c('new_brep', 'new_trep', 'new_frac')
    colnames(RMtagsdata)[3:5] <- c('new_brep', 'new_trep', 'new_frac')
    if(RMbrepsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_brep')], by.x = rawfile_col, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'new_brep')], by.x = cond_spec_col, by.y = 'name')
    }
    
    # Do the same thing for treps
    
    if(RMtrepsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_trep')], by.x = rawfile_col, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'new_trep')], by.x = cond_spec_col, by.y = 'name')
    }
    if(RMbrepsinrawfiles | RMtrepsinrawfiles) {
      
      # Do the same thing for fracs
      
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_frac')], by.x = rawfile_col, by.y = 'name')
    }
    if(!RMbrepsinrawfiles & !RMtrepsinrawfiles) {
      evidence$new_raw_file <- paste0('b', evidence$new_brep, 't', evidence$new_trep)
    } else {
      evidence$new_raw_file <- paste0('b', evidence$new_brep, 't', evidence$new_trep, 'f', evidence$new_frac)
    }
    
    # Now lets refresh the rep_structure array the pseudo raw files we created are descriptive and contain the breps treps and fracs
    
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
    
    # Now erase the old columns for raw files and conditions and replace with the new ones
    
    evidence[, c(rawfile_col, cond_spec_col)] <- list(NULL)
    colnames(evidence)[colnames(evidence) == "new_raw_file"] <- rawfile_col
    colnames(evidence)[colnames(evidence) == "cond"] <- cond_spec_col
    if(length(unique(rep_structure$biorep)) == 1){
      levellog("Error User: Cannot accept dataset with just one biological replicate. Aborting ...")
      return(F)
    }
  }
  
  levellog("read.pgroups: Assigning labels ...")
  levellog("",change=1)
  
  # First define some variables:
  evidence$label_<-NA
  background_species_lbl<-NA
  
  
  if (!PDdata & !LabelFree)
  {
    # In case of Labelled maxquant data each evidence line has some columns that provide us with quantitation values. The interpretation of these quantitation values
    # depend on the way that the evidence line (otherwise called "hit") was detected. This is defined by the "Type" column that can have some specific values:
    #
    # MULTI-MSMS is the easiest to understand meaning that a peptide was detected in MS and quantified in all the labels of the experiment,
    #   this will bring us for example two valid quantitation numbers for a SILAC 2plex
    # MSMS means that no MS1 label state could be detected so the peptide was hit and identified but not quantified so this should be ommited from the analysis
    # MULTI-SECPEP means that the hit is of a peptide that was detected using "second peptide algorithm". This is an algorithm that is used to quantify peptides that are
    #   co-elluted and sent to MSMS together with others. The algorithm is well established and also produces valid quantification values
    # MULTI-MATCH means that the hit was detected using match between runs. In match between runs MQ guesses that a peptide should be present in a elution time point
    #   in an MS run by matching its elution profile with the profiles of other runs. It also provides valid quantitation values
    # ISO-MSMS means that an isotope cluster was detected but with no pair. In this case the labeling state column states the label the hit was assigned
    #   normally all other label intensity columns should be set to 0
    #
    # A safe way to get all necessary quantitation values is to ommit all lines that have MSMS as type and to force to 0 all intensity values in ISO-MSMS that match
    # to labels different from the one stated under Labeling state
    #
    
    evidence <- evidence[!(evidence$Type == "MSMS"),]
    
    # To force 0 intensity to all other labels lets see this example:
    
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | id  | Sequence                        | Type     | Labeling.State | Intensity.L | Intensity.H |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | 62  | AAAAAAAIQAK                     | ISO-MSMS | 1              | 33117000    | 5176300     |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | 71  | AAAAAAAIQAK                     | ISO-MSMS | 0              | 873790      | 17029000    |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | 90  | AAAAGGGGPGTAVGATGSGIAAAAAGIAVYR | ISO-MSMS | 1              | 298290      | 96692       |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | 102 | AAAAVVEFQR                      | ISO-MSMS | 0              | 672120      | 9695200     |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | 104 | AAAAVVEFQR                      | ISO-MSMS | 0              | 4325700     | 105140000   |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    # | 116 | AAAAVVEFQR                      | ISO-MSMS | 0              | 509660      | 8551300     |
    # +-----+---------------------------------+----------+----------------+-------------+-------------+
    
    # Here we have isolated some records in evidence that were identified with ISO-MSMS type and have quantitation values for both H and L
    # The labeling state defines which one is the valid intensity value so the other one should be set to 0
    # The labeling state is encoded so that 0 stands for Light 1 for Medium and 2 for Heavy (or 0 Light 1 Heavy in 2plex). When the isotope can not be assigned safely to a
    # label then labeling state is set to -1
    # First erase all records with labeling state -1 and type ISO-MSMS
    
    evidence <- evidence[!(evidence$Type == "ISO-MSMS" & evidence$Labeling.State == "-1"),]
    
    # Now check if there is an Intensity.M label to decide if we have a 2plex or 3plex experiment
    
    if (any(grepl("Intensity.M", x = colnames(evidence))))
    {
      
      # 3plex example:
      # Set to 0 all ISO-MSMS Intensities with labeling state = 0 except the light one
      evidence[evidence$Type == "ISO-MSMS" & evidence$Labeling.State == "0", grep("Intensity\\.[^L]", colnames(evidence))] <- 0
      
      # Do the same for M and H
      evidence[evidence$Type == "ISO-MSMS" & evidence$Labeling.State == "1", grep("Intensity\\.[^M]", colnames(evidence))] <- 0
      evidence[evidence$Type == "ISO-MSMS" & evidence$Labeling.State == "2", grep("Intensity\\.[^H]", colnames(evidence))] <- 0
    } else {
      
      # 2plex example: (as with the 3plex example above)
      evidence[evidence$Type == "ISO-MSMS" & evidence$Labeling.State == "0", grep("Intensity\\.[^L]", colnames(evidence))] <- 0
      evidence[evidence$Type == "ISO-MSMS" & evidence$Labeling.State == "1", grep("Intensity\\.[^H]", colnames(evidence))] <- 0
    }
    
    # Now each hit has one valid quantitation value per condition as seen in these selected columns of evidence:
    #
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | id | Sequence    | Type         | Labeling.State | Intensity.L | Intensity.H |
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | 58 | AAAAAAAIQAK | MULTI-SECPEP | NA             | 28764       | 29688       |
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | 59 | AAAAAAAIQAK | MULTI-MSMS   | NA             | 345050      | 369960      |
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | 60 | AAAAAAAIQAK | MULTI-MSMS   | NA             | 58373000    | 66849000    |
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | 61 | AAAAAAAIQAK | MULTI-MSMS   | NA             | 7549500     | 9036800     |
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | 62 | AAAAAAAIQAK | ISO-MSMS     | 1              | 0           | 5176300     |
    # +----+-------------+--------------+----------------+-------------+-------------+
    # | 63 | AAAAAAAIQAK | MULTI-MSMS   | NA             | 131360      | 209020      |
    # +----+-------------+--------------+----------------+-------------+-------------+
    #

    # Melt the data frame to create an Intensity column with the intensity and a label_ with the respective label in the long format
    
    evidence <- reshape::melt(data = evidence, id.vars = colnames(evidence)[!grepl("Intensity\\..+", colnames(evidence))], variable_name = "label")
    evidence$label_ <- evidence$label
    evidence$Intensity <- evidence$value
    
    # Drop the unneccessary columns
    evidence <- evidence[,!(colnames(evidence) %in% c("label", "value"))]
    
    # Remove the unneccessary Intensity prefix from the label_ column
    evidence$label_ <- str_replace(evidence$label_, "^Intensity\\.", "")
    
    # Discard all evidence lines with 0 intensity
    evidence <- evidence[evidence$Intensity != 0,]
    
    # Now evidence looks like this (chosen columns and rows displayed):
    # 
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | id    | Protein.Group.IDs | Sequence         | Type       | Labeling.State | Intensity | label_ |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 0     | 419               | AAAAAAAAAAGAAGGR | MULTI-MSMS | NA             | 373520    | L      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 1     | 419               | AAAAAAAAAAGAAGGR | MULTI-MSMS | NA             | 10432000  | L      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 2     | 419               | AAAAAAAAAAGAAGGR | MULTI-MSMS | NA             | 4183800   | L      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 99432 | 7005              | ANIPQSFQVDTSK    | MULTI-MSMS | NA             | 4303800   | H      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 99433 | 7006              | ANIPQSFQVDTSK    | MULTI-MSMS | NA             | 3313400   | H      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 99434 | 7007              | ANIPQSFQVDTSK    | MULTI-MSMS | NA             | 1608800   | H      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 99435 | 7008              | ANIPQSFQVDTSK    | ISO-MSMS   | 1              | 6892100   | H      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    # | 99436 | 7009              | ANIPQSFQVDTSK    | MULTI-MSMS | NA             | 3708700   | H      |
    # +-------+-------------------+------------------+------------+----------------+-----------+--------+
    
  }
  
  
  
  # TODO test PD, Label free and iso labeling
  for(i in 1:length(conditions)){
    if(PDdata){
      if(LabelFree){
        if(!IsobaricLabel)
        {
          mi<-which(grepl(conditions[i], LFQ_conds[, "condition"]))
          mi2<-which(grepl(paste(LFQ_conds[mi,]$raw_file, collapse="|"), evidence[, cond_spec_col]))
          evidence[mi2,]$label_<-conditions[i]
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
          mi<-which(grepl(conditions[i], LFQ_conds[, "condition"]))
          mi2<-which(grepl(paste(LFQ_conds[mi,]$raw_file, collapse="|"), evidence[, cond_spec_col]))
          evidence[mi2,]$label_<-conditions[i]
        } else {
          evidence$label_<-evidence$Labeling.State
        }
        
      } else {
        
       # The case of MQ labelled data is already processed above
        
      }
    }
    levellog(paste0("read.pgroups: Assigned label '", conditions[i],"'."))
  }
  
  # Rename any labels if necessary
  # TODO test label renaming
  if (AllowLabelRename == T)
  {
    if (length(Rename_Array$old_label) != 0)
    {
      for(i in 1:length(Rename_Array$old_label))
      {
        if(Rename_Array$old_label[i] != Rename_Array$new_label[i])
        {
          # The case where old label = new label is common since if the user did not ask for a label rename (merge) for a label
          # this label is sent to be renamed in R to itself
          # in any other case rename the labels that are merged to the same label so that they become indistinguishable
          # and refresh the conditions labels by erasing the old label and adding the new if necessary
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
          mi<-which(conditions == Rename_Array$old_label[i])
          if (length(mi) > 0)
          {
            conditions <- conditions[-mi]
          }
          mi<-which(conditions == Rename_Array$new_label[i])
          if (length(mi) == 0)
          {
            conditions <- c(conditions, Rename_Array$new_label[i])
          }
        }
      } 
      conditions <<- conditions
      nConditions<<-length(conditions)
    }
  }
  
  # Rename again the labels in case of label swapping
  # TODO test label swap
  if (AllowLS == T)
  {
    for(i in 1:length(Ls_array$first_label))
    {
      mi1<-which(evidence[, rawfile_col] == Ls_array$selected_raw_file[i] & evidence$label_ == Ls_array$first_label[i])
      mi2<-which(evidence[, rawfile_col] == Ls_array$selected_raw_file[i] & evidence$label_ == Ls_array$second_label[i])
      evidence$label_[mi1] <- as.character(Ls_array$second_label[i])
      evidence$label_[mi2] <- as.character(Ls_array$first_label[i])
    }
  }
  
  
  levellog("",change=-1)
  
  # TODO test background species
  mi<-which(is.na(evidence$label_))
  if(is.na(background_species_lbl))
  {
    if(length(mi) > 0)
    {
      evidence<-evidence[-mi,]
      levellog(paste("read.pgroups: Discarded PSM records due to unassigned label: ",length(mi),sep=""))
    }
  }else{
    evidence[mi,]$label_<-background_species_lbl
  }

  # At this point each hit (evidence record) is matched against a label and its respective intensity (under label_ and Intensity columns respectively)
  

  
  # PS needs to know which replicate coordinate (biorep - trep - frac) each evidence record corresponds to. Since a coordinate is matched
  # against a specific raw file we will merge rep_structure and evidence. rep_structure looks like this:
  # 
  # +--------------------------------------+--------+---------+----------+----------+
  # | raw_file                             | biorep | techrep | fraction | rep_desc |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF01_01 | 1      | 1       | 1        | b1t1f1   |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF02_01 | 1      | 1       | 2        | b1t1f2   |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF03_01 | 1      | 1       | 3        | b1t1f3   |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF04_01 | 1      | 1       | 4        | b1t1f4   |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF05_01 | 1      | 1       | 5        | b1t1f5   |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF06_02 | 1      | 1       | 6        | b1t1f6   |
  # +--------------------------------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF07_01 | 1      | 1       | 7        | b1t1f7   |
  # +--------------------------------------+--------+---------+----------+----------+ (...)
  
  # We will create a new column in rep_structure to erase any fractionation information - we start preparing the data to be in a format that limma will like!
  
  if(length(unique(.GlobalEnv[["rep_structure"]]$fraction)) > 1){
    .GlobalEnv[["rep_structure"]]$old_rep_desc <- .GlobalEnv[["rep_structure"]]$rep_desc
    .GlobalEnv[["rep_structure"]]$rep_desc <- paste0('b',.GlobalEnv[["rep_structure"]]$biorep,'t',.GlobalEnv[["rep_structure"]]$techrep)
  }
  
  # The next line will merge the table above with evidence by raw_file meaning that evidence will be widened: for each record data biorep techrep fraction and rep_desc will
  # be added to the end of evidence according to the Raw file the record was assigned to
  
  levellog("Adding experimental structure information...")
  
  # Add the experimental structure information (merge all columns in evidence with all columns in rep_structure except the old_rep_desc one)
  evidence<-merge(evidence, .GlobalEnv[["rep_structure"]][,!grepl("old_rep_desc", colnames(.GlobalEnv[["rep_structure"]]))], by.x=c(rawfile_col), by.y=c('raw_file'))
  
  # Refresh rep_structure by dropping the rep_desc column
  .GlobalEnv[["rep_structure"]]$rep_desc <- .GlobalEnv[["rep_structure"]]$old_rep_desc
  .GlobalEnv[["rep_structure"]] <- .GlobalEnv[["rep_structure"]][,!grepl("old_rep_desc", colnames(.GlobalEnv[["rep_structure"]]))]
  
  # evidence now looks like this (dropped unnecessary columns - random rows)
  
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | Raw.file                                    | id    | Sequence    | biorep | techrep | fraction | rep_desc |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2013-08-06_DC-SILAC-BR3Fr2_01   | 72068 | QITVNDIPVGR | 3      | 1       | 2        | b3t1     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2013-01-15_DMC-SILAC_HL_NF-02_3 | 91469 | TVAVITSDGR  | 2      | 1       | 2        | b2t1     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2013-01-15_DMC-SILAC_HL_NF-11_4 | 4863  | AIAIGAIQNIR | 2      | 2       | 11       | b2t2     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF05_02        | 74312 | REPTVSSFYVK | 1      | 2       | 5        | b1t2     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF04_01        | 38925 | IDFVDAIK    | 1      | 1       | 4        | b1t1     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF09_02        | 47140 | IIIIDWPEIK  | 1      | 2       | 9        | b1t2     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF10_01        | 51450 | INQIKPGIQYK | 1      | 1       | 10       | b1t1     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF05_01        | 56339 | ITFCTGIR    | 1      | 1       | 5        | b1t1     |
  # +---------------------------------------------+-------+-------------+--------+---------+----------+----------+
  
  # Now make sure that all conditions are represented in the label_ column. Warn the user if one condition is not found and
  # refresh conditions so that it contains only valid conditions
  
  # TODO test filtering
  new_cond_labels <- NULL
  for (cond_i in conditions)
  {
    if (!(cond_i %in% evidence$label_))
    {
      levellog(paste0("Warn User: ", cond_i, " condition was not found in the selected raw files and was not used in comparisons!"))
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
    conditions <- new_cond_labels
    conditions <<- conditions
    nConditions<<-length(conditions)
  }else{
    levellog(paste0("Error User: Not enough labels left, aborting..."))
    return(F)
  }
  
  # Bring Labeled or Label-free data to the following common format (table headers):
  # rep_desc Protein.IDs UniqueSequences.Intensity.condition_1 ... UniqueSequences.Intensity.condition_N Intensity.condition_1 ... Intensity.condition_N
  
  # Now we need to have a unique ID per peptide to merge all evidence information that point to the same peptide, we will call
  # this column Unique.Sequence.ID and in fact it will be the Peptide.ID column from MQ
  
  levellog("read.pgroups: Standarizing data format ...")
  if(!PDdata)
  {
    colnames(evidence)[grepl('^Peptide.ID$',colnames(evidence))]<-'Unique.Sequence.ID'
  }
  
  # If the Peptide.ID did not exist we will create a unique identifier for each peptide from the peptide sequence
  # TODO test datasets where Unique.Sequence.ID is not present
  if (!'Unique.Sequence.ID' %in% colnames(evidence))
  {
    if ('Annotated.Sequence' %in% colnames(evidence))
    {
      levellog("Unique Sequence ID column was not found!")
      colnames(evidence)[colnames(evidence) == 'Annotated.Sequence'] <- 'Unique.Sequence.ID'
      evidence$Unique.Sequence.ID <- sub(".*?\\.(.*?)\\..*", "\\1", evidence$Unique.Sequence.ID)
    }
  }
  
  # The following part of evidence table shows that Unique Sequence ID matches all evidence to a specific peptide
  
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | id    | Unique.Sequence.ID | Sequence  | Intensity | label_ | rep_desc |
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | 3968  | 349                | AGIQFPVGR | 362720000 | H      | b1t1     |
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | 3973  | 349                | AGIQFPVGR | 150410    | H      | b1t1     |
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | 3974  | 349                | AGIQFPVGR | 222510    | H      | b1t1     |
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | 3975  | 349                | AGIQFPVGR | 166730    | H      | b1t1     |
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | 17377 | 1728               | DVNQQEFVR | 191520    | L      | b1t1     |
  # +-------+--------------------+-----------+-----------+--------+----------+
  # | 17378 | 1728               | DVNQQEFVR | 23279     | L      | b1t1     |
  # +-------+--------------------+-----------+-----------+--------+----------+
  
  
  # Rename Intensity.(condition) column to (condition)
  if (!PDdata & !IsobaricLabel)
  {
    colnames(evidence) <- sub('Intensity\\.(.+)', "\\1", colnames(evidence))
  }
  
  # TODO test for PD data and label free
  if(LabelFree)
  {
    if(PDdata)
    {
      
      # Ismini edit: in PD 2.4 there is no Intensity column, instead there is the "PRECURSOR.ABUNDANCE" column
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
    
    # Retrieve the following information for all evidence: Protein ID, Unique sequence ID (so one column for protein and one
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
      evidence.dt<-data.table(evidence[, c('Quan.Usage','Protein.IDs', 'Unique.Sequence.ID', conditions,'rep_desc', 'label_')])
    }else{
      
      # MQ - labelled data:
      # Retrieve the following information for all evidence: Protein ID, Unique sequence ID its Intensity, the condition it derives from (label_) and the repetition it derived from (rep_desc)
      evidence.dt<-data.table(evidence[, c('Protein.IDs', 'Unique.Sequence.ID','rep_desc', 'Intensity', 'label_')])
    }
    setkey(evidence.dt, rep_desc, Protein.IDs, Unique.Sequence.ID)    
  }
  
  # Here evidence.dt contains hits as rows. Each row contains info that matches the respective hit to a specific peptide, to each parent protein group
  # and to each intensity for a specific repetition and condition.
  #
  # An example could be:
  
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | Protein.IDs                                                                  | Unique.Sequence.ID | rep_desc | Intensity | label_ |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 374510    | L      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 6085700   | H      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 4570100   | H      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 497290    | L      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31368.1;ABV71643.1;Q68400|: Membrane   protein UL148;Orf UL148            | 7649               | b1t1     | 4084800   | H      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31368.1;ABV71643.1;Q68400|: Membrane   protein UL148;Orf UL148            | 7649               | b1t1     | 3270300   | L      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31368.1;ABV71643.1;Q68400|: Membrane   protein UL148;Orf UL148            | 7649               | b1t1     | 1305700   | H      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  # | AAR31368.1;ABV71643.1;Q68400|: Membrane   protein UL148;Orf UL148            | 7649               | b1t1     | 835190    | L      |
  # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
  
  # The data table above is the only one that will be used afterwards so we can erase evidence from memory
  rm(evidence)
  
  levellog("Calculating peptide counts...")
  ## Calculate peptide counts per protein, condition and replicate; This procedure will take 4 steps:
  
  # 1. For each condition (per sequence, protein and replicate), set a corresponding column to TRUE if there are > 0 evidence records, (FALSE otherwise)
  #
  # To analyse the line below this comment first take a look at this line:
  #   evidence.dt[, n := .N, by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)]
  #   will create a data table with a new column called n that will show how many instances of each quad {rep_desc, Protein.IDs, Unique.Sequence.ID, label_} exist
  # the example below shows this for the protein [AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein 1] :
  #
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | id     | Protein.IDs                                                                  | Unique.Sequence.ID | rep_desc | label_ | n |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 58850  | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b1t2     | L      | 1 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 1      | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b1t1     | L      | 2 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 4      | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b1t1     | L      | 2 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 117381 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b2t1     | H      | 2 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 117382 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b2t1     | H      | 2 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 137009 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b2t2     | H      | 2 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  # | 137010 | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | b2t2     | H      | 2 |
  # +--------+------------------------------------------------------------------------------+--------------------+----------+--------+---+
  #
  # Instead, executing:
  #   evidence.dt[, .(n=.N), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)]
  #   will do the same thing but will collapse the rows with the same quad {rep_desc, Protein.IDs, Unique.Sequence.ID, label_} as seen below:
  #
  # +----------+------------------------------------------------------------------------------+--------------------+--------+---+
  # | rep_desc | Protein.IDs                                                                  | Unique.Sequence.ID | label_ | n |
  # +----------+------------------------------------------------------------------------------+--------------------+--------+---+
  # | b2t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | H      | 2 |
  # +----------+------------------------------------------------------------------------------+--------------------+--------+---+
  # | b2t2     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | H      | 2 |
  # +----------+------------------------------------------------------------------------------+--------------------+--------+---+
  # | b1t2     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | L      | 1 |
  # +----------+------------------------------------------------------------------------------+--------------------+--------+---+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | L      | 2 |
  # +----------+------------------------------------------------------------------------------+--------------------+--------+---+
  #
  # Since the only thing we want to do is to merge the quads as seen above, the n column is useless so we will fill it with TRUE values executing instead:
  #   evidence.dt[, .(n=.N > 0), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)]
  #   
  # This data table will then be casted (widened) so that the label_ column will be replaced by as many columns as the conditions. These columns will be filled with T ans F values
  # to match each triad {rep_desc, Protein.IDs, Unique.Sequence.ID} to each condition stating if the triad (a specific peptide matched to a specific protein in a certain experiment replicate)
  # was quantified in this condition (for more technical details the 'value' variable is casting is 'n'). An example is:
  #
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | rep_desc | Protein.IDs                                                                               | Unique.Sequence.ID | H    | L     |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1              | 6633               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 1450               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 1852               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 2148               | TRUE | FALSE |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 2188               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 2576               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 3401               | TRUE | FALSE |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 3679               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 3689               | TRUE | TRUE  |
  # +----------+-------------------------------------------------------------------------------------------+--------------------+------+-------+
  #
  # This result is saved in evidence.dt.seqCounts and all the aforementioned procedures are executed in a single line:
  
  evidence.dt.seqCounts<-dcast.data.table(evidence.dt[, .(n=.N > 0), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)], rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=FALSE, value.var = 'n')
  
  # 2. Add a column flagging the common triads {ProteinID - SequenceID - rep_desc}, between at least two labels (note: SD means SubData in data tables)

  evidence.dt.seqCounts[, 'common' := rowSums(.SD) > 1, .SDcols=conditions]
  
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | rep_desc | Protein.IDs                                                                                         | N | H | L | common |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                        | 1 | 1 | 1 | 1      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82;P16755|: Uncharacterized protein UL13                                  | 4 | 4 | 0 | 0      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82|:                                                                      | 1 | 1 | 0 | 0      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1;D3YRW5|: Phosphoprotein 85;Phosphoprotein   UL25;Tegument protein UL25 | 8 | 5 | 8 | 5      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1|: Phosphoprotein 85;Phosphoprotein UL25                                | 1 | 1 | 1 | 1      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | ABV71555.1;P16761;D3YRW5|: Phosphoprotein 85;Phosphoprotein UL25;Tegument   protein UL25            | 1 | 1 | 1 | 1      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  # | b1t1     | ABV71557.1;ACZ79950.1;D3YRW7;P16763|: Uncharacterized protein UL27                                  | 1 | 1 | 0 | 0      |
  # +----------+-----------------------------------------------------------------------------------------------------+---+---+---+--------+
  #
  # Now the number of rows is the same with the number of protein groups
  # TODO: add the option to merge the protein groups to the master protein
  
  # 3. Collapse the records for each protein (per replicate) and count the TRUEs.
  #
  # The following line can be split into two commands:
  #
  # First, evidence.dt.seqCounts[,c(n.Unique.Sequence.IDs=.N), by=.(rep_desc,Protein.IDs)]
  # can produce a data table where all peptides are merged to pairs of {rep_desc, Protein IDs}. A new column (N) will state the amount of peptides that
  # belonged to each one of the pairs:
  #
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | rep_desc | Protein.IDs                                                                                         | N  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                        | 1  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1;D3YRW5|: Phosphoprotein 85;Phosphoprotein   UL25;Tegument protein UL25 | 8  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1|: Phosphoprotein 85;Phosphoprotein UL25                                | 1  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71555.1;P16761;D3YRW5|: Phosphoprotein 85;Phosphoprotein UL25;Tegument   protein UL25            | 1  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71557.1;ACZ79950.1;D3YRW7;P16763|: Uncharacterized protein UL27                                  | 1  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71559.1;P16764;ACZ79951.1;C0H677;D3YRW8|: Uncharacterized protein UL29                           | 11 |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71559.1|:                                                                                        | 1  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  # | b1t1     | ABV71561.1;ACZ79953.1;D3YRX0;P16848|: Protein UL31;Uncharacterized   protein UL31                   | 3  |
  # +----------+-----------------------------------------------------------------------------------------------------+----+
  #
  # Second, mytemp2 = evidence.dt.seqCounts[,c(lapply(.SD, function(x){return(length(which(x)))})), .SDcols=c(conditions, 'common'), by=.(rep_desc, Protein.IDs)]
  # will produce a data table where each pair {rep_desc - Protein ID} will be matched to how many peptides belonging to the pair where identified in each condition and how many in at least two
  # conditions together (Protein IDs truncated)
  #
  # +----------+--------------------------------+----+----+--------+
  # | rep_desc | Protein.IDs                    | H  | L  | common |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;A | 1  | 1  | 1      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q | 2  | 1  | 1      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;K | 1  | 0  | 0      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400;A | 2  | 2  | 2      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400|: | 1  | 1  | 1      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P | 19 | 14 | 14     |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82;P | 4  | 0  | 0      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82|: | 1  | 0  | 0      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1;D | 5  | 8  | 5      |
  # +----------+--------------------------------+----+----+--------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1|: | 1  | 1  | 1      |
  # +----------+--------------------------------+----+----+--------+
  #
  # Notice that common is always <= than min(conditions) which is rational
  #
  # The wider data table from the combined line below looks like:
  #
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | rep_desc | Protein.IDs                                                                                         | N  | H  | L | common |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                        | 1  | 1  | 1 | 1      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82;P16755|: Uncharacterized protein UL13                                  | 4  | 4  | 0 | 0      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82|:                                                                      | 1  | 1  | 0 | 0      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1;D3YRW5|: Phosphoprotein 85;Phosphoprotein   UL25;Tegument protein UL25 | 8  | 5  | 8 | 5      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1|: Phosphoprotein 85;Phosphoprotein UL25                                | 1  | 1  | 1 | 1      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71555.1;P16761;D3YRW5|: Phosphoprotein 85;Phosphoprotein UL25;Tegument   protein UL25            | 1  | 1  | 1 | 1      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71557.1;ACZ79950.1;D3YRW7;P16763|: Uncharacterized protein UL27                                  | 1  | 1  | 0 | 0      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  # | b1t1     | ABV71559.1;P16764;ACZ79951.1;C0H677;D3YRW8|: Uncharacterized protein UL29                           | 11 | 11 | 9 | 9      |
  # +----------+-----------------------------------------------------------------------------------------------------+----+----+---+--------+
  #
  # Also notice that N >= than max(conditions) and of course N > common that is rational too.
  
  evidence.dt.seqCounts<-evidence.dt.seqCounts[,c(n.Unique.Sequence.IDs=.N, lapply(.SD, function(x){return(length(which(x)))})), by=.(rep_desc,Protein.IDs), .SDcols=c(conditions, 'common')]
  
  # 4. Calculate the percentage columns
  #
  # For each one of the conditions calculate the percentage of the {rep_desc - protein} pair's peptides that were quantified in the specific condition
  #
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | rep_desc | Protein.IDs                                                                                         | N  | UniqueSequences.H | UniqueSequences.L | common | Hp       | Lp       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                        | 1  | 1                 | 1                 | 1      | 50       | 50       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400;ACZ80031.1|: Membrane protein UL148;Orf   UL148                        | 2  | 2                 | 2                 | 2      | 50       | 50       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400|: Membrane protein UL148;Orf UL148                                     | 1  | 1                 | 1                 | 1      | 50       | 50       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22           | 19 | 19                | 14                | 14     | 57.57576 | 42.42424 |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82;P16755|: Uncharacterized protein UL13                                  | 4  | 4                 | 0                 | 0      | 100      | 0        |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71543.1;ACZ79936.1;E2RU82|:                                                                      | 1  | 1                 | 0                 | 0      | 100      | 0        |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1;D3YRW5|: Phosphoprotein 85;Phosphoprotein   UL25;Tegument protein UL25 | 8  | 5                 | 8                 | 5      | 38.46154 | 61.53846 |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71555.1;P16761;ACZ79948.1|: Phosphoprotein 85;Phosphoprotein UL25                                | 1  | 1                 | 1                 | 1      | 50       | 50       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71555.1;P16761;D3YRW5|: Phosphoprotein 85;Phosphoprotein UL25;Tegument   protein UL25            | 1  | 1                 | 1                 | 1      | 50       | 50       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71557.1;ACZ79950.1;D3YRW7;P16763|: Uncharacterized protein UL27                                  | 1  | 1                 | 0                 | 0      | 100      | 0        |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71559.1;P16764;ACZ79951.1;C0H677;D3YRW8|: Uncharacterized protein UL29                           | 11 | 11                | 9                 | 9      | 55       | 45       |
  # +----------+-----------------------------------------------------------------------------------------------------+----+-------------------+-------------------+--------+----------+----------+
  #
  # (the table has renamed columns as produced by the command below)
  
  evidence.dt.seqCounts[, paste0(conditions,'p') := lapply(.SD, function(x){return((x/sum(.SD))*100)}), by=.(rep_desc,Protein.IDs), .SDcols=c(conditions)]
  
  # Rename the peptide counts columns
  setnames(evidence.dt.seqCounts,colnames(evidence.dt.seqCounts)[which(colnames(evidence.dt.seqCounts) %in% conditions)],paste('UniqueSequences',conditions,sep='.')) 
  
  levellog("Calculating protein intensity...")
  
  # Calculate the protein intensity = (sum of unique peptide intensities) for each condition - replicate in the following steps
  # TODO test for label free and PD
  if(LabelFree)
  {
    # Cast the data so that we have columns for each label and intensity separately
    evidence.dt<-dcast.data.table(evidence.dt, rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=0)    
  }else{
    if(PDdata)
    {
      evidence.dt<-evidence.dt[Quan.Usage == 'Used' | Quan.Usage == 'Use', lapply(.SD, max), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID), .SDcols=conditions]    
    }else{
      # MQ labelled data
      # evidence dt may have many rows per {replicate - protein - peptide - label} quad that means that a peptide belonging to a protein may have been quantified more than once in a replicate
      # (possibly in more than one fractions or in more than one modifications or charges) and condition. We assume that the sum of these intensities is the most accurate quanttification value for this peptide
      # aggregating these intensities will do the work:
      #
      # The aggregation in the line below will for example merge these evidence intensities to their sum - peptide intensity:
      #
      # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
      # | Protein.IDs                                                                  | Unique.Sequence.ID | rep_desc | Intensity | label_ |
      # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
      # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 6085700   | H      |
      # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
      # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 4570100   | H      |
      # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
      # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 6633               | b1t1     | 3864500   | H      |
      # +------------------------------------------------------------------------------+--------------------+----------+-----------+--------+
      #
      # +----------+------------------------------------------------------------------------------+--------------------+--------+-----------+
      # | rep_desc | Protein.IDs                                                                  | Unique.Sequence.ID | label_ | Intensity |
      # +----------+------------------------------------------------------------------------------+--------------------+--------+-----------+
      # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1 | 6633               | H      | 14520300  |
      # +----------+------------------------------------------------------------------------------+--------------------+--------+-----------+
      #
      
      
      evidence.dt <- evidence.dt[, .(Intensity = sum(Intensity)), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)]
      
      # We will cast the data to widen them once more:
      #
      # evidence.dt will become:
      #
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | rep_desc | Protein.IDs                                                                               | Unique.Sequence.ID | H        | L       |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1              | 6633               | 14520300 | 871800  |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | AAR31368.1;ABV71643.1;Q68400;ACZ80031.1|: Membrane protein UL148;Orf   UL148              | 2400               | 4408400  | 3374800 |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | AAR31368.1;ABV71643.1;Q68400;ACZ80031.1|: Membrane protein UL148;Orf   UL148              | 9860               | 2110900  | 2635200 |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | AAR31368.1;ABV71643.1;Q68400|: Membrane protein UL148;Orf UL148                           | 7649               | 6010890  | 4459470 |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 12                 | 41383670 | 2524151 |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 1050               | 11698700 | 611150  |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclear protein   HWLF1;Tegument protein US22 | 1450               | 27244455 | 2502488 |
      # +----------+-------------------------------------------------------------------------------------------+--------------------+----------+---------+
      #
      # Notice the sum 14520300 validating the example above.
      
      evidence.dt<-dcast.data.table(evidence.dt, rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=0, value.var = 'Intensity')  
      
      # The following lines are obsolete
      # # Get another column (sumI) with the sum of intensities in all conditions: 
      # 
      # evidence.dt[, sumI := rowSums(.SD, na.rm = T), .SDcols=conditions]
      # 
      # evidence.dt <- evidence.dt[sumI > 0, lapply(.SD, max), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID), .SDcols=conditions]    
      # evidence.dt[, sumI := NULL]
    }
  }
  
  # Here, each evidence.dt row represents a detected peptide in a specific experiment replicate. Let's sort the unique sequecne IDs to have a look
  # at our missing values (Protein IDs truncated):
  #
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | rep_desc | Unique.Sequence.ID | H        | L        | Protein.IDs                                        |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b1t1     | 0                  | 43635927 | 24237101 | H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446| |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b1t2     | 0                  | 40498440 | 20943890 | H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446| |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b2t1     | 0                  | 88589670 | 45784620 | H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446| |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b2t2     | 0                  | 80149090 | 42010380 | H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446| |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b3t1     | 0                  | 2134140  | 986747   | H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446| |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b3t2     | 0                  | 22462057 | 14065705 | H-INV:HIT000035131;Q86U42;ENSEMBL:ENSP00000380446| |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b2t1     | 2                  | 5145100  | 5252100  | REFSEQ:NP_001104262;TREMBL:C6FGV0;Q59FJ6|: Methyl- |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b1t1     | 3                  | 3290000  | 2731700  | REFSEQ:NP_001104262|: Methyl-CpG-binding protein 2 |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b1t2     | 3                  | 2532900  | 2199800  | REFSEQ:NP_001104262|: Methyl-CpG-binding protein 2 |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b1t1     | 7                  | 96692    | 0        | H-INV:HIT000220955;Q92922;ENSEMBL:ENSP00000414266; |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b1t2     | 7                  | 116640   | 56799    | H-INV:HIT000220955;Q92922;ENSEMBL:ENSP00000414266; |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b2t1     | 7                  | 2619680  | 1468100  | H-INV:HIT000220955;Q92922;ENSEMBL:ENSP00000414266; |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b2t2     | 7                  | 2386500  | 1309500  | H-INV:HIT000220955;Q92922;ENSEMBL:ENSP00000414266; |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b3t1     | 7                  | 771920   | 473710   | H-INV:HIT000220955;Q92922;ENSEMBL:ENSP00000414266; |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  # | b3t2     | 7                  | 725750   | 387590   | H-INV:HIT000220955;Q92922;ENSEMBL:ENSP00000414266; |
  # +----------+--------------------+----------+----------+----------------------------------------------------+
  #
  # Since this is a SILAC example most of the peptides are detected in both conditions (an LFQ experiment could have more peptides detected not in all conditions)
  # Indeed, 10.2% of the peptides in evidence.dt are detected in a single condition and 5.1% of the intensity cells are missing values. 
  # Also not all peptides are detected in all repetitions.. In fact aggregating all peptides by times of appearance (length - in this experiment it should be 1-6)
  # and then aggregating the times of appearance by legth again will give us how many peptides were detected in 1 sample how many in 2 etc. and we can compute the respective
  # percentages as follows:
  #
  # +---------------------+------------+
  # | times of appearance | percentage |
  # +---------------------+------------+
  # | 1                   | 0.22       |
  # +---------------------+------------+
  # | 2                   | 0.35       |
  # +---------------------+------------+
  # | 3                   | 0.09       |
  # +---------------------+------------+
  # | 4                   | 0.17       |
  # +---------------------+------------+
  # | 5                   | 0.04       |
  # +---------------------+------------+
  # | 6                   | 0.13       |
  # +---------------------+------------+
  #
  # Only 13% of the peptides are detected in all repetitions of the experiment! This can be due to variance of biological or technical reasons
  # 
  # The following lines perform imputation of the missing values.
  # The intensities detected show a left skewed distribution as expected.. According to "A comparative study of evaluating missing value imputation methods in label???free proteomic" (Jin et al.)
  # seven commonly used imputation methods are: lowest of detection (LOD), random drawing from a left-censored normal distribution (ND), k-nearest neighbors (kNN), local least squares (LLS), random
  # forest (RF), singular value decomposition (SVD) and Bayesian principal component analysis (BPCA)
  #
  # TODO: add more imputation algorithms based on "Missing value imputation" section of the aforementioned paper
  
  # Here we apply the LOD (lowest of detection) imputation algorithm:
  
  # Get a vector of unique peptides intensities
  tmp.I<-sort(unique(evidence.dt[,get(conditions)]))
  
  # If the minimum intensity is zero
  if(tmp.I[1] == 0)
  {
    # Replace 0's with minimum intensity
    minI<-tmp.I[2]
    evidence.dt[, (conditions) := lapply(.SD, function(x){ t<-which(x == 0); if(length(t) > 0){x[t] <- minI}; return(x) }), .SDcols=conditions]
  }else{
    minI<-tmp.I[1]
  }
  
  # If enabled, do filter out peptides where all conditions except the condition 'filterL_lbl' have noise-level intensity (peptide-level filtering)
  # This block is used to exclude background peptides as described in methods such as pSILAC. In pSILAC we usually let the production of aminoacids happen in the light label
  # and then apply different conditions to our samples for a short period of time. The peptides that were not affected by the pulses should have high abundances
  # in the light label and low (noise level) in the other ones and should be excluded by the analysis:
  if(filterL && filterL_lvl)
  {
    # The following line counts how many conditions (except the background condition - filterL_lbl) have noise level intensity (noise level is considered as only
    # the lowest intensity detected) - the result is written in the minIcount column
    # TODO: is the noise calculating algorithm correct? (noise level is considered as only the lowest intensity detected)
    
    evidence.dt[, minIcount := rowSums(.SD == minI), .SDcols=conditions[! conditions %in% filterL_lbl]]
    n1<-nrow(evidence.dt)
    
    # Drop all rows peptides that have noise level intensity in all conditions except the background condition
    evidence.dt<-evidence.dt[minIcount < (nConditions - 1)]
    n2<-nrow(evidence.dt)
    if(n2 < n1)
    {
      levellog(paste0("read.pgroups: Filtered out ", (n1-n2)," peptides having noise-level intensity in all channels except the '", filterL_lbl,"' channel ..."));
    }
    
    # Drop the minIcount column that is no longer used
    evidence.dt[, minIcount := NULL]
  }
  
  # 2. Calculate the protein intensity (sum of peptide intensities) for each condition and replicate
  # Also count the number of quantifiable peptides (those which do not have intensity NA)
  
  # TODO: test fol LFQ
  if(LabelFree)
  {
    # Top three peptides in abundance
    #evidence.dt<-evidence.dt[, lapply(.SD, function(x){x<-x[!is.na(x)]; x<-sort(x, decreasing<-T); if(length(x)<3){return(sum(x))}else{return(sum(x[1:3]))}}), by=.(rep_desc, Protein.IDs), .SDcols=conditions]
    evidence.dt<-evidence.dt[, c(n=.N, nas=length(which(is.na(.SD))) ,lapply(.SD, function(x){x<-x[!is.na(x)]; x<-sort(x, decreasing<-T); if(length(x)<3){return(sum(x))}else{return(sum(x[1:3]))}})), by=.(rep_desc, Protein.IDs), .SDcols=conditions]
  }else{
    # All peptides
    # The following command aggregates the intensity in all condition of all peptides of the same protein in a specifi replicate by adding all the intensities of the respective peptides
    # It also creates a column that counts how many NAs there are in intensity values for each pair {protein - replicate}
    evidence.dt<-evidence.dt[, c(n=.N, nas=length(which(is.na(.SD))) ,lapply(.SD, sum, na.rm = T)), by=.(rep_desc, Protein.IDs), .SDcols=conditions] 
  }
  evidence.dt[, N := NULL]
  # Rename the intensity columns form [condition] to Intensity.[condition]
  setnames(evidence.dt,colnames(evidence.dt)[which(colnames(evidence.dt) %in% conditions)],paste('Intensity',conditions,sep='.'))
  
  # evidence.dt now seems like:
  #
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | rep_desc | Protein.IDs                                        | nas | Intensity.H | Intensity.L |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa i | 0   | 14520300    | 871800      |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 6315898     | 290628.7    |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1;AAR | 0   | 3846360     | 4330.7      |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1|: 5 | 0   | 42086141    | 804846.7    |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31362.1;ABV71635.1;P19893;Q6SWJ2;ACZ80024.1|: 4 | 0   | 4330.7      | 2835300     |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400;ACZ80031.1|: Membrane | 0   | 6519300     | 6010000     |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400|: Membrane protein UL | 0   | 6010890     | 4459470     |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclea | 0   | 2.25E+08    | 12415101    |
  # +----------+----------------------------------------------------+-----+-------------+-------------+
  # 
  # Merge with the evidence.dt.seqCounts table
  # both evidence.dt and evidence.dt.seqCounts have the pair {Protein.ID - rep_desc} as key so per line and for each of these pairs all the other columnns will be simply put together
  
  evidence.dt<-merge(evidence.dt, evidence.dt.seqCounts)
  
  # evidence.dt starts building up...
  #
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | rep_desc | Protein.IDs                                        | nas | Intensity.H | Intensity.L | N  | UniqueSequences.H | UniqueSequences.L | common | Hp       | Lp       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa i | 0   | 14520300    | 871800      | 1  | 1                 | 1                 | 1      | 50       | 50       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 6315898     | 290628.7    | 2  | 2                 | 1                 | 1      | 66.66667 | 33.33333 |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1;AAR | 0   | 3846360     | 4330.7      | 1  | 1                 | 0                 | 0      | 100      | 0        |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1|: 5 | 0   | 42086141    | 804846.7    | 3  | 3                 | 2                 | 2      | 60       | 40       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31362.1;ABV71635.1;P19893;Q6SWJ2;ACZ80024.1|: 4 | 0   | 4330.7      | 2835300     | 1  | 0                 | 1                 | 0      | 0        | 100      |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400;ACZ80031.1|: Membrane | 0   | 6519300     | 6010000     | 2  | 2                 | 2                 | 2      | 50       | 50       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | AAR31368.1;ABV71643.1;Q68400|: Membrane protein UL | 0   | 6010890     | 4459470     | 1  | 1                 | 1                 | 1      | 50       | 50       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  # | b1t1     | ABV71513.1;ACZ80070.1;D3YS51;P09722|: Early nuclea | 0   | 2.25E+08    | 12415101    | 19 | 19                | 14                | 14     | 57.57576 | 42.42424 |
  # +----------+----------------------------------------------------+-----+-------------+-------------+----+-------------------+-------------------+--------+----------+----------+
  
  # Add the experimental structure information to evidence.dt based on rep_desc (raw file at this point has no information and is dropped)
  # First drop the fraction and rawfile columns from rep_structure - also rename repdesc to indicate the brep and trep only
  #
  # +--------+---------+----------+
  # | biorep | techrep | rep_desc |
  # +--------+---------+----------+
  # | 1      | 1       | b1t1     |
  # +--------+---------+----------+
  # | 1      | 2       | b1t2     |
  # +--------+---------+----------+
  # | 2      | 1       | b2t1     |
  # +--------+---------+----------+
  # | 2      | 2       | b2t2     |
  # +--------+---------+----------+
  # | 3      | 1       | b3t1     |
  # +--------+---------+----------+
  # | 3      | 2       | b3t2     |
  # +--------+---------+----------+
  
  tmp.rep_struct<-.GlobalEnv[["rep_structure"]][! duplicated(.GlobalEnv[["rep_structure"]][,c('biorep','techrep')]), !grepl('raw_file', colnames(.GlobalEnv[["rep_structure"]])) & !grepl('fraction', colnames(.GlobalEnv[["rep_structure"]]) )]
  tmp.rep_struct$rep_desc<-paste0('b',tmp.rep_struct$biorep,'t',tmp.rep_struct$techrep)
  
  # Add information for bioreps and treps by merging the evidence.dt and the temporary rep struct table shown above with common key being "rep_desc"
  evidence.dt<-merge(evidence.dt, data.table(tmp.rep_struct), by='rep_desc')
  
  # If enabled, filter out proteins based on percentage labeling for the desired label (protein-level filtering)
  # This is also used in pSILAC to exclude proteins that were not de novo produced during the pulses
  if(filterL && !filterL_lvl)
  {
    n1<-length(unique(evidence.dt[get(paste0(filterL_lbl,"p")) == 100.0]$Protein.IDs))
    
    # The following command simply drops the proteins that were not entirely quantified in the background label
    evidence.dt<-evidence.dt[get(paste0(filterL_lbl,"p")) < 100.0]
    levellog(paste0("read.pgroups: Filtered out ", n1," proteins which where identified solely by '", filterL_lbl, "' peptides ..."));
  }
  
  # Get protein IDs that match the following eligibility criteria: they were quantified with a total of at least 'nRequiredLeastPeps' different peptides
  # accross at least 'nRequiredLeastBioreps' biological replicates.
  #
  # Let's explain the following line:
  # Protein.Biorep.aggregation <- evidence.dt[, .(c1 = sum(N-nas)) , by=.(Protein.IDs, biorep)] is a data table that aggregates all proteins by their biological replicate and computes a new variable called c1
  # which is the sum of all their N (peptides assigned to them) minus their NA values
  #
  # This transforms evidence.dt from:
  #
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | rep_desc | Protein.IDs                                        | nas | Intensity.H | Intensity.L | N | UniqueSequences.H | UniqueSequences.L | common | Hp       | Lp       | biorep | techrep |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa i | 0   | 14520300    | 871800      | 1 | 1                 | 1                 | 1      | 50       | 50       | 1      | 1       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b1t2     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa i | 0   | 10515500    | 300750      | 1 | 1                 | 1                 | 1      | 50       | 50       | 1      | 2       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b2t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa i | 0   | 5686600     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 1       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b2t2     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa i | 0   | 4669800     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 2       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 6315898     | 290628.7    | 2 | 2                 | 1                 | 1      | 66.66667 | 33.33333 | 1      | 1       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b1t2     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 7997380     | 63117.7     | 2 | 2                 | 1                 | 1      | 66.66667 | 33.33333 | 1      | 2       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b2t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 2115170     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 1       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b2t2     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 1054030     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 2       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  # | b3t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV | 0   | 8142900     | 3170900     | 1 | 1                 | 1                 | 1      | 50       | 50       | 3      | 1       |
  # +----------+----------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+
  #
  # to: (Protein ID trunction makes it a little difficult to show the aggregation of the table above to the table below...)
  # +------------------------------------------------------------------------------+--------+----+
  # | Protein.IDs                                                                  | biorep | c1 |
  # +------------------------------------------------------------------------------+--------+----+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1 | 1      | 2  |
  # +------------------------------------------------------------------------------+--------+----+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa i                         | 2      | 2  |
  # +------------------------------------------------------------------------------+--------+----+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV                           | 1      | 4  |
  # +------------------------------------------------------------------------------+--------+----+
  #
  # Notice that now the key is the pair {Protein ID - biorep} so each line represents a protein as detected in whole on a single biological replicate
  #
  # Protein.Least.Required.Stats <- Protein.Biorep.aggregation[, .(nQuantPeps = sum(c1), geqXnRequiredLeastBioreps = .N >= .GlobalEnv[["nRequiredLeastBioreps"]]), by=.(Protein.IDs)]
  # will combine data from all bioreps for each protein and compute: all the peptide detections for the protein by aggregating by sum (storing the result to nQuantPeps) and
  # the number of bioreps where the protein was by computing N (length of aggregation i.e. how many lines aggregated to a single row) a new column (geqXnRequiredLeastBioreps) will
  # be set to false if N >= LeastRequiredBreps
  #
  # An example of the result is:
  #
  # +------------------------------------------------------+------------+---------------------------+
  # | Protein.IDs                                          | nQuantPeps | geqXnRequiredLeastBioreps |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa i | 4          | TRUE                      |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV   | 7          | TRUE                      |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1;AAR   | 2          | FALSE                     |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1|:   5 | 12         | TRUE                      |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31362.1;ABV71635.1;P19893;Q6SWJ2;ACZ80024.1|:   4 | 2          | FALSE                     |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31368.1;ABV71643.1;Q68400;ACZ80031.1|:   Membrane | 9          | TRUE                      |
  # +------------------------------------------------------+------------+---------------------------+
  # | AAR31368.1;ABV71643.1;Q68400|: Membrane   protein UL | 2          | FALSE                     |
  # +------------------------------------------------------+------------+---------------------------+
  # | ABV71513.1;ACZ80070.1;D3YS51;P09722|:   Early nuclea | 77         | TRUE                      |
  # +------------------------------------------------------+------------+---------------------------+
  #
  # Protein.Least.Required.Stats[nQuantPeps >= .GlobalEnv[["nRequiredLeastPeps"]] & geqXnRequiredLeastBioreps == T] will get only the proteins that match the eligibility criteria
  # The following line summarizes all the aforementioned commands to get a vector of the proteins that match all eligibility criteria:
  
  
  Protein.IDs.quant <- evidence.dt[, .(c1 = sum(N-nas)) , by=.(Protein.IDs, biorep)][, .(nQuantPeps = sum(c1), geqXnRequiredLeastBioreps = .N >= .GlobalEnv[["nRequiredLeastBioreps"]]), by=.(Protein.IDs)][nQuantPeps >= .GlobalEnv[["nRequiredLeastPeps"]] & geqXnRequiredLeastBioreps == T]$Protein.IDs
  levellog(paste0("read.pgroups: Filtered out ", (length(unique(evidence.dt$Protein.IDs)) - length(Protein.IDs.quant))," proteins which were not identified in at least ",nRequiredLeastBioreps," biological replicate(s) with at least a total of ",nRequiredLeastPeps," peptide(s)"));
  
  # Add the nQuantPeps in evidence.dt to show hoe many peptides were detected per protein
  evidence.dt[,nQuantPeps := N-nas]
  
  # Keep only the proteins that match the eligibility criteria
  evidence.dt<-evidence.dt[Protein.IDs %in% Protein.IDs.quant]

  # Cast the table to the following format
  # Protein.IDs Intensity.[<rep_desc_X>.<condition_Y> ...] [<rep_desc_X>.Ratio.counts ...] [<rep_desc_X>.uniqueSequences ...] exp_desc [<condition_Y> ...] [<condition_Y>p ...]
  # This will be the final form of our data (data table repdesc)
  
  # Step 1: For each 'rep_desc', add the evidence.dt data to the growing dataframe pgroups
  # Also, calculate the missing columns required by the target format and drop the unnecessary columns
  
  # evidence.dt has now information only for the proteins that meet the eligibility criteria and the key is the triad {rep_desc - ProteinID - Unique.sequence.ID} (each line
  # represents a peptide respective to a protein detected in a specific replication of the experiment)
  #
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | rep_desc | Protein.IDs                                                                                                                                                                                                 | nas | Intensity.H | Intensity.L | N | UniqueSequences.H | UniqueSequences.L | common | Hp       | Lp       | biorep | techrep | nQuantPeps |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                                                                                                                                | 0   | 14520300    | 871800      | 1 | 1                 | 1                 | 1      | 50       | 50       | 1      | 1       | 1          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b1t2     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                                                                                                                                | 0   | 10515500    | 300750      | 1 | 1                 | 1                 | 1      | 50       | 50       | 1      | 2       | 1          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b2t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                                                                                                                                | 0   | 5686600     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 1       | 1          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b2t2     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                                                                                                                                | 0   | 4669800     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 2       | 1          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV71635.1;P19893;Q6SWJ2|:   55 kDa immediate-early protein 1;Regulatory protein IE1;45 kDa   immediate-early protein 2;Protein UL122;Regulatory protein IE2 | 0   | 6315898     | 290628.7    | 2 | 2                 | 1                 | 1      | 66.66667 | 33.33333 | 1      | 1       | 2          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b1t2     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV71635.1;P19893;Q6SWJ2|:   55 kDa immediate-early protein 1;Regulatory protein IE1;45 kDa   immediate-early protein 2;Protein UL122;Regulatory protein IE2 | 0   | 7997380     | 63117.7     | 2 | 2                 | 1                 | 1      | 66.66667 | 33.33333 | 1      | 2       | 2          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  # | b2t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV71635.1;P19893;Q6SWJ2|:   55 kDa immediate-early protein 1;Regulatory protein IE1;45 kDa   immediate-early protein 2;Protein UL122;Regulatory protein IE2 | 0   | 2115170     | 4330.7      | 1 | 1                 | 0                 | 0      | 100      | 0        | 2      | 1       | 1          |
  # +----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+-------------+-------------+---+-------------------+-------------------+--------+----------+----------+--------+---------+------------+
  #
  # Start by getting all the unique Protein IDs
  #
  
  setkey(evidence.dt, Protein.IDs)
  pgroups<-data.frame(Protein.IDs = unique(evidence.dt$Protein.IDs))
  setkey(evidence.dt, rep_desc)
  
  # For each replication of the experiment:
  # The comments below wwill assume that rep_desc_i == "b1t1"
  for(rep_desc_i in unique(evidence.dt$rep_desc))
  {
    # Get all rows of evidence.dt that correspond to b1t1 replicate
    rep_desc_i_pgroups<-data.frame(evidence.dt[rep_desc == rep_desc_i,])
    
    # Rename the columns
    allcols<-colnames(rep_desc_i_pgroups)
    
    # allcols before renaming:  [1] "rep_desc"    "Protein.IDs"  "nas"   "Intensity.H"   "Intensity.L"   "N"    "UniqueSequences.H" "UniqueSequences.L" "common"           
    # "Hp"     "Lp"       "biorep"    "techrep"  "nQuantPeps"  
    
    # Rename Intensity columns
    colsl<-grepl('^Intensity' ,allcols)
    colnames(rep_desc_i_pgroups)[colsl]<-gsub("^Intensity(.+)$",paste("Intensity\\1",rep_desc_i,sep='.'), allcols[colsl])
    
    # Rename UniqueSequences columns
    colsl<-grepl('^UniqueSequences' ,allcols)
    colnames(rep_desc_i_pgroups)[colsl]<-gsub("^UniqueSequences(.+)$",paste(rep_desc_i,"uniqueSequences\\1",sep='.'), allcols[colsl])
    
    
    # Add new column <rep_desc_X>.uniqueSequences
    rep_desc_i_pgroups[, paste(rep_desc_i,'uniqueSequences',sep='.')]<-rowSums(rep_desc_i_pgroups[, colnames(rep_desc_i_pgroups)[colsl]])
    
    # The new column has the total amount of peptides attributed to this prtein in all conditions for the specific replication of the experiment:
    #
    # +----------+------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+------------------------+------------------------+----------------------+--------+----------+----------+--------+---------+------------+
    # | rep_desc | Protein.IDs                                                                                                | nas | Intensity.H.b1t1 | Intensity.L.b1t1 | N | b1t1.uniqueSequences.H | b1t1.uniqueSequences.L | b1t1.uniqueSequences | common | Hp       | Lp       | biorep | techrep | nQuantPeps |
    # +----------+------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+------------------------+------------------------+----------------------+--------+----------+----------+--------+---------+------------+
    # | b1t1     | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                               | 0   | 14520300         | 871800           | 1 | 1                      | 1                      | 2                    | 1      | 50       | 50       | 1      | 1       | 1          |
    # +----------+------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+------------------------+------------------------+----------------------+--------+----------+----------+--------+---------+------------+
    # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV71635.1;P19893;Q6SWJ2|:   55 kDa immediate-early (...)   | 0   | 6315898          | 290628.7         | 2 | 2                      | 1                      | 3                    | 1      | 66.66667 | 33.33333 | 1      | 1       | 2          |
    # +----------+------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+------------------------+------------------------+----------------------+--------+----------+----------+--------+---------+------------+
    # | b1t1     | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1|: 55 kDa immediate-early   protein 1;Regulatory protein IE1 | 0   | 42086141         | 804846.7         | 3 | 3                      | 2                      | 5                    | 2      | 60       | 40       | 1      | 1       | 3          |
    # +----------+------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+------------------------+------------------------+----------------------+--------+----------+----------+--------+---------+------------+
    
    # Rename 'p' (percentage) cols
    colsl<-allcols %in% paste0(conditions,'p')
    colnames(rep_desc_i_pgroups)[colsl]<-gsub("^(.+)$",paste("\\1",rep_desc_i,sep='.'), allcols[colsl])
    
    # Rename the 'nQuantPeps' column to <rep_desc_i>.Ratio.counts
    colsl<-allcols %in% c('nQuantPeps')
    colnames(rep_desc_i_pgroups)[colsl]<-paste(rep_desc_i,'Ratio.counts',sep='.')
    
    # column names after renaming:
    
    #  "rep_desc"               "Protein.IDs"            "nas"                    "Intensity.H.b1t1"       "Intensity.L.b1t1"       "N"                      "b1t1.uniqueSequences.H"
    #  "b1t1.uniqueSequences.L" "common"                 "Hp.b1t1"                "Lp.b1t1"                "biorep"                 "techrep"                "b1t1.Ratio.counts"     
    #  "b1t1.uniqueSequences" 
    
    # Merge with the growing data frame
    
    # First find which columns are present in pgroups and rep_desc_i_pgroups at the same time (in the beggining it should be Protein.IDs only)
    cc<-intersect(names(pgroups), names(rep_desc_i_pgroups))
    
    # Extend proups by merging itself with rep_desc_i_pgroups excluding some columns
    pgroups<-merge(pgroups, rep_desc_i_pgroups[, ! colnames(rep_desc_i_pgroups) %in% c('biorep', 'techrep', 'fraction', 'rep_desc', cc[! grepl('Protein.IDs', cc)] )], all.x = T)
  }
  
  # pgroups now should look like:
  #
  # +-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------+------------------+----------+------------------------+------------------------+--------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+
  # | Protein.IDs | nas                                                                                                                                                                                                         | Intensity.H.b1t1 | Intensity.L.b1t1 | N        | b1t1.uniqueSequences.H | b1t1.uniqueSequences.L | common | Hp.b1t1 | Lp.b1t1  | b1t1.Ratio.counts | b1t1.uniqueSequences | Intensity.H.b1t2 | Intensity.L.b1t2 | b1t2.uniqueSequences.H | b1t2.uniqueSequences.L | Hp.b1t2 | Lp.b1t2  | b1t2.Ratio.counts | b1t2.uniqueSequences | Intensity.H.b2t1 | Intensity.L.b2t1 | b2t1.uniqueSequences.H | b2t1.uniqueSequences.L | Hp.b2t1 | Lp.b2t1 | b2t1.Ratio.counts | b2t1.uniqueSequences | Intensity.H.b2t2 | Intensity.L.b2t2 | b2t2.uniqueSequences.H | b2t2.uniqueSequences.L | Hp.b2t2 | Lp.b2t2 | b2t2.Ratio.counts | b2t2.uniqueSequences | Intensity.H.b3t1 | Intensity.L.b3t1 | b3t1.uniqueSequences.H | b3t1.uniqueSequences.L | Hp.b3t1 | Lp.b3t1 | b3t1.Ratio.counts | b3t1.uniqueSequences | Intensity.H.b3t2 | Intensity.L.b3t2 | b3t2.uniqueSequences.H | b3t2.uniqueSequences.L | Hp.b3t2 | Lp.b3t2 | b3t2.Ratio.counts | b3t2.uniqueSequences |
  # +-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------+------------------+----------+------------------------+------------------------+--------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+-----------+----------+
  # | 1           | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|: 55 kDa immediate-early protein   1                                                                                                                                | 0                | 14520300         | 871800   | 1                      | 1                      | 1      | 1       | 50       | 50                | 1                    | 2                | 10515500         | 300750                 | 1                      | 1       | 50       | 50                | 1                    | 2                | 5686600          | 4330.7                 | 1                      | 0       | 100     | 0                 | 1                    | 1                | 4669800          | 4330.7                 | 1                      | 0       | 100     | 0                 | 1                    | 1                | NA               | NA                     | NA                     | NA      | NA      | NA                | NA                   | NA               | NA               | NA                     | NA                     | NA      | NA      | NA                | NA        | NA       |
  # +-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------+------------------+----------+------------------------+------------------------+--------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+-----------+----------+
  # | 2           | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV71635.1;P19893;Q6SWJ2|:   55 kDa immediate-early protein 1;Regulatory protein IE1;45 kDa   immediate-early protein 2;Protein UL122;Regulatory protein IE2 | 0                | 6315898          | 290628.7 | 2                      | 2                      | 1      | 1       | 66.66667 | 33.33333          | 2                    | 3                | 7997380          | 63117.7                | 2                      | 1       | 66.66667 | 33.33333          | 2                    | 3                | 2115170          | 4330.7                 | 1                      | 0       | 100     | 0                 | 1                    | 1                | 1054030          | 4330.7                 | 1                      | 0       | 100     | 0                 | 1                    | 1                | 8142900          | 3170900                | 1                      | 1       | 50      | 50                | 1                    | 2                | NA               | NA                     | NA                     | NA      | NA      | NA                | NA        | NA       |
  # +-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------+------------------+----------+------------------------+------------------------+--------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+-----------+----------+
  # | 3           | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1|: 55 kDa immediate-early   protein 1;Regulatory protein IE1                                                                                                  | 0                | 42086141         | 804846.7 | 3                      | 3                      | 2      | 2       | 60       | 40                | 3                    | 5                | 32540397         | 1955784                | 3                      | 2       | 60       | 40                | 3                    | 5                | 17497600         | 8661.4                 | 2                      | 0       | 100     | 0                 | 2                    | 2                | 10065000         | 8661.4                 | 2                      | 0       | 100     | 0                 | 2                    | 2                | 15538740         | 4850292                | 1                      | 1       | 50      | 50                | 1                    | 2                | 16142000         | 5161000                | 1                      | 1       | 50      | 50                | 1         | 2        |
  # +-------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------+------------------+----------+------------------------+------------------------+--------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+----------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+----------------------+------------------+------------------+------------------------+------------------------+---------+---------+-------------------+-----------+----------+
  #
  # Representing all data in a very wide format
  
  # Step 2: Calculate the columns [<label/condition_Y> ...] containing the number of peptides found per condition in all replicates - this is a simple sum of the peptide count in all replicates for each condition
  allcols<-colnames(pgroups)
  for(cond_i in conditions)
  {
    colsl<-grepl(paste('uniqueSequences', cond_i,sep='\\.') ,allcols)
    pgroups[, cond_i]<-rowSums(pgroups[, allcols[colsl]])
  }
  
  # Step 3: Calculate the columns [<label/condition_Y>p ...] containing the percentage of unique sequences that were found in a specific condition in all replicates
  # this is also aa simple division of the peptide count per condition over the total peptide count in all conditions
  allcols<-colnames(pgroups)
  for(cond_i in conditions)
  {
    colsl<-allcols %in% conditions & ! allcols %in% cond_i
    pgroups[, paste0(cond_i,'p')]<-(pgroups[, cond_i]/rowSums(pgroups[, c(cond_i, allcols[colsl])]))*100
  }
  
  # Step 4: Add time-point column
  pgroups$exp_desc <- exp_desc
  
  # Step 5: in case there is replicate mismatch between conditions i.e. there is at least one replicate of the experiment that contains quantification values not from all conditions, there are some columns in pgroups at the moment
  # that have the format <replicate_description>.uniqueSequences.<condition> e.g. b1t1.uniqueSequences.WildType that contain 0 values solely. These columns should be emmited and the mismatch must be taken into account while sending the data to limma
  replicate_mismatch <<- F
  
  # We have replicate mismatch in case not all conditions correspond to the same experimental replicates
  for (rep_desc_i in unique(evidence.dt$rep_desc))
  {
    for (cond_i in conditions)
    {
      # For each of the columns of interest check if it contains solely 0 values and if so delete the respective intensity column
      if (all(pgroups[, paste0(rep_desc_i, ".uniqueSequences.", cond_i)] == 0, na.rm = TRUE) || all(is.nan(pgroups[, paste0(rep_desc_i, ".uniqueSequences.", cond_i)])))
      {
        allcols <- colnames(pgroups)
        pgroups <- pgroups[, - which(grepl(paste0("Intensity", ".", cond_i, ".", rep_desc_i), allcols))]
        replicate_mismatch <<- T
      }
    }
  }
  
  # Step 6: Remove unnecessary columns (uniqueSequences per rep_desc and percentage unique peptides per rep_desc)
  allcols<-colnames(pgroups)
  pgroups<-pgroups[,-which(grepl('uniqueSequences\\.', allcols) | grepl('p\\.b[0-9]+t[0-9]+$', allcols) | grepl('^common$', allcols))]
  
  levellog(paste0("read.pgroups: Quantifiable proteins: ", nrow(pgroups)," (",exp_desc,")"))
  levellog("",change=-1)
  
  # pgroups should now look like:
  #
  # +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+-------------------+----------------------+------------------+------------------+-------------------+
  # | Protein.IDs                                                                                                                                                                                                 | nas | Intensity.H.b1t1 | Intensity.L.b1t1 | N | b1t1.Ratio.counts | b1t1.uniqueSequences | Intensity.H.b1t2 | Intensity.L.b1t2 | b1t2.Ratio.counts |
  # +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+-------------------+----------------------+------------------+------------------+-------------------+
  # | AAR31361.1;ABV71636.1;P13202;ACZ80025.1|:   55 kDa immediate-early protein 1                                                                                                                                | 0   | 14520300         | 871800           | 1 | 1                 | 2                    | 10515500         | 300750           | 1                 |
  # +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+-------------------+----------------------+------------------+------------------+-------------------+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;AAR31362.1;ABV71635.1;P19893;Q6SWJ2|:   55 kDa immediate-early protein 1;Regulatory protein IE1;45 kDa   immediate-early protein 2;Protein UL122;Regulatory protein IE2 | 0   | 6315898          | 290628.7         | 2 | 2                 | 3                    | 7997380          | 63117.7          | 2                 |
  # +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+-------------------+----------------------+------------------+------------------+-------------------+
  # | AAR31361.1;ABV71636.1;P13202;Q6SWJ1;ACZ80025.1|:   55 kDa immediate-early protein 1;Regulatory protein IE1                                                                                                  | 0   | 42086141         | 804846.7         | 3 | 3                 | 5                    | 32540397         | 1955784          | 3                 |
  # +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----+------------------+------------------+---+-------------------+----------------------+------------------+------------------+-------------------+
  #
  # (...)
  #
  # +----------------------+------------------+------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+------------------+
  # | b1t2.uniqueSequences | Intensity.H.b2t1 | Intensity.L.b2t1 | b2t1.Ratio.counts | b2t1.uniqueSequences | Intensity.H.b2t2 | Intensity.L.b2t2 | b2t2.Ratio.counts | b2t2.uniqueSequences | Intensity.H.b3t1 |
  # +----------------------+------------------+------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+------------------+
  # | 2                    | 5686600          | 4330.7           | 1                 | 1                    | 4669800          | 4330.7           | 1                 | 1                    | NA               |
  # +----------------------+------------------+------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+------------------+
  # | 3                    | 2115170          | 4330.7           | 1                 | 1                    | 1054030          | 4330.7           | 1                 | 1                    | 8142900          |
  # +----------------------+------------------+------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+------------------+
  # | 5                    | 17497600         | 8661.4           | 2                 | 2                    | 10065000         | 8661.4           | 2                 | 2                    | 15538740         |
  # +----------------------+------------------+------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+------------------+
  #
  # (...)
  #
  # +------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+----+----+----------+----------+----------+
  # | Intensity.L.b3t1 | b3t1.Ratio.counts | b3t1.uniqueSequences | Intensity.H.b3t2 | Intensity.L.b3t2 | b3t2.Ratio.counts | b3t2.uniqueSequences | H  | L  | Hp       | Lp       | exp_desc |
  # +------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+----+----+----------+----------+----------+
  # | NA               | NA                | NA                   | NA               | NA               | NA                | NA                   | NA | NA | NA       | NA       | 0h       |
  # +------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+----+----+----------+----------+----------+
  # | 3170900          | 1                 | 2                    | NA               | NA               | NA                | NA                   | NA | NA | NA       | NA       | 0h       |
  # +------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+----+----+----------+----------+----------+
  # | 4850292          | 1                 | 2                    | 16142000         | 5161000          | 1                 | 2                    | 12 | 6  | 66.66667 | 33.33333 | 0h       |
  # +------------------+-------------------+----------------------+------------------+------------------+-------------------+----------------------+----+----+----------+----------+----------+
  #
  # Return the final data frame
  
  return(pgroups)  
}

prepare_working_pgroups <- function(working_pgroups){
  
  # Prepare protein intensity table for differential expression analysis (the format limma requires)  
  
  rownames(working_pgroups)<-working_pgroups[,paste(quantitated_items_lbl,".IDs",sep="")]
  inten_cols<-c()
  for(cond_i in conditions){
    inten_cols<-c(inten_cols,sort(colnames(working_pgroups)[grep(paste("Intensity.",cond_i,".b",sep=""),colnames(working_pgroups))]))
  }  
  working_pgroups<-working_pgroups[,inten_cols]
  colnames(working_pgroups)<-sub("Intensity\\.","",inten_cols)
  return(working_pgroups)
}

addLabel <- function(lblname, ...){
  
  # This routine adds labels or conditions (the two terms can be used interchangeably in PS) to the experiment
  levellog("", change=1)
  
  # If label name is a number some routines won't work, it has to be converted to some acceptable variable name
  lblname<-make.names(lblname)
  
  # labeltxt is used only for making the log seem nicer
  labeltxt <- "label"
  if(LabelFree)
  {
    labeltxt <- "condition";
  }
  
  # Make sure that there is not another condition with the same name
  lblname_i<-which(grepl(paste("^",lblname,"$",sep=""),conditions))
  if(length(lblname_i) != 0){
    levellog(paste("addLabel: Could not add ",labeltxt," '",lblname,"': An existing ",labeltxt," with name '",lblname,"' already exists. Please try a different name.",sep=""), change=-1)
    return(FALSE)
  }
  
  # Add the condition to the conditions global variable
  conditions<<-c(conditions, lblname)
  
  # Recalculate the conditions length
  nConditions<<-length(conditions)
  levellog("", change=-1)
}

clearConditions <- function(){
  # clearConditions will clear the conditions vector that stores all conditions in the current experiment
  # and will update the nConditions variable that stores how many conditions we have
  
  levellog("", change=1)
  conditions<<-c()
  nConditions<<-length(conditions)
  levellog("", change=-1)
}

my_grep <- function(x) {
  # Uniprot IDs pattern
  grep('^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',x, value=TRUE)
}

get_uniprot_ids <- function(results, cond1, cond2) {
  
  # Get the UNIPROT IDs from the protein names
  
  #get p values associated to DE proteins
  col_desc_<-paste("P-value adjusted ",paste(cond2,"/",cond1,sep=""),sep="")
  col_picker <- grepl(col_desc_, x = colnames(results), ignore.case = T)
  
  # Find DE proteins in results data frame
  ind_diffexp_tmp<-which(results[,col_picker]<pThreshold)
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

run_enrichment_analysis <- function(UniProtList, myFNorganism, cond1, cond2) {
  # Enrichment analysis utilizing gprofiler2 R package
  enrich <- gost(query= UniProtList, organism = myFNorganism, domain_scope = "annotated", significant = T, evcodes = TRUE, sources = c("GO", "KEGG", "REAC", "HPA", "WP"))
  enrich.matrix <- as.matrix(enrich$result[enrich$result$p_value < 0.05,c( "source", "term_name", "term_id", "p_value", "term_size", "query_size", "intersection_size", "intersection")])
  colnames(enrich.matrix) = c("Data Source", "Function", "Term ID", "p-Value", "Term size", "Query size", "Intersection Size", "Intersection")
  
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
  
  # Which describes the Function detected per line, the term ID etc and all information relative to a functional enrichment analysis
  
  write.table(enrich.matrix, paste(outputFigsPrefix,"_enrichment_results_" , cond2, ".", cond1 , ".txt",sep=""), row.names=FALSE, sep = "\t", dec = ".", quote = F)
}

generate_Venn_diagrams <- function(results_intensities, replicate_descs) {
  
  # Venn diagrams can be very useful for visualising the reproducibility of the experiment. In general two replicates
  # of the experiment should not have great differences in the amount of proteins quantified. The Venn diagram might
  # show these possible differences between biological replicates, or technical replicates of the same biorep
  # the following function gets the experimental structure and the columns of the protein intensities across all replicates
  
  # First of all suppress the log file from Venn Diagram package
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  # This function gets the results data frame (only the columns that show the intensities per replicate) and a
  # vector that matches an index to the respective replicate description like: "b1t1" "b1t2" "b2t1" "b2t2" "b3t1" "b3t2"
  # in venn production we only care if a protein was quantified or not so get the results_intensities and wherevere we see
  # an NA value set it to FALSE (not quantified) and set the rest to TRUE
  results_intensities <- !(is.na(results_intensities))
  
  nRep<-length(replicate_descs)
  nCond<-length(conditions)
  
  # The following formula produces a matrix that match a condition and a replicate to a results_intensities column
  mymap <- sapply(c(1:nRep), function(x) (c(0:(nCond-1)) * nRep) + x)
  rownames(mymap) = conditions
  colnames(mymap) = replicate_descs
  
  #example of mymap:
  
  #    b1t1 b1t2 b2t1 b2t2 b3t1 b3t2
  # H    1    2    3    4    5    6
  # L    7    8    9   10   11   12
  
  # since we care if a protein was overally quantified in a replicate but do not care in which condition it was quantified
  # merge the columns of same replicate and different condition
  
  Quant_bool_per_rep <- sapply(c(1:nRep), function(x) apply(results_intensities[,mymap[,x]], 1, any))
  colnames(Quant_bool_per_rep) <- replicate_descs
  
  # this data frame has now nRep columns showing if the proteins are present in each replicate
  # To understand the concept of the Venn diagrams lets assume we want to see the reproducibility between 2 bio reps that
  # have no technical replication. Then Quant_bool_per_rep will have two columns b1t1 and b2t1. To construct the venn diagram the only
  # needed information is how many proteins where quantified in each one of the bioreps and how many in both (only 3 numbers)
  # To generalize this we will search for all the columns in Quant_bool_per_rep that match a specific biorep, the following lines
  # do exactly that
  
  # first get all the biorep indices
  
  all_brep_indices <- unique(sapply(1:nRep, function(x) str_extract(str_extract(colnames(Quant_bool_per_rep)[x], "b\\d+?[t$]"), "\\d+")))
  
  # for example if the descriptions of the replicates are b1t1 b2t1 b4t1 all_brep_indices is "1" "2" "4"
  
  # Now get all the columns that refer to a specific biorep and compute which proteins where quantified in each column
  
  pdf(paste(outputFigsPrefix,"_Venn_for_all_bioreps_separately",exp_desc,".pdf",sep="") ,width=10, height=10, family = "Helvetica", pointsize=8)
  
  for(i in as.numeric(all_brep_indices))
  {
    # First get the columns in Quant_bool_per_rep that are respective to a specific biorep
    col_idxs_of_a_biorep <- grep(paste0("b", i), colnames(Quant_bool_per_rep))
    
    lst_prots_per_rep_quantified <- c()
    # The following gets a list of the names of all proteins quantified in each techrep of this current biorep
    lst_prots_per_rep_quantified <- lapply(col_idxs_of_a_biorep, function(x) rownames(Quant_bool_per_rep)[Quant_bool_per_rep[,x] == T])
    
    # for two treps of a specific biorep the result could be:
    # lst_prots_per_rep_quantified[[1]] (truncated for preview reasons)
    
    # [1] "H-INV:HIT000007749;Q58E" "H-INV:HIT000038585;P220" "H-INV:HIT000013642;Q96M" "H-INV:HIT000029498;Q135" "SWISS-PROT:P55265-4;E7E"
    # [6] "H-INV:HIT000037112;P467" "H-INV:HIT000195837;P010" "H-INV:HIT000196945;Q010" "REFSEQ:NP_002096;P16104" "H-INV:HIT000063977;Q135"
    # [11] "REFSEQ:NP_001008;P62277" "F5H2A4;H-INV:HIT0003245" "REFSEQ:NP_001518;J3KPW7" "H-INV:HIT000038368;P390" "TREMBL:Q59GA1;REFSEQ:NP"
    # [16] "H-INV:HIT000262212;P158" "H-INV:HIT000301385;P300" "H-INV:HIT000032127;Q132" "H0YK46;H0YN88;H-INV:HIT" "H-INV:HIT000032918;P381"
    
    # and for the lst_prots_per_rep_quantified[[2]] (second trep)
    
    # [1] "H-INV:HIT000007749;Q58E" "H-INV:HIT000038585;P220" "H-INV:HIT000013642;Q96M" "H-INV:HIT000029498;Q135" "SWISS-PROT:P55265-4;E7E"
    # [6] "H-INV:HIT000037112;P467" "H-INV:HIT000195837;P010" "REFSEQ:NP_055644;Q7L014" "REFSEQ:NP_002096;P16104" "H-INV:HIT000063977;Q135"
    # [11] "H-INV:HIT000288710;P073" "F5H2A4;H-INV:HIT0003245" "REFSEQ:NP_001518;J3KPW7" "H-INV:HIT000038368;P390" "TREMBL:Q59GA1;REFSEQ:NP"
    # [16] "H-INV:HIT000262212;P158" "H-INV:HIT000301385;P300" "H0YK46;H0YN88;H-INV:HIT" "H-INV:HIT000032918;P381" "REFSEQ:NP_001035807;NP_"
    
    # Notice that some proteins are quantified only in one trep but some are in both
    # this is the only list that VennDiagram package needs to create our Venn Diagram
    
    # In case the experiment has more than 5 bio reps the Venn can not be created since it will be no informative at all
    
    # Abort plot drawing in such a case
    
    if (length(col_idxs_of_a_biorep)>5)
    {
      levellog(paste0("Warn User: Biological Replicate ", i , " Venn diagram failed - too many technical replicates"))
      next
    }
    
    #Create the diagram
    
    result <- tryCatch({
      png(paste("../",outputFigsPrefix,"_Venn_for_bio_rep_", i, "_",exp_desc,".png",sep=""),  width = 1000, height = 1000)
      
      names(lst_prots_per_rep_quantified) = colnames(Quant_bool_per_rep)[col_idxs_of_a_biorep]
      VennPalette <- c("#00a8ff", "#9c88ff", "#fbc531", "#4cd137", "#487eb0")
      g <- grid.newpage()
      venn.plot <- venn.diagram(lst_prots_per_rep_quantified, NULL , fill=VennPalette[1:length(col_idxs_of_a_biorep)], lwd=1, col=VennPalette[1:length(col_idxs_of_a_biorep)], margin = 0.08, cex=2.5, cat.cex=2.5)
      g <- grid.draw(venn.plot)
      g <- grid.text(paste0("Venn diagram for Biological Replicate ", i), x = unit(0.5, "npc"), y = unit(0.97, "npc"), gp = gpar(cex=2.5), draw = TRUE, vp = NULL)
      
      dev.off()
      
      # The followig lines print to the pdf dev
      g <- grid.newpage()
      venn.plot <- venn.diagram(lst_prots_per_rep_quantified, NULL , fill=VennPalette[1:length(col_idxs_of_a_biorep)], lwd=1, col=VennPalette[1:length(col_idxs_of_a_biorep)], margin = 0.08, cex=2.5, cat.cex=2.5)
      g <- grid.draw(venn.plot)
      g <- grid.text(paste0("Venn diagram for Biological Replicate ", i), x = unit(0.5, "npc"), y = unit(0.97, "npc"), gp = gpar(cex=2.5), draw = TRUE, vp = NULL)
      
      
    }, error = function(err){
      levellog(paste0("Warn User: Biological Replicate ", i , " Venn diagram failed"))
    })
    
  }
  
  dev.off() # Close the pdf dev
  
  # Now create the same plot for reproducibility between bioreps
  
  # Merge the columns that refer to the same biorep
  # e.g. if we have the columns "b1t1" "b1t2" "b2t1" "b2t2" "b3t1" "b3t2" merge b1t1 with b1t2 to see if a protein was overally
  # quantified in b1 etc. The final data frame will have nBioRep columns
  
  # First create a matrix where columns correspond to bioreps and its contents are the column idxs in Quant_bool_per_rep that
  # correspond to this biorep
  
  col_idxs <- sapply(all_brep_indices, function(x) grep(paste0("b", x), colnames(Quant_bool_per_rep)))
  
  # If all breps have only 1 trep the col_idxs will be a single dimensional vector
  
  if (is.matrix(col_idxs))
  {
    colnames(col_idxs) <- sapply(all_brep_indices, function(x) paste0("b", x))
  }
  
  # e.g of col_idxs
  #     b1 b2 b3
  # [1,]  1  3  5
  # [2,]  2  4  6
  
  # Now feed all columns of a biorep one at a time to any function to merge the desired columns
  # if the col_idxs is a vector then there is no need to merge anything
  
  if (is.matrix(col_idxs))
  {
    Quant_bool_per_bio_rep <- apply(col_idxs, 2, function(x) apply(Quant_bool_per_rep[, x], 1, any))
  } else {
    Quant_bool_per_bio_rep <- Quant_bool_per_rep
  }
  
  # Quant_bool_per_bio_rep contains one column per biorep and TRUE/FALSE values showing if the protein ws overally
  # quantified in a biorep
  
  #e.g.:
  
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # |                                                                                                                                     | b1   | b2   | b3    |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | H-INV:HIT000007749;Q58EX7;REFSEQ:NP_001123203;  [Pleckstrin homology domain-containing family G mem ...]                            | TRUE | TRUE | FALSE |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | H-INV:HIT000038585;P22087;ENSEMBL:ENSP00000339522;ENSP00000397121;A6NHQ2   [34 kDa nucleolar scleroderma antigen;rRNA 2-O-meth ...] | TRUE | TRUE | TRUE  |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | H-INV:HIT000013642;Q96ME7;ENSEMBL:ENSP00000369040;   [Zinc finger protein 512;cDNA FLJ52441, highly simi ...]                       | TRUE | TRUE | FALSE |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | H-INV:HIT000029498;Q13547;F5GXM1;TREMBL:B4DSK9;   [Histone deacetylase 1;cDNA FLJ51764, highly simila ...]                          | TRUE | TRUE | FALSE |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | SWISS-PROT:P55265-4;E7ENU4;REFSEQ:NP_001102; [136 kDa double-stranded RNA-binding protein;Double ...]                               | TRUE | TRUE | TRUE  |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | H-INV:HIT000037112;P46782 [40S ribosomal   protein S5;40S ribosomal protein S5, ...]                                                | TRUE | TRUE | TRUE  |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # | H-INV:HIT000195837;P01031;TREMBL:Q59GS8   [C3 and PZP-like alpha-2-macroglobulin domain-conta ...]                                  | TRUE | TRUE | FALSE |
  # +-------------------------------------------------------------------------------------------------------------------------------------+------+------+-------+
  # 
  
  # Create a list of all protein names quantified in each biorep
  
  
  lst_prots_per_brep_quantified <- c()
  # The following gets a list of the names of all proteins quantified in each techrep of this current biorep
  lst_prots_per_brep_quantified <- lapply(1:ncol(Quant_bool_per_bio_rep), function(x) rownames(Quant_bool_per_bio_rep)[Quant_bool_per_bio_rep[,x] == T])
  
  
  
  # Now simply create a VennDiagram for these data
  
  if (n_bioreps>5)
  {
    levellog(paste0("Warn User: Overall biological Replicates Venn diagram failed, too many bioreps to plot"))
  } else {
    
    result <- tryCatch({
      png(paste("../",outputFigsPrefix,"_Venn_for_bio_reps_",exp_desc,".png",sep=""),  width = 1000, height = 1000)
      
      
      if (is.matrix(col_idxs))
      {
        names(lst_prots_per_brep_quantified) = colnames(col_idxs)
      } else {
        names(lst_prots_per_brep_quantified) = names(col_idxs)
      }
      
      VennPalette <- c("#00a8ff", "#9c88ff", "#fbc531", "#4cd137", "#487eb0")
      g <- grid.newpage()
      venn.plot <- venn.diagram(lst_prots_per_brep_quantified, NULL , fill=VennPalette[1:n_bioreps], lwd=1, col=VennPalette[1:n_bioreps], margin = 0.08, cex=2.5, cat.cex=2.5)
      g <- grid.draw(venn.plot)
      g <- grid.text("Venn diagram for all Biological Replicates", x = unit(0.5, "npc"), y = unit(0.97, "npc"), gp = gpar(cex=2.5), draw = TRUE, vp = NULL)
      
      dev.off()
      
      # The following lines print to the pdf dev
      
      pdf(paste(outputFigsPrefix,"_Venn_for_bio_reps_",exp_desc,".pdf",sep="") ,width=10, height=10, family = "Helvetica", pointsize=8)
      g <- grid.newpage()
      venn.plot <- venn.diagram(lst_prots_per_brep_quantified, NULL , fill=VennPalette[1:n_bioreps], lwd=1, col=VennPalette[1:n_bioreps], margin = 0.08, cex=2.5, cat.cex=2.5)
      g <- grid.draw(venn.plot)
      g <- grid.text("Venn diagram for all Biological Replicates", x = unit(0.5, "npc"), y = unit(0.97, "npc"), gp = gpar(cex=2.5), draw = TRUE, vp = NULL)
      dev.off()
      
      
    }, error = function(err){
      levellog(paste0("Warn User: Overall biological Replicates Venn diagram failed"))
    })
  }
  # The last Venn to produce would be a Venn between conditions. We will go back to results_intensities to merge all columns
  # refering to the same condition - the procedure is similarto what we did for the overall bre veVn
  col_idxs <- sapply(conditions, function(x) grep(paste0("^", x, " "), colnames(results_intensities)))
  
  
  # e.g of col_idxs
  #      H  L
  # [1,] 1  7
  # [2,] 2  8
  # [3,] 3  9
  # [4,] 4 10
  # [5,] 5 11
  # [6,] 6 12
  
  Quant_bool_per_condition <- apply(col_idxs, 2, function(x) apply(results_intensities[, x], 1, any))
  # Notice that this Venn Diagram may not have any meaning in precursor ion experiments or any experiment where proteins are labelled
  # thus if a protein is detected it is then quantified to all possible conditions. In this case the Venn diagram will show all
  # circles being totally overlapped
  
  lst_prots_per_cond_quantified <- c()
  
  lst_prots_per_cond_quantified <- lapply(1:ncol(Quant_bool_per_condition), function(x) rownames(Quant_bool_per_condition)[Quant_bool_per_condition[,x] == T])
  if (length(colnames(col_idxs))>5)
  {
    levellog(paste0("Warn User: Conditions Venn diagram failed, too many conditions to plot"))
  } else {
    
    result <- tryCatch({
      png(paste("../",outputFigsPrefix,"_Venn_for_conditions_",exp_desc,".png",sep=""),  width = 1000, height = 1000)
      
      names(lst_prots_per_cond_quantified) = colnames(col_idxs)
      VennPalette <- c("#00a8ff", "#9c88ff", "#fbc531", "#4cd137", "#487eb0")
      g <- grid.newpage()
      venn.plot <- venn.diagram(lst_prots_per_cond_quantified, NULL , fill=VennPalette[1:length(colnames(col_idxs))], lwd=1, col=VennPalette[1:length(colnames(col_idxs))], margin = 0.08, cex=2.5, cat.cex=2.5)
      g <- grid.draw(venn.plot)
      g <- grid.text("Venn diagram for all Conditions", x = unit(0.5, "npc"), y = unit(0.97, "npc"), gp = gpar(cex=2.5), draw = TRUE, vp = NULL)
      
      dev.off()
      
      # The following lines print to the pdf dev
      
      pdf(paste(outputFigsPrefix,"_Venn_for_conditions_",exp_desc,".pdf",sep="") ,width=10, height=10, family = "Helvetica", pointsize=8)
      g <- grid.newpage()
      venn.plot <- venn.diagram(lst_prots_per_cond_quantified, NULL , fill=VennPalette[1:length(colnames(col_idxs))], lwd=1, col=VennPalette[1:length(colnames(col_idxs))], margin = 0.08, cex=2.5, cat.cex=2.5)
      g <- grid.draw(venn.plot)
      g <- grid.text("Venn diagram for all Conditions", x = unit(0.5, "npc"), y = unit(0.97, "npc"), gp = gpar(cex=2.5), draw = TRUE, vp = NULL)
      dev.off()
      
      
    }, error = function(err){
      levellog(paste0("Warn User: Conditions Venn diagram failed"))
    })
  }
  
  
}

initProteoSign <- function() {
  
  # Set the duplicate correlation trim for limma
  duplicateCorrelation_trim<<-0.15 # use 0.22 for "bad" datasets (too many missing values)
  
  #Prepare the conditions vector
  clearConditions()
  
  #get the current working directory
  working_directory<<-getwd()
  
  # limma_out_dir is the directory that will be created to output limma results files and more. In order to work properly with the rest of ProteoSign it should be hardcoded to "msdiffexp_out"
  # but it can be changed to preference if ProteoSign is run in command line
  limma_out_dir<<-"msdiffexp_out"
  
  # An experiment that PS works on can be either an Isobaric Labelled one (e.g. TMT) a LabelFree one (known otherwise as LFQ)
  # or a metabolically labelled one. If both of the following boolean vars are false it is a metabolically labelled one otherwise
  # an isobaric or LFQ
  IsobaricLabel<<-F
  LabelFree<<-F
  
  # The following variables define if PS may use the features label rename and label swap (see ur documentation for more info)
  AllowLabelRename<<-F
  AllowLS<<-F
  
  
  if (!DEBUG) {
    cmd_args = commandArgs(T)
    if (length(cmd_args) < 1)
    {
      # if the script is not run through with arguments
      # MSdiffexp_definitions.R is a source R file that is used to load some experiment specific options
      source("MSdiffexp_definitions.R")
    } else {
      # If there is at least one argument we will load the preferences from the arguments
      # The arguments are passed in flags starting with a minus (all arguments here are optional)
      # and are usually followed by another argument that sets a value to the corresponding variable
      
      # The following table shows all valid arguments:
      #
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | ==== Argument ==== |                             ==== Description ====                             |         Default value         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -df <FILENAME>     | Definition filename - an R script to load before PS backend                   | MSdiffexp_definitions.R       |
      # |                    | runs to load any definitions                                                  |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -hc <COLOR>        | The histogram's box colour                                                    | cyan                          |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -lc <COLOR>        | The color of the best fit line in scatter plots                               | red                           |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -PD                | Set Proteome Discoverer as the preprocessing program                          | (PD is considered by default) |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -MQ                | Set MaxQuant as the preprocessing program                                     | (PD is considered by default) |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -rb <NUMBER>       | Set the least required bioreps for a protein to be valid                      | 2                             |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -rp <NUMBER>       | Set the least required peptides quantified for a protein to consider          | 2                             |
      # |                    | it valid                                                                      |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -p <NUMBER>        | Set the p threshold                                                           | 0.05                          |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -exf <FILENAME>    | The experimental structure file                                               | exp_struct.txt                |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -LFf <FILENAME>    | The Label to raw file matching file for LFQ experiments                       | LFQ_conditions.txt            |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RNf <FILENAME>    | The file containing all information to Rename - merge conditions              | Rename_array.txt              |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -LSf <FILENAME>    | The file containing all necessarry information for label swap                 | LS_array.txt                  |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -mf <FILENAME>     | The main file to retrieve data from (MultiConsensus for PD                    | msdiffexp_protein.txt         |
      # |                    | and evidence for MQ)                                                          |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -pgf <FILENAME>    | The protein groups file used for MQ data processing                           | msdiffexp_peptide.txt         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RMrf <FILENAME>   | The file containing all Raw file data for replication multiplexing            | RMrawfiles.txt                |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RMtf <FILENAME>   | The file containing all Tags data for replication multiplexing                | RMtags.txt                    |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -op <PREFIX>       | A prefix to be added to all output files                                      | exp_1                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -fl                | Set background filtering to true                                              | False                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -flL <LABEL>       | The label to be considered as background label                                | Light                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -flp               | Set the filtering level to peptide                                            | False                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -et <TYPE>         | Set the experiment type: ML for metabolic labeling (SILAC etc.) IL for        | ML                            |
      # |                    |  Isobaric Labeling (ITRAQ etc.) and LF for label free                         |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -ed <DESCRIPTION>  | Set the experiment description                                                | df_1                          |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -LR                | Enables label renaming                                                        | False                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -LS                | Enables Label Swapping                                                        | False                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RM                | Enables Replication Multiplexing                                              | False                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RMbir             | For Replication Multiplexing experiments: defines that biological reps        | False                         |
      # |                    | are represented as different raw files (False for being represented as        |                               |
      # |                    | different tags)                                                               |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RMtir             | For Replication Multiplexing experiments: defines that technical reps         | False                         |
      # |                    |  are represented as different raw files (False for being represented as       |                               |
      # |                    |  different tags)                                                              |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -RMcir             | For Replication Multiplexing experiments: defines that conditions are         | False                         |
      # |                    |  represented as different raw files (False for being represented as different |                               |
      # |                    |  tags)                                                                        |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -FE                | Enables funtional enrichment analysis                                         | False                         |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -Fo <ORGANISM>     | Sets the target organism for functional enrichment. Notice that the           | human                         |
      # |                    | only valid organism values are the ones in the second column of the tab       |                               |
      # |                    | separated file 'valid_GO_Organisms.txt' under 'cgi-bin' folder in             |                               |
      # |                    | ProteoSign's distribution                                                     |                               |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      # | -of <FOLDER>       | The output folder's path - name                                               | msdiffexp_out                 |
      # +--------------------+-------------------------------------------------------------------------------+-------------------------------+
      
      # First load all default values:
      levellog("PS is executed with command line arguments...")
      levellog("Setting default arguments")
      
      ratios.hist.colour<<-"cyan"
      reps.scatter.lmline.colour<<-"red"
      PDdata<<-T
      nRequiredLeastBioreps<<-2
      nRequiredLeastPeps<<-2
      pThreshold<<-0.05
      exportFormat<<-"pdf"
      experimental_structure_file<<-"exp_struct.txt"
      LFQ_conditions_file<<-"LFQ_conditions.txt"
      Rename_Array_file<<-"Rename_array.txt"
      LS_Array_file<<-"LS_array.txt"
      pgroups_fname<<-"msdiffexp_protein.txt"
      evidence_fname<<-"msdiffexp_peptide.txt"
      RMrawfilesdata_file<<-"RMrawfiles.txt"
      RMtagsdata_file<<-"RMtags.txt"
      outputFigsPrefix<<-"exp_1"
      filterL<<-F
      filterL_lbl<<-"Light"
      filterL_lvl<<-F
      LabelFree<<-F
      IsobaricLabel<<-F
      All_MQ_Labels<<-c()
      exp_desc<<-"df_1"
      ProteinQuantitation<<-T
      AllowLabelRename<<-F
      AllowLS<<-F
      RMisused<<-F
      RMbrepsinrawfiles<<-F
      RMtrepsinrawfiles<<-F
      RMconditionsinrawfiles<<-F
      FNenrichment<<-F
      FNorganism<<-'human'
      
      
      levellog("Reading arguments", change = 1)
      
      # ProteoSign will always search for the -df arument that defines a new definitions file to be loaded first
      # if it finds one it will try to execute. If it does not it will try to execute the default MSdiffexp_definitions.R
      definitions_file <<- "MSdiffexp_definitions.R"
      
      for (i in 1:length(cmd_args))
      {
        if (cmd_args[i] == "-df")
        {
          definitions_file <<- cmd_args[i + 1]
          levellog(paste0("Changed definitions file to ", definitions_file))
        }
      }
      
      
      if (file.exists(definitions_file))
      {
        source(definitions_file)
        levellog("Loaded definitions file")
      }
      else
      {
        levellog(paste0("Definitions file ", definitions_file, " was not found"))
      }
      
      cur_arg <- ""
      for (i in 1:length(cmd_args))
      {
        if (cur_arg == "")
        {
          # Perform operations for '-' starting arguments
          if (cmd_args[i] == '-hc'){
            cur_arg <- 'hc'
            next
          } else if (cmd_args[i] == '-lc'){
            cur_arg <- 'lc'
            next
          } else if (cmd_args[i] == '-PD'){
            PDdata <<- T
            levellog("Set PD as the preprocessing program")
            next
          } else if (cmd_args[i] == '-MQ'){
            PDdata <<- F
            levellog("Set MQ as the preprocessing program")
            next
          } else if (cmd_args[i] == '-rb'){
            cur_arg <- 'rb'
            next
          } else if (cmd_args[i] == '-rp'){
            cur_arg <- 'rp'
            next
          } else if (cmd_args[i] == '-p'){
            cur_arg <- 'p'
            next
          } else if (cmd_args[i] == '-exf'){
            cur_arg <- 'exf'
            next
          } else if (cmd_args[i] == '-LFf'){
            cur_arg <- 'LFf'
            next
          } else if (cmd_args[i] == '-RNf'){
            cur_arg <- 'RNf'
            next
          } else if (cmd_args[i] == '-LSf'){
            cur_arg <- 'LSf'
            next
          } else if (cmd_args[i] == '-mf'){
            cur_arg <- 'mf'
            next
          } else if (cmd_args[i] == '-pgf'){
            cur_arg <- 'pgf'
            next
          } else if (cmd_args[i] == '-RMrf'){
            cur_arg <- 'RMrf'
            next
          } else if (cmd_args[i] == '-RMtf'){
            cur_arg <- 'RMtf'
            next
          } else if (cmd_args[i] == '-op'){
            cur_arg <- 'op'
            next
          } else if (cmd_args[i] == '-fl'){
            filterL <<- T
            levellog("Set background filtering to True")
            next
          } else if (cmd_args[i] == '-flL'){
            cur_arg <- 'flL'
            next
          } else if (cmd_args[i] == '-flp'){
            filterL_lvl <<- T
            levellog("Set filtering level to peptide")
            next
          } else if (cmd_args[i] == '-et'){
            cur_arg <- 'et'
            next
          } else if (cmd_args[i] == '-ed'){
            cur_arg <- 'ed'
            next
          } else if (cmd_args[i] == '-LR'){
            AllowLabelRename <<- T
            levellog("Enabled Label Renaming")
            next
          } else if (cmd_args[i] == '-LS'){
            AllowLS <<- T
            levellog("Enabled Label Swapping")
            next
          } else if (cmd_args[i] == '-RM'){
            RMisused <<- T
            levellog("Enabled Replication Multiplexing")
            next
          } else if (cmd_args[i] == '-RMbir'){
            RMbrepsinrawfiles <<- T
            levellog("(Only for RM experiments) set biological replicate representation to 'different raw files'")
            next
          } else if (cmd_args[i] == '-RMtir'){
            RMtrepsinrawfiles <<- T
            levellog("(Only for RM experiments) set technical replicate representation to 'different raw files'")
            next
          } else if (cmd_args[i] == '-RMcir'){
            RMconditionsinrawfiles <<- T
            levellog("(Only for RM experiments) set conditions' representation to 'different raw files'")
            next
          } else if (cmd_args[i] == '-FE'){
            FNenrichment <<- T
            levellog("Enabled functional enrichment")
            next
          } else if (cmd_args[i] == '-Fo'){
            cur_arg <- 'Fo'
            next
          } else if (cmd_args[i] == '-of'){
            cur_arg <- 'of'
            next
          }
        } else {
          # Perform operations for arguments not starting with '-' that are dependent to the previous argument
          if (cur_arg == 'hc'){
            cur_arg <- ''
            ratios.hist.colour <<- cmd_args[i]
            levellog(paste0("Changed histogram box color to ", ratios.hist.colour))
            next
          } else if (cur_arg == 'lc'){
            cur_arg <- ''
            reps.scatter.lmline.colour <<- cmd_args[i]
            levellog(paste0("Changed best fit line color to ", reps.scatter.lmline.colour))
            next
          } else if (cur_arg == 'rb'){
            cur_arg <- ''
            nRequiredLeastBioreps <<- as.numeric(cmd_args[i])
            levellog(paste0("Least required bioreps to consider a protein valid set to ", nRequiredLeastBioreps))
            next
          } else if (cur_arg == 'rp'){
            cur_arg <- ''
            nRequiredLeastPeps <<- as.numeric(cmd_args[i])
            levellog(paste0("Least required peptides to consider a protein valid set to ", nRequiredLeastPeps))
            next
          } else if (cur_arg == 'p'){
            cur_arg <- ''
            pThreshold <<- as.numeric(cmd_args[i])
            levellog(paste0("Changed differential analysis p threshold to ", pThreshold))
            next
          } else if (cur_arg == 'exf'){
            cur_arg <- ''
            experimental_structure_file <<- cmd_args[i]
            levellog(paste0("Set experimental structure file path to ", experimental_structure_file))
            next
          } else if (cur_arg == 'LFf'){
            cur_arg <- ''
            LFQ_conditions_file <<- cmd_args[i]
            levellog(paste0("Set LFQ file path to ", LFQ_conditions_file))
            next
          } else if (cur_arg == 'RNf'){
            cur_arg <- ''
            Rename_Array_file <<- cmd_args[i]
            levellog(paste0("Set Rename file path to ", Rename_Array_file))
            next
          } else if (cur_arg == 'LSf'){
            cur_arg <- ''
            LS_Array_file <<- cmd_args[i]
            levellog(paste0("Set Label swap file path to ", LS_Array_file))
            next
          } else if (cur_arg == 'mf'){
            cur_arg <- ''
            evidence_fname <<- cmd_args[i]
            levellog(paste0("Set main data file path to ", evidence_fname))
            next
          } else if (cur_arg == 'pgf'){
            cur_arg <- ''
            pgroups_fname <<- cmd_args[i]
            levellog(paste0("Set protein groups data file for MQ path to ", pgroups_fname))
            next
          } else if (cur_arg == 'RMrf'){
            cur_arg <- ''
            RMrawfilesdata_file <<- cmd_args[i]
            levellog(paste0("Set Replication Multiplexing raw files data file to ", RMrawfilesdata_file))
            next
          } else if (cur_arg == 'RMtf'){
            cur_arg <- ''
            RMtagsdata_file <<- cmd_args[i]
            levellog(paste0("Set Replication Multiplexing tags data file to ", RMtagsdata_file))
            next
          } else if (cur_arg == 'op'){
            cur_arg <- ''
            outputFigsPrefix <<- cmd_args[i]
            levellog(paste0("Set output prefix to ", outputFigsPrefix))
            next
          } else if (cur_arg == 'flL'){
            cur_arg <- ''
            filterL_lbl <<- cmd_args[i]
            levellog(paste0("Set filter label to ", filterL_lbl))
            next
          } else if (cur_arg == 'et'){
            cur_arg <- ''
            temp_et <- cmd_args[i]
            if (temp_et == 'ML')
            {
              LabelFree <<- F
              IsobaricLabel <<- F
              levellog("Set experiment type to metabolic labeling")
            } else if (temp_et == 'IL') {
              LabelFree <<- F
              IsobaricLabel <<- T
              levellog("Set experiment type to isobaric labeling")
            } else if (temp_et == 'LF') {
              LabelFree <<- T
              IsobaricLabel <<- F
              levellog("Set experiment type to label free")
            }
            next
          } else if (cur_arg == 'ed'){
            cur_arg <- ''
            exp_desc <<- cmd_args[i]
            levellog(paste0("Set experiment description to ", exp_desc))
            next
          } else if (cur_arg == 'Fo'){
            cur_arg <- ''
            FNorganism <<- cmd_args[i]
            levellog(paste0("Set Functional enrichment organism to ", FNorganism))
            next
          } else if (cur_arg == 'of'){
            cur_arg <- ''
            limma_out_dir <<- cmd_args[i]
            levellog(paste0("Set output folder to ", limma_out_dir))
            next
          }
        }
      }
      levellog("Loaded command line arguments", change = -1)
    }
  }
}

load_table_param_files <- function() {
  # Much information for the experiment is given using some text files
  # that are tab separated parameter tables. This function loads these files
  # Reading the experimental structure is the first step:
  rep_structure<<-read.table(experimental_structure_file,col.names=c('raw_file','biorep','techrep','fraction'), sep="\t")
  
  # A truncated example of rep_structure could be:
  #
  # +--------------------------------------+--------+---------+----------+
  # | raw_file                             | biorep | techrep | fraction |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF01_01 | 1      | 1       | 1        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF01_02 | 1      | 2       | 1        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF02_01 | 1      | 1       | 2        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF02_02 | 1      | 2       | 2        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF03_01 | 1      | 1       | 3        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF03_02 | 1      | 2       | 3        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF04_01 | 1      | 1       | 4        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF04_02 | 1      | 2       | 4        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF05_01 | 1      | 1       | 5        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF05_02 | 1      | 2       | 5        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF06_02 | 1      | 1       | 6        |
  # +--------------------------------------+--------+---------+----------+
  # | OT2_Terhune_2012-10-09_DMC-HLNF06_03 | 1      | 2       | 6        |
  # +--------------------------------------+--------+---------+----------+
  # ...
  #
  # The table above matches each raw file to a biorep a techrep and a fraction
  # (rep_structure is a global variable)
  
  # Lets sort the table:
  rep_structure<<-rep_structure[order(rep_structure[,2],rep_structure[,3],rep_structure[,4]),]
  
  # Especially for LFQ data each raw file is matched to a condition as well. This kind of information
  # is stored in the LFQ_conditions_file
  
  # LFQ_conds will hold this kind of information:
  LFQ_conds<<-c()
  if(LabelFree)
  {
    # Load the necessary information
    LFQ_conds<<-read.table(LFQ_conditions_file, col.names=c('raw_file', 'condition'), stringsAsFactors = F)
    
    # An example of LFQ_conds could be:
    #
    #   +--------------+-----------+
    #   | raw_file     | condition |
    #   +--------------+-----------+
    #   | 20130426EL04 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL05 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL06 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL07 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL08 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL09 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL10 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL11 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL12 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL13 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL14 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL15 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL16 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL17 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL18 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL19 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL20 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL21 | WT        |
    #   +--------------+-----------+
    #   | 20130426EL24 | pfg377KO  |
    #   +--------------+-----------+
    #   | 20130426EL25 | pfg377KO  |
    #   +--------------+-----------+
    #   | 20130426EL26 | pfg377KO  |
    #   +--------------+-----------+
    #   | 20130426EL27 | pfg377KO  |
    #   +--------------+-----------+
    #
    # Matching each raw file to a condition.
    
  }
  
  
  if (AllowLabelRename == T)
  {
    # Label Renaming is a procedure where the user can rename some of their conditions. If two or more conditions
    # point to the same new name they will be merged.
    # Rename_array is a global variable. In cse we have merging for example of conditions M and L to "wt" the
    # array will seem like:
    
    # +-----------+-----------+
    # | old_label | new_label |
    # +-----------+-----------+
    # | H         | H         |
    # +-----------+-----------+
    # | M         | wt        |
    # +-----------+-----------+
    # | L         | wt        |
    # +-----------+-----------+
    
    # Notice that all conditions that are not renamed will have the same conditions in both columns.
    # In case a ProteoSign run does not need any renaming ther are 2 options. Either AllowLabelRename will be disabled
    # or the Rename_array will have the same labels on both columns.
    # The latter is used when running PS from its web interface, so in our example Rename_array will have the following format:
    
    # +-----------+-----------+
    # | old_label | new_label |
    # +-----------+-----------+
    # | H         | H         |
    # +-----------+-----------+
    # | L         | L         |
    # +-----------+-----------+
    
    
    Rename_Array <<- read.table(Rename_Array_file, col.names=c('old_label', 'new_label'), stringsAsFactors = F)
  }
  
  
  if (AllowLS == T)
  {
    # LS_array keeps all necessary information for label swapping (see our documentation for label swap)
    # Proteosign swaps quantitated data for the swapped labels in all raw files designated by LS array
    
    # An example of LS_array is seen below:
    #
    # +--------------------------------------+-------------+--------------+
    # | selected_raw_file                    | first_label | second_label |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF01_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF02_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF03_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF04_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF05_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF06_03 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF07_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF08_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF09_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF10_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF11_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    # | OT2_Terhune_2012-10-09_DMC-HLNF12_02 | L           | H            |
    # +--------------------------------------+-------------+--------------+
    #
    # That contains exactly the necessary information with all raw files where label swap between H and L will happen
    #
    # To make things easy, when PS is run using the web interface the LS array file is filled with a predtermined
    # table:
    # +-------------------+----------+-------+
    # | This run contains | no label | swaps |
    # +-------------------+----------+-------+
    # in this case LS_array is populated with the table above but since (hopefully) no raw files are ever
    # named "This run contains" PS runs properly with this LS_array but simply overlooks it
    
    
    Ls_array <<- read.table(LS_Array_file, col.names=c('selected_raw_file', 'first_label', 'second_label'), stringsAsFactors = F)
  }
  
  # Replication Multiplexing is quite complicated - it lets the user run PS on experiments
  # where bio and techreps might be represented as different tags and conditions as different raw files (MS runs)
  # PS's frontend sets RMisused to true to enable Replication Multiplexing (RM) if RM is used in the experiment
  # all necessary information for RM is contained in 2 tables that are RMrawfilesdata and RMtagsdata and show what is represented
  # by different raw files and tags respectively and 3 boolean variables:
  #
  # RMbrepsinrawfiles that is true if bioreps are represented by different raw files (as in a typical proteomics experiment)
  # RMtrepsinrawfiles that is true if treps are represented by different raw files (as in a typical proteomics experiment)
  # and RMconditionsinrawfiles that is true if conditions are represented as different raw files (that is an atypical way to represent your conditions, conditions are usually different labels or tags)
  #
  # almost all possible combinations of the booleans above are acceptable (if something is not reasonable it will be rejected from PS's frontend)
  
  
  # The following example shows an experiment where the conditions were distributed across MS runs but the
  # bio and techreps were represented as different tags in a TMT 6-ple example. Lets say that 3 bio reps of mice were chosen
  # and we have two samples (2 tech reps) for each one of these:
  #
  #                       +-------------------+
  #                       |                   |
  #                       |       Mice        |
  #                       |                   |
  #                       +---------+---------+
  #                                 |
  #                                 |
  #                  +--------------+----------------+
  #                  |              |                |
  #                  |              |                |
  #           +------v------+  +----v--------+ +-----v------+
  #           |             |  |             | |            |
  #           |   brep1     |  |   brep2     | |   brep3    |
  #           +-------+-----+  +---+---------+ +--------+---+
  #                   |            |                    |
  #      +------------+            +---------+          +-------------+
  #      |            |            |         |          |             |
  # +----v----+  +----v----  +-----v---+ +---v-----+ +--v------+ +----v----+
  # |         |  |         | |         | |         | |         | |         |
  # |  b1t1   |  |  b1t2   | |  b2t1   | |  b2t2   | |  b3t1   | |  b3t2   |
  # +---------+  +---------+ +---------+ +---------+ +---------+ +---------+
  # 
  # 
  # 
  # We tag each one of the samples with a different TMT tag (0 to 5) and we thus represent the technical replication in different tags
  #
  # in this example the aforementioned booleans will be 
  #
  # RMbrepsinrawfiles       F
  # RMtrepsinrawfiles       F
  # RMconditionsinrawfiles  T
  #
  # and RMtagsdata will be:
  #
  # +----+------+------+------+------+------+------+----------+
  # | id | name | brep | trep | frac | cond | used | selected |
  # +----+------+------+------+------+------+------+----------+
  # | 1  | 0    | 1    | 1    | -    | -    | TRUE | FALSE    |
  # +----+------+------+------+------+------+------+----------+
  # | 2  | 1    | 1    | 2    | -    | -    | TRUE | FALSE    |
  # +----+------+------+------+------+------+------+----------+
  # | 3  | 2    | 2    | 1    | -    | -    | TRUE | FALSE    |
  # +----+------+------+------+------+------+------+----------+
  # | 4  | 3    | 2    | 2    | -    | -    | TRUE | FALSE    |
  # +----+------+------+------+------+------+------+----------+
  # | 5  | 4    | 3    | 1    | -    | -    | TRUE | FALSE    |
  # +----+------+------+------+------+------+------+----------+
  # | 6  | 5    | 3    | 2    | -    | -    | TRUE | FALSE    |
  # +----+------+------+------+------+------+------+----------+
  #
  # Notice that all TMT tags (0-5) are matched against the respective brep-trep pair. Tags are not fractionated
  # so the frac column is empty and conditions are represented by MS runs so this column is also empty
  # The user can alway choose to opt out some tags so if the used field is set to FALSE the tag will not be used in the
  # analysis. The selected field is deprecated and always set to false
  #
  # Now lets assume that we have 4 conditions that we want our samples to be studied on. We can mix all the
  # samples above since they are tagged and divide the mixture into 4 parts and set them under different conditions
  # Afterwards, we can pass the parts through a mass spectrometer, thus the conditions will
  # be represented by different MS runs. Don't forget that many MS runs might correspond to a single condition
  # since MS runs are lmost always fractionated
  # Thus, a possible RMrawfilesdata could be:
  # 
  # +----+----------+------+------+------+-------+------+----------+
  # | id | name     | brep | trep | frac | cond  | used | selected |
  # +----+----------+------+------+------+-------+------+----------+
  # | 1  | M456-B06 | -    | -    | 1    | cond2 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 2  | M456-C03 | -    | -    | 1    | cond3 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 3  | M456-C04 | -    | -    | 2    | cond3 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 4  | M456-C09 | -    | -    | 3    | cond3 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 5  | M456-C10 | -    | -    | 4    | cond3 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 6  | M456-D01 | -    | -    | 1    | cond4 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 7  | M456-D05 | -    | -    | 2    | cond4 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 8  | M456-D10 | -    | -    | 3    | cond4 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 9  | M456-E04 | -    | -    | 1    | cond1 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  # | 10 | M456-E08 | -    | -    | 2    | cond1 | TRUE | FALSE    |
  # +----+----------+------+------+------+-------+------+----------+
  #
  # Each line corresponds to a different MS run - a different raw file - and matched against a condition - fraction pair
  # As before, all unnecessary columns are left blank
  # 
  # Using the 3 aformentioned boolean variables and these two tables all reasonable combinations of representation
  # of breps, treps and conditions by tags or raw files are possible
  
  
  if (RMisused == T)
  {
    # Load the raw files data for RM:
    RMrawfilesdata <<- read.table(RMrawfilesdata_file, col.names=c('id', 'name', 'brep', 'trep', 'frac', 'cond', 'used', 'selected'), stringsAsFactors = F)
  }
  if (RMisused == T)
  {
    # Load the tafs data for RM
    RMtagsdata <<- read.table(RMtagsdata_file, col.names=c('id', 'name', 'brep', 'trep', 'frac', 'cond', 'used', 'selected'), stringsAsFactors = F)
  }
}

edit_rep_structure <- function() {
  # This function edits the rep structure variable that holds the experimental structure
  # and writes the original rep structure variable as seen below:
  
  
  # Since sometimes some breps or treps might be discarded during experimental processes we might have
  # vaid experimental strucctures like the one below:
  #
  #   +--------------------------------------+------+------+------+
  #   | raw_file                             | brep | trep | frac |
  #   +--------------------------------------+------+------+------+
  #   | OT2_Terhune_2012-10-09_DMC-HLNF01_01 | 1    | 1    | 1    |
  #   +--------------------------------------+------+------+------+
  #   | OT2_Terhune_2012-10-09_DMC-HLNF01_02 | 1    | 3    | 1    |
  #   +--------------------------------------+------+------+------+
  #   | OT2_Terhune_2012-10-09_DMC-HLNF02_01 | 1    | 1    | 2    |
  #   +--------------------------------------+------+------+------+
  #   | OT2_Terhune_2012-10-09_DMC-HLNF02_02 | 1    | 4    | 1    |
  #   +--------------------------------------+------+------+------+
  #
  # Notice in the example above that brep 1 has 3 treps named 1-3-4 - three non-consecutive numbers.
  # This will break PS so we will keep a copy of the original structure to a global variable
  # and we will refresh rep_structure so that it contains only consecutive biorep numbers and
  # consecutive trep numbers for each brep
  
  # We will keep a copy of the original rep_structure to display in the graphs
  original_rep_structure <- rep_structure
  
  # We are not sure if the biorep and techrep numbers the user typed are sequential
  
  # Get all breps without duplicates
  unique_reps <- unique(rep_structure$biorep)
  
  counter <- 1
  for(rep_i in unique_reps){
    
    # Get all records in rep structure that contain the current brep
    mi <- which(rep_structure$biorep == rep_i)
    
    # Renumber the brep to a sequential counter
    rep_structure$biorep[mi] <- counter
    counter <- counter + 1
    
    # Get all the treps for the current brep and set them to sequential numbers as we do for all the breps
    unique_techreps <- unique(rep_structure$techrep[mi])
    counter2 <- 1
    for(techrep_i in unique_techreps){
      mi2 <- which(rep_structure$biorep == counter - 1 & rep_structure$techrep == techrep_i)
      rep_structure$techrep[mi2] <- counter2
      counter2 <- counter2 + 1
    }
  }
  
  # The following line aborts PS run in case we have only one biological replicate.
  # This is done because limma needs at least 2 bioreps to perform the differential analysis
  
  if (!RMisused)
  {
    if(length(unique(rep_structure$biorep)) == 1){
      levellog("Error User: Cannot accept dataset with just one biological replicate. Aborting ...")
      return(F)
    }
  }
  
  # The following lines sets a description for each biorep - trep - frac group.
  # Example of original_rep_structure
  #
  # +---+--------------------------------------+--------+---------+----------+----------+
  # |   | raw_file                             | biorep | techrep | fraction | rep_desc |
  # +---+--------------------------------------+--------+---------+----------+----------+
  # | 1 | OT2_Terhune_2012-10-09_DMC-HLNF01_01 | 1      | 1       | 1        | b1t1f1   |
  # +---+--------------------------------------+--------+---------+----------+----------+
  # | 2 | OT2_Terhune_2012-10-09_DMC-HLNF01_02 | 1      | 2       | 1        | b1t2f1   |
  # +---+--------------------------------------+--------+---------+----------+----------+
  # | 3 | OT2_Terhune_2012-10-09_DMC-HLNF02_01 | 1      | 1       | 2        | b1t1f2   |
  # +---+--------------------------------------+--------+---------+----------+----------+
  # | 4 | OT2_Terhune_2012-10-09_DMC-HLNF02_02 | 1      | 2       | 2        | b1t2f2   |
  # +---+--------------------------------------+--------+---------+----------+----------+
  # | 5 | OT2_Terhune_2012-10-09_DMC-HLNF03_01 | 1      | 1       | 3        | b1t1f3   |
  # +---+--------------------------------------+--------+---------+----------+----------+
  # (...)
  #
  #
  # The description depends in whether we have fractionation or not (in the case of fractionation the f<#> suffix is added)
  
  
  
  if(length(unique(rep_structure$fraction)) > 1){
    # we have fractionation
    rep_structure$rep_desc<-paste(paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep=''),'f',rep_structure$fraction,sep='')
    original_rep_structure$rep_desc<-paste(paste(paste('b',original_rep_structure$biorep,sep=''),'t',original_rep_structure$techrep,sep=''),'f',original_rep_structure$fraction,sep='')
  }else{
    # we have no fractionation
    rep_structure$rep_desc<-paste(paste('b',rep_structure$biorep,sep=''),'t',rep_structure$techrep,sep='')
    original_rep_structure$rep_desc<-paste(paste('b',original_rep_structure$biorep,sep=''),'t',original_rep_structure$techrep,sep='')
  }
  
  # Make rep structure and original rep structure globa variables
  .GlobalEnv[["rep_structure"]]<-rep_structure
  .GlobalEnv[["original_rep_structure"]]<-original_rep_structure
  
}

remove_double_quotes <- function() {
  
  # Removes the double quotes from all quantitation input datafiles
  
  # Remove double quotes (if any) from the evidence file
  if(grepl("\"",readLines(evidence_fname, n=1))){
    levellog("Removing double quotes from input data file #1 ...")
    tmpdata<-gsub("\"", "", readLines(evidence_fname))
    evidence_fname_cleaned<-file(evidence_fname, open="w")
    writeLines(tmpdata, con=evidence_fname_cleaned)
    close(evidence_fname_cleaned)
  }
  
  # If we are in an MQ preprocessed experiment remove the double quotes (if any) from the protein groups file as well
  if (!PDdata)
  {
    if(grepl("\"",readLines(pgroups_fname, n=1))){
      levellog("Removing double quotes from input data file #2 ...")
      tmpdata<-gsub("\"", "", readLines(pgroups_fname))
      pgroups_fname_cleaned<-file(pgroups_fname, open="w")
      writeLines(tmpdata, con=pgroups_fname_cleaned)
      close(pgroups_fname_cleaned)
    }
  }
}

write_proteinGroupsDF <- function() {
  # This writes the proteinGroupsDF file that contains information for each protein necessary to proceed to
  # the differential analysis
  
  #Restore the original rep descriptions to add to proteinGroupsDF output file
  if (!RMisused)
  {
    newcolumns <- names(protein_groups)
    oldcolumns <- newcolumns
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
  setwd(limma_out_dir)
  write.table(protein_groups[, -which(names(protein_groups) %in% c("N.x","N.y"))],file=paste(outputFigsPrefix,"_proteinGroupsDF.txt",sep=""),row.names=F,sep="\t")
  setwd(working_directory)
  
  if (!RMisused)
  {
    colnames(protein_groups) <- oldcolumns
  }
}

create_expdesign <- function() {
  # Creates the expdesign table that describes the experimental design in a limma friendly manner
  # It also saves it in a file:
  
  # Create the expdesign table:
  expdesign<-c()
  
  for(cond_i in conditions){
    expdesign<-rbind(expdesign,cbind(paste(sub("Intensity\\.","",sort(colnames(protein_groups)[grep(paste("Intensity.",cond_i,".b",sep=""),colnames(protein_groups))]))),cond_i))  
  }
  
  colnames(expdesign)<-c("Sample","Category")
  
  # The first step is to assign a label to each one of the <condition>.<rep_desc> values. An example of expdesign now
  # can be:
  
  # +--------+--------+
  # | Sample |Category|
  # +--------+--------+
  # | H.b1t1 | H      |
  # +--------+--------+
  # | H.b1t2 | H      |
  # +--------+--------+
  # | H.b2t1 | H      |
  # +--------+--------+
  # | H.b2t2 | H      |
  # +--------+--------+
  # | H.b3t1 | H      |
  # +--------+--------+
  # | H.b3t2 | H      |
  # +--------+--------+
  # | L.b1t1 | L      |
  # +--------+--------+
  # | L.b1t2 | L      |
  # +--------+--------+
  # | L.b2t1 | L      |
  # +--------+--------+
  # | L.b2t2 | L      |
  # +--------+--------+
  # | L.b3t1 | L      |
  # +--------+--------+
  # | L.b3t2 | L      |
  # +--------+--------+
  #
  

  
  if(!RMisused){
    # The following lines deal with restoring the original breps and treps numbers from original_rep_structure
    # temp_vector created below will get only the rep description from the first column of exp design
    temp_vector <- sub("(.*)\\.","", expdesign[,1])
    
    # e.g. : temp_vector: "b1t1" "b1t2" "b2t1" "b2t2" "b3t1" "b3t2" "b1t1" "b1t2" "b2t1" "b2t2" "b3t1" "b3t2"
    
    # The following line executes match to find the indices in rep_structure$rep_desc that
    # contain rep description information for the repdescs in temp_vector
    # for example b1t1 in rep_structure$rep_desc is located in the first element (in fact the first element is b1t1f1)
    # b1t2 is in index 13 etc..
    #
    # This way PS can trace back to the original rep structure$rep_desc vector and fetch the original rep
    # descriptions. All this is done in the line below
    temp_vector <- original_rep_structure$rep_desc[match(temp_vector, sub("f.*", "", rep_structure$rep_desc))]
    
    # The only issue with the method above is that the original rep structure might contain fraactionation information
    # (as seen in the example above) and this should be removed afterwards
    
    # Make sure that expdesign (column Sample) contains the original rep descriptions
    # to do so merget temp_vector with the sample column in expdesign
    
    tmp_counter <- 0
    for (expdesign_i in expdesign[,1]){
      expdesign[tmp_counter + 1,1] <- sub("(.*)\\..*",paste0("\\1.", temp_vector[tmp_counter + 1]), expdesign_i)
      tmp_counter <- tmp_counter + 1
    }
  }
  
  # Remove the fractionation information: (if any) as noted above
  expdesign[,1] <- sub("(.*\\..*)f.*", "\\1", expdesign[,1], perl = TRUE)
  
  # expdesign now seems like:
  
  # +--------+----------+
  # | Sample | Category |
  # +--------+----------+
  # | H.b1t1 | H        |
  # +--------+----------+
  # | H.b1t2 | H        |
  # +--------+----------+
  # | H.b2t1 | H        |
  # +--------+----------+
  # | H.b2t2 | H        |
  # +--------+----------+
  # | H.b3t1 | H        |
  # +--------+----------+
  # | H.b3t2 | H        |
  # +--------+----------+
  # | L.b1t1 | L        |
  # +--------+----------+
  # | L.b1t2 | L        |
  # +--------+----------+
  # | L.b2t1 | L        |
  # +--------+----------+
  # | L.b2t2 | L        |
  # +--------+----------+
  # | L.b3t1 | L        |
  # +--------+----------+
  # | L.b3t2 | L        |
  # +--------+----------+
  
  # that in our example is the same as before but if the original rep descs were different than the
  # ones used in the analysis it would now contain the correct original rep descs
  
  # Write the expdesign table to a file
  write.table(expdesign,file="curr_exp_design.txt",row.names=F,quote=F,sep = "\t")
  exp_design_fname<<-"curr_exp_design.txt"
  
  # Make the variable global
  .GlobalEnv[["expdesign"]]<-expdesign
}

perform_analysis<-function() {
  
  # Perform analysis is the main function of ProteoSign's backend that will properly execute
  # subroutines such as read.pgroups that will read the main file and transform all data to a common format
  # do limma analysis and do results plots that do the analysis and draw the plots respectively more
  # functions such as venn diagram creation and functional enrichment analysis will be
  # executed afterwards
  
  levellog("",change=1)
  
  # First we will navigate to the current directory since all paths to parameters files can be relative
  setwd(working_directory)
  
  # Load all tabe separated parameter files taht usually accompany a PS run:
  load_table_param_files()
  
  # Edit the rep structure and original reo structure variables that hold all necessary information for our current experimental structure
  edit_rep_structure()
  
  # Count how many treps and breps we have and store the information to global variables
  # Note that some experiments are not "orthogonal" that means that some breps might have more treps than others
  # In this case n_techreps contains the maximum amount of treps per brep
  .GlobalEnv[["n_bioreps"]]<-max(rep_structure$biorep)
  .GlobalEnv[["n_techreps"]]<-max(ddply(rep_structure[,c("biorep","techrep")],c("biorep"),function(x){return(max(x$techrep))})$V1)
  
  # Protein quantitation is deprecated - always set to true
  if(ProteinQuantitation){
    quantitated_items_lbl<<-"Protein"
  }else{
    quantitated_items_lbl<<-"Peptide"
  }
  
  # Create the output folder:
  if(file.exists(limma_out_dir)){
    levellog("Warn User: The output working folder already exists...")
  } else {
    dir.create(limma_out_dir)
  }
  
  # Remove double quotes from proteomics quantitation datafiles (if any)
  remove_double_quotes()

  # Read the proteomics data
  levellog("Reading input data ...")
  protein_groups<<-read.pgroups(pgroups_fname,evidence_fname,exp_desc)
  
  # Write the proteinGroupsDF file:
  write_proteinGroupsDF()
  
  # The exp_design variable is a very important one that describes the experimental design
  # in a limma - friendly manner. the following function creates thuis gobal variable and saves a copy of that
  # in a file called "curr_exp_design.txt"
  
  create_expdesign()
 
  
  levellog("Performing the analysis ...")
  
  norm.median.intensities <- do_limma_analysis(prepare_working_pgroups(protein_groups),exp_desc,exp_design_fname,exportFormat="pdf",outputFigsPrefix=outputFigsPrefix)
  
  levellog("Generating analysis plots ...")
  

  
  results<-do_results_plots(norm.median.intensities, exp_desc, exportFormat="pdf", outputFigsPrefix=outputFigsPrefix)
  
  #{
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
  #}
  
  # Perform functional enrichment analysis:
  
  if (FNenrichment)
  {
    # Perform the enrichment analysis using the R package gprofiler
    levellog("Perform enrichment analysis.")
    
    # First get the uniprot IDs of the Differentially expressed proteins for each combination of conditions
    ratio_combs<-combinations(nConditions,2,1:nConditions)
    for(j in 1:nrow(ratio_combs))
    {
      result<-tryCatch({
        # In this line conditions[ratio_combs[i, 1]] and conditions[ratio_combs[i, 2]] have the names of the conditions to compare
        uniprot_ids <- c()
        uniprot_ids <- get_uniprot_ids(results, conditions[ratio_combs[j, 1]], conditions[ratio_combs[j, 2]])
        if(length(uniprot_ids) == 0)
        {
          levellog(paste0("Warn User: Functional enrichment analysis for conditions: ", conditions[ratio_combs[j, 1]], " and ", conditions[ratio_combs[j, 2]], " failed (", FNorganism, " was selected as target organism) because no UNIprot IDs were found for the differentially expressed proteins (if any)"))
          next;
        }
        # Since uniprot_ids contain the IDs of the DE expressed proteins lets run a GO analysis for them
        run_enrichment_analysis(uniprot_ids,FNorganism, conditions[ratio_combs[j, 1]], conditions[ratio_combs[j, 2]])
        
        # run_enrichment_analysis produces a file with all enrichment analysis information - take a look at its comments to find out more
        
      }, error = function(err){
        levellog(paste0("Warn User: Functional enrichment analysis for conditions: ", conditions[ratio_combs[j, 1]], " and ", conditions[ratio_combs[j, 2]], " failed (", FNorganism, " was selected as target organism)"))
      })
    }
  }
  
  # Create the Venn diagrams

  levellog("Create Venn Diagrams.")
  
  # First lets build a matrix with columns the different conditions, rows the replicate descriptions (b1t1, b1t2 etc..)
  # and values the indices of the respective columns in results dataframe that store the intensities for the specific
  # condition - replicate pair
  
  intensity_cols_idxs <- sapply(conditions, function(x) grep(paste0("^", x, " \\d+"), colnames(results)))
  
  # The following lines adds the description of the replicates as rownames:
  tmp_vector <- expdesign[expdesign[,"Category"] == conditions[1],"Sample"]
  tmp_vector <- substr(tmp_vector, str_length(conditions[1]) + 2, str_length(tmp_vector))
  rownames(intensity_cols_idxs) <- tmp_vector
  # An example of intensity_cols_idxs is:
  #
  #       H  L
  # b1t1  8 14
  # b1t2  9 15
  # b2t1 10 16
  # b2t2 11 17
  # b3t1 12 18
  # b3t2 13 19
  
  # The following function will generate the diagrams and save them to image files - see the integrated comments for more information
  
  generate_Venn_diagrams(results[,as.vector(intensity_cols_idxs)], rownames(intensity_cols_idxs))
  
  # Return to the main working directory and log the completion of the procedure
  setwd(working_directory)
  levellog("",change=-1)
  levellog("Data analysis finished.")
  levellog("",change=-1)
  return(T)
}

#================ PRODUCTION ===============

if(DEBUG){  
  # Kept for debugging purposes
}else{
  initProteoSign()
  perform_analysis()
}

