options(warn=1)
# DEPLOYMENT VERSION
# ======================================
# WARNING: Make sure to install the following packages as administrator/root (for them to be available to all users)
# ======================================
DEBUG <- FALSE;
# source("http://www.bioconductor.org/biocLite.R")
# latest limma version might not be compatible!! limma 3.18.13 suits best for PS
# if(!require("limma")){ biocLite("limma") }
# if(!require("statmod")){ biocLite("statmod") }
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
# #if(!require("devtools")){ { install.packages("devtools", repos="http://cran.fhcrc.org") }
# #if(!require("lineprof")){ devtools::install_github("hadley/lineprof") }

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
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex )
}

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
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

# FROM: http://www-personal.umich.edu/~ladamic/presentations/Rtutorial/Rtutorial.R
panel.lmline = function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = reps.scatter.lmline.colour, ...){
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
#Called by do_results_plot (with smooth=TRUE,scale=TRUE,lm=TRUE)
pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE,lm=FALSE){
  if (smooth){
    if (scale) {
      if(lm){
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

lm_eqn<-function(df){
  m = lm(y ~ x, data=df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  return(as.character(as.expression(eq)));
}

# to be used with apply on a limma-results data.frame, called by do_result_plots
# calculate mean, sd and N of ratios between available channels
# conds_cols_idxs the column indexes of quant channels intensities
# x a limma-results data.frame row
calcRowStats<-function(x,conds_cols_idxs,ratio_combs){
  ratios<-c()
  m<-c()
  std<-c()
  N<-c()
  avg.I<-c()
  # what should we do here in case there is replication mismatch?
  for(i in 1:nrow(ratio_combs)){
    tmp1<-as.numeric(x[conds_cols_idxs[ratio_combs[i,2],]]) #columns of "heavier" label
    tmp2<-as.numeric(x[conds_cols_idxs[ratio_combs[i,1],]]) #columns of "lighter" label    
    ratio_i<-(tmp1-tmp2) #ratios of log-2 transformed data (one per replicate)
    ratios<-rbind(ratios, ratio_i)
    m<-rbind(m, mean(ratio_i,na.rm=T)) #mean of all ratios calculated
    std<-rbind(std, sd(ratio_i,na.rm=T)) #standard deviation of all ratios calculated
    N<-rbind(N, length(which(!is.na(ratio_i) & ratio_i != 0))) # number of all ratios calculated
    avg.I<-rbind(avg.I, mean(c(tmp1,tmp2),na.rm=T))  # average of a
  }
  
  return(c(m[,1],std[,1],N[,1],avg.I[,1],as.vector(t(ratios)))) # columns: mean/standard deviation/number of ratios calculated/average of log2 quant data/the rest of the columns are ratios calculated (one per replicate)
}

# Produces Reproducibility, Volcano plot, MA plot and Scatterplot (matrix). Called by do_limma_analysis subroutine.
do_results_plots<-function(norm.median.intensities,time.point,exportFormat="pdf",outputFigsPrefix=""){
  levellog("",change=1)
  
  ratio_combs<-combinations(nConditions,2,1:nConditions)
  levellog("Preparing data ...")
  results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "")
	if (nrow(ratio_combs) == 1) {
	# limma does not indicate the conditions compared in case only two conditions are compared in the column names so if ratio_combs = 1 add the conditions manually to the colnames of results
	  colnames(results)[grep("p\\.value\\.adj",colnames(results))]<-paste("p.value.adj.",conditions.labels[2],".",conditions.labels[1],sep="")
  }
  
  # Due to the mysterious bug (undefined results$ID column), the following does not make sense ...
  #rownames(results)<-results$ID
  #
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
  d<-data.frame(t(apply(results,1, function(x) calcRowStats(x,conds_cols_idxs,ratio_combs))))
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
	#the following columns are named depending on the experimental structure
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep="") # for example Heavy.Light
    colnames_d_<-c(colnames_d_, paste(paste("log2.",ratio_i_str,sep=""),1:nsamples))
  }
  colnames(d)<-colnames_d_
  # in case of replicate mismatch there will be some columns that contain NA ratios since some replicates do not contain all the conditions
  # Due to the mysterious bug, the following was added ...
  results$ID <- rownames(results)
  results<-cbind(results,d)
  #
  results$ID <- factor(results$ID, levels=unique(as.character(results$ID)))
  results$nID<-1:nrow(results)
  
  levellog("Plotting time ...")
  theme_set(theme_bw())
  # customized colorblind-friendly palette from http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2)/
  cbPalette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")
  #Save the results and othere parameters to a RData file:
  save(results, pThreshold, quantitated_items_lbl, nConditions, calcRowStats, time.point, outputFigsPrefix, conditions.labels, IsobaricLabel, PDdata, log.intensities, norm.intensities, fit2.coefficients, file = "Plot_Generator.RData")
  
  #Plot generation:
  for(i in 1:nrow(ratio_combs)){
    levellog(paste("Generating plots for combination #",i," ..."),change=1,after=T)
	ratio_i_str <- paste(conditions.labels[ratio_combs[i, 2]], ".", conditions.labels[ratio_combs[i, 1]], sep = "") # for example Heavy.Light
    
    # 1 - volcano - -log10 P-value vs log ratio
    result <- tryCatch({
      levellog("Making volcano plot ...")
      figsuffix<-paste("_",ratio_i_str,"-volcano","_",sep="")
      if(exportFormat == "pdf"){
        pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      }
      
      ratio_i_<-paste("log2.",ratio_i_str,sep="")
      ratio_i_sd_col<-paste("log2.sd.",ratio_i_str,sep="")
      tmp2<-results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_, " ")),colnames(results))]]+results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_sd_col, "$")),colnames(results))]]
      
      tmp1<-results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_, " ")),colnames(results))]]-results[,colnames(results)[grep(gsub("\\.","\\\\.",paste0(ratio_i_sd_col, "$")),colnames(results))]]
      ratiolim<-ceiling(max(max(range(tmp1,na.rm=T),range(tmp2,na.rm=T)),abs(min(range(tmp1,na.rm=T),range(tmp2,na.rm=T)))))
      #If two conditions contain exactly the same data ratiolim will be equal to 0. In this case add all the intensities to the same block
      if(ratiolim == 0)
      {
        ratiolim <- 5
      }
      panel.hist.breaks<<-(-ratiolim:ratiolim)
      
      ratio_i_p.value.adj<-paste("p.value.adj.",paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep=""),sep="")
      ratio_i_avg_col<-paste("log2.avg.",ratio_i_str,sep="")
      mlog10_ratio_i_p.value.adj<-paste("mlog10_",ratio_i_p.value.adj,sep="")
      diffexp_ratio_i<-paste("diffexp_",ratio_i_str,sep="")
      
      results[,mlog10_ratio_i_p.value.adj]<-(-log10(results[,ratio_i_p.value.adj]))
      
      na_indexes<-which(is.na(results[,ratio_i_p.value.adj]))
      if(length(na_indexes)>0){
        results[na_indexes,ratio_i_p.value.adj]<-1
        results[,diffexp_ratio_i]<-results[,ratio_i_p.value.adj]<pThreshold
        results[na_indexes,ratio_i_p.value.adj]<-NA
      }else{
        results[,diffexp_ratio_i]<-results[,ratio_i_p.value.adj]<pThreshold
      }
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
	  # volcano plot is nothing more than a plot of the proteins ratios (the average in all replicates) vs -log10 of the p-value for this ratio:
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
      
      results<-results[with(results, order(results[,c(ratio_i_avg_col)])),] # sort the data
      results$nID<-1:nrow(results)
      ratio_i_avg_col_ymax<-paste(ratio_i_avg_col,".ymax",sep="")
      ratio_i_avg_col_ymin<-paste(ratio_i_avg_col,".ymin",sep="")
      results[,ratio_i_avg_col_ymax]<-results[,ratio_i_avg_col]+results[,ratio_i_sd_col]
      results[,ratio_i_avg_col_ymin]<-results[,ratio_i_avg_col]-results[,ratio_i_sd_col]
      
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
	  # the MA plot shows the average intensity of each protein (between the replicates) against the intensity ratio
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
      
	  allratios <- results[, colnames(results)[grep(paste0(ratio_i_, " "), colnames(results))]]
	  allratios <- allratios[, colSums(is.na(allratios)) < nrow(allratios)]
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
  
  levellog("Saving plots results to files ...")
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
  
  results<-results[,which(!grepl("^-log10",colnames(results)))]
  results<-results[,which(!grepl("sd$",colnames(results)))]
  results<-results[,which(!grepl("^nID$",colnames(results)))]
  results<-results[,which(!grepl("^diffexp",colnames(results)))]
  
  quant_species<-"proteins"
  if(!ProteinQuantitation){
    quant_species<-"peptides"
  }  
  
  levellog(paste("do_results_plots: Quantified ",quant_species,": ",nrow(results)," (",time.point,")",sep=""))
  
  conditions.labels <- gsub("\\.", " ", conditions.labels)
  
  for(i in 1:nrow(ratio_combs)){
    col_desc_<-paste("P-value adjusted ",paste(conditions.labels[ratio_combs[i,2]],"/",conditions.labels[ratio_combs[i,1]],sep=""),sep="")
    ndiffexp_tmp<-length(which(results[,col_desc_]<pThreshold))
    levellog(paste("do_results_plots: Differentially expressed for ",conditions.labels[ratio_combs[i,2]]," vs ",conditions.labels[ratio_combs[i,1]]," : ",ndiffexp_tmp,sep=""))
  }
  
  if(nrow(ratio_combs) > 1){
    levellog(paste("do_results_plots: Differentially expressed in at least one combination of conditions: ",ndiffexp,sep=""))
  }
  
  diffexp<-results[,c(grep("^P-value adjusted",colnames(results)),grep("avg log2 [^I]+",colnames(results)),grep("sd log2 ",colnames(results)),grep("N log2 ",colnames(results)),grep("avg log2 I ",colnames(results)))]
  diffexp[,quantitated_items_lbl]<-rownames(diffexp)
  
  diffexp<-diffexp[,c(grep(quantitated_items_lbl,colnames(diffexp)),grep("avg log2 [^I]+",colnames(diffexp)),grep("log2 sd ",colnames(diffexp)),grep("P-value ",colnames(diffexp)), grep("log2 N ",colnames(diffexp)),grep("log2 avg I ",colnames(diffexp)))]
  tmp_protein_groups<-protein_groups
  colnames(tmp_protein_groups)[grep(paste(quantitated_items_lbl,".IDs",sep=""),colnames(tmp_protein_groups))]<-quantitated_items_lbl
  diffexp<-merge(diffexp,tmp_protein_groups[,c(quantitated_items_lbl,sort(colnames(tmp_protein_groups)[grep("Ratio\\.counts",colnames(tmp_protein_groups))]))],by=quantitated_items_lbl,all.x=T)
  diffexp$newcol<-rowSums(diffexp[,colnames(diffexp)[grep("Ratio\\.counts",colnames(diffexp))]],na.rm=T)
  colnames(diffexp)[length(colnames(diffexp))]<-"Ratio.counts.total"  
  
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
  #Restore original rep description for output
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
  write.table(diffexp[signTruth,],dec=dec,file=paste(outputFigsPrefix,"_diffexp_",time.point,".txt",sep=""),sep="\t",row.names=F,quote=F)
  colnames(diffexp) <- oldcolumns
  diffexp<-merge(diffexp,results[,-grep("^(avg log2|P-value adjusted)",colnames(results))],by.x=c(quantitated_items_lbl),by.y=c("ID"),all.x=T)
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
    results<-read.table(paste(outputFigsPrefix,"_condition-i_vs_condition-j_",time.point,".txt",sep=""), header = T, sep = "\t",quote='',stringsAsFactors=F,comment.char = "")
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
      pgroups_colname<-'Protein.Accessions'
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
      intensityCol <- 'Intensity'
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
			if (all(pgroups[, paste0(rep_desc_i, ".uniqueSequences.", cond_i)] == 0)) {
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

