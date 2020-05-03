# Dear ProteoSign user,
# Please find below the code that ProteoSign uses to generate the data plots.
# The two main functions are: do_results_plots, which produces the Reproducibility plot, the Volcano plot, the MA plot and the Scatterplot (matrix),
# and do_limma_plots, which produces the replicates' intensities boxplots before and after normalization, as well as the average intensity histogram.

options(warn=1)

source("http://www.bioconductor.org/biocLite.R")
if(!require("ggplot2"))
{
  install.packages("ggplot2", repos="http://cran.fhcrc.org")
  library(ggplot2)
}
if(!require("gtools"))
{
  install.packages("gtools", repos="http://cran.fhcrc.org")
  library(gtools)
}

# do_results_plots produces the Reproducibility plot, the Volcano plot, the MA plot and the Scatterplot (matrix)

do_results_plots<-function(){
  #ratio_combs contains the combinations of the conditions
  ratio_combs<-combinations(nConditions,2,1:nConditions)
  
  #Set the theme in ggplot2:
  theme_set(theme_bw())
  
  # cbPalette will be used in creating the plots
  # the default one is a customized colorblind-friendly palette from http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2)/
  cbPalette <- c("#999999", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7")
    
  #Plot generation:
  for(i in 1:nrow(ratio_combs)){
    #Prepare the combination:
    print(paste("Generating plots for combination #",i," ..."),change=1,after=T)
    result <- tryCatch({
      ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],".",conditions.labels[ratio_combs[i,1]],sep="")
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
    }, error = function(err){
      print(paste0("Warning! ", ratio_i_str, " combination preparation failed!"))
    })
	
    # 1 - volcano - -log10 P-value vs log ratio
    result <- tryCatch({
      print("Making volcano plot ...")
	  #Customize the filename and the plot size by editing the following two lines:
      figsuffix<-paste("_",ratio_i_str,"-volcano","_",sep="")
      pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
	  #Data preparation:
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
      #The following lines optimize the plot's x-label in specific dataset types
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
      
      # p is a plot created by the ggplot library
	  # Change the next command to suit your needs:
      p<-ggplot(data=results, aes_string(x=ratio_i_avg_col, y=mlog10_ratio_i_p.value.adj, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=1.75) +
        theme(legend.position = "none", axis.title.y=element_text(vjust=0.2), axis.title.x=element_text(vjust=0), plot.title = element_text(vjust=1.5, lineheight=.8, face="bold")) +
        xlim(c(-ratiolim, ratiolim)) + ylim(c(0, 6)) + scale_colour_manual(values=cbPalette) +
        xlab(myxlab) + ylab("-log10 P-value") + ggtitle("P-value vs Fold change") +
        geom_hline(aes(yintercept=-log10(pThreshold)), colour="#990000", linetype="dashed") +
        geom_text(size=2.5, hjust=1, vjust=-0.5,aes(x=-4.2, y=-log10(pThreshold)), label=paste0("P-value=", pThreshold),colour="#990000")
      print(p)
      dev.off()
    }, error = function(err){
      print(paste0("Warning! ", ratio_i_str, " volcano plot failed"))
    })
    
    # 2 - value-ordered - log ratio
    result <- tryCatch({
      print("Making value-ordered plot ...")
	  #Customize the filename and the plot size by editing the following two lines:
      figsuffix<-paste("_",ratio_i_str,"-value-ordered-log-ratio","_",sep="")
      pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      #Data preparation:
      results<-results[with(results, order(results[,c(ratio_i_avg_col)])),]
      results$nID<-1:nrow(results)
      ratio_i_avg_col_ymax<-paste(ratio_i_avg_col,".ymax",sep="")
      ratio_i_avg_col_ymin<-paste(ratio_i_avg_col,".ymin",sep="")
      results[,ratio_i_avg_col_ymax]<-results[,ratio_i_avg_col]+results[,ratio_i_sd_col]
      results[,ratio_i_avg_col_ymin]<-results[,ratio_i_avg_col]-results[,ratio_i_sd_col]
      #The following lines optimize the plot's y-label in specific dataset types
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
      
      # p is a plot created by the ggplot library
	  # Change the next command to suit your needs:
      p<-ggplot(data=results, aes_string(x="nID", y=ratio_i_avg_col, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=1.5) +
        geom_errorbar(aes_string(ymin=ratio_i_avg_col_ymin, ymax=ratio_i_avg_col_ymax), width=1.5) +
        theme(legend.position = "none", axis.title.y=element_text(vjust=0.2), axis.title.x=element_text(vjust=0), plot.title = element_text(vjust=1.5, lineheight=.8, face="bold")) +
        ylim(c(-ratiolim, ratiolim)) + scale_colour_manual(values=cbPalette) +
        xlab(paste(quantitated_items_lbl,"ID")) + ylab(myylab) + ggtitle("Value-ordered fold change")
      print(p)
      dev.off()    
    }, error = function(err){
      print(paste0("Warning! ", ratio_i_str, " value-ordered plot failed"))
    })
    
    # 3 - MA plot
    result <- tryCatch({
      print("Making MA plot ...")
	  #Customize the filename and the plot size by editing the following two lines:
      figsuffix<-paste("_",ratio_i_str,"-MA","_",sep="")
      ratio_i_avgI_col<-paste("log2.avg.I.",ratio_i_str,sep="")
      pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      #The following lines optimize the plot's y-label in specific dataset types
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
      
      # p is a plot created by the ggplot library
	  # Change the next command to suit your needs:
      p<-ggplot(data=results, aes_string(x=ratio_i_avgI_col, y=ratio_i_avg_col, colour=diffexp_ratio_i)) +
        geom_point(alpha=0.7, size=1.75) +
        theme(legend.position = "none", axis.title.y=element_text(vjust=0.2), axis.title.x=element_text(vjust=0), plot.title = element_text(vjust=1.5, lineheight=.8, face="bold")) +
        ylim(c(-ratiolim, ratiolim)) + scale_colour_manual(values=cbPalette) +
        xlab("M (average log2 Intensity)") + ylab(myylab) + ggtitle("MA plot")
      print(p)
      dev.off()
    }, error = function(err){
      print(paste0("Warning! ", ratio_i_str, " MA plot failed"))
    })
    
    # 4 - Reproducibility plots & histograms
    result <- tryCatch({
      print("Making reproducibility plot ...")
	  #Customize the filename suffix by editing the following line:
      figsuffix<-paste("_",ratio_i_str,"-reproducibility","_",sep="")
      allratios<-results[,colnames(results)[grep(paste0(ratio_i_, " "),colnames(results))]]
      #The following lines optimize the plot's y-label in specific dataset types
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
	  #Customize the filename and the size of the plot by editing the following line:
      pdf(file=paste(outputFigsPrefix,figsuffix,time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
      pairs.panels(allratios,scale=T,lm=T)
      dev.off()
    }, error = function(err){
      print(paste0("Warning! ", ratio_i_str, " reproducibility plot failed"))
    })
  }
}

# do_limma_plots draws the limma boxplots and the limma histograms in one pdf file:

do_limma_plots<-function()
{
  ratio_combs<-combinations(nConditions,2,1:nConditions)
  pdf(file=paste(outputFigsPrefix,"_limma-graphs_",time.point,".pdf",sep=""),width=10, height=7, family = "Helvetica", pointsize=8)
  # Create the intensities before normalisation boxplot
  print("Making Intensities before normalisation limma boxplot")
  boxplot(log.intensities)
  title(main="Intensities Before Normalisation")
  # Create the intensities after normalisation boxplot
  print("Making Intensities after normalisation limma boxplot")
  boxplot(norm.intensities)
  title(main="Intensities After Normalisation")
  #Create the limma histograms for each combination:
  print("Making limma histograms")
  for(i in 1:nrow(ratio_combs)){
    ratio_i_str<-paste(conditions.labels[ratio_combs[i,2]],"/",conditions.labels[ratio_combs[i,1]],sep="")
    hist(fit2.coefficients[,i],main=paste("Log2 Fold Change ",ratio_i_str,sep=""), xlab="Log2 Fold Change", breaks=50 )
  }
  dev.off()
}
# For more customization options someone can modify the following functions:

# FROM: http://musicroamer.com/blog/2011/01/16/r-tips-and-tricks-modified-pairs-plot/

# The following functions are used to draw the different parts of the Reproducibility plot

# panel.cor.scale displays the correllation coeeficients in the upper right half of the Reproducibility plot
# the size of the text is proportional to the value of the coefficient

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

#panel.cor is not called by default but can replace panel.cor.scale if scaling the text acording to the R value is not desirable

panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex )
}

# panel.hist draws the histograms in the diagonal of the Reproducibility plot

panel.hist <- function(x, ...){
  #ratios.hist.colour is the colour of the histogram columns
  ratios.hist.colour<-"cyan"
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

# panel.lmline creates the scatterplots displayed in the bottom left half of the Reproducibility plot

# FROM: http://www-personal.umich.edu/~ladamic/presentations/Rtutorial/Rtutorial.R
panel.lmline = function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "red", ...){
  #Note: col.smooth is the colour of the linear regression line (by default red)
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
      print("Warning!: panel.lmline: found abline with NA slope, the regression line will not be drawn")
    }
  }
}

#Called by do_results_plot (by default smooth=TRUE,scale=TRUE,lm=TRUE)
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

#MAIN proccess:

#Load the necessary variables: the file Plot_Generator.RData must be contained in the same folder with this script
load("Plot_Generator.RData", .GlobalEnv)
#Draw the basic plots:
do_results_plots()
# Draw the limma plots:
do_limma_plots()
print("Procedure finished")