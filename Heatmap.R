library("DESeq2")
library("tximport")
library("pheatmap")
library("dplyr")
#parse command line arguments
#command line usage 1:

#Rscript Heatmap.R -Z zeroSet zerofile1.genes.results zerofile2.genes.results ... zerofileN.genes.results -S setN setNfile1.genes.results ... setNfile2.genes.results -d path/to/directory/of/files -o /path/to/output/file -n filename 

zeroName <- ""
zeroFiles <- c()
setNames <- c()
setFile <- c()
setFiles <- list()


args <- commandArgs(trailingOnly = TRUE)
zeroSet <- FALSE
otherSets <- FALSE
directory <- FALSE
dirPath <- ""
outName <- FALSE
outPath <- ""
fNamed <- FALSE
fName <- ""
named <- FALSE
zeroVsAll <- FALSE

indexZ <- 0
indexS <- c()
cIndex <- 0
setIndex <- 0

filtFlag <- FALSE
filter <- 0

numFlag <- FALSE
numClusters <- 0

mathArtifact <- FALSE

for (arg in args) {
    if(substring(arg, first = 1, last = 1) == "-")
    {
          if(otherSets)
          {
            #other sets was previously set, finish it up
            indexS[setIndex] = cIndex
            print("adding set file")
            print(setFile)
            setFiles[[setIndex]] <- setFile

            setFile = c()
            cIndex = 0
          }
          if (substring(arg, first = 1, last = 2) == "-Z")
          {
            zeroSet = TRUE
            otherSets = FALSE
            directory = FALSE
            named = FALSE
            outName = FALSE
            fNamed = FALSE
            filtFlag = FALSE
            numFlag = FALSE
            #reset the index
            index = 0
          } else if (substring(arg, first = 1, last = 2) == "-S") {
            #assume it is -2
            #reset the index
            index = 0
            zeroSet = FALSE
            otherSets = TRUE
            directory = FALSE
            named = FALSE
            outName = FALSE
            fNamed = FALSE
            filtFlag = FALSE
            numFlag = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-d") {
            #assume it is -2
            zeroSet = FALSE
            otherSets = FALSE
            directory = TRUE
            outName = FALSE
            fNamed = FALSE
            filtFlag = FALSE
            numFlag = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-o") {
            #assume it is -2
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = TRUE
            fNamed = FALSE
            filtFlag = FALSE
            numFlag = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-n") {
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = FALSE
            fNamed = TRUE
            filtFlag = FALSE
            numFlag = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-c") {
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = FALSE
            fNamed = FALSE
            zeroVsAll = TRUE
            filtFlag = FALSE
            numFlag = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-f") {
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = FALSE
            fNamed = FALSE
            filtFlag = TRUE
            numFlag = FALSE
          }else if (substring(arg, first = 1, last = 2) == "-G") {
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = FALSE
            fNamed = FALSE
            filtFlag = FALSE
            numFlag = TRUE
          } else if (substring(arg, first = 1, last = 2) == "-m") {
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = FALSE
            fNamed = FALSE
            mathArtifact = TRUE
            filtFlag = FALSE
          }
    } else {
          if ( zeroSet )
          {
             if ( !named )
             {
                  zeroName = arg
                  named = TRUE
             } else {
                  indexZ = indexZ + 1
                  zeroFiles[indexZ] = arg
             }
          } else if (otherSets) {
            #assume second
            if ( !named )
            {
                  #add the last setFile to setFiles and reset it so that we can start a new one
                  
                  setIndex = setIndex + 1
                  setNames[setIndex] = arg
                  named = TRUE
             } else {
                  print(cat("adding ",arg," to set file"))
                  cIndex = cIndex + 1
                  setFile[cIndex] = arg
             }
          } else if (directory)
          {
              dirPath = arg
          } else if (outName)
          {
              outPath = arg
          } else if (fNamed)
          {
              fName = arg
          } else if (filtFlag)
          {
              filter = arg
          } else if (numFlag)
          {
              numClusters = arg
          }
     
    }
 }
#the last set file is not added yet
#indexS[setIndex] = cIndex
#setFiles[[setIndex]] <- setFile
#print("set files:")
#print(setFiles)

loopIndex <- 0

rlogSet <- c()
#resultsSet <- c()
foldChangesShared <- data.frame()
foldChangesAll <- data.frame()
foldChangeInit <- FALSE
#make true if we actually want to compare the zeroset to itself to make it come out zero in the heatmap, DESEQ does not permit this. to accomplish this use the following code
#foldChanges[,"control_vs_control"] <- 0
#print("outfilename:")
#print(paste(outPath, paste(fName, "csv", sep="."), sep="/"))

#for the zscored heatmap we need a dataframe for counts

#this must include the zeroset
zfiles <- c()
for (file in zeroFiles) {
      fname <- paste(dirPath,file,sep="/")
      zfiles <- append(zfiles,fname)
}
ztxi <- tximport(zfiles, type = "rsem")
#head(ztxi$counts)
countsShared<-data.frame(col1=rowMeans(ztxi$counts,na.rm=TRUE))

row.names(countsShared)<-row.names(ztxi$counts)
colnames(countsShared)[1] <- zeroName

#head(countsShared)

#now add all the other sets
for (set in setFiles) {
    loopIndex = loopIndex + 1
    files <- c()
    for (file in set) {
      fname <- paste(dirPath,file,sep="/")
      files <- append(files,fname)
    }
    txi <- tximport(files, type = "rsem")
    tempFrame<-data.frame(col1=rowMeans(txi$counts,na.rm=TRUE))
    row.names(tempFrame)<-row.names(txi$counts)
    colnames(tempFrame)[1] <- setNames[loopIndex]
    
    #merge the dataframes
    countsShared <- merge(countsShared, tempFrame, by = 0, all = FALSE)
    row.names(countsShared) = countsShared[,"Row.names"]
    countsShared <- countsShared[,-1]
}

cnname <- paste(outPath,"counts",sep="/")
cname <- paste(cnname,"csv",sep=".")
write.csv(countsShared, file=cname)
head(countsShared)

#if this dataframe is zscored, meaningful data may be obscured by the pressence of rows where only one sample has counts, which produces the highest zscore. These will be considered mathematical artifacts with the -m flag and removed
if(mathArtifact)
{
  #countsShared <- countsShared[rowSums(countsShared == 0) < ncol(countsShared)-1,]
  #the above line of code removes lines where only one of the rows is not 0, however it appears to truly remove artifacts, all rows that are not zero must be removed
  countsShared <- countsShared[rowSums(countsShared == 0) ==0,]
}

#create a dataframe with the averages and standard deviations of each row for counts
stats<-data.frame(avg=rowMeans(countsShared,na.rm=TRUE),stdev=apply(countsShared, 1, sd, na.rm=TRUE))
#zscore the counts

for(i in 1:nrow(countsShared)) {
    countsShared[i,] <- (countsShared[i,] - stats[i,'avg'])/stats[i,'stdev']
}
head(countsShared)

#create a heatmap for all genes, with no gene names
df_all = as.matrix(countsShared)
pheatmap(df_all,cluster_rows=FALSE,cluster_cols=FALSE,legend=TRUE,show_rownames=FALSE,show_colnames=TRUE,fontsize_row=1,color=colorRampPalette(c("navy", "white", "red"))(50),filename=paste(outPath, paste(paste(fName,"zscore_all",sep="_"), "pdf", sep="."), sep="/"))

znname <- paste(outPath,"zscoredcounts",sep="/")
zname <- paste(znname,"csv",sep=".")
write.csv(countsShared, file=zname)
#create a filtered heatmap for filtered genes
countsFiltered <- countsShared %>% filter_all(any_vars(.>as.double(filter)|-(as.double(filter))>.))
df_filt = as.matrix(countsFiltered)

zfnname <- paste(outPath,"zscoredcounts_filtered",sep="/")
zfname <- paste(zfnname,"csv",sep=".")
write.csv(countsFiltered, file=zfname)

heat <- pheatmap(df_filt,cluster_rows=TRUE,cluster_cols=FALSE,legend=TRUE,show_rownames=TRUE,show_colnames=TRUE,fontsize_row=1,color=colorRampPalette(c("navy", "white", "red"))(50),filename=paste(outPath, paste(paste(fName,"zscore_filtered",sep="_"), "pdf", sep="."), sep="/"))
if(numClusters > 1)
{
    #get all clusters
    hclusters <-heat$tree_row
    #turn into the # of clusters desired
    nclusters <- cutree(hclusters, numClusters)
    #print the clusters to a txt
    write.table(nclusters,file=paste(outPath, paste(paste(fName,"zscore_filtered_clusters",sep="_"), "txt", sep="."), sep="/"),sep=",",col.names=TRUE,row.names=TRUE);
    #separate by cluster
    #convert to dataframe
    dfc = data.frame(nclusters)
    #print(rownames(dfc)[dfc[,1] == 1])
    for (x in 1:numClusters) {
        write.table(rownames(dfc)[dfc[,1] == x],file=paste(outPath, paste(paste(fName,paste("zscore_filtered_clusters",as.character(x),sep="_"),sep="_"), "txt", sep="."), sep="/"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
    write.table(rownames(dfc),file=paste(outPath, paste(paste(fName,paste("zscore_filtered_clusters","background",sep="_"),sep="_"), "txt", sep="."), sep="/"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

#the following code produces a heatmap comparing the zero set to all other sets and taking a heatmap of the log2fold.
if(zeroVsAll)
{
    for (set in setFiles) {
      #create conditions
      loopIndex = loopIndex + 1
      conditions <- c(rep(zeroName,indexZ),rep(setNames[loopIndex],indexS[loopIndex]))

      #resume converting from below
      files <- c()
      snames <- c()
      #zero set first
      for (file in zeroFiles) {
          snames <- append(snames, substring(file, first = 0, last = nchar(file) - 14))
          fname <- paste(dirPath,file,sep="/")
          files <- append(files,fname)
      }

      for (file in set) {
          print(file)
          snames <- append(snames, substring(file, first = 0, last = nchar(file) - 14))
          fname <- paste(dirPath,file,sep="/")
          files <- append(files,fname)
      }
      print(cat("set ",loopIndex," ",setNames[loopIndex]))
      samples <- data.frame("run"=snames,"condition"=conditions)
      names(files) = samples$run
      #convert RSEM results
      txi <- tximport(files, type = "rsem")
      txi$length[txi$length == 0] <- 1
      ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples, design = ~ condition)
      #filtering, filter low counts to ignore them
      keep <- rowSums(counts(ddsTxi)) >= 10
      ddsTxi <- ddsTxi[keep,]

      ##### Perform deseq2 #####
      ddsTxi <- DESeq(ddsTxi)
      res <- results(ddsTxi)
      #resultsSet <- append(resultsSet,res)
      #res contains the results for this one
      #write results
      ofnnname <- paste(zeroName,setNames[loopIndex],sep="_vs_")
      ofnname <- paste(outPath,ofnnname,sep="/")
      ofname <- paste(ofnname,"csv",sep=".")
      write.csv(as.data.frame(res), file=ofname)

      #create normalized counts for heatmap
      rlog_out <- assay(rlog(ddsTxi, blind=FALSE)) #normalized count data from the DESeq object

      rlogSet <- append(rlogSet,rlog_out)  

      nomnnnname <- paste(zeroName,setNames[loopIndex],sep="_vs_")
      nomnnname <- paste(nomnnnname,"normalizedCounts",sep="_")
      nomnname <- paste(outPath,nomnnname,sep="/")
      nomname <- paste(nomnname,"csv",sep=".")
      write.csv(as.data.frame(rlog_out), file=nomname)


      #rlog out will contain data for comparing each replicant. we do not want this as this is whats produced by heatmap.py
      #res contains log2foldchange for the entire sample against the zero.

      #make a temporary dataframe
      tdf <- data.frame(matrix(0, ncol = 1, nrow = length(row.names(res))))
      #print(length(row.names(tdf)))
      row.names(tdf) <- row.names(res)
      #we may want to filter by pvalue before this step
      tdf[,1] <- res[,"log2FoldChange"]
      #tdf.cbind(res[,c('log2FoldChange')])
      #names(tdf)[names(tdf) == 'log2FoldChange'] <- setNames[loopIndex]
      #print(setNames[loopIndex])
      colnames(tdf)[1] <- setNames[loopIndex]

      if(foldChangeInit)
      { 
        #fuse new data with previous dataframe, keeping all genes even if they do not appear in each sample. this is for the csv.
        foldChangesAll <- merge(foldChangesAll, tdf, by = 0, all = TRUE)
        #the gene names have now become a column, make the row names the gene names
        row.names(foldChangesAll) = foldChangesAll[,"Row.names"]
        #remove the column containing row names
        foldChangesAll <- foldChangesAll[,-1]

        #repeat for the shared dataframe, but removing genes only existing in one file. this is for the heatmap.
        #fuse new data with previous dataframe, keeping all genes even if they do not appear in each sampl
        foldChangesShared <- merge(foldChangesShared, tdf, by = 0, all = FALSE)
        row.names(foldChangesShared) = foldChangesShared[,"Row.names"]
        foldChangesShared <- foldChangesShared[,-1]
      }else{
        #the dataframes are empty. they should just be identical to tdf
        foldChangesAll <- tdf
        foldChangesShared <- tdf
        foldChangeInit <- TRUE
      }
      #print("debug fold changes ")
      #print(head(foldChangesAll,10))
      #print(tail(foldChangesAll,10))
      #print("fold changes all^")
      #print(head(foldChangesShared,10))
      #print(tail(foldChangesShared,10))
      #print("fold changes shared^")
    }
    write.csv(foldChangesAll, file=paste(outPath, paste(fName, "csv", sep="."), sep="/"))
    #use pheatmap to create the heatmap. we want to cluster by genes which are rows
    #we cannot use a dataframe, foldChanges must be turned into a matrix
    df_num = as.matrix(foldChangesShared)

    pheatmap(df_num,cluster_rows=TRUE,legend=TRUE,show_rownames=TRUE,show_colnames=TRUE,fontsize_row=1,color=colorRampPalette(c("navy", "white", "red"))(50),filename=paste(outPath, paste(fName, "pdf", sep="."), sep="/"))
    #notes for the future of pheatmap, annotation row will take a dataframe that combines rows into larger groups which will be displayed with an annotation
    #annotation_col does the same for columns. annotation col should be used to group samples. annotation_names_col will display the names
    #main can give a name to the entire plot, fontsize_row and fontsize_col can change the font size which may be helpful

    #the following code is for a filtered heatmap
    foldChangesFiltered <- foldChangesShared %>% filter_all(any_vars(.>6|-6>.))
    df_num2 = as.matrix(foldChangesFiltered)
    pheatmap(df_num2,cluster_rows=TRUE,legend=TRUE,show_rownames=TRUE,show_colnames=TRUE,fontsize_row=3,color=colorRampPalette(c("navy", "white", "red"))(50),filename=paste(outPath, paste(paste(fName,"filtered",sep="_"), "pdf", sep="."), sep="/"))
    
}
