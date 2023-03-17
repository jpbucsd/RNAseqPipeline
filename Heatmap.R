library("DESeq2")
library("tximport")
library("pheatmap")

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

indexZ <- 0
indexS <- c()
cIndex <- 0
setIndex <- 0


for (arg in args) {
    if(substring(arg, first = 1, last = 1) == "-")
    {
          if(otherSets)
          {
            #other sets was previously set, finish it up
            indexS[setIndex] = cIndex
            setFiles <- append(setFiles,setFile)
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
          } else if (substring(arg, first = 1, last = 2) == "-d") {
            #assume it is -2
            zeroSet = FALSE
            otherSets = FALSE
            directory = TRUE
            outName = FALSE
            fNamed = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-o") {
            #assume it is -2
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = TRUE
            fNamed = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-o") {
            zeroSet = FALSE
            otherSets = FALSE
            directory = FALSE
            outName = FALSE
            fNamed = TRUE
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
          }
     
    }
 }
 #the last set file is not added yet
 setFiles <- append(setFiles,setFile)
 indexS[setIndex] = cIndex

 
loopIndex <- 0

rlogSet <- c()
resultsSet <- c()
foldChanges <- data.frame()

#we actually want to compare the zeroset to itself to make it come out zero in the heatmap
if(TRUE){
  #create conditions
  loopIndex = loopIndex + 1
  conditions <- c(rep(zeroName,indexZ),rep(zeroName,indexZ))
  #resume converting from below
  files <- c()
  snames <- c()
    
  
    
  #zero set first
  for (file in zeroFiles) {
      snames <- append(snames, substring(file, first = 0, last = nchar(file) - 14))
      fname <- paste(dirPath,file,sep="/")
      files <- append(files,fname)
  }
  for (file in zeroFiles) {
      snames <- append(snames, substring(file, first = 0, last = nchar(file) - 14))
      fname <- paste(dirPath,file,sep="/")
      files <- append(files,fname)
  }
    
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
    
  resultsSet <- append(resultsSet,res)
  
  #res contains the results for this one
  #write results
  ofnnnname <- paste(zeroName,zeroName,sep="_vs_")
  ofnnname <- paste("control",ofnnnname,sep="/")
  ofnname <- paste(dirPath,ofnnname,sep="/")
  ofname <- paste(ofnname,"csv",sep=".")
  write.csv(as.data.frame(res), file=ofname)

  #create normalized counts for heatmap
  rlog_out <- assay(rlog(ddsTxi, blind=FALSE)) #normalized count data from the DESeq object
    
  rlogSet <- append(rlogSet,rlog_out)  
  
  nomnnnnname <- paste(zeroName,zeroName,sep="_vs_")
  nomnnnname <- paste(nomnnnnname,"normalizedCounts",sep="_")
  nomnnname <- paste("control",nomnnnname,sep="/")
  nomnname <- paste(dirPath,nomnnname,sep="/")
  nomname <- paste(nomnname,"csv",sep=".")
  write.csv(as.data.frame(rlog_out), file=nomname)   
  
  #copy row names from rlog_out
  row.names(foldChanges) <- row.names(rlog_out)
  
  #make a temporary dataframe
  tdf <- data.frame()
  row.names(tdf) <- row.names(rlog_out)
  #we may want to filter by pvalue before this step
  tdf.cbind(select(rlog_out, c('log2FoldChange'))
  names(tdf)[names(tdf) == 'log2FoldChange'] <- zeroName
  #fuse dataframes
  foldChanges <- merge(foldChanges, tdf, by = 'row.names', all = TRUE)
}




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
      snames <- append(snames, substring(file, first = 0, last = nchar(file) - 14))
      fname <- paste(dirPath,file,sep="/")
      files <- append(files,fname)
  }
    
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
    
  resultsSet <- append(resultsSet,res)
  
  #res contains the results for this one
  #write results
  ofnnnname <- paste(zeroName,zeroName,sep="_vs_")
  ofnnname <- paste("control",ofnnnname,sep="/")
  ofnname <- paste(dirPath,ofnnname,sep="/")
  ofname <- paste(ofnname,"csv",sep=".")
  write.csv(as.data.frame(res), file=ofname)

  #create normalized counts for heatmap
  rlog_out <- assay(rlog(ddsTxi, blind=FALSE)) #normalized count data from the DESeq object
  
  rlogSet <- append(rlogSet,rlog_out)  
  
  nomnnnnname <- paste(zeroName,zeroName,sep="_vs_")
  nomnnnname <- paste(nomnnnnname,"normalizedCounts",sep="_")
  nomnnname <- paste("control",nomnnnname,sep="/")
  nomnname <- paste(dirPath,nomnnname,sep="/")
  nomname <- paste(nomnname,"csv",sep=".")
  write.csv(as.data.frame(rlog_out), file=nomname)
    
  #make a temporary dataframe
  tdf <- data.frame()
  row.names(tdf) <- row.names(rlog_out)
  #we may want to filter by pvalue before this step
  tdf.cbind(select(rlog_out, c('log2FoldChange'))
  names(tdf)[names(tdf) == 'log2FoldChange'] <- zeroName
  #fuse dataframes
  foldChanges <- merge(foldChanges, tdf, by = 'row.names', all = TRUE)
}
write.csv(foldChanges, file=paste(outName, paste(outName, paste(fName, "csv", sep="."), sep="/"))
#use pheatmap to create the heatmap. we want to cluster by genes which are rows
#we cannot use a dataframe, foldChanges must be turned into a matrix
df_num = as.matrix(foldChanges)

pheatmap(df_num,cluster_rows=TRUE,legend=TRUE,show_rownames=TRUE,show_colnames=TRUE,filename=paste(outName, paste(fName, "pdf", sep="."), sep="/"))
#notes for the future of pheatmap, annotation row will take a dataframe that combines rows into larger groups which will be displayed with an annotation
#annotation_col does the same for columns. annotation col should be used to group samples. annotation_names_col will display the names
#main can give a name to the entire plot, fontsize_row and fontsize_col can change the font size which may be helpful
