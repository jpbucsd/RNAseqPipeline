library("DESeq2")
library("tximport")

#parse command line arguments
#command line usage 1:
#Rscript Rtest.R -1 firstSet set1file1.genes.results set1file2.genes.results ... set1fileN.genes.results -2 secondSet set2file1.genes.results set2file2.genes.results ... set2fileN.genes.results -d path/to/directory
#command line usage 2:
#Rscript Rtest.R -1 set1file1.genes.results set1file2.genes.results ... set1fileN.genes.results -2 set2file1.genes.results set2file2.genes.results ... set2fileN.genes.results -d path/to/directory

args <- commandArgs(trailingOnly = TRUE)
first <- FALSE
second <- FALSE
directory <- FALSE
dirPath <- ""
named <- FALSE
firstName <- ""
secondName <- ""
files1 <- c()
files2 <- c()
files1ind <- 1
files2ind <- 1
for (arg in args) {
    if(substring(arg, first = 1, last = 1) == "-")
    {
          if (substring(arg, first = 1, last = 2) == "-1")
          {
            first = TRUE
            second = FALSE
            directory = FALSE
            named = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-2") {
            #assume it is -2
            second = TRUE
            first = FALSE
            directory = FALSE
            named = FALSE
          } else if (substring(arg, first = 1, last = 2) == "-d") {
            #assume it is -2
            second = FALSE
            first = FALSE
            directory = TRUE
          }
    } else {
          if ( first )
          {
             if ( !named )
             {
                  if (substring(arg, first = nchar(arg) - 13, last = nchar(arg)) == ".genes.results" )
                  {
                        #this is the name of the first file and no default name for the group was provided so we will steal this name
                        firstName = substring(arg, first = 0, last = nchar(arg) - 14)
                        #we also must declare the first variable
                        files1[1] = arg
                        files1ind = files1ind + 1
                        named = TRUE
                  } else {
                        firstName = arg
                        named = TRUE
                  }
             } else {
                  files1[files1ind] = arg
                  files1ind = files1ind + 1
             }
          } else if (second) {
            #assume second
            if ( !named )
             {
                  if (substring(arg, first = nchar(arg) - 13, last = nchar(arg)) == ".genes.results" )
                  {
                        #this is the name of the first file and no default name for the group was provided so we will steal this name
                        secondName = substring(arg, first = 0, last = nchar(arg) - 14)
                        #we also must declare the first variable
                        files2[1] = arg
                        #files2ind = files2ind + 1
                        named = TRUE
                  } else {
                        secondName = arg
                        named = TRUE
                  }
             } else {
                  files2ind = files2ind + 1
                  files2[files2ind] = arg
             }
          } else if (directory)
          {
              dirPath = arg
          }
     
    }
 }

conditions <- c(rep(firstName,files1ind),rep(secondName,files2ind))
files <- c()
snames <- c()
for (file in files1) {
    append(snames, substring(file, first = 0, last = nchar(file) - 14))
    fname <- paste(dirPath,file,sep="/")
    append(files,fname)
}
for (file in files2) {
    append(snames, substring(file, first = 0, last = nchar(file) - 14))
    fname <- paste(dirPath,file,sep="/")
    append(files,fname)
}
samples <- data.frame("run"=snames,"condition"=conditions)
names(files) = samples$run

#convert RSEM results
txi <- tximport(files, type = "rsem")
txi$length[txi$length == 0] <- 1
ddsTxi <- DESeqDataSetFromTximport(txi,coData = samples, design = ~ condition)

#filtering, filter low counts to ignore them
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep,]

#write results
ofnname <- paste(firstName,secondName,sep="_vs_")
write.csv(as.data.frame(res), file=ofname)
                