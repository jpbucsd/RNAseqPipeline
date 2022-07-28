library("DESeq2")
library("tximport")

#parse command line arguments
#command line usage 1:
#Rscript Rtest.R -1 firstSet set1file1.genes.results set1file2.genes.results ... set1fileN.genes.results -2 secondSet set2file1.genes.results set2file2.genes.results ... set2fileN.genes.results
#command line usage 2:
Rscript Rtest.R -1 set1file1.genes.results set1file2.genes.results ... set1fileN.genes.results -2 set2file1.genes.results set2file2.genes.results ... set2fileN.genes.results

args <- commandArgs(trailingOnly = TRUE)
first <- FALSE
second <- FALSE
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
            named = FALSE
          } else {
            #assume it is -2
            second = TRUE
            first = FALSE
            named = FALSE
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
          } else {
            #assume second
            if ( !named )
             {
                  if (substring(arg, first = nchar(arg) - 13, last = nchar(arg)) == ".genes.results" )
                  {
                        #this is the name of the first file and no default name for the group was provided so we will steal this name
                        secondName = substring(arg, first = 0, last = nchar(arg) - 14)
                        #we also must declare the first variable
                        files2[1] = arg
                        files2ind = files2ind + 1
                        named = TRUE
                  } else {
                        secondName = arg
                        named = TRUE
                  }
             } else {
                  files2[files2ind] = arg
                  files2ind = files2ind + 1
             }
          }
     
    }
 }
 #check arguments
 print(firstName)
 print(files1)
 print(secondName)
 print(files2)
