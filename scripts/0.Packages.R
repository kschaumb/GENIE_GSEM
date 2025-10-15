

#check if devtools is installed
if(!require(devtools)){
  install.packages("devtools")
}

#check if data.table is installed
if(!require(data.table)){
  install.packages("data.table")
}

# check if dplyr is installed
if(!require(dplyr)){
  install.packages("dplyr")
}

if(!require(data.table)){
  install.packages("data.table")
}

# install gdata
if(!require(gdata)){
  install.packages("gdata")
}

# install mgsub
if(!require(mgsub)){
  install.packages("mgsub")
}

#check if GenomicSEM is installed

if(!require(GenomicSEM)){
  install.packages("GenomicSEM", repos = NULL, type="source")
}

library(GenomicSEM) #Genomic SEM
library(data.table) #Data management and fast file reading/writing
library(dplyr)
library(data.table)

