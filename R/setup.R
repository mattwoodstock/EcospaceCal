rm(list=ls())
library('dplyr'); library('reshape2'); library('tidyr'); library('stringr'); library('foreach')
library('PBSmodelling'); library('snowfall'); library('parallel'); library('snow'); 
library('doSNOW'); library('openxlsx'); library('tcltk')
library('plyr')
#source(paste0(getwd(),"/R_EwE_console_functions.R"))
source("./R_EwE_console_functions.R")

file.console <- "C:/Program Files/EcoSpace Console 1.4.3.1 64-bit/EwEClientConsole.exe" #You'll need to insert the location of your Command-line executable

file.cmdbase <- paste0("./GIST CommandFile_template.txt")
file.obsts <- paste0("./updated_TS.csv")
file.predprey <- paste0("./predprey_pairs.csv")
file.basevul <- paste0("./GIST Vulnerabilities.csv")

#Years to test
startyear    <- 1980
endyear      <- 2016
endyear_sens <- 2016
n_years_sens <- endyear_sens - startyear + 1; n_years_sens

#get base command file and predator/prey list--------------------------------------------
cmd_base = readLines(file.cmdbase)
predprey = read.csv(file.predprey)

#LOAD REFERENCE DATA--------------------------------------------------------------------------------
#read observed timeseries csv file input to Ecosim
obs.ts = fn.read_ecosim_timeseries(file.obsts)


#Load group names
vuls.base <- read.csv(file.basevul,row.names=1)
group.names <- gsub(" ","_",vuls.base[,1]); group.names
df.names <- data.frame(group.names = group.names, num = 1:length(group.names)) ## 2 column df with standardized name and pool code
write.csv(df.names, "./FG_names_standardized.csv", row.names = FALSE)




