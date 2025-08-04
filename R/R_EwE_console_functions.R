#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.runEwE-----------------------------------------------------------------------------------------
fn.runEwE <-  function(runlist,do.obj=1,i=1){
  #runlist=runlist_sens; i=3; do.obj=T
  #create command line and run Ecospace
  dir.cmdfile <- runlist$cmd_file[i]
  cmd = paste(paste0('"', file.console, '"'),
              paste0('"', dir.cmdfile, '"'))
  system(cmd,intern = F)
  
  #calculate objective function
  if(do.obj==1){
    objout <- fn.objfxn1(i,runlist)
    return(objout)
  }
  if(do.obj==2){
    objout <- fn.objfxn2(i,runlist)
    return(objout)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#OBJECTIVE FUNCTIONS---------------------------------------------------------------------------------  
#Loop to run objfxn for sensitivity and fill in runlist ////////////////////////
#print(paste("Sensitivity run obj. fxn.", start, "to", stop))
##Objective fxn 1-----
#fits to ecosim timeseries data
fn.objfxn1 <- function(i,runlist){
  #i = 1
  #runlist = readRDS(file="C:\\SREM\\SREM Console\\multistanza fit\\out\\Runlists\\i15 - Step 2 - Alternate test value runs.rds")
  #i=2; runlist = runlist_sens
  #runlist$dir.out[i] <- paste0(dir.ewe,"/output/sp01_s88 base")
  #get Ecospace predicted annual average biomass for this run
  fnm.pred = paste0(runlist$dir.out[i], "/Ecospace_Annual_Average_Biomass.csv")  ## Read directory from run list
  runname = basename(runlist$dir.out[i])
  metadat = read.csv(fnm.pred, row.names = NULL,header=F, blank.lines.skip = F)
  metadat = metadat[1:(which(substr(metadat$V1,1,4)=='Year')-1),]
  run.timestamp = metadat[which(metadat[,1]=="Date"),2] ## Extract date and time metadata for when the console ran
  #tmp.predB = read.table(fnm.pred,header=F,sep=":")
  
  idx.skip = nrow(metadat) #which(metadat[,1]=="<HEADER end/>")+1
  predB = read.csv(fnm.pred, skip=idx.skip, header=T,check.names=F)
  predB$Year = startyear:endyear 
  rownames(predB) = predB$Year
  predB$Year = NULL
  names(predB) = group.names
  #predB = predB[1:nrow(obsB), ]  #remove extra years in ecospace run
  #get ecospace predicted annual catch for this run
  fnm.predC = paste0(runlist$dir.out[i], "/Ecospace_Annual_Average_Catch.csv")  ## Read directory from run list
  predC = read.csv(fnm.predC,skip=idx.skip,header=T,check.names = F)
  predC = predC[,-1]
  idx.nm = unlist(gregexpr('\\|', names(predC)))+1
  predC.grps = substr(names(predC),idx.nm,100)
  names(predC) = predC.grps
  predC.grpnum = as.numeric(df.names$num[match(predC.grps,gsub("/","_",df.names$group.names))])
  predC2 = t(rowsum(t(predC), predC.grps))
  obs.ts.head = obs.ts$obsB.head
  lk.ts.biomass <- cbind(obs.ts.head, dattype='timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=48
    #grp = group.names[j]
    grpnum = obs.ts.head$pool_code[j]
    grp = group.names[grpnum]
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1e-6
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    #if observed timeseries is biomass
    this_biom = data.frame(Year = startyear:endyear_sens ,pred = predB[,grpnum])
    this_obs = data.frame(Year = startyear:endyear_sens ,obs = predB[1:n_years_sens,grpnum])
    lk.dat = merge(this_biom,this_obs,by=c("Year"))[,2:3]
    rownames(lk.dat) <- rownames(predB)
    if(obs.ts.head$type[j]==0){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$obs = lk.dat$obs*q
    lk.dat$ll = 0.5*(log(lk.dat$pred/lk.dat$obs))^2/obs.sd^2
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.biomass$obs.cv[j] = obs.cv
    lk.ts.biomass$obs.sd[j] = obs.sd
    lk.ts.biomass$loglik[j] = lk.sum
    
  }
  
  obs.ts.head = obs.ts$obsC.head
  lk.ts.catch <- cbind(obs.ts.head, dattype='timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=48
    #grp = group.names[j]
    grpnum = obs.ts.head$pool_code[j]
    grp = gsub("_", " ", group.names[grpnum])
    idx_group = match(grp,colnames(predC2))
    #grp.wt = 1
    obs.wt = obs.ts.head$weight[j]
    if (obs.wt == 0){
      obs.wt = 1e-6
    }
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    this_catch = data.frame(Year = startyear:endyear_sens ,pred = predC2[,idx_group])
    this_obs = data.frame(Year = startyear:endyear_sens ,obs = obs.ts$obsC[1:n_years_sens,j])
    lk.dat = merge(this_catch,this_obs,by=c("Year"))[,2:3]
    
    lk.dat$obs[lk.dat$obs <=0] <- NA
    rownames(lk.dat) <- rownames(predC2)
    if(obs.ts.head$type[j]==61){
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    } else {
      q = 1
    }
    lk.dat$obs = lk.dat$obs*q
    lk.dat$ll = 0.5*(log(lk.dat$pred/lk.dat$obs))^2/obs.sd^2
    lk.sum = sum(lk.dat$ll, na.rm = T)
    #lk.sum = lk.sum * obs.wt
    lk.ts.catch$obs.cv[j] = obs.cv
    lk.ts.catch$obs.sd[j] = obs.sd
    lk.ts.catch$loglik[j] = lk.sum
  } 
  
  lk.ts = rbind(lk.ts.biomass,lk.ts.catch)
  
  ## Make output vector................................
  ##
  lk.vec <- rep(0,length(group.names))
  lk.agg <- aggregate(loglik~pool_code, data=lk.ts, sum,na.rm=T)
  lk.vec[lk.agg$pool_code] <- lk.agg$loglik
  
  outvec = round(c(sum(lk.vec), lk.vec),2)
  write.csv(outvec,paste0(runlist$dir.out[i],"/objvals.csv"))
  return(outvec)
}
##Objective fxn 2-----
#includes biomass and catch (scaled or unscaled) timeseries, plus avg cpue rasters
  
#Loop to run objfxn for sensitivity and fill in runlist ////////////////////////
#print(paste("Sensitivity run obj. fxn.", start, "to", stop))
fn.objfxn2 <- function(i,runlist){
  #runlist = readRDS(file="C:\\SREM\\SREM Console\\multistanza fit\\out\\Runlists\\i15 - Step 2 - Alternate test value runs.rds")
  #i=2; runlist = runlist_sens
  #runlist$dir.out[i] <- paste0(dir.ewe,"/output/sp01_s88 base")
  #get Ecospace predicted annual average biomass for this run
  fnm.pred = paste0(runlist$dir.out[i], "/Ecospace_Annual_Average_Biomass.csv")  ## Read directory from run list
  runname = basename(runlist$dir.out[i])
  metadat = read.csv(fnm.pred, row.names = NULL,header=F, blank.lines.skip = F)
  metadat = metadat[1:(which(substr(metadat$V1,1,4)=='Year')-1),]
  run.timestamp = metadat[which(metadat[,1]=="Date"),2] ## Extract date and time metadata for when the console ran
  #tmp.predB = read.table(fnm.pred,header=F,sep=":")
  idx.skip = nrow(metadat) #which(metadat[,1]=="<HEADER end/>")+1
  predB = read.csv(fnm.pred, skip=idx.skip, header=T)
  predB$Year = startyear:endyear_sens ## Changed to start 1997 for SREM
  rownames(predB) = predB$Year
  predB$Year = NULL
  names(predB) = group.names
  #predB = predB[1:nrow(obsB), ]  #remove extra years in ecospace run
  
  #get ecospace predicted annual catch for this run - need to add a character flag in Ecopath fleet names to be able to parse out groups 
  fnm.predC = paste0(runlist$dir.out[i], "/Ecospace_Annual_Average_Catch.csv")  ## Read directory from run list
  predC = read.csv(fnm.predC,skip=idx.skip,header=T)
  predC = predC[,-1]
  #idx.nm = unlist(gregexpr('\\.', names(predC)))+4
  #predC.grps = substr(names(predC),idx.nm,100)
  #predC.grps = gsub("1_","1+",gsub("2_","2+",gsub("4_","4+",gsub("1_4","1-4",gsub("2_5_","2.5+",gsub("1_2_5","1-2.5",gsub("\\.","_",predC.grps)))))))
  #names(predC) = predC.grps
  #predC.grpnum = as.numeric(df.names$num[match(predC.grps,gsub("/","_",df.names$group.names))])
  #predC2 = t(rowsum(t(predC), predC.grps))
  

  lk.ts <- cbind(obs.ts.head, dattype='timeseries',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.ts.head)){
    #j=48
    #grp = group.names[j]
    grpnum = obs.ts.head$Pool_code[j]
    grp = group.names[grpnum]
    #grp.wt = 1
    obs.wt = obs.ts.head$Weight[j]
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    #if observed timeseries is biomass
    if(obs.ts.head$Type[j]==0){
      lk.dat = data.frame(pred=predB[,grpnum], obs=obs.ts[,j])
      rownames(lk.dat) <- rownames(predB)
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
      lk.dat$obs = lk.dat$obs*q
      matplot(1985:2022,lk.dat,type='l')
      lk.dat$ll = 0.5*(log(lk.dat$pred/lk.dat$obs))^2/obs.sd^2
      lk.sum = sum(lk.dat$ll, na.rm = T)
      lk.sum
      #lk.sum = lk.sum * obs.wt
      lk.ts$obs.cv[j] = obs.cv
      lk.ts$obs.sd[j] = obs.sd
      lk.ts$loglik[j] = lk.sum
    }
   
    # #if observed timeseries is catch
    # if(obs.ts.head$Type[j]==6){
    #   idx.pred = which(colnames(predC2)==grp)
    #   idx.obs = which(obsC.head$Pool_code==j)
    #   
    #   grpC.cv = 1/ifelse(is.na(idx.obs),1,obsC.head$Weight[idx.obs])
    #   grpC.sd = sqrt(log(1+(grpC.cv/grp.wt)^2))
    #   
    #   #get pred and observed
    #   lk.dat = data.frame(obs=obsC[,idx.obs], predC2[,idx.pred])
    #   lk.dat$obs[lk.dat$obs==0] = NA
    #   lk.dat = lk.dat[complete.cases(lk.dat),] ## Discard NAs
    #   
    #   #rescale obsB similar to Ecosim
    #   #q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
    #   #lk.dat$obs = lk.dat$obs*q
    #   
    #   #calculate weighted log-likelihood
    #   lk.dat$ll = 0.5*(log(lk.dat$pred/lk.dat$obs))^2/grpC.sd^2
    #   lk.sum = sum(lk.dat$ll, na.rm = T)
    #   lk.sum = lk.sum * obsC.head$Weight[idx.obs]
    #   lk.C = lk.sum
    # }  
  } #end obs timeseries loop 

  #map data - single average maps (no time component)
  #get list of ascii output files from Ecospace run
  files.asc1 <- list.files(paste0(runlist$dir.out[i],'/asc'), pattern=".asc$",full.names=T)
  if(length(files.asc1)==0)  files.asc1 <- list.files(paste0(runlist$dir.out[i]), pattern=".asc$",full.names=T)
  files.asc1 <- gsub("red snapper 1-2","red snapper 1", files.asc1)
  
 
  files.asc <- files.asc1[grepl('EcospaceMapBiomass-red grouper',files.asc1)]  #/////// for testing only, remove when ready to fit all maps
  files.asc.splt <- as.data.frame(matrix(unlist(strsplit(basename(files.asc),"-")),nrow=length(files.asc),ncol=3,byrow=T))
  files.asc.del <- files.asc1[!files.asc1%in%files.asc] 
  unlink(gsub("red snapper 1","red snapper 1-2", files.asc.del))
  files.keep <- character()
  
  dims.asc <- dim(raster(files.asc[1]))
  
  lk.map <- cbind(obs.map.head, dattype='raster',obs.cv=NA, obs.sd=NA, loglik=NA)
  for(j in 1:nrow(obs.map.head)){
    obs.map.head
    #j=7
    grpnum = obs.map.head$Pool_code[j]
    grp = group.names[grpnum]
    print(paste(obs.map.head$datasource[j],grp));flush.console()
    obs.wt = obs.map.head$Weight[j]
    obs.cv = 1/obs.wt#,obsB.head$Weight[which(obsB.head$Pool_code==j)])
    obs.sd = sqrt(log(1+(obs.cv/obs.wt)^2))
    
    #get predicted biomass maps and averaged over all years for functional group j
    files.j <- files.asc[which(files.asc.splt$V2==gsub("_"," ",grp) & files.asc.splt$V1=='EcospaceMapBiomass')]
    files.keep <- c(files.keep,files.j)
    
    if(length(files.j)>0){
      #stack.j <- stack(files.j)
      #mean.j <- calc(stack.j,mean)

      #converting to array is faster than raster::calc
      vals.j <- unlist(sapply(files.j,FUN=function(x) read.table(file=x, skip=6),simplify=T))
      vals.j[vals.j<0] <- NA
      arr.j <- array(vals.j,dim=c(dims.asc[1:2],length(files.j)))
      
      mean.j <- raster(files.j[1])
      mean.j[] <- NA
      meanarr.j <- apply(arr.j,c(1,2),mean,na.rm=T)
      mean.j[] <- meanarr.j

      obs.j <- obs.map[[j]]
      obs.j <- resample(obs.j, mean.j)
      # mask.j <- mask(mean.j,obs.j, updatevalue=NA, maskvalue=NA)
      # par(mfrow=c(2,2))
      #  plot(mean.j,colNA='black'); plot(obs.j,colNA='black')
      #  plot(mask.j, colNA='black')
      lk.dat = data.frame(pred=getValues(mean.j), obs=getValues(obs.j))
      lk.dat = lk.dat[complete.cases(lk.dat),]
      q = mean(lk.dat$pred,na.rm=T)/mean(lk.dat$obs,na.rm=T)
      lk.dat$obs = lk.dat$obs*q
  
      #matplot(1985:2022,lk.dat,type='l')
      lk.dat$ll = 0.5*(log((lk.dat$pred+1)/(lk.dat$obs+1))^2)/obs.sd^2
      lk.sum = sum(lk.dat$ll, na.rm = T)
      lk.sum
      #lk.sum = lk.sum * obs.wt
      lk.map$obs.cv[j] = obs.cv
      lk.map$obs.sd[j] = obs.sd
      lk.map$loglik[j] = lk.sum
      rm(mean.j,obs.j,mask.j,lk.dat)
    }
  }
  
  files.mapcsv <- list.files(paste0(runlist$dir.out[i]), pattern=".csv$",full.names=T)
  files.mapcsv <- files.mapcsv[which(substr(basename(files.mapcsv),1,11)=='EcospaceMap')]
  unlink(files.mapcsv)
  
  #files.del <- files.asc[which(!files.asc%in%files.keep)]
  
  #make lk output table
  lk.map2 = lk.map
  lk.map2$indexname = paste0(gsub(" ","_",lk.map2$datasource),lk.map2$layername)
  lk.map2$Type=0
  lk.out = rbind(lk.ts,lk.map2[,-c(1,2)])
  
  write.csv(lk.out,file=paste0(runlist$dir.out[i],"/objvals_",runname,".csv"), row.names=F)
  
  ## Make output vector................................
  ##
  lk.vec <- rep(0,length(group.names))
  lk.agg <- aggregate(loglik~Pool_code, data=lk.out, sum,na.rm=T)
  lk.vec[lk.agg$Pool_code] <- lk.agg$loglik
  
  outvec = round(c(sum(lk.vec), lk.vec),2)

  #write.csv(outvec,paste0(runlist$dir.out,"/objvals.csv"))
  return(outvec)
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#fn.longvuls ---------------------------------------------------------------------------------------
fn.longvuls <- function(vuls=vuls){
  
  longvuls <- reshape2::melt(vuls, id.vars <- 1) ## Melt to long
  names(longvuls)[2:3] <- c('pred','baseval')  ## Rename
  longvuls$prey <- as.factor(longvuls$prey)    
  longvuls$grp.pred <- as.integer(longvuls$pred)  ## Add group number for pred. Note this doesn't work for prey
  longvuls <- merge(longvuls, df.names, by.x="prey", by.y ="group.names") ## Merge in prey group number from df.names
  names(longvuls)[which(names(longvuls)=='num')] = 'grp.prey'
  #longvuls <- rename(longvuls, c(grp.prey = "num"))
  longvuls <- longvuls[!is.na(longvuls$baseval),]
  return(longvuls)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#fn.test_pred_col_vuls------------------------------------------------------------------------------
fn.pred_vuls_sensitivity_tests <- function(testvals=c(vul.min,vul.max)){
  #testvals <- c(vul.min,vul.max)
  basevals <- aggregate(baseval~grp.pred, data=longvuls, mean)
  basevals <- rep(basevals,each=2)
  predvuls <- data.frame(prey=NA,
                         pred=rep(group.names[1:n_consumers],each=length(testvals)),
                         baseval=rep(basevals$baseval, each=2),
                         grp.pred=rep(1:n_consumers,each=length(testvals)),
                         grp.prey=NA)
  predvuls$testval <- rep(testvals,n_consumers)
  predvuls$test.par <- 'pred_col'
  #param = "<ECOSIM_VULNERABILITIES_BY_PRED>"
  param <- "<ECOSIM_VULNERABILITIES_INDEXED>"
  predvuls$tag <- paste0(param,"(",predvuls$grp.pred,"),",predvuls$testval,",Indexed.Single[]")
  return(predvuls)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#fn.make_cmd_files----------------------------------------------------------------------------------
fn.make_cmd_files = function(runlist,iter,nyrs=n_years_sens){
  #runlist = runlist_sens; iter=1
  for(i in 1:nrow(runlist)){
    #i=2
    cmd_i = cmd_base
    
    #set output directory
    param = "<ECOSPACE_OUTPUT_DIR>"
    update = paste0(runlist$dir.out[i])
    n = which(substr(cmd_i[,1],1,nchar(param))==param)
    cmd_i[n,1]= paste(param, update, "System.String, Updated", sep = ", ")
    
    #set run length
    param = "<N_ECOSPACE_YEARS>"
    update = nyrs
    n = which(substr(cmd_i[,1],1,nchar(param))==param)
    cmd_i[n,1]= paste(param, update, "System.Int32, Updated", sep = ", ")
    
    #add taglines
    cmd_i = rbind(cmd_i,runlist$tag[i])
    
    
    #save command file to run folder
    write.table(cmd_i, runlist$cmd_file[i],  row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fn.make_cmd_files_old = function(runlist,iter){
    #runlist = runlist_atv.k
    for(i in 1:nrow(runlist)){
      #i=1
      cmd_i = cmd_template
      #set output directory
      param = "<ECOSPACE_OUTPUT_DIR>"
      update = paste0(runlist$dir.out[i])
      n = which(substr(cmd_i[,1],1,nchar(param))==param)
      cmd_i[n,1]= paste(param, update, "System.String, Updated", sep = ", ")
      
      if(iter==1){
        #add tagline to end
        cmd_i = rbind(cmd_i,runlist$tag[i])
      } else {
        param = unlist(strsplit(runlist$tag[i],split=","))[1]
        n = which(substr(cmd_i[,1],1,nchar(param))==param)
        
        #if tagline already exists, add new testvals
        if(length(n)>0){
          param_split = trimws(unlist(strsplit(cmd_i[n,1],split=','))[2])
          if(param=="<ECOSIM_VULNERABILITIES_BY_PREY_PRED>"){
            update.tag = paste0(param,", ",paste(param_split,runlist$grp.prey[i],runlist$grp.pred[i],runlist$testval[i]),", System.Single[], updated")
          }
          if(param=="<ECOSIM_VULNERABILITIES_BY_PRED>"){
            update.tag = paste0(param,", ",paste(param_split,runlist$grp.pred[i],runlist$testval[i]),", System.Single[], updated")
          }
          cmd_i[n,1] = update.tag
          
        } else{
          #if tagline is not in command file, append it to bottom 
          cmd_i = rbind(cmd_i,runlist$tag[i])
        }
      }
      #save command file in input folder
      write.table(cmd_i,paste0(runlist$dir.out[i],"/",runlist$cmd_file[i]),  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.read ecosim timeseries---------------------------------------------------------------------------------------
fn.read_ecosim_timeseries = function(filename){
  fnm.obs_ts = filename
  obs.ts.head = as.data.frame(t(read.csv(fnm.obs_ts,header=F,nrows=4)))
  names(obs.ts.head) = obs.ts.head[1,]; obs.ts.head = obs.ts.head[-1,]
  obs.ts.head[,2:4] = as.numeric(as.matrix(obs.ts.head[,2:4]))
  obs.ts = read.csv(fnm.obs_ts,header=F,skip=4)
  rownames(obs.ts) = obs.ts[,1]; obs.ts[,1] = NULL
  if(nrow(obs.ts.head) != ncol(obs.ts)) print('HEADER AND TIMESERIES DIMENSION DO NOT MATCH!!!')
  
  obsB.head = obs.ts.head[obs.ts.head[,4] %in% c(0,1),]
  obsB = obs.ts[,which(obs.ts.head[,4] %in% c(0,1))]
  obsC.head = obs.ts.head[obs.ts.head[,4] %in% c(6,61,-6),]
  obsC = obs.ts[,which(obs.ts.head[,4] %in% c(6,61,-6))]
  names(obsB.head) = names(obsC.head) = gsub(" ","_",names(obsB.head))
  ret = list(obsB.head,obsB,obsC.head,obsC)
  names(ret) = c('obsB.head','obsB','obsC.head','obsC')
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
#fn.pred_prey_testvals-----
fn.pred_prey_testvals = function(longvuls, maxvul=10000, minvul=1.01){
  
  #maxvul = 10000
  #minvul = 1.01
  longvuls$testval = maxvul
  longvuls$dir.test = "Higher"
  
  longvuls2 = longvuls 
  longvuls2$testval = minvul
  longvuls2$dir.test = "Lower"
  ## Combine lower testvalues back into vuls data frame; Constrain values
  longvuls = rbind(longvuls, longvuls2); rm(longvuls2)
  longvuls = subset(longvuls, longvuls$baseval != longvuls$testval)
  
  
  #runlist$prey = 5
  #runlist$pred = 3
  #runlist$vul = c(1.01,1e6,round(runif(nruns-2,min=1.01,max=1000),2))
  longvuls$tag = paste0("<ECOSIM_VULNERABILITIES_INDEXED>(",longvuls$grp.pred," ",longvuls$grp.prey,"), ",longvuls$testval,", Indexed.Single")
  longvuls = longvuls[order(longvuls$grp.pred, longvuls$grp.prey, longvuls$testval), ]
  return(longvuls)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Ecospace output to arrays---------------------------------------------------------------------------------  
fn.ecospace_predB_ts2array = function(dir.out=dir.out, timestep='annual'){
  if(timestep=='annual'){
    files.bio = list.files(dir.out,pattern="Ecospace_Annual_Average_Biomass.csv",recursive=T,full.names = T)
    nskip = which(substr(readLines(files.bio[1]),1,4)=='Year')-1
    bio = lapply(files.bio,read.csv,as.is=T,skip=nskip, row.names=1)
    bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],length(files.bio)),
                      dimnames=list(rownames(bio[[1]]),names(bio[[1]]),basename(dirname(files.bio))))
    for(r in 1:length(bio)){
      tmp = as.matrix(bio[[r]])
      bio.array[,,r] <- tmp
    }
  }
  
  if(timestep=='monthly'){
    files.bio = list.files(dir.out,pattern="Ecospace_Average_Biomass.csv",recursive=T,full.names = T)
    nskip = which(substr(readLines(files.bio[1]),1,8)=='TimeStep')-1
    bio = lapply(files.bio,read.csv,as.is=T,skip=nskip, row.names=1)
    bio.array = array(dim=c(dim(bio[[1]])[1],dim(bio[[1]])[2],length(files.bio)),
                      dimnames=list(rownames(bio[[1]]),names(bio[[1]]),basename(dirname(files.bio))))
    for(r in 1:length(bio)){
      tmp = as.matrix(bio[[r]])
      bio.array[,,r] <- tmp
    }
  }
  
  return(bio.array)
}

fn.ecospace_predC_ts2array = function(dir.out=dir.out, timestep='annual'){
  if(timestep=='annual'){
    files.cat = list.files(dir.out,pattern="Ecospace_Annual_Average_Catch.csv",recursive=T,full.names = T)
    nskip = which(substr(readLines(files.cat[1]),1,4)=='Year')-1
    cat = lapply(files.cat,read.csv,as.is=T,skip=nskip, check.names=F, row.names=1)
    names.split = strsplit(names(cat[[1]]),split="\\|")
    cat.grp = sapply(names.split, function(x) x[2])
    cat = lapply(cat,FUN=function(x) as.data.frame(t(rowsum(t(x), group=cat.grp))))
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],length(files.cat)),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),basename(dirname(files.cat))))
    for(r in 1:length(cat)){
      tmp = as.matrix(cat[[r]])
      cat.array[,,r] <- tmp
    }
  }
  
  if(timestep=='monthly'){
    files.cat = list.files(dir.out,pattern="Ecospace_Average_Catch.csv",recursive=T,full.names = T)
    nskip = which(substr(readLines(files.cat[1]),1,8)=='TimeStep')-1
    cat = lapply(files.cat,read.csv,as.is=T,skip=nskip, check.names=F, row.names=1)
    names.split = strsplit(names(cat[[1]]),split="\\|")
    cat.grp = sapply(names.split, function(x) x[2])
    cat = lapply(cat,FUN=function(x) as.data.frame(t(rowsum(t(x), group=cat.grp))))
    cat.array = array(dim=c(dim(cat[[1]])[1],dim(cat[[1]])[2],length(files.cat)),
                      dimnames=list(rownames(cat[[1]]),names(cat[[1]]),basename(dirname(files.cat))))
    for(r in 1:length(cat)){
      tmp = as.matrix(cat[[r]])
      cat.array[,,r] <- tmp
    }
  }
  
  return(cat.array)
}


