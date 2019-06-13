#Same as GeneDrive in main, except that this time, Res virus grow with the  the same fitness cost as GD virus

GeneDriveRes <- function(N,numGen,alpha,f,startfreq,MOI=1){
  start_time <- Sys.time()
  Percentdf <- matrix(c(1-startfreq,startfreq,0,0), ncol = 4, nrow = 1)
  colnames(Percentdf)= c('WT','GD','GDo','R')
  VirusNb <- matrix(0, ncol = 4, nrow = 1)
  colnames(VirusNb)= c('WT','GD','GDo','R')
  
  #each row represent a cell, which can be infected by a number of the different viruses
  Celldf_empty<- matrix(0, ncol = 4, nrow = N)
  colnames(Celldf_empty)= c('WT','GD','GDo','R')
  
  InfectedCell <- Celldf_empty
  
  for (gen in 0:(numGen-1)){
    
    #create the Nviruses that are going to infect the cells, starting from the observed proportion of virus
    if (gen==0){ 
      Viruslist <- sample( colnames(Percentdf), MOI*N, replace=TRUE, prob=Percentdf)
    } else {
      Viruslist <- sample( colnames(Percentdf), MOI*N, replace=TRUE, prob=Percentdf[gen,] )
    }
    
    #infect the N cells with the MOI*N virus, sending randomly each virus to a cell
    InfectedCell <- Celldf_empty
    for (virus in 1:(MOI*N)) {
      ChooseCell  <- sample(1:N, 1)
      InfectedCell[ChooseCell,Viruslist[virus]]  <- InfectedCell[ChooseCell,Viruslist[virus]] + 1
    }
    
    #create the second generation of viruses
    Virus_sg <- c()
    for (Cell in 1:N) {
      NumberofVirus <- sum(InfectedCell[Cell,])
      C <- InfectedCell[Cell,]
      P <- prop.table(C)
      if (gen==0 & NumberofVirus >0){ 
        Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=100, sd=30)), replace=TRUE, prob=P ))
      }
      if (NumberofVirus >0 & gen>0){
        #if no GD or R virus, make new virus
        if ( (C['GD']+C['GDo']+C['R'])==0 ){ 
          Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=100, sd=30)), replace=TRUE, prob=P ))
        }
        #if no WT , make new virus with a fitness cost of f
        #This is where it differs from the main function
        if ( (C['WT'])==0 ){ 
          Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=f*100, sd=f*30)), replace=TRUE, prob=P ))
        }
        #if  WT and GD virus present,WT are converted to GDo and R
        if ( (C['GD']+C['GDo']) > 0 & (C['WT']) > 0){ 
          Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=100, sd=30)), replace=TRUE, 
                                         prob=  c(0,P['GD'],P['GDo']+(1-alpha)*P['WT'],P['R']+alpha*P['WT'])))
        }
      }
    }
    #calculate the new proportion of viruses
    if (gen > 0) {
      Percentdf <- rbind(Percentdf,0)
      VirusNb <- rbind(VirusNb,0)
    }
    Percentdf[gen+1,'WT'] = prop.table(table(Virus_sg))["WT"]
    Percentdf[gen+1,"GD"] <- prop.table(table(Virus_sg))["GD"]
    Percentdf[gen+1,'GDo'] <- prop.table(table(Virus_sg))["GDo"]
    Percentdf[gen+1,'R'] <- prop.table(table(Virus_sg))["R"]
    Percentdf[is.na(Percentdf)] <- 0
    VirusNb[gen+1,'WT'] <- table(Virus_sg)["WT"]
    VirusNb[gen+1,"GD"] <- table(Virus_sg)["GD"]
    VirusNb[gen+1,'GDo'] <- table(Virus_sg)["GDo"]
    VirusNb[gen+1,'R'] <- table(Virus_sg)["R"]
    VirusNb[is.na(VirusNb)] <- 0
    #table(Virus_sg)
    #print(Percentdf[gen+1,])
  }
  
  return(VirusNb)
}

GenerateMultipleDriveRes <-function(Simul=10,N=1000,numGen=50,alpha=0.1,f=0.9,startfreq=0.1,MOI=1){
  listdf= list()
  time <-0
  for (i in 1:Simul){
    start_time <- Sys.time()
    listdf[[i]] <- GeneDriveRes(N=N,numGen=numGen,alpha=alpha,f=f,startfreq=startfreq,MOI=MOI)
    end_time <- Sys.time()
    time <- time + as.numeric(end_time-start_time)
    toaffich <-as.numeric(time/60)
    message(cat("simulation",i,"took",toaffich,"minutes"))
  }
  
  return(listdf)   
}  

Xres <-list(
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=1)
)
mclapply(Xres,
         FUN=function(X) ggsave(paste0(paste(names(Xres),Xres,sep = "",collapse="_"), ".pdf"),
                                PlotMultipleDrive(do.call(GenerateMultipleDriveRes, as.list(Xres))), 
                                height =80 , width = 210, units = "mm"),
         mc.cores = 30)
