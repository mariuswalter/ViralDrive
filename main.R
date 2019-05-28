library(ggpubr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(dplyr)
library(abind)
library(parallel)

theme_set(theme_classic()) + theme_replace(
  #plot.margin = unit(c(0.5,0.5,0,0.5), "cm"), 
  aspect.ratio=1,
  legend.position="none")        

N <-1000   #number of cells in a simulation
alpha <- 0.1 #efficiency of the drive: proportion of resistant virus created when coinfecting with a WT and Drive virus
Simul <- 50   #number of simulations to be run in parralel
numGen <- 50  #number of viral generations in a simulation
f <-0.9        #fitness cost of gene drive viruses
startfreq <- 0.1 #starting proportion of gene drive viruses
MOI=1

#Main viral gene drive function. All the other functions are just to plot/run simulations in parallele
GeneDrive <- function(N,numGen,alpha,f,startfreq,MOI=1){
  #  start_time <- Sys.time()
  Percent <- matrix(c(1-startfreq,startfreq,0,0), ncol = 4, nrow = 1) #at each generation, Percent is a matrix that give the percentage of the different viruses in the viral population. It is initated with a starting population of WT and GD viruses
  colnames(Percent)= c('WT','GD','GDo','R') #4 different virus types
  VirusNb <- matrix(0, ncol = 4, nrow = 1) #virus Nb give the actual number of viruses.
  colnames(VirusNb)= c('WT','GD','GDo','R')
  
  #each row represent a cell, which can be infected and coinfected by the different viruses. Each column give the number of the different viruses infected a cell.
  Cells_empty<- matrix(0, ncol = 4, nrow = N)
  colnames(Cells_empty)= c('WT','GD','GDo','R')
  
  InfectedCells <- Cells_empty #Matrix of the N cells, each of them is going to be randomly infected by viruses
  
  for (gen in 0:(numGen-1)){
    
    #create the N viruses that are going to infect the cells, starting from the observed proportion of virus 
    if (gen==0){ 
      Viruslist <- sample( colnames(Percent), MOI*N, replace=TRUE, prob=Percent)
    } else {
      Viruslist <- sample( colnames(Percent), MOI*N, replace=TRUE, prob=Percent[gen,] )
    }
    
    #infect the N cells with the MOI*N virus. The algorithm go through the virus list, sending randomly each virus to one cell, and updating the InfectedCells matrix
    InfectedCells <- Cells_empty
    for (virus in 1:(MOI*N)) {        
      ChooseCell  <- sample(1:N, 1)   #choose a random cell
      InfectedCells[ChooseCell,Viruslist[virus]]  <- InfectedCells[ChooseCell,Viruslist[virus]] + 1
    }
    
    #create the second generation of viruses (virus_sg), an average of 100 virus is created per cells. The algorithm goes through the matrix of cells, and compute the number and the type of new viruses produced by every cell
    for (Cell in 1:N) {
      NumberofVirus <- sum(InfectedCells[Cell,])  #total number of viruses that infected the cell
      C <- InfectedCells[Cell,]                   #table giving the number of of virus of the different type ('WT','GD','GDo','R')
      P <- prop.table(C)                          #Proportion of the different viruses that infected that specific cell, computed from C
      if (gen==0 & NumberofVirus >0){ 
        Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=100, sd=30)), replace=TRUE, prob=P ))
      }
      if (NumberofVirus >0 & gen>0){
        #if no GD virus, make new WT and R viruses with no fitness cost (average 100, sd 30), with the proportions given by P
        if ( (C['GD']+C['GDo'])==0 ){ 
          Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=100, sd=30)), replace=TRUE, prob=P ))
        }
        #if no WT or Res virus, make new viruses with a fitness cost of f
        if ( (C['WT']+C['R'])==0 ){ 
          Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=f*100, sd=f*30)), replace=TRUE, prob=P ))
        }
        #if  WT and GD (or GDo) viruses coinfected the cell, WT are converted to GDo and R. No fitness cost
        if ( (C['GD']+C['GDo']) > 0 & (C['WT']+C['R']) > 0){ 
          Virus_sg <- c(Virus_sg,sample( names(C), abs(rnorm(1, mean=100, sd=30)), replace=TRUE, 
                                         prob=  c(0,P['GD'],P['GDo']+(1-alpha)*P['WT'],P['R']+alpha*P['WT'])))
        }
      }
    }
    #calculate the new proportion of viruses
    if (gen > 0) {
      Percent <- rbind(Percent,0)
      VirusNb <- rbind(VirusNb,0)
    }
    Percent[gen+1,'WT'] = prop.table(table(Virus_sg))["WT"]
    Percent[gen+1,"GD"] <- prop.table(table(Virus_sg))["GD"]
    Percent[gen+1,'GDo'] <- prop.table(table(Virus_sg))["GDo"]
    Percent[gen+1,'R'] <- prop.table(table(Virus_sg))["R"]
    Percent[is.na(Percent)] <- 0
    VirusNb[gen+1,'WT'] <- table(Virus_sg)["WT"]
    VirusNb[gen+1,"GD"] <- table(Virus_sg)["GD"]
    VirusNb[gen+1,'GDo'] <- table(Virus_sg)["GDo"]
    VirusNb[gen+1,'R'] <- table(Virus_sg)["R"]
    VirusNb[is.na(VirusNb)] <- 0

  }
  
  #  end_time <- Sys.time()
  #  print(end_time-start_time)
  return(VirusNb)
}

#Plot one drive - Useless
PlotDrive <-function(VirusNb){
  VirusNbforplot <- data.frame(cbind(1:length(VirusNb[,1]),VirusNb))
  Percentforplot <- data.frame(cbind(1:numGen,prop.table(data.matrix(VirusNb),margin=1)))
  names(VirusNbforplot) <- c('gen','WT','GD','GDo','R')
  names(Percentforplot) <- c('gen','WT','GD','GDo','R')
  PlotV <- ggplot(VirusNbforplot, aes(x=gen))   + xlim(1,numGen)+
    geom_line(aes(y=GD),color='red')+
    geom_line(aes(y=GDo),color='orange') +
    geom_line(aes(y=WT),color='green') +
    geom_line(aes(y=R),color='black') 
  PlotP <- ggplot(Percentforplot, aes(x=gen))  + ylim(0,1) + xlim(0,numGen)+
    geom_ribbon(aes(ymin=GDo+WT+R,ymax=1),fill='red') +
    geom_ribbon(aes(ymin=WT+R,ymax=GDo+WT+R),fill='orange') +
    geom_ribbon(aes(ymin=R,ymax=WT+R),fill='green') +
    geom_ribbon(aes(ymin=0,ymax=R),fill='black') 
  
  ggarrange(PlotP, PlotV, ncol=2)
}

#Generate multiple drive simulations, return a list of dataframe, each one of them being the result of one simulation
GenerateMultipleDrive <-function(Simul=50,N=1000,numGen=50,alpha=0.1,f=0.9,startfreq=0.05,MOI=1){
  listdf= list()
  time <-0
  for (i in 1:Simul){
    start_time <- Sys.time()
    listdf[[i]] <- GeneDrive(N=N,numGen=numGen,alpha=alpha,f=f,startfreq=startfreq,MOI=MOI)
    end_time <- Sys.time()
    time <- time + (end_time-start_time)
    toaffich <- (time)
    print(toaffich)
    message(cat("for simulation",i))
  }
  names(listdf) <- rep(paste("Simulation=",Simul,"N=",N,", f=",f,", Startfreq=",startfreq,", alpha=",alpha,", MOI=",MOI))
  return(listdf)   
}  

#Plot multiple drives
PlotMultipleDrive <-function(listdf) {
  numGen <- length(listdf[[1]][,1])
  StoreNames <- names(listdf[1])
  listdf <- unname(listdf)
  listforplot <- lapply(listdf, function(x) cbind(x, gen=1:numGen))
  listforplot <- lapply(listforplot, function(x) data.frame(x))
  
  PlotAll <- ggplot(bind_rows(listforplot, .id="df"), aes(x=gen)) +
    geom_line(aes(y=GDo,group=df),alpha=0.5,color='orange',size=1) +
    geom_line(aes(y=WT,group=df),alpha=0.5,color='green',size=1) +
    geom_line(aes(y=R,group=df),alpha=0.5,color='darkgreen',size=1) +
    geom_line(aes(y=GD,group=df),alpha=0.5, color='red' ,size=1) +
    labs(y = "Titer",x= "Generation")  + ggtitle("All simulation")+ 
    scale_x_continuous(expand = c(0, 0), breaks = seq(0,numGen,10)) +
    scale_y_continuous(expand = c(0, 0))
  
  meanlist <- data.frame(apply(abind(listdf, along=3), c(1,2), mean))
  sdlist <- data.frame(apply(abind(listdf, along=3), c(1,2), sd))
  
  VirusNbforplot <- cbind(1:numGen,meanlist,sdlist)
  names(VirusNbforplot) <- c('gen','WT','GD','GDo','R','sdWT','sdGD','sdGDo','sdR')
  
  PlotSD <- ggplot(VirusNbforplot, aes(x=gen))   + 
    geom_ribbon(aes(ymin=(GD-sdGD),ymax=(GD+sdGD)),fill='red',alpha=0.3) +
    geom_ribbon(aes(ymin=(GDo-sdGDo),ymax=(GDo+sdGDo)),fill='orange',alpha=0.3) +
    geom_ribbon(aes(ymin=(WT-sdWT),ymax=(WT+sdWT)),fill='green',alpha=0.3) +
    geom_ribbon(aes(ymin=(R-sdR),ymax=(R+sdR)),fill='darkgreen',alpha=0.3) +
    geom_line(aes(y=GD),color='red',size=1)+
    geom_line(aes(y=GDo),color='orange',size=1) +
    geom_line(aes(y=WT),color='green',size=1) +
    geom_line(aes(y=R),color='darkgreen',size=1) +
    labs(y = "Titer",x= "Generation") + ggtitle("Mean + SD") +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0,numGen,10)) +
    scale_y_continuous(expand = c(0, 0))
  

  MeanPercentforplot <- data.frame(cbind(1:numGen,prop.table(data.matrix(meanlist),margin=1)))
  names(MeanPercentforplot) <- c('gen','WT','GD','GDo','R')

  PlotP <- ggplot(MeanPercentforplot, aes(x=gen))  + 
    geom_ribbon(aes(ymin=GDo+WT+R,ymax=1),fill='red') +
    geom_ribbon(aes(ymin=WT+R,ymax=GDo+WT+R),fill='orange') +
    geom_ribbon(aes(ymin=R,ymax=WT+R),fill='green') +
    geom_ribbon(aes(ymin=0,ymax=R),fill='darkgreen') +
    labs(y = "% Virus",x= "Generation") + ggtitle("Percentages") + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(0,numGen,10)) +
    scale_y_continuous(expand = c(0, 0))
  
  Fig <- ggarrange(PlotAll,PlotSD,PlotP, ncol=3,nrow=1)
  Fig <- annotate_figure(Fig,
                         #                 top = text_grob("Visualizing Tooth Growth", color = "red", face = "bold", size = 14),
                         bottom = text_grob(StoreNames, color = "blue",
                                            hjust = 1, x = 1, face = "bold", size = 10),
  )
  return(Fig)
}

#Generate a list of condition to be plotted
Conditions <-list(
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.9,startfreq=0.1,MOI=1),
  c(alpha=0.1,f=0.9,startfreq=0.01,MOI=1),
  c(alpha=0.1,f=0.9,startfreq=0.2,MOI=1),
  c(alpha=0.1,f=0.9,startfreq=0.5,MOI=1),
  
  c(alpha=0.1,f=0.99,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.95,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.8,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.7,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.6,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.5,startfreq=0.05,MOI=1),
  
  c(alpha=0,f=0.9,startfreq=0.05,MOI=1),
  c(alpha=0.001,f=0.9,startfreq=0.05,MOI=1),
  c(alpha=0.01,f=0.9,startfreq=0.05,MOI=1),
  c(alpha=0.2,f=0.9,startfreq=0.05,MOI=1),
  c(alpha=0.3,f=0.9,startfreq=0.05,MOI=1),
  
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=1),
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=0.8),
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=0.5),
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=0.4),
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=0.3),
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=0.2),
  c(alpha=0.1,f=0.9,startfreq=0.05,MOI=0.1)
)

#Run all the condition (slow)
for (X in Conditions){
  ggsave(paste0(paste(names(X),X,sep = "",collapse="_"), ".pdf"),
         PlotMultipleDrive(do.call(GenerateMultipleDriveCl, as.list(X))), 
         height =80 , width = 210, units = "mm")
}

#So you can do in parralel in a cluster, otherwise it takes forever
mclapply(X,
         FUN=function(X) ggsave(paste0(paste(names(X),X,sep = "",collapse="_"), ".pdf"),
                                PlotMultipleDrive(do.call(GenerateMultipleDrive, as.list(X))), 
                                height =80 , width = 210, units = "mm"),
         mc.cores = 30)




