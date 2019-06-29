#' Code to fit Dirichelet
#'  
#' @author Michelle Masi & Joe Tarnecky
#' @date 2014
#' PDF distributions were calculated in FitDietDirichelet7.r
#' Modified by Hem Nalini Morzaria Luna January 2017, modified to save PDF
#' Also added paralellization to speed up code


#-------------------------------INITIALIZATION-----------------------------------

# List of packages for session

install.packages(c("XMLSchema", "SSOAP"), repos = c("http://packages.ropensci.org", "http://cran.rstudio.com"))

.packages = c("data.table","dplyr","ggplot2","stringi","tidyr","stringr",
              "R.utils","data.table","ggExtra", "ggthemes","Hmisc","VGAM",
              "future","parallel","doSNOW","magrittr")


# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#' clean up the space

rm(list=ls())
graphics.off()

#workpath <- "E:/Archivos/1Archivos/Articulos/En preparacion/GOM/Analysis"
#savepath <- "E:/Archivos/1Archivos/Articulos/En preparacion/GOM/Analysis/Beta"

#this should be changed to a google drive folder eventually because cam doesnt use dropbox anymore
workpath <- "~/Dropbox/GOM_Diets"
savepath <- "~/Dropbox/GOM_Diets"

setwd(workpath)

#read data 
x <- fread("Final_CSV_For_Dirchelet.csv") %>% mutate(Predator=as.factor(Predator))

npredspp  <- nlevels(x$Predator)
nprey     <- ncol(x)-2

predator.names <- x %>% tbl_df %>% select(Predator) %>% distinct() %>% .$Predator %>% as.character()


options(digits=20)

whatplot=2    # tells which plot to use
#NOTE: We want to calculate all graphs at the same time, otherwise the bootstrap wont match
#whatplot = 2 RedSnapperMarginalBetaA and B and also figs 1-4 

prey.names <- fread("prey_names.csv", header=TRUE) %>% select(prey) %>% 
  filter(prey!="Sediment bacteria") %>% 
  filter(prey!="Sediment bacteria") %>% 
  filter(prey!="Carrion detritus sediment") %>% 
  filter(prey!="Labile detritus sediment") %>% 
  filter(prey!="Refractory detritus sediment")

x <- x %>% as.data.frame()

#dirichlet fit will fail if you feed it 0 or 1, so squish the values a bit on both sides
dat2 <- cbind(x[,1:2],(0.999*(0.00001 + x[,3:(nprey+2)]/rowSums(x[,3:(nprey+2)]))))

#checks to see if number of prey remain after normalization
# dat2 %>% filter(Predator == "Benthic Feeding Sharks") %>% 
#   select(-Predator, -wgt) %>% 
#   summarise_each(funs(sum)) %>% 
#   Filter(function(x) any(x!=0.00039960000000000006),.)

#set up matrices to hold results
mat_names   <- list(c(colnames(dat2[,3:(nprey+2)])),c(levels(dat2$Predator)))
raw_wt_avg  <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
mle         <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
betavar     <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
lower95     <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
upper95     <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
mode_       <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
beta_a      <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
beta_b      <- matrix(nrow=nprey,ncol=npredspp,dimnames=mat_names)
simpleaverage = matrix(nrow=npredspp,ncol=nprey)
the_mode = matrix(nrow=nprey,ncol=npredspp)
bootaverage = matrix(nrow=nprey,ncol=npredspp)
bootmode = matrix(nrow=nprey,ncol=npredspp)



#-------------------------------BOOTSTRAPPING-----------------------------------
#for(i in 1:nlevels(dat2$Predator)){

for(i in 1:nlevels(x$Predator)){         
  
  this.predator <- predator.names[i]
  print(paste("This predator",this.predator))
  dat3 <- subset(dat2[dat2$Predator==levels(dat2$Predator)[i],])
  simpleaverage[i,]=colMeans(dat3[3:(nprey+2)])     #verified correct
  
  #checks to see if number of prey remain after normalization
    dat3 %>% 
     select(-Predator, -wgt) %>% 
     summarise_each(funs(sum)) %>% 
      Filter(function(x) any(x!=0.00039960000000000006),.)
  
  
  nstomachs   <- nrow(dat3)
  nfishavg    <- 50
  nbootstraps <- 1000
  stor        <- matrix(nrow=nbootstraps,ncol=nprey)
  tempstor <- matrix(nrow=nbootstraps,ncol=nprey)
  
  for(j in 1:nbootstraps){
    draw      <- dat3[sample(nstomachs,nfishavg,replace=T),]
    
    #-----
    tempstor[j,] <- apply(draw[,3:(nprey+2)],2,weighted.mean,w=draw$wgt)
    stor[j,] = tempstor[j,] / sum(tempstor[j,])
    
  }
  
  for (k in 1:nprey){
    bootaverage[k,i]=mean(stor[,k])
    bootmode[k,i]=as.numeric(names(sort(-table(stor[,k])))[1])
  }
  #---------------------------FIT TO DIRICHLET DIST-------------------------------
  #fit to distribution
  #stor[stor>1]=0.999999999
  #stor[stor<=0]=0.000000001
  
  print(i)
  print("trying to fit")
  fit <- vglm(stor~1,dirichlet)
  print("fit worked")
  k   <- Coef(fit)
  
  #obtain marginal beta distribution parameters
  betaEst <- matrix(nrow=nprey,ncol=2)
  for(m in 1:nprey){
    betaEst[m,1] <- k[m]
    betaEst[m,2] <- sum(k)-k[m]
  }
  #find pdf of beta for each prey item looping over a range of diet proportions
  
  
  prop = exp(seq(log(0.0000000001), log(0.999), length.out = 10000)) 
  #prop=c(seq(0,0.009999, 0.000001),seq(0.01,0.0999,0.0001),seq(0.1,1,0.001))
  #prop=c(seq(0,0.09999,0.00001),seq(0.1,1,0.001))
  #prop=seq(0,1,0.0001)
  
  bdist   <- data.frame(matrix(nrow=length(prop),ncol=nprey))
  
  print("trying to populate bdist")
  for(n in 1:length(prop)){      
    bdist[n,] <- dbeta(prop[n],betaEst[,1],betaEst[,2])
  }
  print("population successful")
  
  prey.names2 <- prey.names %>% .$prey %>% as.character()
  
  pdf_dist_maxval <- bdist %>% tbl_df %>% 
    setNames(prey.names2) %>% 
      summarise_each(funs(max)) %>% 
       max(.,na.rm = TRUE) %>% 
    floor(.)
    
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  pdf_dist <- bdist %>% tbl_df %>% 
    setNames(prey.names2) %>% 
    summarise_each(funs(max)) %>% 
    Filter(function(x) any(x<pdf_dist_maxval),.) %>% 
    names()
    
  filter_pdf <- . %>% tbl_df %>% 
    setNames(prey.names2) %>% 
    select(one_of(pdf_dist)) %>% 
    mutate(predator=this.predator) %>% 
    select(predator, everything()) %>% 
    group_by(predator)
   
  pdf_dist_sum <- bdist %>% filter_pdf %>% 
    mutate_each(funs(sum)) 
  
  pdf_dist_mode <- bdist %>% filter_pdf %>% 
    summarise_each(funs(getmode)) %>% 
    gather(`prey name`,mode,2:ncol(.))
  
  # PDF cumulative value 
  pdf_dist_cumsum <- bdist %>% filter_pdf %>% 
    mutate_each(funs(cumsum))
  
  pdf_dist_cont <- pdf_dist_cumsum
  
  for(eachcol in 2:ncol(pdf_dist_sum)){
    
    cont.value <- pdf_dist_cumsum[,eachcol] / pdf_dist_sum[,eachcol]
  
    pdf_dist_cont[,eachcol] <- pdf_dist_cumsum[,eachcol] / pdf_dist_sum[,eachcol]
  }
  
  pdf_dist_fin <-   pdf_dist_cont %>% 
    gather(`prey name`,availability,2:ncol(.))%>% 
    inner_join(pdf_dist_mode, by=c("predator","prey name"))
  
  setwd(savepath)
  write.csv(pdf_dist_fin,paste(this.predator,"_PDF.csv",sep=""),row.names = FALSE)
  
  #-----------------------COMPUTE AND STORE THE STATISTICS----------------------------------
  #mean of raw data
  print("raw_wt_avg[,i] <- apply(dat3[,3:(nprey+2)],2,weighted.mean,w=dat3$wgt)")
  raw_wt_avg[,i] <- apply(dat3[,3:(nprey+2)],2,weighted.mean,w=dat3$wgt)
  #the_mode[,i] = apply(dat3[,3:(nprey+2)],2,mode)  
  
  
  #maximum likelihood estimate
  print("Y <- apply(bdist,2,max)")
  Y <- apply(bdist,2,max)
  print("I <- apply(bdist,2,which.max)")  
  I <- apply(bdist,2,which.max)
  
  print("maxes <- prop[c(I)]")  
  maxes <- prop[c(I)]
  
  print("mle[,i] <- maxes")
  mle[,i] <- maxes         
  #variance
  print("betavar[,i] <- betaEst[,1]*betaEst[,2]/((betaEst[,1]+betaEst[,2])^2 * (betaEst[,1]+betaEst[,2]+1))")
  betavar[,i] <- betaEst[,1]*betaEst[,2]/((betaEst[,1]+betaEst[,2])^2 * (betaEst[,1]+betaEst[,2]+1))
  #95% confidence intervals
  print("for(p in 1:nprey){")
  for(p in 1:nprey){
    lower95[p,i] <- qbeta(.025,betaEst[p,1],betaEst[p,2])
    upper95[p,i] <- qbeta(.975,betaEst[p,1],betaEst[p,2])
    
    #sensitivity of the MLE - change this range
    prop2=seq(0.000001,0.999999,0.000001) 
    density_=-dbeta(prop2,betaEst[p,1],betaEst[p,2],log=TRUE)
    mode_[p,i] = prop2[which(density_ == min(density_))]
    
    
  }
  #Beta Parameters
  print("  beta_a[,i] <- t(betaEst)[1,]")
  beta_a[,i] <- t(betaEst)[1,]
  print("beta_b[,i] <- t(betaEst)[2,] ")
  beta_b[,i] <- t(betaEst)[2,] 
  
  
  
} #iloopends

npredlist  <- unique(x$Predator)


#open file

#write column names
columnnames<-c("PredatorName","PreyName","lower95","upper95", "mode", "simpleaverages", "bootaverage", "bootmode")
#write to file here

preynames=names(x)[3:length(names(x))]

writeDF=columnnames
for(i in FirstPredToDo:NumPredsToDo)  {
  
  mlesum = 0
  for(j in 1:nprey){
    mlesum = mlesum + mle[j,i]
  }
  
  for(j in 1:nprey){
    newrow = c(toString(npredlist[i]),preynames[j],lower95[j,i],upper95[j,i],(mle[j,i]/mlesum),raw_wt_avg[j,i],bootaverage[j,i],bootmode[j,i])
    writeDF=rbind(writeDF,newrow)
    
  }
  
  
}

setwd(savepath)
write.csv(writeDF, file = paste(getwd(),"/","DietResults.csv",sep=""),row.names=FALSE)

write.csv(mode_,"Mode_pda",row.names = FALSE)
#---------------------------------OUTPUT----------------------------------------
#output <- list(raw_wt_avg,mle,betavar,lower95,upper95,beta_a,beta_b)
#names(output) <- c("raw_wt_avg","mle","betavar","lower95","upper95","beta_a","beta_b")
#return(output)
#mydata=read.csv("C:/Users/Michelle/Desktop/Atlantis/Diet/morethan10.csv")
#FitDietDirichlet(mydata)
#}
#}

setwd(workpath)
