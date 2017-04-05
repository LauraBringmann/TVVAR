



GAM_TVVAR <- function(Data,
                    np,
                    N,
                    K
                    )



{
  Results_GAM<-array(NA,c(c(np+1,np),N,3))


  #N=100 # 100 500
  tt=1:N

  #%##########################################%###
  ####  Part 1: creating the data #############
  #%##########################################%##
  #The functions in this part are created in the file (source code for simulation article 3.R)



  y<-Data


  #%##########################################%###
  ####  Part 2: ESTIMATING GAM##### #############
  #%##########################################%##

  for (ii in 1:np){

    mod<-Estimate_gam(y,N,np,K)[[ii]]

    mat_dat<-matrix(c(tt,rep(rep(1,N),np)),length(tt),np+1)
    coln_data<-paste("y",1:np,"L",sep="")
    coln_data_full<-c("tt",coln_data)
    colnames(mat_dat)<-coln_data_full
    newd<-as.data.frame(mat_dat)
    Xp=predict(mod,newd,type="lpmatrix",seWithMean = TRUE)
    kdim=dim(Xp)[2]/c(np+1)
    newpre=predict(mod,new.data=newd,type="terms",se=TRUE)
    Results_GAM[1,ii,1:N,2]<-Xp[,1:kdim]%*% coef(mod)[1:kdim]#basis functions intercept!

    Numbrep=1000
    modr<-mvrnorm(Numbrep,coef(mod),mod$Vp+diag((np+1)*kdim)*10^(-14))

    #The confidence intervals
    int.ci<-matrix(NA,N,Numbrep)
    for (m in 1:Numbrep){
      int.ci[,m]<-  Xp[,1:kdim]%*%modr[m,1:kdim]
    }

    Results_GAM[1,ii,1:N,1]<-apply(int.ci,1,quantile,c(.975))
    Results_GAM[1,ii,1:N,3]<-apply(int.ci,1,quantile,c(.025))


    for (j in 1:np){
      Results_GAM[j+1,ii,1:N,2]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*% coef(mod)[(j*kdim+1):((j+1)*kdim)]

      #The confidence intervals
      phi.ci<-matrix(NA,N,Numbrep)
      for (m in 1:Numbrep){
        phi.ci[,m]<-Xp[,(j*kdim+1):((j+1)*kdim)]%*%modr[m,(j*kdim+1):((j+1)*kdim)]
      }
      Results_GAM[j+1,ii,1:N,1]<-apply(phi.ci,1,quantile,c(.975))
      Results_GAM[j+1,ii,1:N,3]<-apply(phi.ci,1,quantile,c(.025))
    }




  }
  return(Results_GAM=Results_GAM)
}

