


Estimate_gam <- function(y,N,np,K){ #y is the y variable,N is the number of time points, np=number of variables and k is the number of knots
  tt=1:N

  data1=matrix(0,N,(np*2),byrow=T)
  for (h in (np+1):(np*2)){
    data1[,(h-np)]=y[,(h-np)]# data wordt hier gelagged

    data1[,h]=c(NA,y[1:(N-1),(h-np)])# data wordt hier gelagged
  }

  data1=as.data.frame(data1)
  colnames(data1)=c(paste("y",1:np,sep=""),paste("y",1:np,"L",sep=""))

  coln=colnames(data1)[1:np]
  colnL=colnames(data1)[(np+1):(np*2)]
  allcol2=c()
  for(i in 1:np){allcol2[i]=paste("s(tt,by=",colnL[i],",k=K",")",sep="")}
  allcol3=paste(allcol2,collapse="+")


  model=list()
  for (j in 1:np){
    ff <- as.formula(paste(coln[j]," ~ ","s(tt,k=K)","+",allcol3))
    model[[j]]<-gam(ff,data=data1,seWithMean=TRUE)

  }

  return(model)

}
