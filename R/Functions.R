transforma=function(mat,lim){

  n=length(mat)

  res=NULL

  for(i in 1:n){

    if((mat[i]>=lim)){

      temp=i
      res=c(res,temp)
    }



  }

  res
}
######################################################
logverouni = function(xx,Tao,nnt,datos){
  theta=xx[1]
  alpha=xx[2]
  res=-(-theta*Tao^alpha+nnt*log(theta)+nnt*log(alpha)+alpha*sum(log(datos),na.rm = T))
  res
}
######################################################
njfunction=function(data,a,b){

  nnj=NULL

  for(i in 1:ncol(data)){

    res=sum(ifelse((data[,i]>a)&(data[,i]<=b),1,0),na.rm=T)
    nnj=c(nnj,res)
  }

  Mdata=array(NA,dim=c(max(nnj),ncol(data)))

  for(i in 1:ncol(data)){
    Mdata[,i]=c(data[which((data[,i]>a)&(data[,i]<=b)),i],rep(NA,(max(nnj)-length(which((data[,i]>a)&(data[,i]<=b))))))
  }

  res=list(Mdata,nnj)
  res

}
######################################################
gCorr<-function(b,def){
  n=nrow(def)
  R=exp(-b*(as.matrix(dist(def))))
  mat=R
  mat


}
######################################################
gSigma<-function(b,v,def){
  n<-nrow(def)
  R<-exp(-b*(as.matrix(dist(def))))
  mat<-v*R
  mat


}
######################################################
sintonizar=function(bar,taxa,tau,mat,i){

  mat=as.matrix(mat)



  mater=(1/50)*sum(mat[(i-49):i,1])

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }



}
######################################################
sintonizarN=function(taxa,tau,mat,i){

  mat=as.matrix(mat)



  mater=(1/50)*sum(mat[(i-49):i,1])

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }


}
#sintonizarN(0.25,SU2,Wtest,i)
#####################################################
amostrarW1=function(WW,MM,S,XX,PPs,bb,vv,nn,TT,ff){

  n=nrow(WW)
  WWprop=as.matrix(mvrnorm(1,WW,ff*diag(1,n)))


  SSig=gSigma(bb,vv,S)

  postWW=sum(WW*nn)-sum(exp(WW)*(TT)^exp(MM))-0.5*t(WW-t(XX)%*%PPs)%*%solve(SSig)%*%(WW-t(XX)%*%PPs)
  postWWprop=sum(WWprop*nn)-sum(exp(WWprop)*(TT)^exp(MM))-0.5*t(WWprop-t(XX)%*%PPs)%*%solve(SSig)%*%(WWprop-t(XX)%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=WWprop

    rejei=1


  }else{

    Wprox=WW
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res
}
#amostrarW1(as.matrix(W1),as.matrix(rMm[1,]),rD1,Ftw,as.matrix(rPsi[1,]),rbw,rvw,as.matrix(nnj[1,]),as.matrix(Tao),SU1[1])
#####################################################
amostrarW=function(WW,MM,S,XX,PPs,bb,vv,nn,taoini,taofim,ff){

  n=nrow(WW)
  WWprop=as.matrix(mvrnorm(1,WW,ff*diag(1,n)))


  SSig=gSigma(bb,vv,S)


  postWW=sum(WW*nn)+sum( exp(WW)*(taoini)^exp(MM)-exp(WW)*(taofim)^exp(MM) )-0.5*t(WW-XX%*%PPs)%*%solve(SSig)%*%(WW-XX%*%PPs)
  postWWprop=sum(WWprop*nn)+sum( exp(WWprop)*(taoini)^exp(MM)-exp(WWprop)*(taofim)^exp(MM) )-0.5*t(WWprop-XX%*%PPs)%*%solve(SSig)%*%(WWprop-XX%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=WWprop

    rejei=1


  }else{

    Wprox=WW
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res
}
######################################################
amostrarM=function(WW,MM,S,ZZ,BBta,bb,vv,data,nn,taoini,taofim,ff){

  MM=as.matrix(MM)
  n=nrow(MM)
  MMprop=as.matrix(mvrnorm(1,MM,ff*diag(1,n)))
  logdata=log(data)
  SSig=gSigma(bb,vv,S)

  sum1=0
  sum2=0

  for(j in 1:ncol(data)){
    res1=sum(exp(MM[j,])*logdata[,j],na.rm=T)
    res2=sum(exp(MMprop[j,])*logdata[,j],na.rm=T)
    sum1=sum1+res1
    sum2=sum2+res2
  }

  postMM=sum(MM*nn)+sum(exp(WW)*(taoini)^exp(MM))-sum(exp(WW)*(taofim)^exp(MM))+sum1-0.5*t(MM-ZZ%*%BBta)%*%solve(SSig)%*%(MM-ZZ%*%BBta)
  postMMprop=sum(MMprop*nn)+sum(exp(WW)*(taoini)^exp(MMprop))-sum(exp(WW)*(taofim)^exp(MMprop))+sum2-0.5*t(MMprop-ZZ%*%BBta)%*%solve(SSig)%*%(MMprop-ZZ%*%BBta)

  prob=min(exp((postMMprop)-(postMM)),1)


  u=runif(1,0,1)

  if(u<prob){

    MMprox=MMprop

    rejei=1


  }else{

    MMprox=MM
    rejei=0
  }





  res=as.matrix(MMprox)
  res=list(MMprox,rejei)
  res



}
#amostrarM(Wm[,1],Mm[,1],sites,Zm[1,,],as.matrix(Beta[1,]),bm,vm,listaD[1][[1]],nnj[1,],ini[1],fim[1],SU2[1])


######################################################
amostrarb=function(W,v,b,loca,ab,bb,X,PPsi,u1){

  bprop=rgamma(1,shape=b*u1, rate = u1)

  SSigprop=gSigma(bprop,v,loca)

  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b,0))
  }

  SSig=gSigma(b,v,loca)
  Corr=gCorr(b,loca)
  SSigprop=gSigma(bprop,v,loca)
  Corrprop=gCorr(bprop,loca)

  acu=0
  acuprop=0

  for(i in 1:nrow(W)){

    acu=acu+t(as.matrix(W[i,])-t(matrix(X[i,,],ncol(PPsi),nrow(loca)))%*%as.matrix(PPsi[i,]))%*%solve(SSig)%*%(as.matrix(W[i,])-t(matrix(X[i,,],ncol(PPsi),nrow(loca)))%*%as.matrix(PPsi[i,]))
    acuprop=acuprop+t(as.matrix(W[i,])-t(matrix(X[i,,],ncol(PPsi),nrow(loca)))%*%as.matrix(PPsi[i,]))%*%solve(SSigprop)%*%(as.matrix(W[i,])-t(matrix(X[i,,],ncol(PPsi),nrow(loca)))%*%as.matrix(PPsi[i,]))

  }

  logp=-0.5*nrow(W)*log(det(Corr))-0.5*acu-bb*b+( ab-1)*log(b)

  logpprop=-0.5*nrow(W)*log(det(Corrprop))-0.5*acuprop-bb*bprop+( ab-1)*log(bprop)

  logprob=logpprop+log(dgamma(b,shape=bprop*u1,rate=u1))-(logp+log(dgamma(bprop,shape=b*u1,rate=u1)))

  prob=min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=bprop

    rejei=1


  }else{

    bprox=b

    rejei=0;

  }

  res=list(bprox,rejei)
  res
}
####################################################################################################
amostrard=function(WW,psi,D,S,phi,sigma,sigmad,Rde,Xl,u1){

  n=nrow(D)
  Do=complex(real=D[,1],imaginary=D[,2]) # Coordenadas complexas
  ztemp=confh(Do) ## Removendo translação e escala
  z=ztemp[[1]] ## Pre-forma de D actual
  esc=ztemp[[2]] ## valor da escala removida
  ## pre-forma de referencia para remover rotação
  marcoc1=complex(real=S[,1],imaginary=S[,2])
  Zotemp=confh(marcoc1)
  Zo=Zotemp[[1]] # Pre-forma de referencia
  theta=Arg(t(Conj(z))%*%Zo)
  ww=as.vector(exp(theta*(0+1i)))*z
  zd=t(H(n))%*%ww ## Configuração Pre-forma de D actual
  vectD=c(Re(zd),Im(zd)) # Configuração pre-forma de D actual vetorizada
  ssigma=u1*diag(rep(1,2))%x%Rde # Matriz de covariancia proposta função de transição
  vectdprop=mvrnorm(1,vectD,ssigma)# Geração de configuração Pre-forma proposta vetorizada
  Dp=cbind(as.matrix(vectdprop[1:n]),as.matrix(vectdprop[(n+1):(2*n)])) # Configuração proposta pre-forma
  dprop=esc*Dp # D proposto escala original

  SSig=gSigma(phi,sigma,D)
  SSigprop=gSigma(phi,sigma,dprop)
  acum1=0
  acum2=0

  for(i in 1:nrow(WW)){
    resu=-0.5*t(as.matrix(WW[i,])-t(matrix(Xl[i,,],ncol(psi),nrow(S)))%*%as.matrix(psi[i,]))%*%solve(SSig)%*%(as.matrix(WW[i,])-t(matrix(Xl[i,,],ncol(psi),nrow(S)))%*%as.matrix(psi[i,]))
    acum1=acum1+resu
    resuprop=-0.5*t(as.matrix(WW[i,])-t(matrix(Xl[i,,],ncol(psi),nrow(S)))%*%as.matrix(psi[i,]))%*%solve(SSigprop)%*%(as.matrix(WW[i,])-t(matrix(Xl[i,,],ncol(psi),nrow(S)))%*%as.matrix(psi[i,]))
    acum2=acum2+resuprop
  }

  logp=-(nrow(psi)/2)*log(det(SSig))+acum1-0.5*sum(diag((D-S)%*%solve(sigmad)%*%t(D-S)%*%solve(Rde)))
  logpprop=-(nrow(psi)/2)*log(det(SSigprop))+acum2-0.5*sum(diag((dprop-S)%*%solve(sigmad)%*%t(dprop-S)%*%solve(Rde)))

  prob=min(exp((logpprop)-(logp)),1)

  u=runif(1,0,1)

  if(u<prob){

    dprox=dprop

    rejei=1

  }else{

    dprox=D
    rejei=0;
  }






  res=list(dprox,rejei)
  res

}
#amostrard(Wm,Psi,sites,sites,bw,vw,diag(1,2),Rd,Xm,0.001)
####################################################################################
H<-function(k){

  matriz<-array(0,dim=c((k-1),k))

  for(i in 1:(k-1)){

    for(j in 1:k){

      if(j<=i){
        matriz[i,j]<--(i*(i+1))^(-1/2)

      }else{
        if(j==(i+1)){
          matriz[i,j]<--i*(-(i*(i+1))^(-1/2))
        }else{

        }
      }

    }


  }

  matriz
}
####################################################################################
confh<-function(vetor){

  n=length(vetor)

  w=as.matrix(H(n)%*%vetor)
  escala=as.numeric(sqrt(t(Conj(w))%*%w))
  z=(1/escala)*w
  res=list(z,escala)
  res
}
####################################################################################
FFBS=function(FFt,GGt,VVt,WWt,mm0,CC0,yy,ttheta){

  tempA=matrix(NA,nrow(ttheta),ncol(ttheta))
  Tt=nrow(yy)
  mm=matrix(0,nrow(yy),nrow(mm0))
  CC=matrix(0,nrow(yy),(nrow(CC0)*ncol(CC0)))

  mmant=mm0
  CCant=CC0

  for(i in 1:nrow(yy)){
    aat=GGt[i,,]%*%mmant
    RRt=GGt[i,,]%*%CCant%*%t(GGt[i,,])+WWt
    fft=t(FFt[i,,])%*%aat
    QQt=t(FFt[i,,])%*%RRt%*%FFt[i,,]+VVt
    AAt=RRt%*%FFt[i,,]%*%solve(QQt)
    eet=as.matrix(yy[i,])-fft
    mmt=aat+AAt%*%eet
    mm[i,]=mmt
    CCt=RRt-AAt%*%QQt%*%t(AAt)
    CC[i,]=as.vector(CCt)
    mmant=mmt
    CCant=CCt
  }

  TTt=nrow(yy)
  mmm=as.matrix(mm[TTt,])
  CCC=matrix(CC[TTt,],ncol(mm),ncol(mm))
  tempA[TTt,]=as.matrix(mvrnorm(1,mmm,CCC))

  for(j in (TTt-1):1){
    mmm=as.matrix(mm[j,])
    CCC=matrix(CC[j,],ncol(mm),ncol(mm))
    varA=solve(t(GGt[j,,])%*%solve(WWt)%*%(GGt[j,,])+solve(CCC))
    mediaA=varA%*%(t(GGt[j,,])%*%solve(WWt)%*%as.matrix(tempA[(j+1),])+solve(CCC)%*%mmm)
    tempA[j,]=as.matrix(mvrnorm(1,mediaA,varA))
  }

  tempA
}

####################################################################################
vero=function(Data,nnn,WW,MMm,tao){

  logvero=array(NA,dim=c(nrow(nnn),1))

  for(x in 1:nrow(MMm)){

    temp1=sum(as.matrix(nnn[x,])*(as.matrix(WW[x,])+as.matrix(MMm[x,])))

    if(x==1){
      temp2=sum(-exp(as.matrix(WW[x,]))*tao[x]^(exp(as.matrix(MMm[x,]))))
    }else{
      temp2=sum( -exp(as.matrix(WW[x,]))*tao[x]^(exp(as.matrix(MMm[x,])))+exp(as.matrix(WW[x,]))*tao[x-1]^(exp(as.matrix(MMm[x,]))) )

    }

    temp3=0

    for(j in 1:ncol(MMm)){

      aux=sum(exp(MMm[x,j])*log( Data[[x]][,j] ),na.rm = T)
      temp3=temp3+aux
    }
    logvero[x,1]=temp1+temp2+temp3



  }

  res=sum(logvero)
  res

}


