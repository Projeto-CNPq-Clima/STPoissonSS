#' Bayesian Estimation for Spatiotemporal Nonhomogeneous Poisson Model with Weibull Intensity and State-Space Prior
#'
#' This function performs Bayesian inference for a class of nonhomogeneous Poisson spatiotemporal models. The intensity function follows a Weibull form with time-varying scale and shape parameters modeled through a state-space prior structure. Spatial dependence is incorporated via anisotropic deformation, and model parameters are estimated using a Markov chain Monte Carlo (MCMC) algorithm.
#'
#' @param data A list or array of observed counts following a nonhomogeneous Poisson process structure, as used in other spatiotemporal Poisson models.
#' @param sites A matrix or data structure specifying the spatial coordinates or structure of observation sites, as in standard spatial Poisson models.
#' @param tau_breaks A numeric vector containing the internal breakpoints \eqn{(\tau_1, \ldots, \tau_{L})} that partition the time interval \eqn{(0, T]} into \eqn{L} disjoint subintervals. The vector should be sorted in increasing order and must not include the initial time 0 or the endpoint T.
#' @param Ftw An array of dimension `L × p × n`. The l-th slice `Xm[l,,]` corresponds to the transposed covariate matrix `X'_l` used to explain the latent vector `W_l`.
#' @param Ftm An array of dimension `L × q × n`. The l-th slice `Zm[l,,]` corresponds to the transposed covariate matrix `X'_l` used to explain the latent vector `M_l`.
#' @param G An array of dimension `L × p × p`. Each slice `G[l,,]` corresponds to the known evolution matrix `G_l` of size `p × p`, which governs the temporal evolution of the latent process associated with `W_l`.
#' @param Q An array of dimension `L × q × q`. Each slice `Q[l,,]` corresponds to the known evolution matrix `Q_l` of size `q × q`, which governs the temporal evolution of the latent process associated with `M_l`.
#' @param prior A list containing the hyperparameters for the prior distributions:
#' \describe{
#'   \item{M0, C0}{Prior mean and covariance for the initial state of `Psi`, i.e., `Psi_0 ~ N(M0, C0)`}
#'   \item{S0, n0}{Scale matrix and degrees of freedom for the inverse-Wishart prior on `V_Psi`, i.e., `V_Psi ~ Inv-Wishart(S0, n0)`}
#'   \item{M1, C1}{Prior mean and covariance for the initial state of `beta`, i.e., `beta_0 ~ N(M1, C1)`}
#'   \item{S1, n1}{Scale matrix and degrees of freedom for the inverse-Wishart prior on `V_beta`, i.e., `V_beta ~ Inv-Wishart(S1, n1)`}
#'   \item{aa1, bb1}{Shape and scale parameters for the inverse-gamma prior on `sigma^2_w`, i.e., `sigma^2_w ~ IG(aa1, bb1)`}
#'   \item{aa2, bb2}{Shape and scale parameters for the inverse-gamma prior on `sigma^2_m`, i.e., `sigma^2_m ~ IG(aa2, bb2)`}
#'   \item{c3, d3}{Shape and scale parameters for the gamma prior on `phi`, i.e., `phi ~ Gamma(c3, d3)`}
#' }
#' @param iteration Total number of MCMC iterations.
#' @param burnin Number of initial iterations to discard as burn-in.
#' @param jump Thinning interval: number of iterations to skip between two consecutive stored MCMC samples after burn-in.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{MWT}{Posterior samples of the time-varying scale parameter (W).}
#'   \item{MMT}{Posterior samples of the time-varying shape parameter (M).}
#'   \item{MbmT}{Posterior samples of the intercept term for `M`.}
#'   \item{MbwT}{Posterior samples of the intercept term for `W`.}
#'   \item{MPsi}{Posterior samples of the latent state vectors for the scale component.}
#'   \item{MBeta}{Posterior samples of the latent state vectors for the shape component.}
#'   \item{MPsi0}{Posterior samples of the initial state of `Psi`.}
#'   \item{MBeta0}{Posterior samples of the initial state of `Beta`.}
#'   \item{Mvw}{Posterior samples of the variance parameter for `W`.}
#'   \item{Mvm}{Posterior samples of the variance parameter for `M`.}
#'   \item{Mbw}{Posterior samples of regression coefficients for `W`.}
#'   \item{Mbm}{Posterior samples of regression coefficients for `M`.}
#'   \item{MVpsi}{Posterior samples of the covariance matrix for `Psi`.}
#'   \item{MVbta}{Posterior samples of the covariance matrix for `Beta`.}
#'   \item{MW}{Posterior mean of the W process (scale).}
#'   \item{MM}{Posterior mean of the M process (shape).}
#' }
#'
#' @export


bayesSTWeibull<-function(data,sites,tau_breaks,Ftw,Ftm,G,Q,
                         prior=list(M0=as.matrix(rep(0,nrow(Ftw[1,,]))),
                                    C0=diag(100,nrow(Ftw[1,,])),
                                    S0=diag(0.001,nrow(Ftw[1,,])),
                                    n0=10,M1=as.matrix(rep(0,nrow(Ftw[1,,]))),
                                    M1=as.matrix(rep(0,nrow(Ftm[1,,]))),
                                    C1=diag(100,nrow(Ftw[1,,])),S1=diag(0.001,nrow(Ftw[1,,])),
                                    S1=diag(0.001,nrow(Ftm[1,,])),
                                    n1=10,aa1=2.1,
                                    bb1=1.1,
                                    aa2=2.1,
                                    bb2=1.1,
                                    c3= 5.197782e-08,
                                    d3= 2.279865e-05,
                                    c2=5.197782e-08,
                                    d2=2.279865e-05),iteration,burnin,jump){




  end_indices=tau_breaks
  start_indices=c(0,end_indices[1:(length(end_indices)-1)])

  tau_breaks-tau_breaks[1]
  start_indices<-c(0,end_indices[1:(length(end_indices)-1)])

  L=length(start_indices)
  nnj=array(NA, dim=c(L,ncol(data)))
  listaD=vector("list", length = length(start_indices))


  for(j in 1:L){
    temp1= njfunction(data,start_indices[j],end_indices[j])
    nnj[j,]=temp1[[2]]
    listaD[[j ]]=temp1[[1]]
  }



  ######################################

  #Psi_0~N(M0,C0)
  M0=prior$M0
  C0=prior$C0
  S0=prior$S0
  n0=prior$n0
  M1=prior$M1
  C1=prior$C1
  #V_{beta} ~ W_{n1}^{-1}(S1)
  S1=prior$S1
  n1=prior$n1
  #sigma^2_w~GI(aa1,bb1)
  aa1=prior$aa1
  bb1=prior$bb1
  #sigma^2_m~GI(aa2,bb2)
  aa2=prior$aa2
  bb2=prior$bb2
  #phi~G(c3,d3), E(phi)=rbd, V(phi)=100

  c3= prior$c3
  d3= prior$d3

  c2=prior$c2
  d2=prior$d2
  ####################

  S=sites
  S=complex(real=S[,1],imaginary=S[,2])
  stemp=H(nrow(sites))%*%S ## Removendo translação
  Sz=t(H(nrow(sites)))%*%stemp ## configuração com translação removida
  Sz=cbind(as.matrix(Re(Sz)),as.matrix(Im(Sz)))

  vw=1
  bw=1
  vm=1
  bm=1
  Psi0=rep(0.1,nrow(Ftw[1,,]))
  Psi=array(0.1,dim=c(L,nrow(Ftw[1,,])))
  Beta=array(0.1,dim=c(L,nrow(Ftw[1,,])))

  Wm=matrix(NA,L,ncol(data))
  Mm=matrix(NA,L,ncol(data))

  temp1= njfunction(data,0,end_indices[length(end_indices)])
  nnjj=temp1[[2]]
  Wmtemp=NULL
  Mmtemp=NULL

  for(i in 1:ncol(data)){
    mmmtemp=optim(c(1,1),logverouni,c(1,1),end_indices[length(end_indices)],nnjj[i],data[,i])
    Wmtemp=c(Wmtemp,log(mmmtemp$par[1]))
    Mmtemp=c(Mmtemp,log(mmmtemp$par[2]))

  }

  for(x in 1:L){
    Wm[x,]=Wmtemp
    Mm[x,]=Mmtemp
  }

  Beta0=rep(0.1,nrow(Ftw[1,,]))
  Vpsi=diag(0.1,nrow(Ftw[1,,]))
  Vbta=diag(0.1,nrow(Ftw[1,,]))
  #########################
  # Hiperparametros
  n0est=n0+L
  n1est=n1+L

  SU1=rep(0.0001,L)
  SU2=rep(0.0001,L)
  SU3=5000
  SU4=5000
  SU5=0.0001
  SU6=0.0001

  MWT=NULL
  MMT=NULL
  MbmT=NULL
  MbwT=NULL
  MPsi=NULL
  MBeta=NULL
  MPsi0=NULL
  MBeta0=NULL

  Mvw=NULL
  Mvm=NULL
  Mbw=NULL
  Mbm=NULL
  MVpsi=NULL
  MVbta=NULL
  MW=NULL
  MM=NULL



  for(i in 1:iteration){



    if(i<=burnin){

      gamma0=solve(solve(C0)+t(G[1,,])%*%solve(Vpsi)%*%G[1,,])
      eta0=gamma0%*%(solve(C0)%*%M0+t(G[1,,])%*%solve(Vpsi)%*%as.matrix(Psi[1,]))
      Psi0=as.matrix(mvrnorm(1,eta0,gamma0))

      Vt=gSigma(bw,vw,Sz)

      tem=tryCatch(FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi),error = function(e) e)

      if(is.matrix(tem)==T){
        Psi=FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi)
      }else{

      }

      B0=solve(solve(C1)+t(Q[1,,])%*%solve(Vbta)%*%Q[1,,])
      A0=B0%*%(solve(C1)%*%M1+t(Q[1,,])%*%solve(Vbta)%*%as.matrix(Beta[1,]))
      Beta0=as.matrix(mvrnorm(1,A0,B0))

      Vmt=gSigma(bm,vm,Sz)

      tem1=tryCatch(FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta),error = function(e) e)

      if(is.matrix(tem1)==T){
        Beta=FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta)
      }else{

      }



      auxw=NULL
      auxm=NULL
      acum0=0
      acum1=0
      acum0w=0
      acum0m=0


      for(x in 1:L){

        temp=amostrarW(as.matrix(Wm[x,]),as.matrix(Mm[x,]),Sz,t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz))),as.matrix(Psi[x,]),bw,vw,as.matrix(nnj[x,]),start_indices[x],end_indices[x],SU1[x])
        Wm[x,]=temp[[1]]
        auxw=c(auxw,temp[[2]])
        aux1=temp[[2]]

        temp=amostrarM(Wm[x,],Mm[x,],Sz,t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz))),as.matrix(Beta[x,]),bm,vm,listaD[x][[1]],nnj[x,],start_indices[x],end_indices[x],SU2[x])
        Mm[x,]=temp[[1]]
        auxm=c(auxm,temp[[2]])

        if(x==1){
          acum0=acum0+( as.matrix(Psi[x,])-G[1,,]%*%Psi0 )%*% t( as.matrix(Psi[x,])-G[1,,]%*%Psi0 )
          acum1=acum1+( as.matrix(Beta[x,])-Q[1,,]%*%Beta0 )%*% t( as.matrix(Beta[x,])-Q[1,,]%*%Beta0 )
          acum0w=acum0w+t(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))%*%solve(gCorr(bw,Sz))%*%(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))
          acum0m=acum0m+t(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))%*%solve(gCorr(bm,Sz))%*%(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))

        }else{

          acum0=acum0+( as.matrix(Psi[x,])-G[1,,]%*%as.matrix(Psi[(x-1),]) )%*% t( as.matrix(Psi[x,])-G[1,,]%*%as.matrix(Psi[(x-1),]) )
          acum1=acum1+( as.matrix(Beta[x,])-Q[1,,]%*%as.matrix(Beta[(x-1),]) )%*% t( as.matrix(Beta[x,])-Q[1,,]%*%as.matrix(Beta[(x-1),]) )
          acum0w=acum0w+t(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))%*%solve(gCorr(bw,Sz))%*%(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))
          acum0m=acum0m+t(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))%*%solve(gCorr(bm,Sz))%*%(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))

        }


        if((i%%50)==0){
          Wtest=c(MWT[,x],aux1)
          Mtest=c(MMT[,x],temp[[2]])
          SU1[x]=sintonizarN(0.25,SU1[x],Wtest,i)
          SU2[x]=sintonizarN(0.25,SU2[x],Mtest,i)
        }else{

        }



      }
      MWT=rbind(MWT,t(auxw))
      MMT=rbind(MMT,t(auxm))

      aa=(nrow(Sz)*L/2)+aa2
      bb=0.5*acum0w+bb2
      vw=1/rgamma(1,shape=aa, rate = bb)

      aaa=(nrow(Sz)*L/2)+aa2
      bbb=0.5*acum0m+bb2
      vm=1/rgamma(1,shape=aaa, rate = bbb)

      temp=amostrarb(Wm,vw,bw,Sz,c3,d3,Ftw,Psi,SU3)
      bw=temp[[1]]
      MbwT=c(MbwT,temp[[2]])

      temp=amostrarb(Mm,vm,bm,Sz,c3,d3,Ftm,Beta,SU4)
      bm=temp[[1]]
      MbmT=c(MbmT,temp[[2]])

      S0est=solve(S0+acum0)
      temp=solve(array(as.vector(rWishart(1, n0est, S0est)),dim=c(ncol(G[1,,]),ncol(G[1,,]))))
      Vpsi=temp

      S1est=solve(S1+acum1)
      temp=solve(array(as.vector(rWishart(1, n1est, S1est)),dim=c(ncol(Q[1,,]),ncol(Q[1,,]))))
      Vbta=temp

      if((i%%50)==0){
        SU3=sintonizar(burnin,0.44,SU3,MbwT,i)
        SU4=sintonizar(burnin,0.44,SU4,MbmT,i)

      }else{

      }


    }else{

      if((i%%150)==0){
        MWT=NULL
        MMT=NULL
        MbmT=NULL
        MbwT=NULL

      }else{

      }

      if((i%%jump)==0){

        gamma0=solve(solve(C0)+t(G[1,,])%*%solve(Vpsi)%*%G[1,,])
        eta0=gamma0%*%(solve(C0)%*%M0+t(G[1,,])%*%solve(Vpsi)%*%as.matrix(Psi[1,]))
        Psi0=as.matrix(mvrnorm(1,eta0,gamma0))
        MPsi0=rbind(MPsi0,t(Psi0))

        Vt=gSigma(bw,vw,Sz)

        tem=tryCatch(FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi),error = function(e) e)

        if(is.matrix(tem)==T){
          Psi=FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi)
        }else{

        }

        MPsi=rbind(MPsi,t(as.matrix(as.vector(Psi))))


        B0=solve(solve(C1)+t(Q[1,,])%*%solve(Vbta)%*%Q[1,,])
        A0=B0%*%(solve(C1)%*%M1+t(Q[1,,])%*%solve(Vbta)%*%as.matrix(Beta[1,]))
        Beta0=as.matrix(mvrnorm(1,A0,B0))
        MBeta0=rbind(MBeta0,t(Beta0))

        Vmt=gSigma(bm,vm,Sz)

        tem1=tryCatch(FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta),error = function(e) e)

        if(is.matrix(tem1)==T){
          Beta=FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta)
        }else{

        }

        MBeta=rbind(MBeta,t(as.matrix(as.vector(Beta))))


        auxw=NULL
        auxm=NULL
        acum0=0
        acum1=0
        acum0w=0
        acum0m=0


        gamma0=solve(solve(C0)+t(G[1,,])%*%solve(Vpsi)%*%G[1,,])
        eta0=gamma0%*%(solve(C0)%*%M0+t(G[1,,])%*%solve(Vpsi)%*%as.matrix(Psi[1,]))
        Psi0=as.matrix(mvrnorm(1,eta0,gamma0))

        Vt=gSigma(bw,vw,Sz)

        tem=tryCatch(FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi),error = function(e) e)

        if(is.matrix(tem)==T){
          Psi=FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi)
        }else{

        }

        B0=solve(solve(C1)+t(Q[1,,])%*%solve(Vbta)%*%Q[1,,])
        A0=B0%*%(solve(C1)%*%M1+t(Q[1,,])%*%solve(Vbta)%*%as.matrix(Beta[1,]))
        Beta0=as.matrix(mvrnorm(1,A0,B0))

        Vmt=gSigma(bm,vm,Sz)

        tem1=tryCatch(FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta),error = function(e) e)

        if(is.matrix(tem1)==T){
          Beta=FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta)
        }else{

        }



        auxw=NULL
        auxm=NULL
        acum0=0
        acum1=0
        acum0w=0
        acum0m=0



        for(x in 1:L){

          temp=amostrarW(as.matrix(Wm[x,]),as.matrix(Mm[x,]),Sz,t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz))),as.matrix(Psi[x,]),bw,vw,as.matrix(nnj[x,]),start_indices[x],end_indices[x],SU1[x])
          Wm[x,]=temp[[1]]
          auxw=c(auxw,temp[[2]])

          temp=amostrarM(Wm[x,],Mm[x,],Sz,t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz))),as.matrix(Beta[x,]),bm,vm,listaD[x][[1]],nnj[x,],start_indices[x],end_indices[x],SU2[x])
          Mm[x,]=temp[[1]]
          auxm=c(auxm,temp[[2]])

          if(x==1){
            acum0=acum0+( as.matrix(Psi[x,])-G[1,,]%*%Psi0 )%*% t( as.matrix(Psi[x,])-G[1,,]%*%Psi0 )
            acum1=acum1+( as.matrix(Beta[x,])-Q[1,,]%*%Beta0 )%*% t( as.matrix(Beta[x,])-Q[1,,]%*%Beta0 )
            acum0w=acum0w+t(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))%*%solve(gCorr(bw,Sz))%*%(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))
            acum0m=acum0m+t(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))%*%solve(gCorr(bm,Sz))%*%(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))

          }else{

            acum0=acum0+( as.matrix(Psi[x,])-G[1,,]%*%as.matrix(Psi[(x-1),]) )%*% t( as.matrix(Psi[x,])-G[1,,]%*%as.matrix(Psi[(x-1),]) )
            acum1=acum1+( as.matrix(Beta[x,])-Q[1,,]%*%as.matrix(Beta[(x-1),]) )%*% t( as.matrix(Beta[x,])-Q[1,,]%*%as.matrix(Beta[(x-1),]) )
            acum0w=acum0w+t(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))%*%solve(gCorr(bw,Sz))%*%(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))
            acum0m=acum0m+t(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))%*%solve(gCorr(bm,Sz))%*%(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))

          }


        }

        MWT=rbind(MWT,t(auxw))
        MW=rbind(MW,t(as.matrix(as.vector(Wm))))
        MMT=rbind(MMT,t(auxm))
        MM=rbind(MM,t(as.matrix(as.vector(Mm))))

        aa=(nrow(Sz)*L/2)+aa2
        bb=0.5*acum0w+bb2
        vw=1/rgamma(1,shape=aa, rate = bb)
        Mvw=c(Mvw,vw)

        aaa=(nrow(Sz)*L/2)+aa2
        bbb=0.5*acum0m+bb2
        vm=1/rgamma(1,shape=aaa, rate = bbb)
        Mvm=c(Mvm,vm)

        temp=amostrarb(Wm,vw,bw,Sz,c3,d3,Ftw,Psi,SU3)
        bw=temp[[1]]
        MbwT=c(MbwT,temp[[2]])
        Mbw=c(Mbw,bw)

        temp=amostrarb(Mm,vm,bm,Sz,c3,d3,Ftm,Beta,SU4)
        bm=temp[[1]]
        MbmT=c(MbmT,temp[[2]])
        Mbm=c(Mbm,bm)

        S0est=solve(S0+acum0)
        temp=solve(array(as.vector(rWishart(1, n0est, S0est)),dim=c(ncol(G[1,,]),ncol(G[1,,]))))
        Vpsi=temp
        MVpsi=rbind(MVpsi,t(as.matrix(as.vector(Vpsi))))

        S1est=solve(S1+acum1)
        temp=solve(array(as.vector(rWishart(1, n1est, S1est)),dim=c(ncol(Q[1,,]),ncol(Q[1,,]))))
        Vbta=temp
        MVbta=rbind(MVbta,t(as.matrix(as.vector(Vbta))))


      }else{
        gamma0=solve(solve(C0)+t(G[1,,])%*%solve(Vpsi)%*%G[1,,])
        eta0=gamma0%*%(solve(C0)%*%M0+t(G[1,,])%*%solve(Vpsi)%*%as.matrix(Psi[1,]))
        Psi0=as.matrix(mvrnorm(1,eta0,gamma0))

        Vt=gSigma(bw,vw,Sz)

        tem=tryCatch(FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi),error = function(e) e)

        if(is.matrix(tem)==T){
          Psi=FFBS(Ftw,G,Vt,Vpsi,M0,C0,Wm,Psi)
        }else{

        }

        B0=solve(solve(C1)+t(Q[1,,])%*%solve(Vbta)%*%Q[1,,])
        A0=B0%*%(solve(C1)%*%M1+t(Q[1,,])%*%solve(Vbta)%*%as.matrix(Beta[1,]))
        Beta0=as.matrix(mvrnorm(1,A0,B0))

        Vmt=gSigma(bm,vm,Sz)

        tem1=tryCatch(FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta),error = function(e) e)

        if(is.matrix(tem1)==T){
          Beta=FFBS(Ftm,Q,Vmt,Vbta,M1,C1,Mm,Beta)
        }else{

        }



        auxw=NULL
        auxm=NULL
        acum0=0
        acum1=0
        acum0w=0
        acum0m=0



        for(x in 1:L){

          temp=amostrarW(as.matrix(Wm[x,]),as.matrix(Mm[x,]),Sz,t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz))),as.matrix(Psi[x,]),bw,vw,as.matrix(nnj[x,]),start_indices[x],end_indices[x],SU1[x])
          Wm[x,]=temp[[1]]
          auxw=c(auxw,temp[[2]])

          temp=amostrarM(Wm[x,],Mm[x,],Sz,t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz))),as.matrix(Beta[x,]),bm,vm,listaD[x][[1]],nnj[x,],start_indices[x],end_indices[x],SU2[x])
          Mm[x,]=temp[[1]]
          auxm=c(auxm,temp[[2]])

          if(x==1){
            acum0=acum0+( as.matrix(Psi[x,])-G[1,,]%*%Psi0 )%*% t( as.matrix(Psi[x,])-G[1,,]%*%Psi0 )
            acum1=acum1+( as.matrix(Beta[x,])-Q[1,,]%*%Beta0 )%*% t( as.matrix(Beta[x,])-Q[1,,]%*%Beta0 )
            acum0w=acum0w+t(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))%*%solve(gCorr(bw,Sz))%*%(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))
            acum0m=acum0m+t(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))%*%solve(gCorr(bm,Sz))%*%(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))

          }else{

            acum0=acum0+( as.matrix(Psi[x,])-G[1,,]%*%as.matrix(Psi[(x-1),]) )%*% t( as.matrix(Psi[x,])-G[1,,]%*%as.matrix(Psi[(x-1),]) )
            acum1=acum1+( as.matrix(Beta[x,])-Q[1,,]%*%as.matrix(Beta[(x-1),]) )%*% t( as.matrix(Beta[x,])-Q[1,,]%*%as.matrix(Beta[(x-1),]) )
            acum0w=acum0w+t(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))%*%solve(gCorr(bw,Sz))%*%(as.matrix(Wm[x,])-t(matrix(Ftw[x,,],ncol(Psi),nrow(Sz)))%*%as.matrix(Psi[x,]))
            acum0m=acum0m+t(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))%*%solve(gCorr(bm,Sz))%*%(as.matrix(Mm[x,])-t(matrix(Ftm[x,,],ncol(Beta),nrow(Sz)))%*%as.matrix(Beta[x,]))

          }


        }

        MWT=rbind(MWT,t(auxw))
        MMT=rbind(MMT,t(auxm))

        aa=(nrow(Sz)*L/2)+aa2
        bb=0.5*acum0w+bb2
        vw=1/rgamma(1,shape=aa, rate = bb)

        aaa=(nrow(Sz)*L/2)+aa2
        bbb=0.5*acum0m+bb2
        vm=1/rgamma(1,shape=aaa, rate = bbb)

        temp=amostrarb(Wm,vw,bw,Sz,c3,d3,Ftw,Psi,SU3)
        bw=temp[[1]]
        MbwT=c(MbwT,temp[[2]])

        temp=amostrarb(Mm,vm,bm,Sz,c3,d3,Ftm,Beta,SU4)
        bm=temp[[1]]
        MbmT=c(MbmT,temp[[2]])

        S0est=solve(S0+acum0)
        temp=solve(array(as.vector(rWishart(1, n0est, S0est)),dim=c(ncol(G[1,,]),ncol(G[1,,]))))
        Vpsi=temp

        S1est=solve(S1+acum1)
        temp=solve(array(as.vector(rWishart(1, n1est, S1est)),dim=c(ncol(Q[1,,]),ncol(Q[1,,]))))
        Vbta=temp



      }



    }

    print(i)


  }
  resul<-list(MWT,MMT,MbmT,MbwT,MPsi,MBeta,MPsi0,MBeta0,Mvw,Mvm,Mbw,Mbm,MVpsi,MVbta,MW,MM)
  names(resul)<-c("MWT","MMT","MbmT","MbwT","MPsi","MBeta","MPsi0","MBeta0","Mvw","Mvm","Mbw","Mbm","MVpsi","MVbta","MW","MM")
  return(resul)

}


