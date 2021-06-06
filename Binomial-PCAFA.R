####################################################################################################
#Police data specific function 
#This assumes structured patterns are not equal across ethnicity.

#This is code for running a single chain Gibbs sampler to fit the principal component (PG-PCA) 
#and factor analysis (PG-FA) model when observed data is multivariate binomial with identification
#constraints as described in Polya-gamma augmentation and latent variable models for multivariate binomial data.

#To run these functions, the R packages BayesLogit and rARPACK need to be installed.

#The function Binom.PGPCA fits the principal component model.
#Inputs required are: 
#y:        	The matrix of binomial successes. Please store as class matrix.
#f:        	The number of components to be extracted.
#maxtrial:	The matrix of binomial trial sizes. Please store as class matrix.
#iter:     	The number of iterations after burnin.
#burnin:   	The number of iterations to be discarded.
#intercept 	TRUE/FALSE statement to specify whether an intercept should be included.
#alphaL, betaL: Parameter values for inverse-gamma prior on loading variances.    

#Outputs are:
#Loading:	A list with entry i corresponding to the identified loading matrix estimate at iteration i.
#Score:		A list with entry i corresponding to the identified score matrix estimate at iteration i.
#var.loading:	A list with entry i corresponding to the vector of estimated loading variances at iteration i.
#Deviance:	A list containing the estimated deviance of the fit at iteration i.
#link:		A list with entry i corresponding to the matrix of estimated realisation of the logit link at iteration i.
#main.effect	If intercept was set to TRUE, a list with entry i corresponding to the estimated of observed variable intercepts at iteration i.


#Version 1: This assumes structured patterns (factor loadings) are not equal across ethnicity.
Binom.PGPCAtie<-function(y,f,maxtrial,iter,burnin,intercept,alphaL,betaL){
  library(BayesLogit)
  library(rARPACK)
  n<-dim(y)[1]
  k<-dim(y)[2]
  
  Kappa<-y-0.5*maxtrial	     				  #corrected y
  keep<-svds(scale(Kappa,center=TRUE,scale=FALSE),k=f,nv=f,nu=f) #Step to help ensure sign invariance.
  
  F<-matrix(rnorm(n*f),n,f)    				  #Initial Factors
  L<-matrix(0,f,k) 					  #Initial Loadings
  
  #Initial estimate of PG rv.
  omega<-matrix(rpg(num=n*k, h=as.numeric(maxtrial), z=0.0),n,k) 
  
  if(intercept == TRUE) {muind<-1} else {muind<-0}
  mu<-muind*rep(0,k)
  
  Devnew<-0					#Initial deviance.
  
  Loading<-list()
  Score<-list()
  eigenval<-list()
  main.effect<-list()
  Link.fun<-list()
  d<-rgamma(f,alphaL,betaL)			#Initial precisions for factor loadings.
  d<-sort(d)
  alphad<-alphaL+0.5*k			        #Posterior for parameter.
  
  iter <- iter + burnin
  for(i in 1:iter){
    
    partmeanL<-t(F)%*%Kappa
    for(j in 1:k){
      partmean2L<-t(F)%*%omega[,j]*mu[j]
      Fp<-diag(sqrt(omega[,j]))%*%F
      Lerr<-t(Fp)%*%rnorm(n)+rnorm(f)*sqrt(d)
      FFT<-crossprod(Fp)
      p2<-partmeanL[,j]-partmean2L+Lerr
      pLinv<-diag(d)+FFT
      pL<-solve(pLinv,p2)  
      L[,j]<-pL          		#Update L by column.
    }
    
    #Update d.
    betad<-betaL+0.5*rowSums(L^2)
    d<-rgamma(f,alphad,betad)
    
    partmeanF<-Kappa%*%t(L)
    for(l in 1:n){
      partmean2F<-t(omega[l,]%*%(mu*t(L)))
      Lp <-diag(sqrt(omega[l,]))%*%t(L)
      p1F<-diag(f)+crossprod(Lp)
      Ferr<-t(Lp)%*%rnorm(k)+rnorm(f)
      p2F<-partmeanF[l,]-partmean2F+Ferr
      pF <-solve(p1F,p2F)
      F[l,]<-pF				#Update F by row.
    }
    
    #Combine F and L
    FL<- F%*%L
    
    muerr  <-rnorm(k)/sqrt(colSums(omega))
    part1mu<-colSums(Kappa-omega*FL)
    mu     <-muind*(part1mu/colSums(omega)+muerr)		#Update mu
    
    
    #new estimate of link function.
    link<-FL+t(matrix(mu,k,n))
    muuse<-colMeans(link)
    if(intercept == TRUE) {decomp = scale(link,center=TRUE,scale=FALSE)} else {decomp = link}
    
    #pg random variable.	
    omega<-matrix(rpg(n*k,h=as.numeric(maxtrial),z=link),n,k)
    
    #Rotation corrected model. 
    use<-svds(decomp,k=f,nu=f,nv=f)
    mcor<-cor(use$v,keep$v)              #This is included to minimise shuffling between factors with similar proportions of variance explained.
    if(n < k) {mcor<-cor(use$u,keep$u) }    #Fix for if permutation results in same factor picked twice or more.
    #To ensure the same factor is not picked twice, do permutation fix as a for loop.
    permuteind<-0
    for(fo in 1:(f-1)){
    permuteind[fo] <- which(abs(mcor[fo,]) == max(abs(mcor[fo,])))
    mcor[,permuteind[fo]]<-0
    }
    permuteind[f]<-(1:f)[-which(1:f %in% permuteind[1:(f-1)])]

    signL<-sign(diag(mcor[permuteind,]))   #This is included to stop bi-modality due to sign invariance.

    Loadingp1	    <-use$v[,permuteind]%*%diag(signL)	   #Fixing sign invariance.	
    
    if(i > burnin){
      Score[[i-burnin]]          <-use$u[,permuteind]%*%diag(signL)*sqrt(n-1)
      eigenvalp1	               <-use$d[permuteind]/sqrt((n-1))
      Loading[[i-burnin]]        <-Loadingp1%*%diag(eigenvalp1)
      eigenval[[i-burnin]]       <-eigenvalp1^2/(k-1) #Note in simulation are constrained should be divide by k -1.
      main.effect[[i-burnin]]    <-muind*muuse
      Link.fun[[i-burnin]]       <-link
      Devnew[i-burnin]	   <-Binom.Deviance.sat(y=y,link=link,maxtrial=as.numeric(maxtrial))
    }
    
  }
  
  x<-list(Loading,Score,eigenval,Devnew,Link.fun)
  names(x)<-c('Loading','Score','var.loading','Deviance','link')
  if(intercept == TRUE){ 
    x  <- c(x,list(main.effect))
    names(x)[length(x)]<-'main.effect'
  }
  return(x)
}

###########################
#Police data specific function 

#Version 2: PG PCA assuming structured patterns (factor loadings) are equal across ethnicity.
Binom.PGpolicePCA<-function(y,f,maxtrial,iter,burnin,intercept,alphaL,betaL){
  library(BayesLogit)
  library(rARPACK)
  n<-dim(y)[1]
  k<-dim(y)[2]
  
  Kappa<-y-0.5*maxtrial	     				  #corrected y
  keep<-svds(scale(0.5*Kappa[,0.5*k+1:(0.5*k)]+0.5*Kappa[,0.5*k+1:(0.5*k)],center=TRUE,scale=FALSE),k=f,nv=f,nu=f)
  
  F<-matrix(rnorm(n*f),n,f)    				  #Initial Factors
  L<-matrix(0,f,k) 					  #Initial Loadings
  
  #Initial estimate of PG rv.
  omega<-matrix(rpg(num=n*k, h=as.numeric(maxtrial), z=0.0),n,k) 
  
  if(intercept == TRUE) {muind<-1} else {muind<-0}
  mu<-muind*rep(0,k)
  
  Devnew<-0					#Initial deviance.
  
  Loading<-list()
  Score<-list()
  eigenval<-list()
  main.effect<-list()
  Link.fun<-list()
  mumean        <-list()
  mu.var        <-list()
  d<-rgamma(f,alphaL,betaL)			#Initial precisions for factor loadings.
  d<-sort(d)
  alphad<-alphaL+0.5*k			        #Posterior for parameter.
  tau.mu  <-rgamma(1,0.1,0.1) #Initial precision for inverse of sigma^2_u
  mubar   <-rep(0,2)          #mean of mean within groups.
  
  
  
  iter <- iter + burnin
  for(i in 1:iter){
    
    partmeanL<-t(F)%*%Kappa
    for(j in 1:(0.5*k)){
      partmean2L<-t(F)%*%omega[,j]*mu[j]+t(F)%*%omega[,j+0.5*k]*mu[j+0.5*k]
      Fp1<-diag(sqrt(omega[,j]))%*%F 
      Fp2<-diag(sqrt(omega[,j+0.5*k]))%*%F 
      Lerr<-t(Fp1)%*%rnorm(n)+t(Fp2)%*%rnorm(n)+rnorm(f)*sqrt(d)
      FFT<-crossprod(Fp1) + crossprod(Fp2)
      p2<-partmeanL[,j]+partmeanL[,j+0.5*k]-partmean2L+Lerr
      pLinv<-diag(d)+FFT
      pL<-solve(pLinv,p2)  
      L[,j]<-pL          		#Update L by column.
    }
    L[,0.5*k + 1:(0.5*k)]<-L[,1:(0.5*k)]
    
    #Update d.
    betad<-betaL+0.25*rowSums(L^2) #Half because L is doubled.
    d<-rgamma(f,alphad,betad)
    
    partmeanF<-Kappa%*%t(L)
    for(l in 1:n){
      partmean2F<-t(omega[l,]%*%(mu*t(L)))
      Lp <-diag(sqrt(omega[l,]))%*%t(L)
      p1F<-diag(f)+crossprod(Lp)
      Ferr<-t(Lp)%*%rnorm(k)+rnorm(f)
      p2F<-partmeanF[l,]-partmean2F+Ferr
      pF <-solve(p1F,p2F)
      F[l,]<-pF				#Update F by row.
    }
    
    #Combine F and L
    FL<- F%*%L
    
    mubaruse<-rep(mubar,each=0.5*k)   
    muerr  <-rnorm(k)/sqrt(colSums(omega)+tau.mu)
    part1mu<-colSums(Kappa-omega*FL)+mubaruse*tau.mu
    mu     <-muind*(part1mu/(colSums(omega)+tau.mu)+muerr)		#Update mu
    
    
    tau.mu<-rgamma(1,0.5*k,0.5*sum( (mu-mubaruse)^2) )    #Update tau.mu
    
    #Update mean of means
    mubar[1] <- mean(mu[1:(0.5*k)]) + rnorm(1)/sqrt(0.5*k*tau.mu)
    mubar[2] <- mean(mu[0.5*k+1:(0.5*k)]) + rnorm(1)/sqrt(0.5*k*tau.mu)
    
    
    #new estimate of link function.
    link<-FL+t(matrix(mu,k,n))
    muuse<-colMeans(link)
    mumean.use <-c(mean(muuse[1:(0.5*k)]),mean(muuse[0.5*k+1:(0.5*k)])) 
    if(intercept == TRUE) {decomp = scale(link,center=TRUE,scale=FALSE)} else {decomp = link}
    #decomp merging assuming equivalence.
    decomp  <-0.5*(decomp[,1:(0.5*k)]+decomp[,0.5*k+1:(0.5*k)])
    
    
    #pg random variable.	
    omega<-matrix(rpg(n*k,h=as.numeric(maxtrial),z=link),n,k)
    

    #Rotation corrected model. 
    use<-svds(decomp,k=f,nu=f,nv=f)
    mcor<-cor(use$v,keep$v)              #This is included to minimise shuffling between factors with similar proportions of variance explained.
    if(n < k) {mcor<-cor(use$u,keep$u) }    #Fix for if permutation results in same factor picked twice or more.
    #To ensure the same factor is not picked twice, do permutation fix as a for loop.
    permuteind<-0
    for(fo in 1:(f-1)){
      permuteind[fo] <- which(abs(mcor[fo,]) == max(abs(mcor[fo,])))
      mcor[,permuteind[fo]]<-0
    }
    permuteind[f]<-(1:f)[-which(1:f %in% permuteind[1:(f-1)])]
    
    signL<-sign(diag(mcor[permuteind,]))   #This is included to stop bi-modality due to sign invariance.
    
    Loadingp1	    <-use$v[,permuteind]%*%diag(signL)	   #Fixing sign invariance.    
    
    if(i > burnin){
      Score[[i-burnin]]          <-use$u[,permuteind]%*%diag(signL)*sqrt(n-1)
      eigenvalp1	               <-use$d[permuteind]/sqrt((n-1))
      Loading[[i-burnin]]        <-Loadingp1%*%diag(eigenvalp1)
      eigenval[[i-burnin]]       <-eigenvalp1^2/(0.5*k-1)
      main.effect[[i-burnin]]    <-muind*muuse
      Link.fun[[i-burnin]]       <-link
      mumean[[i-burnin]]        <-mumean.use
      mu.var[[i-burnin]]        <-sum( (muuse - rep(mumean.use,each=0.5*k))^2)/(k-1)
      Devnew[i-burnin]	   <-Binom.Deviance.sat(y=y,link=link,maxtrial=as.numeric(maxtrial))
    }
    
  }
  
  x<-list(Loading,Score,eigenval,Devnew,Link.fun,mumean,mu.var)
  names(x)<-c('Loading','Score','var.loading','Deviance','link','mu.mean','mu.var')
  if(intercept == TRUE){ 
    x  <- c(x,list(main.effect))
    names(x)[length(x)]<-'main.effect'
  }
  return(x)
}


####################################################################################
##################################################################################
#The function Binom.PGFA fits the factor analysis model.

#Inputs required are: 
#y:        	The matrix of binomial successes. Please store as class matrix.
#f:        	The number of components to be extracted.
#maxtrial:	The matrix of binomial trial sizes. Please store as class matrix.
#iter:     	The number of iterations after burnin.
#burnin:   	The number of iterations to be discarded.
#intercept 	TRUE/FALSE statement to specify whether an intercept should be included.
#alphaL, betaL  Parameter values for inverse-gamma priors on loading variances.
#alpha, beta	Parameter values for inverse-gamma priors on error variances.    


#Outputs are:
#Loading:	A list with entry i corresponding to the identified loading matrix estimate at iteration i.
#Score:		A list with entry i corresponding to the identified score matrix estimate at iteration i.	
#var.unique:	A list with entry i corresponding to the vector of estimated error variances at iteration i.
#var.loading:	A list with entry i corresponding to the vector of estimated loading variances at iteration i.
#Deviance:	A list containing the estimated deviance of the fit at iteration i.
#link:		A list with entry i corresponding to the matrix of estimated realisation of the logit link at iteration i.
#main.effect:	If intercept was set to TRUE, a list with entry i corresponding to the estimated of observed variable intercepts at iteration i.
#mu.mean: 
#mu.var:

#Version 1: PG FA assuming structured patterns (factor loadings) are not equal across ethnicity.
Binom.PGFA<-function(y,f,maxtrial,iter,burnin,intercept,alphaL,betaL,alpha,beta){
  library(BayesLogit)
  library(rARPACK)
  n<-dim(y)[1]
  k<-dim(y)[2]
  
  Kappa<-y-0.5*maxtrial	     		  #corrected y
  keep<-svds(scale(Kappa,center=TRUE,scale=FALSE),k=f,nv=f,nu=f)
  
  F<-matrix(rnorm(n*f),n,f)    		  #Initial Factors   
  L<-matrix(0,f,k)    		  	  #Initial L	   (Set to zero matrix.)	
  
  tau<-rgamma(k,alpha,beta)        	  #Initial error precision
  FL <-F%*%L
  
  if(intercept == TRUE) {thetaind<-1} else {thetaind<-0}
  theta<-rep(0,k)
  
  #Initial estimate of PG rv.
  omega<-matrix(rpg(num=n*k, h=as.numeric(maxtrial), z=0.0),n,k)
  
  #Initial storage of estimates.
  Loading       <-list()
  Score         <-list()
  main.effect   <-list()
  Link.fun      <-list()
  sigma2.unique <-list()        
  sigma2.loading<-list()
  
  d<-rgamma(f,alphaL,betaL)		#Initial precisions for factor loadings.
  alphad<-alphaL+0.5*k				#alpha parameters.
  alpha <-alpha+0.5*n				#alpha parameter for error variance.
  thetamat<-t(matrix(theta,k,n))
  
  Devnew<-0					#Initial deviance.
  
  
  iter <- iter + burnin
  for(i in 1:iter){
    
    taumat<-t(matrix(tau,k,n))
    
    p1<-matrix(rnorm(n*k),n,k)
    varmat<-1/(omega+taumat)
    Psimean<-varmat*(Kappa+taumat*(FL+thetamat))
    Psi<-Psimean+p1*sqrt(varmat)			#Update Psi
    
    
    
    FFT<-crossprod(F)
    p1<-t(F)%*%t(matrix(rnorm(n*k),k,n)*sqrt(tau))+sqrt(d)*matrix(rnorm(f*k),f,k)
    p2<-t(F)%*%(Psi-thetamat)%*%diag(tau)
    pL<-p1+p2
    for(j in 1:k){
      varLinv<-diag(d)+FFT*tau[j]
      L[,j]<-solve(varLinv,pL[,j])         		#Update L.
    }
    
    #Update d.
    betad<-betaL+0.5*rowSums(L^2)
    d<-rgamma(f,alphad,betad)
    
    
    p1F<-diag(f)+L%*%diag(tau)%*%t(L)
    errF<-matrix(rnorm(f*n),f,n)+L%*%diag(sqrt(tau))%*%matrix(rnorm(k*n),k,n)
    p2F<-L%*%diag(tau)%*%t(Psi-thetamat)+errF
    F<-t(solve(p1F,p2F))				#Update F
    
    
    
    #-FL
    theta<-(colMeans(Psi)+rnorm(k)/sqrt(tau*n))*thetaind				#update theta (really mu)
    FL<- F%*%L
    
    #New estimate of PG rv.
    omega<-matrix(rpg(num=n*k, h=as.numeric(maxtrial), z=Psi),n,k)			#Update omega (Set to Psi here)
    
    thetamat<-t(matrix(theta,k,n))
    err<-colSums((Psi-FL-thetamat)^2)
    tau<-rgamma(k,alpha,beta+0.5*err) 						#Update tau 
    
    #Apply decomposition to mu to FL
    
    decomp <- FL + thetamat
    if(intercept == TRUE) {
      mu.use <- colMeans(decomp)
      decomp <-scale(decomp,center=TRUE,scale=FALSE)
    }
    
    #Rotation corrected model. 
    use<-svds(decomp,k=f,nu=f,nv=f)
    mcor<-cor(use$v,keep$v)              #This is included to minimise shuffling between factors with similar proportions of variance explained.
    if(n < k) {mcor<-cor(use$u,keep$u) }    #Fix for if permutation results in same factor picked twice or more.
    #To ensure the same factor is not picked twice, do permutation fix as a for loop.
    permuteind<-0
    for(fo in 1:(f-1)){
      permuteind[fo] <- which(abs(mcor[fo,]) == max(abs(mcor[fo,])))
      mcor[,permuteind[fo]]<-0
    }
    permuteind[f]<-(1:f)[-which(1:f %in% permuteind[1:(f-1)])]
    
    signL<-sign(diag(mcor[permuteind,]))   #This is included to stop bi-modality due to sign invariance.
    
    Loadingp1	    <-use$v[,permuteind]%*%diag(signL)	   #Fixing sign invariance.    
  
    if(i > burnin){
      Score[[i-burnin]]          <-use$u[,permuteind]%*%diag(signL)*sqrt(n-1)
      eigenvalp1	               <-use$d[permuteind]/sqrt((n-1))
      Loading[[i-burnin]]        <-Loadingp1%*%diag(eigenvalp1)
      sigma2.loading[[i-burnin]]<-eigenvalp1^2/(k-1)
      main.effect[[i-burnin]]		<-thetaind*mu.use
      Link.fun[[i-burnin]]   		<-Psi
      sigma2.unique[[i-burnin]]	<-1/tau   
      Devnew[i]<-Binom.Deviance.sat(y=y,link=Psi,maxtrial=as.numeric(maxtrial))
    }
    
  }
  
  x<-list(Loading,Score,Link.fun,sigma2.unique,sigma2.loading)
  names(x)<-c('Loading','Score','link','var.unique','var.loading')
  if(intercept == TRUE){ 
    x  <- c(x,list(main.effect))
    names(x)[length(x)]<-'main.effect'
  }
  return(x)
}


################################################################################################################
#Version 2: PG FA assuming structured patterns (factor loadings) are equal across ethnicity.

Binom.PGpoliceFA<-function(y,f,maxtrial,iter,burnin,intercept,alphaL,betaL,alpha,beta){
  library(BayesLogit)
  library(rARPACK)
  n<-dim(y)[1]
  k<-dim(y)[2]  #dataset has repeated variables.
  
  Kappa<-y-0.5*maxtrial	     		  #corrected y
  keep<-svds(scale(0.5*Kappa[,0.5*k+1:(0.5*k)]+0.5*Kappa[,0.5*k+1:(0.5*k)],center=TRUE,scale=FALSE),k=f,nv=f,nu=f)
  
  F<-matrix(rnorm(n*f),n,f)    		  #Initial Factors   
  L<-matrix(0,f,k)    		  	  #Initial L	   (Set to zero matrix.)	
  
  tau<-rgamma(k,alpha,beta)        	  #Initial error precision
  FL <-F%*%L
  
  if(intercept == TRUE) {thetaind<-1} else {thetaind<-0}
  theta<-rep(0,k)
  
  #Initial estimate of PG rv.
  omega<-matrix(rpg(num=n*k, h=as.numeric(maxtrial), z=0.0),n,k)
  
  #Initial storage of estimates.
  Loading       <-list()
  Score         <-list()
  main.effect   <-list()
  Link.fun      <-list()
  sigma2.unique <-list()        
  sigma2.loading<-list()
  mumean        <-list()
  mu.var        <-list()
  
  d<-rgamma(f,alphaL,betaL)		#Initial precisions for factor loadings.
  alphad<-alphaL+0.5*k				#alpha parameters.
  alpha <-alpha+0.5*n				#alpha parameter for error variance.
  thetamat<-t(matrix(theta,k,n))
  tau.mu  <-rgamma(1,0.1,0.1) #Initial precision for sigma^2_u
  mubar   <-rep(0,2)          #mean of mean within groups.
  
  Devnew<-0					#Initial deviance.
  
  
  iter <- iter + burnin
  for(i in 1:iter){
    
    taumat<-t(matrix(tau,k,n))
    
    p1<-matrix(rnorm(n*k),n,k)
    varmat<-1/(omega+taumat)
    Psimean<-varmat*(Kappa+taumat*(FL+thetamat))
    Psi<-Psimean+p1*sqrt(varmat)			#Update Psi (The link function)
    
    
    
    FFT<-crossprod(F)
    p1<-t(F)%*%t(matrix(rnorm(n*k),k,n)*sqrt(tau))
    p1b<-sqrt(d)*matrix(rnorm(f*k),f,0.5*k)
    p2<-t(F)%*%(Psi-thetamat)%*%diag(tau)
    pL<-p1+p2
    for(j in 1:(0.5*k)){
      varLinv<-diag(d)+FFT*(tau[j]+tau[j+0.5*k])
      L[,j]<-solve(varLinv,pL[,j]+pL[,j+0.5*k]+p1b[,j])         		#Update L.
    }
    L[,0.5*k + 1:(0.5*k)]<-L[,1:(0.5*k)]
    
    #Update d.
    betad<-betaL+0.25*rowSums(L^2)    #0.25 not 0.5 due to double up.
    d<-rgamma(f,alphad,betad)
    
    
    p1F<-diag(f)+L%*%diag(tau)%*%t(L)
    errF<-matrix(rnorm(f*n),f,n)+L%*%diag(sqrt(tau))%*%matrix(rnorm(k*n),k,n)
    p2F<-L%*%diag(tau)%*%t(Psi-thetamat)+errF
    F<-t(solve(p1F,p2F))				#Update F
    
    
    
    #-FL
    mubaruse<-rep(mubar,each=0.5*k)
    theta<-( (colMeans(Psi)*n*tau+mubaruse*tau.mu)/(tau.mu+n*tau)+rnorm(k)/sqrt(tau*n+tau.mu))*thetaind				#update theta (really mu)
    FL<- F%*%L
    
    tau.mu<-rgamma(1,0.5*k,0.5*sum( (theta-mubaruse)^2))    #Update tau.mu
    
    #Update mean of means
    mubar[1] <- mean(theta[1:(0.5*k)]) + rnorm(1)/sqrt(0.5*k*tau.mu)
    mubar[2] <- mean(theta[0.5*k+1:(0.5*k)]) + rnorm(1)/sqrt(0.5*k*tau.mu)
    
    
    #New estimate of PG rv.
    omega<-matrix(rpg(num=n*k, h=as.numeric(maxtrial), z=Psi),n,k)			#Update omega (Set to Psi here)
    
    thetamat<-t(matrix(theta,k,n))
    err<-colSums((Psi-FL-thetamat)^2)
    tau<-rgamma(k,alpha,beta+0.5*err) 						#Update tau 
    
    
    #Apply decomposition to mu + FL
    
    decomp <- FL + thetamat
    if(intercept == TRUE) {
      mu.use <- colMeans(decomp)
      mumean.use <-c(mean(mu.use[1:(0.5*k)]),mean(mu.use[0.5*k+1:(0.5*k)])) 
      decomp <-scale(decomp,center=TRUE,scale=FALSE)
    }
    #decomp merging assuming equivalence.
    decomp  <-0.5*(decomp[,1:(0.5*k)]+decomp[,0.5*k+1:(0.5*k)])
    
    
    #Rotation corrected model. 
    use<-svds(decomp,k=f,nu=f,nv=f)
    mcor<-cor(use$v,keep$v)              #This is included to minimise shuffling between factors with similar proportions of variance explained.
    if(n < k) {mcor<-cor(use$u,keep$u) }    #Fix for if permutation results in same factor picked twice or more.
    #To ensure the same factor is not picked twice, do permutation fix as a for loop.
    permuteind<-0
    for(fo in 1:(f-1)){
      permuteind[fo] <- which(abs(mcor[fo,]) == max(abs(mcor[fo,])))
      mcor[,permuteind[fo]]<-0
    }
    permuteind[f]<-(1:f)[-which(1:f %in% permuteind[1:(f-1)])]
    
    signL<-sign(diag(mcor[permuteind,]))   #This is included to stop bi-modality due to sign invariance.
    
    Loadingp1	    <-use$v[,permuteind]%*%diag(signL)	   #Fixing sign invariance.
    
    
    #Storing results.	
    if(i > burnin){
      Score[[i-burnin]]          <-use$u[,permuteind]%*%diag(signL)*sqrt(n-1)
      eigenvalp1	               <-use$d[permuteind]/sqrt((n-1))
      Loading[[i-burnin]]        <-Loadingp1%*%diag(eigenvalp1)
      sigma2.loading[[i-burnin]]<-eigenvalp1^2/(0.5*k-1)
      main.effect[[i-burnin]]		<-thetaind*mu.use
      Link.fun[[i-burnin]]   		<-Psi
      sigma2.unique[[i-burnin]]	<-1/tau   
      mumean[[i-burnin]]        <-mumean.use
      mu.var[[i-burnin]]        <-sum( (mu.use - rep(mumean.use,each=0.5*k))^2)/(k-1)
      Devnew[i-burnin]<-Binom.Deviance.sat(y=y,link=Psi,maxtrial=as.numeric(maxtrial))
    }
    
  }
  
  x<-list(Loading,Score,Link.fun,sigma2.unique,sigma2.loading,mumean,mu.var,Devnew)
  names(x)<-c('Loading','Score','link','var.unique','var.loading','mu.mean','mu.var','Deviance')
  if(intercept == TRUE){ 
    x  <- c(x,list(main.effect))
    names(x)[length(x)]<-'main.effect'
  }
  return(x)
}

#################################################################################################################
#Deviance calculations.
Binom.Deviance.int<-function(y,link,maxtrial){
  yhat<-maxtrial*(1+exp(-link))^-1
  yMeans<-colMeans(y)
  k<-dim(y)[2];n<-dim(y)[1]
  yMmat<-t(matrix(yMeans,k,n))
  
  y1<-y[y>0]
  yMmat1<-yMmat[y>0]
  yhat1<-yhat[y>0]
  y2<-y[y<maxtrial]
  yhat2<-yhat[y<maxtrial]
  yMmat2<-yMmat[y<maxtrial]
  
  p1<-y1*log(yMmat1/yhat1)
  p2<-(maxtrial[y<maxtrial]-y2)*log((maxtrial[y<maxtrial]-yMmat2)/(maxtrial[y<maxtrial]-yhat2))
  
  D<- -2*(sum(p1)+sum(p2))
  return(D)
}


#Compared to saturated model.
Binom.Deviance.sat<-function(y,link,maxtrial){
  yhat<-maxtrial*(1+exp(-link))^-1
  y1<-y[y>0]
  yhat1<-yhat[y>0]
  y2<-y[y<maxtrial]
  yhat2<-yhat[y<maxtrial]
  
  p1<-y1*log(yhat1/y1)
  p2<-(maxtrial[y<maxtrial]-y2)*log((maxtrial[y<maxtrial]-yhat2)/(maxtrial[y<maxtrial]-y2))
  
  D<- -2*(sum(p1)+sum(p2))
  return(D)
}

