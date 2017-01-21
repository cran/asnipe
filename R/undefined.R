	
infer_graph_from_datastream_mmVB = function(DATA,verbose=T){
	DATA<-as.matrix(DATA)
	prior_on_K=length(DATA[,2])-1
	
	
	gmmresults = gmmvar(DATA[,1],prior_on_K,gmmoptions=list(cyc=50,tol=1e-5,
		cov='diag',display=as.numeric(verbose),testcovmat=1));
	
	gmm<-gmmresults[[1]][[1]]
	
	centroids=c()
	
	for (i in c(1:length(gmm$post))){
		centroids=c(centroids,gmm$post[[i]]$Norm_Mu)
	}
	
	centroids=mlmatrix(centroids)
	Y = t(gmm$pjgx)
	
	
	Y_hard=get_hard_incidence_matrix(Y)
	if(!is.matrix(Y_hard)){
		Y_hard=mlmatrix(Y_hard)
	}
    	Y_hard = t(Y_hard)
	

	active_clusters = colSums(Y_hard)>0
	Y_hard = Y_hard[,active_clusters]
	if(!is.matrix(Y_hard)){
		Y_hard<-as.matrix(Y_hard)
	}
	centroids = centroids[,active_clusters]
	
	EVENT_TIMES = (get_T_matrix(Y_hard,DATA))
	
	B_hard_incidence_matrix = build_B_from_Y(DATA,t(Y_hard))
	A_hard_cooccurences = get_coocurences_in_bipartite_graph(B_hard_incidence_matrix)
	
	
	output = list('centroids'=centroids,'A_hard_cooccurences'=A_hard_cooccurences,'B_hard_incidence_matrix'=B_hard_incidence_matrix,
	'EVENT_TIMES'=EVENT_TIMES,'Y_hard'=Y_hard)
	
	igfdresult=list(output,gmm)
	return(igfdresult)
}
	
gmmvar = function(data,K,gmmoptions=list(cyc=50,tol=1e-5,
		cov='diag',display=1,testcovmat=1)){
	
	data=matrix(data)

	eps=2^-52
	
	oldFrEn=1
	MIN_COV=eps
	
	N=dim(data)[1]
	ndim=dim(data)[2]
	
	A=mlmatrix(length(K));
	FrEn=t(matrix(rep(0,A)))
	
	U=array(NA, dim=c(gmmoptions[[1]],A,20))
	
	gmm=as.list(rep(NA,A))
	gmminit=as.list(rep(NA,A))
	
	# initialise
	for (a in 1:A){
		gmm[[a]]=gmmvarinit(data,K[a],gmmoptions)
		# save for resetting during convergence problems
		gmminit[[a]]$post=gmm[[a]]$post
	
	}
	
	
	FrEn=rep(NA,A)
	for(a in c(1:A)){# iterating over q-mixture components
		
		for (cyc in c(1:gmmoptions$cyc)){
		
			# The E-Step, i.e. estimating Q(hidden variables)
			gmm[[a]]$pjgx=estep(gmm[[a]],data,N,ndim,K[a]);
			
			# The M-Step, i.e. estimating Q(model parameters)
			gmm[[a]]=mstep(gmm[[a]],data,N,ndim,K[a]);
			
			# computation of free energy 
			feeenerresult=freeener(gmm[[a]],data,N,ndim,K[a]);
			FrEn[a]=feeenerresult[[1]]
			U[cyc,a,1:3]=feeenerresult[[2]]
			
			# check change of Free Energy
			if (abs((FrEn[a] - oldFrEn)/oldFrEn*100) < gmmoptions$tol){ 
				break() 
			}else{
				oldFrEn=FrEn[a];
			}
			
			 if (gmmoptions$display==1){
				print(sprintf('Iteration %d ; Free-Energy = %f',cyc,FrEn[a])); 
			}
		
		}					
		if (gmmoptions$display==1){
			print(sprintf('Model %d: %d kernels, %d dimensions, %d data samples',
			a,K[a],ndim,N));
			print(sprintf('Final Free-Energy (after %d iterations)  = %f',
			cyc,FrEn[a])); 
		}
	
	}				
	
	pFrEn=matrix(NA,nrow=1,ncol=A)
	
	for (a in c(1:A)){
		gmm[[a]]$pa=1/sum(exp(-FrEn+FrEn[a]))
		pFrEn[a]=gmm[[a]]$pa
	}
	
	pFrEn=c(pFrEn[]);
	FrEn=cbind(FrEn,pFrEn);
	
	gmmresults=list(gmm,FrEn,U)
	return(gmmresults)
}

gmmvarinit=function(data,K,options){
	# initialises the gaussian mixture model for Variational GMM algorithm
	data=matrix(data)
	eps=eps=2^-52
	MIN_COV=eps
	
	gmm=list()
	gmm$K=K
	N=dim(data)[1]
	ndim=dim(data)[2]
	
	midscale=t(mlmatrix(median(data)))
	drange=mlmatrix(max(data)-min(data))
	# educated guess with scaling
	
	# define P-priors
	defgmmpriors<-list()
	
	defgmmpriors$Dir_alpha=mlmatrix(rep(1,K));
	defgmmpriors$Norm_Mu=(midscale);
	defgmmpriors$Norm_Cov=matrix(diag(drange^2))
	defgmmpriors$Norm_Prec=solve(defgmmpriors$Norm_Cov);
	defgmmpriors$Wish_B=matrix(diag(drange))
	defgmmpriors$Wish_iB=solve(defgmmpriors$Wish_B);
	defgmmpriors$Wish_alpha=mlmatrix(ndim+1);
	defgmmpriors$Wish_k=mlmatrix(ndim);
	
	gmm$priors=defgmmpriors
	
	
	# initialise posteriors
	
	# sample mean from data
	ndx=(floor(mlmatrix(runif(K))*N+1))
	Mu=t(matrix(data[ndx,]))
	
	# sample precision
	alpha=gmm$priors$Wish_alpha*2
	Prec=wishartrnd(alpha,gmm$priors$Wish_iB,K);
	# sample weights
	kappa=gmm$priors$Dir_alpha;
	# assign values
	
	gmm$post<-as.list(rep(NA,K))
	for (k in c(1:K)){
		gmm$post[[k]]<-list()
		# mean posterior
		gmm$post[[k]]$Norm_Prec=gmm$priors$Norm_Prec
		gmm$post[[k]]$Norm_Cov=gmm$priors$Norm_Cov
		gmm$post[[k]]$Norm_Mu=mlmatrix(Mu[,k])
		# covariance posterior
		gmm$post[[k]]$Wish_alpha=alpha
		gmm$post[[k]]$Wish_iB=drop(Prec[,,k])/alpha #gmm.priors.Wish_B ?
		gmm$post[[k]]$Wish_B=solve(gmm$post[[k]]$Wish_iB)
		# weights posterior
		gmm$post[[k]]$Dir_alpha=kappa[k]
		
	}
	
		


	return(gmm)
}
	
estep=function(gmm,data,N,ndim,K){
	
	Dir_alpha=mlmatrix(rep(NA,K))
	for (i in c(1:K)){
		Dir_alpha[,i]=gmm$post[[i]]$Dir_alpha
	}
	PsiDiralphasum=digamma(sum(Dir_alpha))
	gmm$pjgx=matrix(NA,nrow=N,ncol=K)
	for (k in c(1:K)){
		qp=gmm$post[[k]];	# for ease of referencing
		
		ldetWishB=0.5*log(det(qp$Wish_B))
		PsiDiralpha=digamma(qp$Dir_alpha)
		PsiWish_alphasum=0
		for (d in c(1:ndim)){
			PsiWish_alphasum=PsiWish_alphasum+
			digamma(qp$Wish_alpha+0.5-d/2)
		}
		PsiWish_alphasum=mlmatrix(PsiWish_alphasum*0.5)
		
		dist=mdist(data,qp$Norm_Mu,qp$Wish_iB%*%qp$Wish_alpha)
		
		NormWishtrace=0.5*sum((qp$Wish_alpha%*%qp$Wish_iB%*%qp$Norm_Cov))
		
		gmm$pjgx[,k]=exp(PsiDiralpha-PsiDiralphasum+as.numeric(PsiWish_alphasum)-ldetWishB+
		dist-NormWishtrace-ndim/2*log(2*pi))
	
	}
	
	eps=2^-52
	# normalise posteriors of hidden variables.
	gmm$pjgx=gmm$pjgx # +eps
	col_sum=matrix(rowSums(gmm$pjgx))
	gmm$pjgx=gmm$pjgx/(col_sum%*%matrix(1,nrow=1,ncol=K))
	
	pjgx=gmm$pjgx
	
	return(pjgx)
}
	
	
mstep=function(gmm,data,N,ndim,K){
	
	pr=gmm$priors			# model priors
	
	gammasum=colSums(gmm$pjgx)
	
	for (k in c(1:K)){
		qp=gmm$post[[k]]		# temporary structure (q-distrib);
		
		# Update posterior Normals
		postprec=gammasum[k]%*%qp$Wish_alpha%*%qp$Wish_iB+pr$Norm_Prec
		postvar=solve(postprec)
		weidata=t(data)%*%gmm$pjgx[,k] #unnormalised sample mean
		
		Norm_Mu=mlmatrix(postvar%*%(qp$Wish_alpha%*%qp$Wish_iB%*%weidata+ 
		pr$Norm_Prec%*%pr$Norm_Mu))
		
		Norm_Prec=postprec
		Norm_Cov=postvar
		
		#Update posterior Wisharts
		Wish_alpha=0.5*gammasum[k]+pr$Wish_alpha
		dist=data-matrix(1,nrow=N,ncol=1)%*%t(Norm_Mu)
		
		sampvar=matrix(0,nrow=ndim,ncol=ndim)
		
		for (n in c(1:ndim)){
			sampvar[n,]=colSums(matrix(gmm$pjgx[,k]*dist[,n])%*%matrix(1,nrow=1,ncol=ndim)*dist)
		}
		
		Wish_B=0.5%*%(sampvar+gammasum[k]%*%Norm_Cov)+pr$Wish_B
		Wish_iB=solve(Wish_B)
		
		# Update posterior Dirichlet
		Dir_alpha=gammasum[k]+pr$Dir_alpha[k]
		
		gmm$post[[k]]$Norm_Mu=Norm_Mu
		gmm$post[[k]]$Norm_Prec=Norm_Prec
		gmm$post[[k]]$Norm_Cov=Norm_Cov
		gmm$post[[k]]$Wish_alpha=Wish_alpha
		gmm$post[[k]]$Wish_B=Wish_B
		gmm$post[[k]]$Wish_iB=Wish_iB
		gmm$post[[k]]$Dir_alpha=Dir_alpha
	}
	
	return(gmm)
	
}
	
freeener=function(gmm,data,N,ndim,K){
	
	KLdiv=0
	avLL=0
	pr=gmm$priors
	
	Dir_alpha=matrix(NA,nrow=1,ncol=K)
	for (i in c(1:K)){
		Dir_alpha[,i]=gmm$post[[i]]$Dir_alpha
	}
	
	Dir_alphasum=sum(Dir_alpha)
	PsiDir_alphasum=digamma(Dir_alphasum)
	ltpi=ndim/2*log(2*pi)
	gammasum=mlmatrix(colSums(gmm$pjgx))
	
	# entropy of hidden variables, which are not zero
	pjgxndx=gmm$pjgx[gmm$pjgx!=0]
	
	
	Entr=sum(sum(pjgxndx*log(pjgxndx)));
	
	for (k in c(1:K)){
		# average log-likelihood
		qp=gmm$post[[k]];		# for ease of referencing
		
		PsiDiralpha=digamma(qp$Dir_alpha)
		dist=mdist(data,qp$Norm_Mu,qp$Wish_iB%*%qp$Wish_alpha)
		NormWishtrace=0.5*sum(diag(qp$Wish_alpha%*%qp$Wish_iB%*%qp$Norm_Cov))
		
		ldetWishB=0.5*log(det(qp$Wish_B))
		PsiWish_alphasum=mlmatrix(0)
		for (d in 1:ndim){
			PsiWish_alphasum=PsiWish_alphasum+
			digamma(qp$Wish_alpha+0.5-d/2)
		}
		PsiWish_alphasum=0.5*PsiWish_alphasum
		
		avLL=avLL+gammasum[k]%*%(PsiDiralpha-PsiDir_alphasum-ldetWishB+
		PsiWish_alphasum-NormWishtrace-ltpi)+
		sum(gmm$pjgx[,k]%*%dist)
		
		# KL divergences of Normals and Wishart
		VarDiv=wishart_kl(qp$Wish_B,pr$Wish_B,qp$Wish_alpha,pr$Wish_alpha)
		MeanDiv=gauss_kl(qp$Norm_Mu,pr$Norm_Mu,qp$Norm_Cov,pr$Norm_Cov)
		KLdiv=KLdiv+VarDiv+MeanDiv
	}
	
	# KL divergence of Dirichlet
	KLdiv=KLdiv+dirichlet_kl(Dir_alpha,pr$Dir_alpha)
	
	FrEn=Entr-avLL+KLdiv
	U=c(Entr,-avLL,+KLdiv)
	
	freeenerresult=list(FrEn,U)
	
	return(freeenerresult)
}					
	
	
mdist = function (x,mu,C){
	
	d=dim(C)[1]
	if (dim(x)[1]!=d)  {x=t(x)}
	if (dim(mu)[1]!=d)  {mu=t(mu)}
	
	
	
	
	ndim=dim(x)[1]
	N=dim(x)[2]
	
	d=x-mu%*%matrix(1,nrow=1,ncol=N)
	
	Cd=C%*%d
	
	dist=matrix(0,nrow=1,ncol=N);
	
	for (l in c(1:ndim)){
		dist=dist+d[l,]*Cd[l,]
	}
	
	dist=-0.5*t(dist)
	
	return (dist)
}

wishartrnd=function(alpha,C,N){
	
	alpha=ceiling(alpha)
	ndim=dim(C)[1]
	
	m=matrix(0,nrow=ndim,ncol=1)
	y=sampgauss(m,C,alpha)
	
	x2=y%*%t(y)
	x=array(0,c(dim(matrix(x2))[1],dim(matrix(x2))[2],N))
	x[,,1]=x2;
	if(length(x)>1){
		for (i in c(2:N)){
			y=sampgauss(rep(0,ndim),C,alpha);
			x[,,i]=y%*%t(y)
		}
		}else{
		y=sampgauss(rep(0,ndim),C,alpha);
		x[,,1]=y%*%t(y)
	}
	return(x)
}

sampgauss=function(m,C,N){
	m=m[]
	ndim=dim(C)[1]
	
	if (dim(C)[2]!=ndim){
		x=c();
		print("Wrong specification calling sampgauss")
	}
	
	if (ndim==1){
		x=(m+as.numeric(C))%*%matrix(rnorm(N),nrow=1,ncol=N);
		return (x)
	}
	
	# check determinant of covariance matrix
	#if det(C)>1 | det(C)<=0, error('Covariance matrix determinant must be
	#0< |C| <= 1'); end;
	
	# generate zero mean/unit variance samples
	#e=randn(ndim,N);
	e=matrix(rnorm(ndim*N), ncol=ndim,nrow=N)
	# make sure they are unit variance;
	s=sd(t(e))
	
	for (i in c(1:ndim)){
		e[i,]=e[i,]/s[i]
	}
	
	# decompose cov-matrix
	S=chol(solve(C))
	x=solve(S)%*%e+m[,rep(1,N)]
	return (x)
}

get_hard_incidence_matrix=function(P){
	
	B <- apply(P,2,function(x) { as.integer(x == max(x))})
	
	return(B)
}

get_T_matrix=function(Y,DATA){
	
	K = dim(Y)[2];
	
	T = matrix(0,nrow=K,ncol=2)
	
	for (k in c(1:K)){
		first_obs = which(Y[,k]!=0)[1]
		last_obs = which(Y[,k]!=0)
		last_obs=last_obs[length(last_obs)]
		
		T[k,1] = DATA[first_obs,1]
		T[k,2] = DATA[last_obs,1]
	}
	sparseMatrix(i=which(t(T)!=0,arr.ind=T)[,1],j=which(t(T)!=0,arr.ind=T)[,2],x=t(T)[t(T!=0)])
	return(T)
	
}

get_coocurences_in_bipartite_graph = function(B){
	# assumes [K x N] incidence matrix - K clusters N individuals
	N = dim(B)[2]
	K = dim(B)[1]
	
	W = matrix(0,nrow=N,ncol=N)
	
	for (i in c(1:N-1)){
		for (j in c((i+1):N)){
			
			visitation_profile_i = B[,i]
			visitation_profile_j = B[,j]
			
			# to vectorise
			for (k in c(1:K)){
				W[i,j] = W[i,j] + min(c(visitation_profile_i[k],visitation_profile_j[k]))
				
			}
		}
	}
	
	W = W + t(W)
	return(W)
}

build_B_from_Y= function(DATA,Y){
	
	individuals = unique(DATA[,2])
	N_individuals = max(individuals)
	
	K = dim(Y)[1]
	
	
	B = matrix(0,nrow=K,ncol=N_individuals)
	
	for (i in c(1:length(individuals))){
		i_indices = DATA[,2] == individuals[i]
		
		if (sum(i_indices)!=1){
			if(is.vector(Y[,i_indices])){
				
				if(K==1){
					B[,individuals[i]] = sum(Y[,i_indices])
				}
			} else {
					B[,individuals[i]] = rowSums(Y[,i_indices])
			}
		} else {
			B[,individuals[i]] = Y[,i_indices]
		}
	}
	
	return(B)
}


wishart_kl=function(B_q,B_p,alpha_q,alpha_p){
		
	K=dim(B_p)[1]
	
	DBq=det(B_q)
	DBp=det(B_p)
	
	D=alpha_p%*%log(DBq%/%DBp)+alpha_q%*%(sum(diag((B_p%*%solve(B_q))-K)))
	for (k in c(1:K)){
		D=D+lgamma(alpha_p+0.5-0.5*k)-lgamma(alpha_q+0.5-0.5*k)
		D=D+(alpha_q-alpha_p)*digamma(alpha_q+0.5-0.5*k)
	}
	
	return(D)
}


gauss_kl=function(mu_q,mu_p,sigma_q,sigma_p){

	mu_q=mu_q[]
	mu_p=mu_p[]
	
	
	DSq=det(sigma_q)
	DSp=det(sigma_p)
	K=dim(sigma_q)[1]
	
	D=log(DSp/DSq)-K+sum(diag(sigma_q%*%solve(sigma_p)))+t(mu_q-mu_p)%*%solve(sigma_p)%*%(mu_q-mu_p)
	D=D*0.5;
	return(D)
}


dirichlet_kl=function(alpha_q,alpha_p){	
	
	K=length(alpha_q)
	
	aqtot=sum(alpha_q)
	aptot=sum(alpha_p)
	Psiqtot=digamma(aqtot)
	
	D=lgamma(aqtot)-lgamma(aptot)
	for (k in c(1:K)){
		D=D+lgamma(alpha_p[k])-lgamma(alpha_q[k])+
		+(alpha_q[k]-alpha_p[k])%*%(digamma(alpha_q[k])-Psiqtot)
		
	}
return(D)
}

mlmatrix<-function(x){
	x<-t(matrix(x))
	return(x)
}
