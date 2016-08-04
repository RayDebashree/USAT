###------------- Code for USAT for testing association of multiple continuous phenotypes with a single SNP -------------------
#
# Cite: "USAT: A Unified Score-based Association Test for Multiple Phenotype-Genotype Analysis". Genetic Epidemiology, 40(1):20-34, 2016.
#
# USAT: A linear combination of SSU and MANOVA score statistics
# T_MANOVA ~ chisq(k), where k = no of traits
# T_SSU ~ a*chisq(d) + b
# T_rho = rho*T_MANOVA + (1-rho)*T_SSU
# optimal rho found using grid search over [0,1]
# T_USAT = min_rho p_rho, where p_rho is pvalue of T_rho
# calculates the appropriate p-value of T_USAT using a one-dimensional numerical integration

##--------------------------------------------- Version 1.1 (dated April 7, 2016) --------------------------------------------------
# Corresponding Authors: Debashree Ray <debashr@umich.edu> ; Saonli Basu <saonli@umn.edu> 

############################################
library(CompQuadForm)
library(pracma)
library(survey)
library(minqa)

message("==============================")
message("     USAT v1.1 is loaded")
message("==============================")
message("For citation, please refer to:")
message("Ray et al.(2016) USAT: A Unified Score-based Association Test for Multiple")
message("    Phenotype-Genotype Analysis. Genetic Epidemiology, 40(1):20-34")

############################################
eps<-1e-13      # a threshold for checking non-positive values
# function for Moore-Penrose Inverse
mpinv <- function(A, eps = 1e-13, power=-1) {
        s<-eigen(A)
        e<-s$values
        V<-s$vectors
        e[e>eps] <- 1/(e[e > eps])^abs(power)
        return(V%*%diag(e)%*%t(V))
}

##----------------------- Begin: Functions borrowed from Dr. Baolin Wu's website -----------------------
## functions required for calculating higher order cumulants in order to do higher order moment matching (Wu and Pankow 2015)
#####
cum2mnc = function(kappa){
### convert cumulants to non-central moments
###    recursive formula produces as many cumulants as moments
###    References: Kenneth Lange: Numerical Analysis for Statisticians, 2nd ed. Page 15
  N = length(kappa)+1
  mc = rep(0, N); mc[1] = 1
  for(k in 1:(N-1)){
    mc[k+1] = sum(choose(k-1, 0:(k-1))*kappa[k:1]*mc[1:k])
  }
  return(mc[-1])
}
mnc2mc = function(mnc){
### convert non-central to central moments, uses recursive formula
  N = length(mnc)
  mc = rep(0,N); mc[1] = 0
  s1 = rep(c(1,-1), N)
  mnc = c(1,mnc)
  for(k in 1:(N-1)){
    mc[k+1] = sum( choose(k+1, 0:(k+1))*s1[(k+2):1]*mnc[1:(k+2)]*mnc[2]^((k+1):0) )
  }
  return(mc)
}
#### non-central chi-square cumulants
chisq.cum = function(k, lam, N){
### k: DF; lam: ncp
  ik = 1:N
  2^(ik-1)*gamma(ik)*(k+ik*lam)
}
## 1-DF chisq mix cumulants
chi1sqm.cum = function(lam, N){
### lam: weight coef
  ik = 1:N
  a1 = 2^(ik-1)*gamma(ik)
  cl = rep(0, N)
  for(i in 1:N) cl[i] = a1[i]*sum(lam^i)
  cl
}

## match higher moments
wu.lambda = function(lam, N=12){        # Nth cumulant (Wu Pankow suggests using N=12)
  cl = chi1sqm.cum(lam, N)
  muQ = cl[1]; sigmaQ = sqrt(cl[2])
  a1 = mnc2mc(cum2mnc(cl))
  a1 = a1/sqrt(a1[2])^(1:N)
  f1 = function(xpar){
    k = exp(xpar[1])
    v = xpar[2]
    a2 = mnc2mc(cum2mnc(chisq.cum(k,v,N)))
    a2 = a2/sqrt(a2[2])^(1:N)
    (a1[N-1]-a2[N-1])^2 + (a1[N]-a2[N])^2
  }
  tmp = bobyqa(c(0,1), f1, lower=c(-Inf,0),upper=c(Inf,Inf))
  xpar = tmp$par
  l = exp(xpar[1])
  d = xpar[2]
 if(f1(c(xpar[1],0))<=tmp$fval){
    d=0
    f.1 = function(xpar) f1(c(xpar,0))
    l = exp(bobyqa(xpar[1], f.1)$par)
  }
  muX = l+d; sigmaX = sqrt(chisq.cum(l,d,N=2)[2])
  list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
}
wu.pval = function(Q.all, lambda, N=12){
  param = wu.lambda(lambda,N)
  Q.Norm = (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 = Q.Norm*param$sigmaX + param$muX
  pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
}
##-----------------------End: Functions borrowed from Dr. Baolin Wu's website -----------------------




############################################
usat<-function(Y, X, COV=NULL, na.check=TRUE){		
				# Y: continuous traits
				# X: genotype data; each row for a subject, and each column for an SNP;
				#    the value of each element is the # of the copies for an allele;
				# COV: matrix of covariates other than SNP

  #### check for NA's and remove subjects (rows with at least one NA)
   # lm uses na.omit as default na.action (excludes all NA observations altogether, the residuals do not have the NA observations while they are present when na.action=na.exclude)
   # need to do the same for SSU score calculation else the score vector becomes NA
     k<-ncol(Y)                                 # no. of traits
     q<-ncol(X)                                 # no. of SNPs/markers
     if(q!=1) stop("This USAT code is meant for testing association for only 1 SNP. USAT for multiple SNPs is on its way.")
     if(na.check) 
     {
        input.mat<-cbind(Y,X,COV)
        if(sum(is.na(input.mat))>0) {
	   message("Excluding missing observations.")
           input.mat<-na.omit(input.mat)
           Y<-as.matrix(input.mat[,1:k])
           X<-as.matrix(input.mat[,(k+1):(k+q)])
           if(is.null(COV)==FALSE) COV<-as.matrix(input.mat[,(k+q+1):ncol(input.mat)])
	   message(paste("Current sample size:",nrow(Y)))
	}
     }
     n<-nrow(Y)                                 # no. of individuals

   ############SSU#####################
   Xg <- X
   Xbar<-apply(Xg, 2, mean, na.rm=TRUE)
   Xgb<-Xg
   for(i in 1:nrow(Xg))
      Xgb[i,]<-Xg[i,]-Xbar      # centered X
   Yc<-scale(Y, scale=F)        # centered Y

   ##### Estimate covariate parameter under H0
   if(is.null(COV)==FALSE) {
     COV<-scale(COV, scale=F)
     if(ncol(COV)>1) {
       Z<-t(COV)%*%COV
       phi.hat<-chol2inv(chol(Z))%*%t(COV)%*%matrix(rowMeans(Yc),nrow=n)	# phi is a vector of parameters
       one<-matrix(rep(1,k),nrow=1)
       Yc<-Yc - one %x% (COV%*%phi.hat)
     }else {
       phi.hat<-((t(COV)%*%matrix(rowMeans(Yc),nrow=n))/(sum(COV^2)))[1,]		# phi is a scalar
       one<-matrix(rep(1,k),nrow=1)
       Yc<-Yc - one%x%(COV*phi.hat)
     }
   }
   #########SumSq Scores where U is score vector
        yc<-as.vector(t(Yc))            # one vector of all the centered Y-obs/residuals
        sig2<-var(yc)                   # estimate of the common variance of the Y's
        U<-( t(Yc) %*% Xgb )/sig2
        ##cov of the score stats:
        CovU<-(t(Xgb)%*%Xgb)[1,1]*(t(Yc)%*%Yc)/(nrow(Yc)*sig2^2)
   #SumSqU:
    if (is.null(dim(CovU)))
    { # only one-dim:
        Tssu<- sum(U^2 /CovU)
        if (is.na(Tssu) || is.infinite(Tssu) || is.nan(Tssu)) Tssu<-0
        pssu<-as.double(pchisq(Tssu, 1, lower.tail=FALSE))
    }else
    { #it's possible U=0 and Cov(U)=0:
        if (all(abs(U)<1e-20)) pssu<-1
        else{
            Tssu<- t(U) %*% U
            ##distr of Tssu is sum of c_r Chisq_1 (c_r is rth eigen value of CovU):
            eigC<-eigen(CovU, only.values=TRUE)$values
	    m<-sum(eigC>0)		# number of positive eigen values
            ##approximate the distri by alpha Chisq_d + beta:
            alpha1<-as.double( sum(eigC^3)/sum(eigC^2) )
            beta1<-as.double( sum(eigC) - (sum(eigC^2))^2/(sum(eigC^3)) )
            d1<-as.double( (sum(eigC^2))^3/(sum(eigC^3))^2 )
	    pssu<-as.double(pchisq((Tssu-beta1)/alpha1, d1, lower.tail=FALSE))    #p-value for SSU test
        }
    }
 
   ############ MANOVA #####################
   res<-total<-NULL
   if(is.null(COV)==FALSE) {	# if COV is not NULL, the residual for SSU was stored in Yc; so re-store centered Y in Yc
      Yc<-scale(Y, scale=F)        # centered Y
   }
   #### regress each trait on all SNPs
   for (t in 1:k){
     if(is.null(COV)){
   	ff0<-lm(Yc[,t]~0)
	ff1<-lm(Yc[,t]~0+Xgb)
     }else {
   	ff0<-lm(Yc[,t]~0+COV)
	ff1<-lm(Yc[,t]~0+COV+Xgb)
     }
   	total0<-resid(ff0)
   	res0<-resid(ff1)
   	total<-cbind(total,total0)
   	res<-cbind(res,res0)
   }
   E=cov(res)
   TOT=cov(total)
   if (det(E)==0|det(TOT)==0) print("Alert! Singular covariance matrix")

   if(det(TOT)==0){ 
      print("Inverse of Singular covariance matrix not possible")
      # Neither USAT nor MANOVA statistics can be calculated 
      Tc<-dfman<-pman<-pval.opt<-rho.opt<-pval.T<-NA
   }else {
      ## Wilk's Lambda
      if(k>1) lambda<-det(E%*%mpinv(TOT)) else lambda<-det(E)/det(TOT)

      #################################### USAT: combined statistic  ###################################
      Ts = Tssu
      Tc = -n*log(lambda)
	dfman<-k
	pman<-as.double(pchisq(Tc, dfman, lower.tail=FALSE))
      pval.rho<-q.rho<-To.rho<-NULL
      rho.set<-seq(0,1,0.1)	# set of rho values (weights) to consider in the linear combination
      for(rho in rho.set)
      {
	To<-(1-rho)*Ts + rho*Tc
	To.rho<-c(To.rho,To)		# value of unified statistics for different values of rho
	eig<- (1-rho)*eigC+rho
	me <- sum(eig>0)
        dav.qf<-davies(q=To, lambda=eig[1:me], h=rep(1,me),acc=1e-12,lim=1e+4)
        if(dav.qf$ifault!=0 | dav.qf$Qq<0 | dav.qf$Qq==0) pval.qf<-pchisqsum(x=To, df=rep(1,me), a=eig[1:me], lower.tail=FALSE, method="saddlepoint") else pval.qf<-dav.qf$Qq
        pval.rho<-c(pval.rho, pval.qf)
      }
      pval.opt<-min(pval.rho)
      rho.opt<-0+(which(pval.rho==pval.opt)-1)*0.1

      # calculate the p-value using a one-dimensional numerical integration
      # pval.opt is our test statistic T, whose p-value needs to be ccomputed
	 # first find the (1-T)th quantiles of the linear combination of chi-sq distributions for every rho
	 T<-pval.opt
	 n.r<-length(rho.set)
	 qminrho<-rep(0,n.r)
	 c1 <- rep(0, 4)
    	 for (i in 1:n.r) {
        	rho <- rho.set[i]
		eig<- (1-rho)*eigC+rho
		if(T>1e-4) {
        	   c1[1] <- sum(eig)
        	   c1[2] <- sum(eig^2)
	           c1[3] <- sum(eig^3)
	           c1[4] <- sum(eig^4)
	           muQ <- c1[1]
	           varQ <- 2*c1[2]
		   s1 <- c1[3]/c1[2]^1.5
		   s2 <- c1[4]/c1[2]^2
		   if(s1^2>s2){
		      a=1/(s1-sqrt(s1^2-s2))
		      ncp.d=s1*a^3-a^2
		      df.l=a^2-2*ncp.d
		   }else{
		     a=1/s1 ; ncp.d=0 ; df.l=c1[2]^3/c1[3]^2
		   }
			mu.chi<-df.l+ncp.d
			sig.chi<-sqrt(2)*a
			q.org <- qchisq(T, df=df.l, ncp=ncp.d, lower.tail=FALSE)
			q.q <- ( (q.org-mu.chi)/sig.chi )*sqrt(varQ) + muQ
		}else {
		   wu.out<-wu.lambda(lam=eig, N=6)
                        mu.chi<-wu.out$muX
                        sig.chi<-wu.out$sigmaX
                        q.org <- qchisq(T, df=wu.out$l, ncp=wu.out$d, lower.tail=FALSE)
                        q.q <- ( (q.org-mu.chi)/sig.chi )*(wu.out$sigmaQ) + wu.out$muQ
		}
	        qminrho[i] <- q.q
	 }

      # the required pvalue is 1-P(chi2_k<quant|Ts)*f(Ts), where quant=min_{rho}((qminrho-(1-rho)*Ts)/rho) and
      # f(Ts) is the density of Ts, where Ts~alpha1*chi2_(d1)+beta1
      # i.e., pvalue is P(chi2_k>=quant|Ts)*f(Ts)
	f.Liu.chi<-function (x) 
	{
	    temp <- (qminrho - rho.set*x)/(1 - rho.set)
	    temp.min <- min(temp)
	    re <- pchisq((temp.min-beta1)/alpha1, df=d1, ncp=0) * dchisq(x, df=dfman)	# using the scaled and shifted chi-sq distn of SSU by Zhang (2005)
	    return(re)
	}
        ## for very small USAT statistic, pracma library is used for computing p-values for better accuracy
        ## however quadinf from pracma has tolerance parameter limit at 1e-15 and p-values <1e-16 are returned as 0; in this situation integrate is used on a slightly modified function
        if(pval.opt>1e-5){
		 pval.T<-as.double(1-integrate(Vectorize(f.Liu.chi),0,Inf)$value)
	}else if(pval.opt<=5e-16){
                f.Liu.chi.upp<-function (x)
                {
                    temp <- (qminrho - rho.set*x)/(1 - rho.set)
                    temp.min <- min(temp)
                    re <- pchisq((temp.min-beta1)/alpha1, df=d1, ncp=0, lower.tail=FALSE) * dchisq(x, df=dfman)   # using the scaled and shifted chi-sq distn of SSU by Zhang (2005)
                    return(re)
                }
                pval.T<-as.double(integrate(Vectorize(f.Liu.chi.upp),0,Inf)$value)
        }else {   # default tolerance paramater for quadinf is 1e-12
                  # however tolerance can not be decreased below 1e-15
        if(pval.opt<1e-12) pval.T<-as.double(1-quadinf(f.Liu.chi,0,Inf,tol=1e-15)$Q) else pval.T<-as.double(1-quadinf(f.Liu.chi,0,Inf)$Q)
        }
 }# end of else loop for det(TOT)=0
  
 return(list(T.usat=pval.opt, omg.opt=rho.opt[1], p.usat=pval.T, n.obs=n))

}

