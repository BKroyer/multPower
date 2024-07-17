#' @importFrom stats cov2cor pf pt qf qnorm qt rWishart uniroot
#' @importFrom mvtnorm pmvnorm pmvt qmvnorm qmvt rmvnorm


#' @title Simulate Mean Vector and Covariance Estimates
#'
#' @description Simulates the results of sampling from two independent normal distributions and estimating the mean difference and the pooled covariance matrix.
#'   This is done by drawing mean vectors from multivariate normal distributions and covariance matrices from Wishart distributions and subsequently calculating the difference and
#'   the pooled covariance matrix estimate. The function is intended to support simulations to validate sample size and power calculations performed with
#'   \code{\link{n_F_test}} and \code{\link{power_F_test}}.
#'
#' @param n0 sample size in group 0.
#' @param n1 sample size in group 1.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#'
#' @return A list with the simulated vector of mean differences (\code{est}),
#'    covariance matrix (\code{V}), the input sample sizes \code{n0} and \code{n1} and the total sample size \code{N=n0+n1}.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hotelling, Harold. "The Generalization of Student's Ratio" Ann. Math. Statist. 2(3) (1931): 360-378.
#' @seealso \code{\link{n_F_test}}, \code{\link{power_F_test}}
#'
#' @export
sim_1<-function(n0=100,n1=100,delta,K) {
	N<-n0+n1
	k<-length(delta)
	DF0<-n0-1
	DF1<-n1-1
	W0<-rWishart(1,df=DF0,Sigma=K)[,,1]
	W1<-rWishart(1,df=DF1,Sigma=K)[,,1]
	V<-(W0+W1)/(N-2) *(1/n0+1/n1)
	d<-rmvnorm(1,mean=delta,K/n0+K/n1)
	est<-as.numeric(d)
	list(est=est,V=V,n0=n0,n1=n1,N=N)
}


#' @title Create Covariance Matrix with Autoregressive Correlation Structure
#'
#' @description Creates a covariance matrix of structure AR(1) (autoregressive of order 1) based on given time points, correlation between adjacent time points and standard deviations of each time point.
#'
#' @param time vector of time-points.
#' @param sigma vector of the standard deviations assumed at the different time point.
#' @param rho correlation between adjacent time-points of distance 1, see details.
#'
#' @details The correlation between two time points \eqn{a} and \eqn{b} and AR(1) correlation parameter \eqn{\rho} is \eqn{\rho^{|a-b|}}.
#'
#' @return A \eqn{k} by \eqn{k} covariance matrix of AR(1) structure, with \eqn{k} being the number of the time-points.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#'
#' @examples
#' ARmat(1:3,rep(1,3),0.5)
#' K<-ARmat(1:3,c(0.5,0.7,1),0.5)
#' cov2cor(K)
#'
#' @export
ARmat<-function(time,sigma,rho) {
  k<-length(time)
  i<-rep(time,each=k)
  j<-rep(time,k)
  mat<-data.frame(i,j,corr=NA)
  mat$delta<-abs(time[i]-time[j])
  mat$cor<-rho^mat$delta
  Korr<-matrix(mat$cor,k,k)
  K<-diag(sigma)%*%Korr%*%diag(sigma)
  K
}


#' @title Create Covariance Matrix with Exchangeable Correlation Structure
#'
#' @description Creates a covariance matrix with a common correlation between all time-points. The variances may be heterogeneous.
#'
#' @param time vector of time points.
#' @param sigma vector of the standard deviations for each time point.
#' @param rho correlation coefficient between adjacent time points.
#'
#' @return A \eqn{k} by \eqn{k} covariance matrix with exchangeable correlation structure, with \eqn{k} being the number of the time-points.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#'
#' @examples
#' EXmat(1:2,rep(1,2),0.5)
#' K<-EXmat(1:4,c(0.75,0.75,1.5,1.5),0.5)
#' cov2cor(K)
#'
#' @export
EXmat<-function(time,sigma,rho) {
  k<-length(time)
  Korr<-matrix(rho,k,k)
  diag(Korr)<-1
  K<-diag(sigma)%*%Korr%*%diag(sigma)
  K
}



#' @title Multivariate F-test (Hotelling's T2 Test) for Comparing Two Groups.
#'
#' @description Tests the multivariate null hypothesis of a zero mean difference between two multivariate (normal) distributions.
#'
#' @param est vector of estimated mean differences.
#' @param V estimate for the covariance matrix.
#' @param N number of independent observations.
#'
#'
#' @details The function calculates Hotelling's T2 test for the comparison of two independent mean vectors.
#'    It is assumed that \code{est} is the estimated vector of mean differences and \code{V} is the pooled sample
#'    covariance matrix estimate and \code{N} is the total number of subjects from which these estimates where obtained.
#'    The function is mainly included to allow for validation of the power and sample size
#'    calculations via simulation. See examples.
#'
#' @return A data frame with the test statistic, numerator and denominator degrees of freedom and the p-value.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hotelling, Harold. "The Generalization of Student's Ratio" Ann. Math. Statist. 2(3) (1931): 360-378.
#' @seealso \code{\link{n_F_test}}, \code{\link{power_F_test}}
#'
#' @examples
#' library(mvtnorm)
#' #Assume true difference and some correlation
#' delta<-c(0.2,0.6)
#' K<-EXmat(1:2,rep(1,2),0.5)
#' #Assume 1:2 allocation between the two groups
#' r<-2/3
#' #Calculate sample size for 80% power with significance level 0.05
#' alpha<-0.05
#' nn<-n_F_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_F_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-ceiling(nn$n0)
#' n1<-ceiling(nn$n1)
#' R<-1000
#' p<-rep(NA,R)
#' library(mvtnorm)
#' for(i in 1:R) {
#' 	dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#' 	p[i]<-F_test(dat$est,dat$V,dat$N)$p
#' }
#' mean(p<=alpha)
#'
#' @export
F_test<-function(est,V,N) {
	k<-length(est)
	#Hotelling:
	stat<- (N-k-1)/(N-2)/k*est%*%solve(V,est)
	df1<-k
	df2<-N-1-k
	p<-1-pf(stat,df1,df2)
	data.frame(stat,df1,df2,p)
}


#' @title Sample Size Calculation for the Multivariate F-Test (Hotelling's T2 Test)
#'
#' @description Calculates the sample size for Hotelling's T2 test for the comparison of two independent mean vectors.
#'
#' @param power the aimed for power of the test. The default is power=0.8.
#' @param r fraction of subjects in group 1. The default value r=0.5 means equally sized groups.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#' @param alpha significance level, default is 0.05.
#' @param interval vector of length two with the interval to search for the sample size, passed to \code{\link[stats]{uniroot}}.
#     The default is c(5,10000).
#'
#' @return A data frame with the total sample size \code{N}, the sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated actual power
#'   (which should match closely with the aimed for power).
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hotelling, Harold. "The Generalization of Student's Ratio" Ann. Math. Statist. 2(3) (1931): 360-378.
#' @seealso \code{\link{F_test}}, \code{\link{power_F_test}}
#'
#' @examples
#' library(mvtnorm)
#' #Assume true difference and some correlation
#' delta<-c(0.2,0.6)
#' K<-EXmat(1:2,rep(1,2),0.5)
#' #Assume 1:2 allocation between the two groups
#' r<-2/3
#' #Calculate sample size for 80% power with significance level 0.05
#' alpha<-0.05
#' nn<-n_F_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_F_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-ceiling(nn$n0)
#' n1<-ceiling(nn$n1)
#' R<-1000
#' p<-rep(NA,R)
#' library(mvtnorm)
#' for(i in 1:R) {
#' 	dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#' 	p[i]<-F_test(dat$est,dat$V,dat$N)$p
#' }
#' mean(p<=alpha)
#'
#' @export
n_F_test<-function(power=0.8,r=0.5,delta,K,alpha=0.05,interval=c(5,10000)) {
	#r=n1/N
	powfun<-function(N) {
		n1<-r*N
		n0<-(1-r)*N
		pow<-power_F_test(n0,n1,delta,K,alpha)$power
		pow-power
	}
	rt<-uniroot(powfun,interval)
	N<-rt$root
	n1<-r*N
	n0<-(1-r)*N
	data.frame(N,n0,n1,power=rt$f.root+power)
}



#' @title Power Calculation for the Multivariate F-Test (Hotelling's T2 Test)
#'
#' @description Calculates the power for Hotelling's T2 test for the comparison of two independent mean vectors.
#'
#' @param n0 sample size in group 0.
#' @param n1 sample size in group 1. Defaults to \code{n0} if not specified otherwise.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#' @param alpha significance level, default is 0.05.
#'
#' @return A data frame with the input sample sizes in group 0 and group 1, \code{n0} and \code{n1}, and the calculated power.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hotelling, Harold. "The Generalization of Student's Ratio" Ann. Math. Statist. 2(3) (1931): 360-378.
#' @seealso \code{\link{F_test}}, \code{\link{n_F_test}}
#'
#' @examples
#' library(mvtnorm)
#' #Assume true difference and some correlation
#' delta<-c(0.2,0.6)
#' K<-EXmat(1:2,rep(1,2),0.5)
#' #Assume 1:2 allocation between the two groups
#' r<-2/3
#' #Calculate sample size for 80% power with significance level 0.05
#' alpha<-0.05
#' nn<-n_F_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_F_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-ceiling(nn$n0)
#' n1<-ceiling(nn$n1)
#' R<-1000
#' p<-rep(NA,R)
#' library(mvtnorm)
#' for(i in 1:R) {
#' 	dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#' 	p[i]<-F_test(dat$est,dat$V,dat$N)$p
#' }
#' mean(p<=alpha)
#'
#' @export
power_F_test<-function(n0,n1=n0,delta,K,alpha=0.05) {
	VC<-K/n0+K/n1
	ncp<-delta%*%solve(VC,delta)
	k<-length(delta)
	N<-n0+n1
	df1<-k
	df2<-N-1-k
	Q<-qf(1-alpha,df1=df1,df2=df2)
	power<-1-pf(Q,df1=df1,df2=df2,ncp=ncp)
	data.frame(n0,n1,power)
}

