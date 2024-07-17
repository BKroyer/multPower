#' @title Multivariate Maximum t-test for comparing two groups.
#'
#' @description Tests the multivariate null hypothesis of a zero mean difference between two multivariate (normal) distributions.
#'
#' @param est vector of estimated mean differences.
#' @param V estimate for the covariance matrix.
#' @param N number of independent observations. 
#' @param sides defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2.
#' 
#'
#' @details The function calculates Maximum t-test for the comparison of two independent mean vectors.
#'    It is assumed that \code{est} is the estimated vector of mean differences and \code{V} is the pooled sample
#'    covariance matrix estimate and \code{N} is the total number of subjects from which these estimates where obtained.
#'    The function is mainly included to allow for validation of the power and sample size 
#'    calculations via simulation. See examples.
#'
#' @return A data frame with the test statistic, numerator and denominator degrees of freedom, the p-value, and an indicator if the p-value is one-sided or two-sided.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hothorn, Torsten, Frank Bretz, and Peter Westfall. "Simultaneous inference in general parametric models." Biometrical Journal: Journal of Mathematical Methods in Biosciences 50.3 (2008): 346-363.
#' @seealso \code{\link{n_max_t_test}}, \code{\link{power_max_t_test}}
#'
#' @examples
#' library(mvtnorm)
#' #Assume true difference and some correlation
#' delta<-c(0.2,0.6)
#' K<-ARmat(1:2,rep(1,2),0.5)
#' #Assume 1:3 allocation between the two groups
#' r<-1/4
#' #Calculate sample size for 80% power with significance level 0.05
#' alpha<-0.05
#' nn<-n_max_t_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_max_t_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-nn$n0
#' n1<-nn$n1
#' R<-1000
#' p<-rep(NA,R)
#' for(i in 1:R) {
#'  dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#'  p[i]<-max_t_test(dat$est,dat$V,dat$N)$p
#'}
#'mean(p<=alpha)
#' 
#' 
#'
#' @export
max_t_test<-function(est,V,N,sides=2) {
  k<-length(est)
  SE<-sqrt(diag(V))
  z<-est/SE
  df<-N-2
  if(sides==2) {
    stat<-rep(max(abs(z)),k)
    p<-as.numeric(1-pmvt(lower=-stat,upper=stat,df=df,corr=cov2cor(V)))
  }
  if(sides==1) {
    stat<-rep(max(z),k)
    p<-as.numeric(1-pmvt(lower=rep(-Inf,k),upper=stat,df=df,corr=cov2cor(V)))
  }
  
  data.frame(stat=max(z),df,p,sides)
}


#needed to find starting value in n_max_t_test
n_max_z_test<-function(power=0.8,r=0.5,delta,K,alpha=0.05,sides=2,interval=c(5,10000)) { 
  #r ist n1/N
  #also n1=r*N
  powfun<-function(N) {	
    n1<-r*N
    n0<-(1-r)*N
    pow<-power_max_z_test(n0,n1,delta,K,alpha,sides)$power
    pow-power
  }
  rt<-uniroot(powfun,interval)
  N<-rt$root
  n1<-r*N
  n0<-(1-r)*N
  data.frame(N,n0,n1,power=rt$f.root+power)
}

power_max_z_test<-function(n0,n1,delta,K,alpha=0.05,sides=2) {
  VC<-K/n0+K/n1
  SE<-sqrt(diag(VC))
  ncp<-delta/SE
  korrmat<-cov2cor(VC)
  k<-length(delta)
  N<-n0+n1
  if(sides==2) {
    Q<-qmvnorm(1-alpha,mean=rep(0,k),corr=korrmat,tail="both.tails")$quantile
    krit<-rep(Q,k)
    power<-1-pmvnorm(lower=-krit,upper=krit,mean=ncp,corr=korrmat)
  }
  if(sides==1) {
    Q<-qmvnorm(1-alpha,mean=rep(0,k),corr=korrmat,tail="lower.tail")$quantile
    krit<-rep(Q,k)
    power<-1-pmvnorm(lower=-Inf,upper=krit,mean=ncp,corr=korrmat)
  }
  
  data.frame(N,n0,n1,power)
}



#' @title Sample size calculation for the Maximum t-test
#'
#' @description Calculates the sample size for Maximum t-test for the comparison of two independent mean vectors.
#'
#' @param power the aimed for power of the test. The default is power=0.8.
#' @param r fraction of subjects in group 1. The default value r=0.5 means equally sized groups.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#' @param alpha significance level, default is 0.05.
#' @param sides  defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2.
#' @param interval vector of length two with the interval to search for the sample size, passed to \code{\link[stats]{uniroot}}.
#     The default is c(5,10000).
#'
#' @return A data frame with the total sample size \code{N}, the sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated actual power
#'   (which should match closely with the aimed for power), and the indication if a one-sided or two-sided null-hypothesis is tested.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hothorn, Torsten, Frank Bretz, and Peter Westfall. "Simultaneous inference in general parametric models." Biometrical Journal: Journal of Mathematical Methods in Biosciences 50.3 (2008): 346-363.
#' @seealso \code{\link{max_t_test}}, \code{\link{power_max_t_test}}
#'
#' @examples
#' library(mvtnorm)
#' #Assume true difference and some correlation
#' delta<-c(0.2,0.6)
#' K<-ARmat(1:2,rep(1,2),0.5)
#' #Assume 1:3 allocation between the two groups
#' r<-1/4
#' #Calculate sample size for 80% power with significance level 0.05
#' alpha<-0.05
#' nn<-n_max_t_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_max_t_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-nn$n0
#' n1<-nn$n1
#' R<-1000
#' p<-rep(NA,R)
#' for(i in 1:R) {
#'  dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#'  p[i]<-max_t_test(dat$est,dat$V,dat$N)$p
#'}
#'mean(p<=alpha)
#'
#' @export
n_max_t_test<-function(power=0.8,r=0.5,delta,K,alpha=0.05,sides=2,interval=c(5,10000)) { 
  #r ist n1/N
  #also n1=r*N
  #Startwerte
  n_start<-n_max_z_test(power,r,delta,K,alpha,sides,interval)
  N<-ceiling(n_start$N)
  n1<-ceiling(r*N)
  n0<-ceiling((1-r)*N)
  pow<-power_max_t_test(n0,n1,delta,K,alpha,sides)$power
  while(pow<power) {
    N<-N+1
    n1<-ceiling(r*N)
    n0<-ceiling((1-r)*N)
    pow<-power_max_t_test(n0,n1,delta,K,alpha,sides)$power
  }
  data.frame(N=n0+n1,n0,n1,power=pow,sides=sides)
}



#' @title Power calculation for the Maximum t-test
#'
#' @description Calculates the power for Maximum t-test for the comparison of two independent mean vectors.
#'
#' @param n0 sample size in group 0.
#' @param n1 sample size in group 1. Defaults to \code{n0} if not specified otherwise.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#' @param alpha significance level, default is 0.05.
#' @param sides  defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2.
#'
#' @return A data frame with the input sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated power, and the indication if a one-sided or two-sided null-hypothesis is tested.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Hothorn, Torsten, Frank Bretz, and Peter Westfall. "Simultaneous inference in general parametric models." Biometrical Journal: Journal of Mathematical Methods in Biosciences 50.3 (2008): 346-363.
#' @seealso \code{\link{max_t_test}}, \code{\link{n_max_t_test}}
#'
#' @examples
#' library(mvtnorm)
#' #Assume true difference and some correlation
#' delta<-c(0.2,0.6)
#' K<-ARmat(1:2,rep(1,2),0.5)
#' #Assume 1:3 allocation between the two groups
#' r<-1/4
#' #Calculate sample size for 80% power with significance level 0.05
#' alpha<-0.05
#' nn<-n_max_t_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_max_t_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-nn$n0
#' n1<-nn$n1
#' R<-1000
#' p<-rep(NA,R)
#' for(i in 1:R) {
#'  dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#'  p[i]<-max_t_test(dat$est,dat$V,dat$N)$p
#'}
#'mean(p<=alpha)
#'
#' @export
power_max_t_test<-function(n0,n1,delta,K,alpha=0.05,sides=2) {
  VC<-K/n0+K/n1
  SE<-sqrt(diag(VC))
  ncp<-delta/SE
  korrmat<-cov2cor(VC)
  k<-length(delta)
  N<-n0+n1
  df<-N-2
  if(sides==2) {
    Q<-qmvt(1-alpha,mean=rep(0,k),corr=korrmat,df=df,tail="both.tails")$quantile
    krit<-rep(Q,k)
    power<-1-pmvt(lower=-krit,upper=krit,delta=ncp,corr=korrmat,df=df)
  }
  if(sides==1) {
    Q<-qmvt(1-alpha,mean=rep(0,k),corr=korrmat,tail="lower.tail",df=df)$quantile
    krit<-rep(Q,k)
    power<-1-pmvt(lower=-Inf,upper=krit,delta=ncp,corr=korrmat,df=df)
  }
  
  data.frame(N,n0,n1,power,sides)
}

