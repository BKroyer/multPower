#' @title Multivariate Sum Test (O'Brien's Test) for Comparing Two Groups.
#'
#' @description Tests the multivariate null hypothesis of a zero mean difference between two multivariate (normal) distributions.
#'
#' @param est vector of estimated mean differences.
#' @param V estimate for the covariance matrix.
#' @param sides defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2. See details.
#' @param N number of independent observations. 
#' @param dftype type of degrees of freedom. Can be either "OBrien"  (df=N-2k) or "Logan" (df=0.5(N-2)(1-1/k^2)) with k=length(est). Default is dftype="OBrien".
#' 
#'
#' @details Calculates O'Briens OLS test for the comparison of two independent mean vectors. 
#'    For the two-sided test, the test statistic of this test is the absolute value of the standardized sum
#'    of elementary Wald test statistics. For the one-sided case as currently implemented,
#'    the test statistics is standardized sum
#'    of elementary Wald test statistics. Hence the one-sided test is sensitive for deviations from the null hypothesis
#'    corresponding to a positive sum of standardized mean differences.
#'    It is assumed that \code{est} is the estimated vector of mean differences and \code{V} is the pooled sample
#'    covariance matrix estimate and \code{N} is the total number of subjects from which these estimates where obtained.
#'    The function is mainly included to allow for validation of the power and sample size 
#'    calculations via simulation. See examples.
#'
#' @return A data frame with the test statistic, numerator and denominator degrees of freedom, the p-value, and an indicator if the p-value is one-sided or two-sided.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references O'Brien, Peter C. "Procedures for comparing samples with multiple endpoints." Biometrics (1984): 1079-1087.
#' Logan, B. R., & Tamhane, A. C. (2004). On O'Brien's OLS and GLS tests for multiple endpoints. Lecture Notes-Monograph Series, 76-88.
#'
#' @seealso \code{\link{n_sum_t_test}}, \code{\link{power_sum_t_test}}
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
#' nn<-n_sum_t_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_sum_t_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-nn$n0
#' n1<-nn$n1
#' R<-1000
#' p<-rep(NA,R)
#' for(i in 1:R) {
#'  dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#'  p[i]<-sum_t_test(dat$est,dat$V,dat$N)$p
#'}
#'mean(p<=alpha)
#' 
#'
#' @export
sum_t_test<-function(est,V,N,dftype=c("OBrien","Logan")[1],sides=2) {
  SE<-sqrt(diag(V))
  z<-est/SE
  stat<-sum(z)/sqrt(sum(cov2cor(V)))
  k<-length(est)
  if(dftype=="OBrien") df=N-2*k
  if(dftype=="Logan") df=0.5*(N-2)*(1-1/k^2)
  if(sides==2) p<-2*(1-pt(abs(stat),df))
  #one-sided
  if(sides==1) p<-1-pt(stat,df)
  data.frame(stat,df=df,p,sides)
}


#needed to find starting value in n_sum_t_test
n_sum_z_test<-function(power=0.8,r=0.5,delta,K,alpha=0.05,sides=2) { 
  #r=n1/N
  VC_1<-K*(1/(1-r)+1/r)
  SE1<-sqrt(diag(VC_1))
  korrmat<-cov2cor(VC_1)
  ncp1<-abs(sum(delta/SE1)/sqrt(sum(korrmat)))
  if(sides==1) krit<-qnorm(1-alpha)
  if(sides==2) krit<-qnorm(1-alpha/2)
  N<-(krit+qnorm(power))^2/ncp1^2
  n1<-r*N
  n0<-(1-r)*N
  data.frame(N,n0,n1,power=power)
}



#' @title Sample Size Calculation for O'Brien's OLS Sum Test
#'
#' @description Calculates the sample size for O'Brien's OLS sum test for the comparison of two independent mean vectors.
#'
#' @param power the aimed for power of the test. The default is power=0.8.
#' @param r fraction of subjects in group 1. The default value r=0.5 means equally sized groups.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#' @param alpha significance level, default is 0.05.
#' @param sides  defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2. See details.
#' @param dftype type of degrees of freedom. Can be either "OBrien"  (df=N-2k) or "Logan" (df=0.5(N-2)(1-1/k^2)) with k=length(est). Default is dftype="OBrien".
#'
#' @details  For the two-sided test, the test statistic of this test is the absolute value of the standardized sum
#'    of elementary Wald test statistics. For the one-sided case as currently implemented,
#'    the test statistics is the standardized sum
#'    of elementary Wald test statistics. Hence the one-sided test is sensitive for deviations from the null hypothesis
#'    corresponding to a positive sum of standardized mean differences. When performing power or sample size calculations for the one-sided test,
#'    the sign of the entries in \code{est} must be set accordingly.
#'
#' @return A data frame with the total sample size \code{N}, the sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated actual power
#'   (which should match closely with the aimed for power), and the indication if a one-sided or two-sided null-hypothesis is tested.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references O'Brien, Peter C. "Procedures for comparing samples with multiple endpoints." Biometrics (1984): 1079-1087.
#' Logan, B. R., & Tamhane, A. C. (2004). On O'Brien's OLS and GLS tests for multiple endpoints. Lecture Notes-Monograph Series, 76-88.
#'
#' @seealso \code{\link{sum_t_test}}, \code{\link{power_sum_t_test}}
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
#' nn<-n_sum_t_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_sum_t_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-nn$n0
#' n1<-nn$n1
#' R<-1000
#' p<-rep(NA,R)
#' for(i in 1:R) {
#'  dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#'  p[i]<-sum_t_test(dat$est,dat$V,dat$N)$p
#'}
#'mean(p<=alpha)
#'
#' @export
n_sum_t_test<-function(power=0.8,r=0.5,delta,K,alpha=0.05,sides=2,dftype=c("OBrien","Logan")[1]) { 
  #r ist n1/N
  #also n1=r*N
  n_start<-n_sum_z_test(power,r,delta,K,alpha,sides)
  N<-ceiling(n_start$N)
  n1<-ceiling(r*N)
  n0<-ceiling((1-r)*N)
  pow<-power_sum_t_test(n0,n1,delta,K,alpha,sides)$power
  while(pow<power) {
    N<-N+1
    n1<-ceiling(r*N)
    n0<-ceiling((1-r)*N)
    pow<-power_sum_t_test(n0,n1,delta,K,alpha,sides,dftype)$power
  }
  data.frame(N=n0+n1,n0,n1,power=pow,sides=sides)
}



#' @title Power Calculation for O'Brien's OLS Sum Test 
#'
#' @description Calculates the power for O'Brien's OLS sum test for the comparison of two independent mean vectors.
#'
#' @param n0 sample size in group 0.
#' @param n1 sample size in group 1. Defaults to \code{n0} if not specified otherwise.
#' @param delta vector of assumed true mean differences.
#' @param K assumed true covariance matrix (common to both groups).
#' @param alpha significance level, default is 0.05.
#' @param sides  defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2. See details.
#' @param dftype type of degrees of freedom. Can be either "OBrien"  (df=N-2k) or "Logan" (df=0.5(N-2)(1-1/k^2)) with k=length(est). Default is dftype="OBrien".
#'
#' @details  For the two-sided test, the test statistic of this test is the absolute value of the standardized sum
#'    of elementary Wald test statistics. For the one-sided case as currently implemented,
#'    the test statistics is the standardized sum
#'    of elementary Wald test statistics. Hence the one-sided test is sensitive for deviations from the null hypothesis
#'    corresponding to a positive sum of standardized mean differences. When performing power or sample size calculations for the one-sided test,
#'    the sign of the entries in \code{est} must be set accordingly.
#'
#' @return A data frame with the input sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated power, and the indication if a one-sided or two-sided null-hypothesis is tested.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references O'Brien, Peter C. "Procedures for comparing samples with multiple endpoints." Biometrics (1984): 1079-1087.
#' Logan, B. R., & Tamhane, A. C. (2004). On O'Brien's OLS and GLS tests for multiple endpoints. Lecture Notes-Monograph Series, 76-88.
#'
#' @seealso \code{\link{sum_t_test}}, \code{\link{n_sum_t_test}}
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
#' nn<-n_sum_t_test(power=0.8,r=r,delta=delta,K=K,alpha=alpha)
#' #Check power
#' power_sum_t_test(n0=nn$n0,n1=nn$n1,delta=delta,K=K,alpha=alpha)
#' #Check power by simulation
#' n0<-nn$n0
#' n1<-nn$n1
#' R<-1000
#' p<-rep(NA,R)
#' for(i in 1:R) {
#'  dat<-sim_1(n0=n0,n1=n1,delta=delta,K=K)
#'  p[i]<-sum_t_test(dat$est,dat$V,dat$N)$p
#'}
#'mean(p<=alpha)
#'
#' @export
power_sum_t_test<-function(n0,n1,delta,K,alpha=0.05,sides=2,dftype=c("OBrien","Logan")[1]) {
  VC<-K/n0+K/n1
  SE<-sqrt(diag(VC))
  z<-delta/SE
  korrmat<-cov2cor(VC)
  ncp<-abs(sum(z)/sqrt(sum(korrmat)))
  k<-length(delta)
  N<-n0+n1
  if(dftype=="OBrien") df=N-2*k
  if(dftype=="Logan") df=0.5*(N-2)*(1-1/k^2)
  
  if(sides==2) {
    krit<-qt(1-alpha/2,df=df)
    power<-1-pt(krit,ncp=ncp,df=df)
  }
  if(sides==1) {
    krit<-qt(1-alpha,df=df)
    power<-1-pt(krit,ncp=ncp,df=df)
    
  }
  
  data.frame(n0,n1,power,sides)
}

