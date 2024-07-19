# needed for finding the starting value
n_z_test<-function(power=0.8,r=0.5,delta,sd_per_group,alpha=0.05,sides=2) {
  if (length(sd_per_group) == 1){
    K0 <- K1 <- sd_per_group^2
  } else if (length(sd_per_group) == 2){
    K0 <- sd_per_group[1]^2
    K1 <- sd_per_group[2]^2
  } else{
    stop("more standard deviations than groups given")
  }
  v_1<-K0/(1-r)+K1/r
  SE1<-sqrt(v_1)
  ncp1<-abs(delta/SE1)
  if(!(sides%in%c(1,2))) {
    stop("sides must be eiher 1 or 2")
  } else {
    krit<-qnorm(1-alpha/sides)
  }
  N<-(krit+qnorm(power))^2/ncp1^2
  n1<-r*N
  n0<-(1-r)*N
  data.frame(N,n0,n1,power=power,sides=sides)
}

#' @title Sample Size Calculation for the Two-Sample t-test
#'
#' @description Calculates the sample size for the Two-Sample t-test for the comparison of two independent mean vectors.
#'
#' @param power the aimed for power of the test. The default is power=0.8.
#' @param r fraction of subjects in group 1 (\code{n1}). The default value r=0.5 means equally sized groups. For example, a r of 1/3 would mean an allocation of 2:1 for n0:n1.
#' @param delta vector of assumed true mean differences.
#' @param sd_per_group assumed standard deviation per group. Can be either a number if the same standard deviation is assumed for every group, or a vector with individual standard deviations.
#' @param alpha significance level, default is 0.05.
#' @param sides defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2.
#'
#' @return A data frame with the total sample size \code{N}, the sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated actual power
#'   (which should match closely with the aimed for power), and the indication if a one-sided or two-sided null-hypothesis was assumed.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @author Bettina Kroyer \email{bettina.kroyer@@meduniwien.ac.at}
#' @seealso \code{\link{power_t_test}}
#'
#' @examples
#' n_t_test(0.8,0.5,3,c(5,4))
#' n_t_test(power = 0.8, r = 1/3, delta = 3, sd_per_group = c(5,4))
#' # linking sample size and power calculations
#' delta = 3
#' sd_per_group = 5
#' nn <- n_t_test(0.8,1/3,delta,sd_per_group)
#' power_t_test(nn$n0,nn$n1,delta,sd_per_group) # matches the power initiated in n_t_test
#' @export
n_t_test<-function(power=0.8,r=0.5,delta,sd_per_group,alpha=0.05,sides=2) {
  n_start<-n_z_test(power,r,delta,sd_per_group,alpha,sides)
  N<-ceiling(n_start$N)
  n1<-ceiling(r*N)
  n0<-ceiling((1-r)*N)
  pow<-power_t_test(n0,n1,delta,sd_per_group,alpha,sides)$power
  while(pow<power) {
    N<-N+1
    n1<-ceiling(r*N)
    n0<-ceiling((1-r)*N)
    pow<-power_t_test(n0,n1,delta,sd_per_group,alpha,sides)$power
  }
  data.frame(N=n0+n1,n0,n1,power=pow,sides=sides)
}



#' @title Power Calculation for the Two-sample t-test
#'
#' @description Calculates the power for the Two-sample t-test for the comparison of two independent mean vectors.
#'
#' @param n0 sample size in group 0.
#' @param n1 sample size in group 1. Defaults to \code{n0} if not specified otherwise.
#' @param delta assumed true mean difference between groups.
#' @param sd_per_group assumed standard deviation per group. Can be either a number if the same standard deviation is assumed for every group, or a vector with individual standard deviations.
#' @param alpha significance level, default is 0.05.
#' @param sides defines whether the one-sided (sides=1) or the two-sided (sides=2) p-value should be calculated, defaults to sides=2.
#'
#' @return A data frame with the input sample sizes in group 0 and group 1, \code{n0} and \code{n1}, the calculated power, and an indicator if a one-sided or two-sided null-hypothesis was assumed.
#'.
#'
#' @author Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @author Bettina Kroyer \email{bettina.kroyer@@meduniwien.ac.at}
#' @seealso \code{\link{n_t_test}}
#'
#' @examples
#' power_t_test(20,20,3,5,0.05) #sides=2 per default
#' power_t_test(20,20,3,5,0.05,1)
#' power_t_test(20,delta=3,sd_per_group=c(2,3))
#' @export
power_t_test<-function(n0,n1=n0,delta,sd_per_group,alpha=0.05,sides=2) {
  if (length(sd_per_group) == 1){
    K0 <- K1 <- sd_per_group^2
  } else if(length(sd_per_group) == 2){
    K0 <- sd_per_group[1]^2
    K1 <- sd_per_group[2]^2
  } else{
    stop("more standard deviations than groups given")
  }
  v <- K0/n0 + K1/n1
  SE <- sqrt(v)
  z <- delta / SE
  ncp <- abs(z)
  df <- (v^2) / (((K0/n0)^2 / (n0-1)) + ((K1/n1)^2 / (n1-1)))

  if(!(sides%in%c(1,2))) {
    stop("sides must be eiher 1 or 2")
  } else {
    krit<-qt(1-alpha/sides,df=df)
    power<-1-pt(krit,ncp=ncp,df=df)
  }

  data.frame(n0,n1,power,sides)
}

