#' Simulate backscattering coefficient using WCM model
#' @description This function can be used to simulate the backscattering coefficient using WCM.
#' This function can be called in nls function for generation of model coefficients (A,B,C,D).
#' @param X In situ LAI or vegetation descriptor
#' @param Y In situ SM soil moisture
#' @param wcm_sim is simulated backscattering coefficient
#' @param theta incident angle of Satellite sensor
#' @param A fitted coefficient for WCM using non linear least squre using in situ data
#' @param B fitted coefficient for WCM using non linear least squre using in situ data
#' @param C fitted coefficient for WCM using non linear least squre using in situ data
#' @param D fitted coefficient for WCM using non linear least squre using in situ data
#' @import pracma
#' @import stats
#' @return simulated backscattering coefficient
#' @export
#'
#' @examples
#' # For single value.
#'  n <- wcm_sim(4,.3,48.9,-9.596695,-0.005331,-11.758309,0.011344)
#'
#' #For list of value
#' X<-c(5.34, 4.34, 4.32, 4.12, 4.17, 3.58, 5.39, 5.66, 5.47, 5.73, 5.76, 5.93, 4.91, 5.36, 6.15,
#'      4.56, 5.44, 6.54, 6.20, 6.34, 5.56, 5.88, 7.34, 5.74, 4.81, 5.73, 3.63, 4.61, 4.76, 4.02)
#' Y<-c(35.0, 26.0, 18.0, 13.0, 18.0, 22.0, 19.0, 16.5, 20.0, 24.0, 24.0, 21.0, 13.0, 22.0, 25.0,
#'      24.0, 30.0, 23.0, 18.0, 17.6, 15.0, 17.0, 27.0, 22.0, 21.0, 15.0, 15.0, 18.0, 31.0, 10.0)
#'
#' w<-c(-9.604, -11.648, -11.556, -11.556, -11.090, -10.444, -10.444, -10.042,  -9.200,  -9.750,
#'        -9.200,  -9.200,  -9.812,  -9.972,  -8.938,  -9.200,  -8.198,  -7.722,  -7.348,  -7.348,
#'        -8.198, -10.082,  -6.870,  -8.104,  -8.732,  -7.830, -10.686, -10.964, -10.976, -10.976)
#'
#'theta<-48.9
#'\donttest{nlc<-nls.control(maxiter = 50000, tol = 1e-05, minFactor = 1/100000000000,
#' printEval = FALSE, warnOnly = FALSE)}
#' \donttest{k<-nls(w~wcm_sim(X,Y,theta,A,B,C,D),control=nlc,
#'  start=list(A= 0.01,B=0.01,C=-21,D= 0.00014),trace = T)}
#'\donttest{y<-predict(k)}
#'n <- wcm_sim(X,Y,theta,-9.596695,-0.005331,-11.758309,0.011344)

wcm_sim<-function(X,Y,theta,A, B, C, D)
{
  theta<- theta*pi/180.
  return(A*X*cos(theta)*(1- exp(- 2*B*X/cos(theta))) + exp( - 2*B*X/cos(theta))*(C + D*Y))
}

