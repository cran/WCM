#' Simulate backscattering coeficient using WCM model
#' @description This function can be used to simulate the backscattering coeficient using WCM.
#' This function can be called in nls function for generationg model coeficients (A,B,C,D).
#' @param X In situ LAI or vegetation descriptor
#' @param Y In situ SM soil moisture
#' @param wcm_sim is simulate backscattering coeficient
#' @param theta incident angle of sattelite sensor
#' @param A fitted coeficinet for WCM using non linear least squre using in siu data
#' @param B fitted coeficinet for WCM using non linear least squre using in siu data
#' @param C fitted coeficinet for WCM using non linear least squre using in siu data
#' @param D fitted coeficinet for WCM using non linear least squre using in siu data
#' @import pracma
#' @import stats
#' @return simulated backscattering coeficient
#' @export
#'
#' @examples
#' # For single value.
#'  n <- wcm_sim(4,.3,48.9,-9.596695,-0.005331,-11.758309,0.011344)
#'
#' #For list of value
#' X<-c(2.6, 2.7, 2.8, 2.9, 3.1, 3.2, 2.9, 3.6, 4.2, 5, 5.1, 5.2, 3, 4.2, 5,
#'     4, 3.5, 4.1, 5.2, 6, 5, 5, 6, 5, 3.5, 4.1, 2.4, 3, 3.2, 3.5)
#' Y<-c(31.0, 25.0, 21.0, 15.0, 21.0, 28.0, 29.0, 28.5, 24.0, 27.0, 27.0, 28.0, 29, 28.0, 27.0,
#'      29.0, 28.0, 29.0, 19.0, 21.6, 19.0, 21.0, 31.0, 26.0, 29.0, 31.0, 38.0, 37.0, 36.0, 15.0)
#'
#' w<-c(-9.9, -11.2, -11.1, -11.9, -11.4, -10.9, -10.2, -10.1,  -7.200,  -6.750,
#'        -10.200,  -8.200,  -11.812,  -9.972,  -8.938,  -9.200,  -8.198,  -7.722,  -7.348,  -7.348,
#'        -8.198, -10.082,  -6.870,  -8.104,  -8.732,  -7.830, -11.686, -10.964, -12.976, -9.976)
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

