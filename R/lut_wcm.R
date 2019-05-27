
#' Look up table of WCM
#'
#' @param LAI one dimensional row vector or a range of LAI value
#' @param SM one dimensionalrow vector or a range of SM value
#' @param coeff Generated A, B, C, D fitted coefficient for WCM using non linear least square using in situ data
#' @import pracma
#' @return look up table for WCM for given range of LAI and SM
#' @export
#' @examples
#' A= -9.596695
#' B=-0.005331
#' C=-11.758309
#' D=0.011344
#' lookuptable <- lut_wcm(LAI=seq(1,6,0.1), SM=seq(0,.6,.01),coeff=c(A,B,C,D))

lut_wcm <- function(LAI, SM, coeff){
  sig <- array(vector(), c(NROW(LAI), NROW(SM)))
  dfAll <- data.frame(matrix(vector(), nrow=0, ncol=3))
  for (n in 1:NROW(LAI)){
    for (m in 1:NROW(SM)){
      sig[n,m] <- coeff[1]*LAI[n]*cosd(48.9)*(1- exp(- 2*coeff[2]*LAI[n]*secd(48.9))) + exp( - 2*coeff[2]*LAI[n]* secd(48.9))*(coeff[3] + coeff[4]*SM[m])
      df <- cbind(LAI[n], SM[m],sig[n,m])
      dfAll <- rbind(dfAll,df)

    }
  }
  names(dfAll)<- c('LAI', 'SM','sig')
  return(dfAll)
}
