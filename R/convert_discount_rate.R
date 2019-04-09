
#' Title
#' This function converts an annual discount rate to the proper discount rate
#' based on the cycle length of the model
#' @param annual_discount_rate annual discount rate in percentages
#' @param cycle_length cycle length in years
#'
#' @return new_discount_rate
#' @export
#'
#' @examples
#' @ annual discount rate in percentages
#' annual_discount_rate<-3
#' @ our model's cycle lenth in years
#' @ here we use monthly cycles
#' cycle_length<-(1/12)
#' covert_discount_rate(annual_discount_rate,cycle_length)
covert_discount_rate<-function(annual_discount_rate,cycle_length){
  new_discount_rate <- (1+(annual_discount_rate/100))^(cycle_length)-1
  return(new_discount_rate)
}
