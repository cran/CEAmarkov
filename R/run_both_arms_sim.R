#' This is data to be included in my package
#'
#' @name sample_data
#' @docType data
#' @author
#' @references http://apps.who.int/gho/data/?theme=main&vid=60550
#' @keywords WHO, all-cause probability mortality

#' Title
#' This function runs the simulation for both arms
#' @param trans_prob_mat_drug transition matrix for the drug arm
#' @param trans_prob_mat_no_drug transition matrix for the no drug arm
#' @param util_drug utility vector for the drug arm
#' @param util_no_drug utility vector for the no drug arm
#' @param cost_drug cost vector for the drug arm
#' @param cost_no_drug cost vector for the no drug arm
#' @param time_horizon null for lifetime horizon or number of years
#' @param inital_state_vec_drug initial state vector for the drug arm
#' @param inital_state_vec_no_drug initial state vector for the no drug arm
#' @param dead_state_num_drug number of dead states in the drug arm
#' @param dead_state_num_no_drug number of dead states in the no drug arm
#' @param discount_rate discount rate in percentages
#' @param cycle_length cycle length in years
#'
#' @return list
#' @export
#'
#' @examples
#  get the relevent probabilities according to:
#  inital age, country, time horizon gender, year to extract data from, cycle length, API data
#' prob_vec<-extract_prob(initial_age=45,selected_country='Ethiopia',
#' time_horizon=NULL,gender = 'Female',selected_year=2015,cycle_length = (1/12),
#' prob_data = sample_data)
#'
# set up transition matrices
#  yearly tranistion probabilities without all cause mortality
#' drug_arm_5_years<-matrix(data = c(0,0.029,0.087*0.59,0,0,0.1,0,0.261*0.59,
#' 0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#' both_arms_after_20_years<-matrix(data = c(0,0,0,0,0,0.1,0,0.261,0,0,0,0,
#' 0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#' after_5_years_drug_arm_no_drug_arm_5<-matrix(data = c(0,0.029,0.087,0,0,
#' 0.1,0,0.261,0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#'
#  if needed convert the probabilities according to the model's cycle length
#' drug_arm_5_years<-change_current_prob_to_cycle_prob(drug_arm_5_years,1,(1/12))
#' both_arms_after_20_years<-change_current_prob_to_cycle_prob(both_arms_after_20_years,1,(1/12))
#' after_5_years_drug_arm_no_drug_arm_5<-change_current_prob_to_cycle_prob(
#' after_5_years_drug_arm_no_drug_arm_5,1,(1/12))
#'
#  construct each arm matrices and replicate them according to the necessary cycles (each cycle has a matrix)
#' drug_arm<-array(c(replicate(n=5*12,drug_arm_5_years),replicate(n=15*12,
#' after_5_years_drug_arm_no_drug_arm_5),replicate(n=40*12,
#' both_arms_after_20_years)),dim = c(5,5,60*12))
#' no_drug_arm<-array(c(replicate(n=20*12,after_5_years_drug_arm_no_drug_arm_5),
#' replicate(n=40*12,both_arms_after_20_years)),dim = c(5,5,60*12))
#'
#  embed the probability of dying in the transition matrices
#' up_mat<-prepare_trans_matrix(drug_arm,prob_vec,(1/12))
#' up_mat_no_drug<-prepare_trans_matrix(no_drug_arm,prob_vec,(1/12))
#'
#  prepare the costs and utilities matrix
#  costs per cycle
#' cost_vec_ethiopia<-c(0,1323/12,1841/12,0,0)
#' cost_vec_12_months<-c(20000,1323/12,1841/12,0,0)
#'
#  utilities per cycle
#' util_vec<-c(0.94,0.82,0.58,0,0)
#' total_cost_vec_no_drug_arm<-t(replicate(n=60*12,cost_vec_ethiopia))
#' total_util_vec_drug<-t(replicate(n=60*12,util_vec))
#' total_util_vec_no_drug<-total_util_vec_drug
#'
#  replicate for each cycle
#' total_cost_drug_arm<-t(cbind(replicate(n=1,cost_vec_12_months),
#' replicate(n=59*12+11,cost_vec_ethiopia)))
#'
#  define the state the model starts from
#' ini_vec<-c(1,0,0,0,0)
#'
#  run simulation for both arms
#' run_both_arms_sim(up_mat,up_mat_no_drug,total_util_vec_drug,total_util_vec_no_drug,
#' total_cost_drug_arm,total_cost_vec_no_drug_arm,NULL,ini_vec,ini_vec,2,2,
#' covert_discount_rate(3,1/12),(1/12))


run_both_arms_sim <- function(trans_prob_mat_drug, trans_prob_mat_no_drug, util_drug, util_no_drug, cost_drug, cost_no_drug, time_horizon, inital_state_vec_drug,inital_state_vec_no_drug, dead_state_num_drug,dead_state_num_no_drug, discount_rate, cycle_length){
  check_variables_compatibility(trans_prob_mat_drug, util_drug, cost_drug, inital_state_vec_drug, dead_state_num_drug)
  check_variables_compatibility(trans_prob_mat_no_drug, util_no_drug, cost_no_drug, inital_state_vec_no_drug, dead_state_num_no_drug)

  res_drug <- sim_function(trans_prob_mat_drug, util_drug, cost_drug, time_horizon, inital_state_vec_drug, dead_state_num_drug, discount_rate, cycle_length)
  res_no_drug <- sim_function(trans_prob_mat_no_drug, util_no_drug, cost_no_drug, time_horizon, inital_state_vec_no_drug, dead_state_num_no_drug, discount_rate, cycle_length)
  res_table_discounted <- matrix(data = c(res_drug$discounted_costs, res_drug$discounted_LYs, res_drug$discounted_QALYs, res_no_drug$discounted_costs, res_no_drug$discounted_LYs, res_no_drug$discounted_QALYs), nrow = 2, ncol = 3, byrow = TRUE)
  res_table_non_discounted <- matrix(data = c(res_drug$costs, res_drug$LYs, res_drug$QALYs, res_no_drug$costs, res_no_drug$LYs, res_no_drug$QALYs), nrow = 2, ncol = 3, byrow = TRUE)
  colnames(res_table_discounted) <- colnames(res_table_non_discounted) <- c("Costs", "LYs", "QALYs")
  rownames(res_table_discounted) <- rownames(res_table_non_discounted) <- c("Drug", "No Drug")
  discounted_ICER <- (res_drug$discounted_costs - res_no_drug$discounted_costs) / (res_drug$discounted_QALYs - res_no_drug$discounted_QALYs)
  non_discounted_ICER <- (res_drug$costs - res_no_drug$costs) / (res_drug$QALYs - res_no_drug$QALYs)
  return(list(res_table_discounted = res_table_discounted, discounted_ICER = discounted_ICER,res_table_non_discounted = res_table_non_discounted, non_discounted_ICER = non_discounted_ICER))
}


#' Title
#' This function ensures the model's input are compatible
#' @param trans_prob_mat transition matrix
#' @param util utility vector
#' @param cost cost vector
#' @param inital_state_vec initial stste vector
#' @param dead_state_num number of dead states
#'
#' @return warning message
#'
#' @examples
check_variables_compatibility<-function(trans_prob_mat, util, cost, inital_state_vec, dead_state_num){
  if(dim(trans_prob_mat[, , 1])[1] != length(inital_state_vec)){
    stop("Transition Matrix dimensions does not match the initial state vector dimension")
  }
  if(dim(trans_prob_mat[, , 1])[1] != ncol(util)){
    stop("Transition Matrix dimensions does not match the utility dimensions")
  }
  if(dim(trans_prob_mat[, , 1])[1] != ncol(cost)){
    stop("Transition Matrix dimensions does not match the cost dimensions")
  }
  if(dim(trans_prob_mat[, , 1])[1] <= dead_state_num){
    stop("The number of Dead states cannot be equal or larger than the number of states")
  }
}
