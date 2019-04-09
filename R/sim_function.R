
#' Title
#' This function runs the simulation using the given parameters
#' @param trans_prob_mat transition matrix
#' @param util utility vector
#' @param cost cost vector
#' @param time_horizon null for lifetime horizon or number of years
#' @param inital_state_vec initial state vector
#' @param dead_state_num number of dead states
#' @param discount_rate discount rate in percentages
#' @param cycle_length cycle length in years
#'
#' @return list
#' @export
#'
#' @examples
#  get the relevent probabilities according to:
#  inital age, country, time horizon gender, year to extract data from, cycle length, API data
#' prob_vec<-extract_prob(initial_age=45,selected_country='Ethiopia',time_horizon=NULL,
#' gender = 'Female',selected_year=2015,cycle_length = (1/12),prob_data = sample_data)
#'
#  set up transition matrices
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
#'
#  embed the probability of dying in the transition matrices
#' up_mat<-prepare_trans_matrix(drug_arm,prob_vec,(1/12))
#'
#  prepare the costs and utilities matrix
#  costs per cycle
#' cost_vec_ethiopia<-c(0,1323/12,1841/12,0,0)
#' cost_vec_12_months<-c(20000,1323/12,1841/12,0,0)
#'
#  utilities per cycle
#' util_vec<-c(0.94,0.82,0.58,0,0)
#'
#  replicate for each cycle
#' total_util_vec_drug<-t(replicate(n=60*12,util_vec))
#' total_cost_drug_arm<-t(cbind(replicate(n=1,cost_vec_12_months),
#' replicate(n=59*12+11,cost_vec_ethiopia)))
#'
#  define the state the model starts from
#' ini_vec<-c(1,0,0,0,0)
#'
# run simulation
#' sim_function(up_mat,total_util_vec_drug,total_cost_drug_arm,NULL,ini_vec,2,(1+0.03)^(1/12)-1,(1/12))
#prepare trans_matrix
sim_function <- function(trans_prob_mat, util, cost, time_horizon,
                         inital_state_vec, dead_state_num, discount_rate,
                         cycle_length){

  cyc <- 0
  num_of_states <- dim(trans_prob_mat[, , 1])[1]
  #creates the results matrix according to the number of states
  res_matrix <- matrix (data = 0, nrow = 0, ncol = num_of_states)
  #current cycle discount rate
  res_discount_rate <- c ( (1 + discount_rate ) ^ cyc)
  #adds the inital state vector to the results matrix
  res_matrix <- rbind(res_matrix, inital_state_vec)
  cyc <- cyc + 1
  flag <- TRUE
  #runs while the sum of dead states isn't 1 or until cycle = time horizon
  while (sum(inital_state_vec[(length(inital_state_vec) - dead_state_num + 1):length(inital_state_vec)]) <= 0.99 & flag){
    if (!is.null(time_horizon)){
      time_horizon<-time_horizon/(cycle_length)
      # if time horizon is defined and it was reached change flag to false
      if (cyc >= time_horizon){
        flag <- FALSE
      }
    }
    #multiply the state vector with the transition matrix
    inital_state_vec <- (inital_state_vec %*% trans_prob_mat[, , cyc])
    #add the result to the result matrix
    res_matrix <- rbind(res_matrix, inital_state_vec)
    cyc <- cyc + 1
    #calculate the cycle total discount rate
    res_discount_rate <- c(res_discount_rate, (1 + discount_rate) ^ cyc)

  }
  #calculate final results
  res_mat_discounted <- (res_matrix) / (res_discount_rate)
  sum_res_mat <- rowSums(res_matrix[, 1:(length(inital_state_vec) - dead_state_num)]) / ( (cycle_length) ^ -1)
  LYs <- sum( (sum_res_mat))
  QALYs <- sum(res_matrix * util[1:nrow(res_matrix), ]) / ( (cycle_length) ^ -1)
  costs <- sum(res_matrix * cost[1:nrow(res_matrix), ])
  discounted_LYs <- sum( (sum_res_mat) / res_discount_rate)
  discounted_QALYs <- sum(res_matrix * (util[1:nrow(res_matrix), ]) / res_discount_rate) / ( (cycle_length) ^ -1)
  discounted_costs <- sum(res_matrix * (cost[1:nrow(res_matrix), ]) / res_discount_rate)
  return(list(LYs = LYs, QALYs = QALYs, costs = costs, discounted_LYs = discounted_LYs,
              discounted_QALYs = discounted_QALYs, discounted_costs = discounted_costs,
              res_matrix = res_matrix, res_discount_rate = res_discount_rate,
              res_mat_discounted = res_mat_discounted))
}


