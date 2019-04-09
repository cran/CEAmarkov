
#' Title
#' This function embeds the all-cause mortality probabilities in the transition matrices
#' @param trans_matrix_array transition matrices without the all-cause mortality probabilities
#' @param prob_of_dying_vector_array table containing the all-cause mortality probabilities
#' @param cycle_length cycle length in years
#'
#' @return trans_matrix_array
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
#  drug arm - 5 years
#' drug_arm_5_years<-matrix(data = c(0,0.029,0.087*0.59,0,0,0.1,0,0.261*0.59,
#' 0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#'
#  both arms after 20 years
#' both_arms_after_20_years<-matrix(data = c(0,0,0,0,0,0.1,0,0.261,0,0,0,0,
#' 0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#'
#  both arms after 5 years and before 20
#' after_5_years_drug_arm_no_drug_arm_5<-matrix(data = c(0,0.029,0.087,0,0,
#' 0.1,0,0.261,0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#'
#  if needed convert the probabilities according to the model's cycle length
#' drug_arm_5_years<-change_current_prob_to_cycle_prob(drug_arm_5_years,1,(1/12))
#' both_arms_after_20_years<-change_current_prob_to_cycle_prob(both_arms_after_20_years,1,(1/12))
#'
#  construct each arm matrices and replicate them according to the necessary cycles (each cycle has a matrix)
#' drug_arm<-array(c(replicate(n=5*12,drug_arm_5_years),replicate(n=15*12,
#' after_5_years_drug_arm_no_drug_arm_5),replicate(n=40*12,
#' both_arms_after_20_years)),dim = c(5,5,60*12))
#'
#  embed the probability of dying in the transition matrices
#' up_mat<-prepare_trans_matrix(drug_arm,prob_vec,(1/12))

prepare_trans_matrix <- function(trans_matrix_array, prob_of_dying_vector_array, cycle_length){
  num_of_arrays <- dim(trans_matrix_array)[3]
  prob_vector_length <- length(prob_of_dying_vector_array[, "cycle_prob"])
  #replicates the final prob to match the length of the transition matrices
  if (prob_of_dying_vector_array[prob_vector_length, "value"] == 0.999){
    final_prob <- prob_of_dying_vector_array[prob_vector_length, "cycle_prob"]
    rep_prob_vector <- rep(prob_of_dying_vector_array[, "cycle_prob"], each = 5 * cycle_length ^ (-1))
    rep_prob_vector_length <- length(rep_prob_vector)
    add_to_prob_vector <- rep(final_prob, each = num_of_arrays - rep_prob_vector_length)
    prob_of_dying_vector <- c(rep_prob_vector, add_to_prob_vector)
  }
  for (mat_num in 1:num_of_arrays) {
    #fills in the probability of dying and the probability of staying in the same state
    trans_matrix_array[,, mat_num] <- fill_trans_matrix_val(trans_matrix_array[,, mat_num], prob_of_dying_vector[mat_num])
  }
  return(trans_matrix_array)
}




#' Title
#' This function embeds the all-cause mortality probabilities in the transition matrices
#' @param trans_matrix transition matrix
#' @param prob_of_dying all-cause mortality probability
#'
#' @return trans_matrix
#'
#' @examples
fill_trans_matrix_val <- function(trans_matrix, prob_of_dying){
  mat_size <- nrow(trans_matrix)
  #fill the probability of dying in the right place
  for (cell in 1:mat_size) {
    #checks if the sum of probabilities in a row are 1
    if (sum(trans_matrix[cell, ]) != 1){
      #probability of staying in the same state
      trans_matrix[cell, cell] <- 1 - sum(trans_matrix[cell, ]) - prob_of_dying
      #probability of dying
      trans_matrix[cell, mat_size] <- prob_of_dying
    }
  }
  return(trans_matrix)
}
