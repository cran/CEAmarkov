

#' Title
#' This function calculates the standard deviation of the model's outputs
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
#' \donttest{
#' prob_vec<-extract_prob(initial_age=45,selected_country='Ethiopia',
#' time_horizon=NULL,gender = 'Female',selected_year=2015,cycle_length = (1/12),
#' prob_data = sample_data)
#' drug_arm_5_years<-matrix(data = c(0,0.029,0.087*0.59,0,0,0.1,0,0.261*0.59,
#' 0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#' both_arms_after_20_years<-matrix(data = c(0,0,0,0,0,0.1,0,0.261,0,0,0,0,
#' 0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#' after_5_years_drug_arm_no_drug_arm_5<-matrix(data = c(0,0.029,0.087,0,0,
#' 0.1,0,0.261,0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#' drug_arm_5_years<-change_current_prob_to_cycle_prob(drug_arm_5_years,1,(1/12))
#' both_arms_after_20_years<-change_current_prob_to_cycle_prob(both_arms_after_20_years,1,(1/12))
#' after_5_years_drug_arm_no_drug_arm_5<-change_current_prob_to_cycle_prob(
#' after_5_years_drug_arm_no_drug_arm_5,1,(1/12))
#' drug_arm<-array(c(replicate(n=5*12,drug_arm_5_years),replicate(n=15*12,
#' after_5_years_drug_arm_no_drug_arm_5),replicate(n=40*12,
#' both_arms_after_20_years)),dim = c(5,5,60*12))
#' no_drug_arm<-array(c(replicate(n=20*12,after_5_years_drug_arm_no_drug_arm_5),
#' replicate(n=40*12,both_arms_after_20_years)),dim = c(5,5,60*12))
#' up_mat<-prepare_trans_matrix(drug_arm,prob_vec,(1/12))
#' up_mat_no_drug<-prepare_trans_matrix(no_drug_arm,prob_vec,(1/12))
#' cost_vec_ethiopia<-c(0,1323/12,1841/12,0,0)
#' cost_vec_12_months<-c(20000,1323/12,1841/12,0,0)
#' util_vec<-c(0.94,0.82,0.58,0,0)
#' total_cost_vec_no_drug_arm<-t(replicate(n=60*12,cost_vec_ethiopia))
#' total_util_vec_drug<-t(replicate(n=60*12,util_vec))
#' total_util_vec_no_drug<-total_util_vec_drug
#' total_cost_drug_arm<-t(cbind(replicate(n=1,cost_vec_12_months),
#' replicate(n=59*12+11,cost_vec_ethiopia)))
#' ini_vec<-c(1,0,0,0,0)
#' run_both_arms_sim(up_mat,up_mat_no_drug,total_util_vec_drug,total_util_vec_no_drug,
#' total_cost_drug_arm,total_cost_vec_no_drug_arm,NULL,ini_vec,ini_vec,2,2,
#' both_arms_res<-covert_discount_rate(3,1/12),(1/12))
#' SD_main_func(up_mat, up_mat_no_drug, total_util_vec_drug, total_util_vec_no_drug,
#' total_cost_drug_arm,total_cost_vec_no_drug_arm, NULL, ini_vec,ini_vec, 2, 2,
#' covert_discount_rate(3,1/12), (1/12))}
SD_main_func <- function(trans_prob_mat_drug, trans_prob_mat_no_drug, util_drug, util_no_drug, cost_drug, cost_no_drug, time_horizon, inital_state_vec_drug,inital_state_vec_no_drug, dead_state_num_drug,dead_state_num_no_drug, discount_rate, cycle_length){

  check_variables_compatibility(trans_prob_mat_drug, util_drug, cost_drug, inital_state_vec_drug, dead_state_num_drug)
  check_variables_compatibility(trans_prob_mat_no_drug, util_no_drug, cost_no_drug, inital_state_vec_no_drug, dead_state_num_no_drug)

  res_drug <- sim_function(trans_prob_mat_drug, util_drug, cost_drug, time_horizon, inital_state_vec_drug, dead_state_num_drug, discount_rate, cycle_length)
  res_drug_SD<-main_sd_function(res_drug$res_matrix,cost_drug,util_drug,trans_prob_mat_drug,dead_state_num_drug,cycle_length,res_drug$res_discount_rate)

  res_no_drug <- sim_function(trans_prob_mat_no_drug, util_no_drug, cost_no_drug, time_horizon, inital_state_vec_no_drug, dead_state_num_no_drug, discount_rate, cycle_length)
  res_no_drug_SD<-main_sd_function(res_no_drug$res_matrix,cost_no_drug,util_no_drug,trans_prob_mat_no_drug,dead_state_num_no_drug,cycle_length,res_no_drug$res_discount_rate)

  discounted_results <- matrix(data = c(as.numeric(res_no_drug$discounted_costs), res_no_drug$discounted_QALYs,res_no_drug$discounted_LYs,as.numeric(res_drug$discounted_costs),res_drug$discounted_QALYs, res_drug$discounted_LYs, sqrt(res_no_drug_SD$cost_sum_d-res_no_drug$discounted_costs^2),sqrt(res_no_drug_SD$QALY_sum_d-res_no_drug$discounted_QALYs^2), sqrt(res_no_drug_SD$LY_sum_d-res_no_drug$discounted_LYs^2), sqrt(res_drug_SD$cost_sum_d-res_drug$discounted_costs^2), sqrt(res_drug_SD$QALY_sum_d-res_drug$discounted_QALYs^2), sqrt(res_drug_SD$LY_sum_d-res_drug$discounted_LYs^2)), nrow = 2, ncol = 6, byrow = TRUE)

  non_discounted_results <- matrix(data = c(as.numeric(res_no_drug$costs), res_no_drug$QALYs,res_no_drug$LYs, as.numeric(res_drug$costs), res_drug$QALYs, res_drug$LYs,sqrt(res_no_drug_SD$cost_sum_no_d-res_no_drug$costs^2),sqrt(res_no_drug_SD$QALY_sum_no_d-res_no_drug$QALYs^2), sqrt(res_no_drug_SD$LY_sum_no_d-res_no_drug$LYs^2),sqrt(res_drug_SD$cost_sum_no_d-res_drug$costs^2),sqrt(res_drug_SD$QALY_sum_no_d-res_drug$QALYs^2), sqrt(res_drug_SD$LY_sum_no_d-res_drug$LYs^2)), nrow = 2, ncol = 6, byrow = TRUE)


  colnames(discounted_results) <-colnames(non_discounted_results) <- c("No Drug Costs","No Drug QALYs", "No Drug LYs","Drug Costs","Drug QALYs", "Drug LYs")
  rownames(discounted_results) <-rownames(non_discounted_results) <- c("Mean", "SD")

  discounted_ICER <- (res_drug$discounted_costs - res_no_drug$discounted_costs) / (res_drug$discounted_QALYs - res_no_drug$discounted_QALYs)
  non_discounted_ICER <- (res_drug$costs - res_no_drug$costs) / (res_drug$QALYs - res_no_drug$QALYs)

  return(list(discounted_results = discounted_results, discounted_ICER = discounted_ICER,non_discounted_results ,non_discounted_ICER=non_discounted_ICER))
}


#' Title
#' part one of the sd calculation
#' @param prob_matrix state probabilities for each time
#' @param discount_rate discount rate in percentage
#' @param currentT current time
#' @param variable LYs/QALYs/Costs
#'
#' @return list
#'
#' @examples
f1_sum<-function(prob_matrix,discount_rate,currentT,variable){
  f1_sum_no_discount<-sum(prob_matrix[currentT,]*variable[currentT,]^2)
  f1_sum_discount<-sum(prob_matrix[currentT,]*(variable[currentT,]/discount_rate[currentT])^2)
  return(list(f1_sum_no_discount=f1_sum_no_discount,f1_sum_discount=f1_sum_discount))
}



#' Title
#' part two of the sd calculation
#' @param prob_matrix state probabilities for each time
#' @param cost cost vector
#' @param util utility vector
#' @param trans_prob_mat probability matrix
#' @param dead_state_num number of dead states
#' @param cycle_length cycle length in years
#' @param discount_rate discount rate in percentages
#'
#' @return list
#'
#' @examples
main_sd_function<-function(prob_matrix,cost,util,trans_prob_mat,dead_state_num,cycle_length,discount_rate){
  T_final<-nrow(prob_matrix)
  currentT<-2

  ly_val_mat<-matrix(data = 0, ncol = ncol(prob_matrix), nrow = nrow(prob_matrix))
  ly_val_mat[,1:(ncol(prob_matrix)-dead_state_num)]<-(1/cycle_length^-1)

  cost_sum_f_1<-f1_sum(prob_matrix,discount_rate,currentT-1,cost)
  util_sum_f_1<-f1_sum(prob_matrix,discount_rate,currentT-1,util/cycle_length^-1)
  ly_sum_f_1<-f1_sum(prob_matrix,discount_rate,currentT-1,ly_val_mat)

  cost_sum_f_1_no_d<-cost_sum_f_1$f1_sum_no_discount
  cost_sum_f_1_d<-cost_sum_f_1$f1_sum_discount

  util_sum_f_1_no_d<-util_sum_f_1$f1_sum_no_discount
  util_sum_f_1_d<-util_sum_f_1$f1_sum_discount

  ly_sum_f_1_no_d<-ly_sum_f_1$f1_sum_no_discount
  ly_sum_f_1_d<-ly_sum_f_1$f1_sum_discount

  while (currentT<=T_final) {
    #a(n)
    cost_sum_a<-f1_sum(prob_matrix,discount_rate,currentT,cost)
    QALY_sum_a<-f1_sum(prob_matrix,discount_rate,currentT,util/cycle_length^-1)
    LY_sum_a<-f1_sum(prob_matrix,discount_rate,currentT,ly_val_mat)

    sum_b<-cond_prob_calc(prob_matrix,discount_rate,cost,util,currentT,trans_prob_mat,dead_state_num,cycle_length,ly_val_mat)

    cost_sum_f_1_no_d<-rbind(cost_sum_f_1_no_d,sum_b$final_res_cost_no_d+cost_sum_a$f1_sum_no_discount)
    cost_sum_f_1_d<-rbind(cost_sum_f_1_d,sum_b$final_res_cost_d+cost_sum_a$f1_sum_discount)

    util_sum_f_1_no_d<-rbind(util_sum_f_1_no_d,sum_b$final_res_util_no_d+QALY_sum_a$f1_sum_no_discount)
    util_sum_f_1_d<-rbind(util_sum_f_1_d,sum_b$final_res_util_d+QALY_sum_a$f1_sum_discount)

    ly_sum_f_1_no_d<-rbind(ly_sum_f_1_no_d,sum_b$final_res_ly_no_d+LY_sum_a$f1_sum_no_discount)
    ly_sum_f_1_d<-rbind(ly_sum_f_1_d,sum_b$final_res_ly_d+LY_sum_a$f1_sum_discount)

    currentT<-currentT+1
  }
  return(list(cost_sum_no_d = (sum(cost_sum_f_1_no_d)), cost_sum_d = (sum(cost_sum_f_1_d)), QALY_sum_no_d = (sum(util_sum_f_1_no_d)),QALY_sum_d = (sum(util_sum_f_1_d)), LY_sum_no_d = (sum(ly_sum_f_1_no_d)),LY_sum_d = (sum(ly_sum_f_1_d))))
}



#' Title
#' conditional probability sd calculation
#' @param prob_matrix state probabilities for each time
#' @param discount_rate discount rate in percentages
#' @param cost cost vector
#' @param util utility vector
#' @param currentT current time
#' @param trans_prob_mat transition matrix
#' @param dead_state_num number of dead states
#' @param cycle_length cycle length in years
#' @param ly_val_mat LYs matrix
#'
#' @return list
#'
#' @examples
cond_prob_calc<-function(prob_matrix,discount_rate,cost,util,currentT,trans_prob_mat,dead_state_num,cycle_length,ly_val_mat){
  #E[Ct,i/currentT-1]
  cond_prob_cost_1_d<-cond_prob_util_1_d<-cond_prob_ly_1_d<-matrix(data = 0,nrow = 0,ncol = ncol(prob_matrix))
  cond_prob_cost_1_no_d<-cond_prob_util_1_no_d<-cond_prob_ly_1_no_d<-matrix(data = 0,nrow = 0,ncol = ncol(prob_matrix))

  for (t in 1:(currentT-1)) {
    for (i in 1:ncol(prob_matrix)) {

      ini_vector<-vector(mode = 'numeric', length = ncol(prob_matrix))
      ini_vector[i]<-1
      cond_prob_1_temp<-sd_sim_function(trans_prob_mat,t+currentT-t-1,ini_vector,dead_state_num,cycle_length,t)
      cond_prob_temp<-(cond_prob_1_temp[nrow(cond_prob_1_temp)-1,])

      cond_prob_calc_cost_1<-cond_prob_calc_1(cost,discount_rate,prob_matrix,currentT,t,i,cond_prob_temp)
      cond_prob_cost_1_no_d<-rbind(cond_prob_cost_1_no_d,cond_prob_calc_cost_1$cond_prob_calc_1_no_d)
      cond_prob_cost_1_d<-rbind(cond_prob_cost_1_d,cond_prob_calc_cost_1$cond_prob_calc_1_d)


      cond_prob_calc_util_1<-cond_prob_calc_1(util/cycle_length^-1,discount_rate,prob_matrix,currentT,t,i,cond_prob_temp)
      cond_prob_util_1_no_d<-rbind(cond_prob_util_1_no_d,cond_prob_calc_util_1$cond_prob_calc_1_no_d)
      cond_prob_util_1_d<-rbind(cond_prob_util_1_d,cond_prob_calc_util_1$cond_prob_calc_1_d)

      cond_prob_calc_ly_1<-cond_prob_calc_1(ly_val_mat, discount_rate, prob_matrix, currentT, t, i, cond_prob_temp)
      cond_prob_ly_1_no_d<-rbind(cond_prob_ly_1_no_d,cond_prob_calc_ly_1$cond_prob_calc_1_no_d)
      cond_prob_ly_1_d<-rbind(cond_prob_ly_1_d,cond_prob_calc_ly_1$cond_prob_calc_1_d)


    }
  }
  cond_prob_cost_2_no_d<-cond_prob_util_2_no_d<-cond_prob_ly_2_no_d<-matrix(data = 0,nrow = 0,ncol = ncol(prob_matrix))
  cond_prob_cost_2_d<-cond_prob_util_2_d<-cond_prob_ly_2_d<-matrix(data = 0,nrow = 0,ncol = ncol(prob_matrix))

  for (i in 1:ncol(prob_matrix)) {
    ini_vector<-vector(mode = 'numeric', length = ncol(prob_matrix))
    ini_vector[i]<-1
    cond_prob_2_temp<-sd_sim_function(trans_prob_mat,currentT+1,ini_vector,dead_state_num,cycle_length,currentT)

    cond_prob_calc_cost_2<-cond_prob_calc_2(cost,discount_rate,prob_matrix,currentT,cond_prob_2_temp)
    cond_prob_cost_2_no_d<-rbind(cond_prob_cost_2_no_d,cond_prob_calc_cost_2$cond_prob_calc_2_no_d)
    cond_prob_cost_2_d<-rbind(cond_prob_cost_2_d,cond_prob_calc_cost_2$cond_prob_calc_2_d)

    cond_prob_calc_util_2<-cond_prob_calc_2(util/cycle_length^-1,discount_rate,prob_matrix,currentT,cond_prob_2_temp)
    cond_prob_util_2_no_d<-rbind(cond_prob_util_2_no_d,cond_prob_calc_util_2$cond_prob_calc_2_no_d)
    cond_prob_util_2_d<-rbind(cond_prob_util_2_d,cond_prob_calc_util_2$cond_prob_calc_2_d)

    cond_prob_calc_ly_2<-cond_prob_calc_2(ly_val_mat,discount_rate,prob_matrix,currentT,cond_prob_2_temp)
    cond_prob_ly_2_no_d<-rbind(cond_prob_ly_2_no_d,cond_prob_calc_ly_2$cond_prob_calc_2_no_d)
    cond_prob_ly_2_d<-rbind(cond_prob_ly_2_d,cond_prob_calc_ly_2$cond_prob_calc_2_d)

  }


  final_res_cost_no_d<-final_res(prob_matrix,currentT,cond_prob_cost_1_no_d,cond_prob_cost_2_no_d)
  final_res_cost_d<-final_res(prob_matrix,currentT,cond_prob_cost_1_d,cond_prob_cost_2_d)

  final_res_util_no_d<-final_res(prob_matrix,currentT,cond_prob_util_1_no_d,cond_prob_util_2_no_d)
  final_res_util_d<-final_res(prob_matrix,currentT,cond_prob_util_1_d,cond_prob_util_2_d)


  final_res_ly_no_d<-final_res(prob_matrix,currentT,cond_prob_ly_1_no_d,cond_prob_ly_2_no_d)
  final_res_ly_d<-final_res(prob_matrix,currentT,cond_prob_ly_1_d,cond_prob_ly_2_d)


  return(list(final_res_cost_no_d=final_res_cost_no_d,final_res_cost_d=final_res_cost_d,final_res_util_no_d=final_res_util_no_d,final_res_util_d=final_res_util_d,final_res_ly_no_d=final_res_ly_no_d,final_res_ly_d=final_res_ly_d))
}



#' Title
#' part one conditional probability sd calculation
#' @param variable LYs/QALYs/Costs
#' @param discount_rate discount rate in percentages
#' @param prob_matrix state probabilities for each time
#' @param currentT current time
#' @param t time t iterator
#' @param i state i iterator
#' @param cond_prob_temp temp conditional probability
#'
#' @return list
#'
#' @examples
cond_prob_calc_1<-function(variable,discount_rate,prob_matrix,currentT,t,i,cond_prob_temp){
  cond_prob_calc_1_no_d<-variable[t,i]*(1/prob_matrix[currentT-1,])*prob_matrix[t,i]*cond_prob_temp
  cond_prob_calc_1_no_d[is.infinite(cond_prob_calc_1_no_d) | is.nan(cond_prob_calc_1_no_d)]<-0

  cond_prob_calc_1_d<-(variable[t,i]/discount_rate[t])*(1/prob_matrix[currentT-1,])*prob_matrix[t,i]*cond_prob_temp
  cond_prob_calc_1_d[is.infinite(cond_prob_calc_1_d) | is.nan(cond_prob_calc_1_d)]<-0

  return(list(cond_prob_calc_1_no_d=cond_prob_calc_1_no_d,cond_prob_calc_1_d=cond_prob_calc_1_d))
}



#' Title
#' part two conditional probability sd calculation
#' @param variable LYs/QALYs/Costs
#' @param discount_rate discount rate in percentages
#' @param prob_matrix state probabilities for each time
#' @param currentT current time
#' @param cond_prob_temp temp conditional probability
#'
#' @return list
#'
#' @examples
cond_prob_calc_2<-function(variable,discount_rate,prob_matrix,currentT,cond_prob_temp){
  cond_prob_calc_2_no_d<-variable[currentT,]*cond_prob_temp[nrow(cond_prob_temp),]
  cond_prob_calc_2_d<-(variable[currentT,]/discount_rate[currentT])*cond_prob_temp[nrow(cond_prob_temp),]

  return(list(cond_prob_calc_2_no_d=cond_prob_calc_2_no_d,cond_prob_calc_2_d=cond_prob_calc_2_d))
}


#' Title
#' sum of conditional probability one and two
#' @param prob_matrix state probabilities for each time
#' @param currentT current time
#' @param cond_prob_1 temp conditional probability 1
#' @param cond_prob_2 temp conditional probability 2
#'
#' @return number
#'
#' @examples
final_res<-function(prob_matrix,currentT,cond_prob_1,cond_prob_2){
  final_res<-2*sum(prob_matrix[currentT,]*colSums(cond_prob_1)*colSums(cond_prob_2))
  return(final_res)
}



#' Title
#' simulation function used in the sd calculation
#' @param trans_prob_mat transition matrix
#' @param time_horizon null for lifetime horizon or number of years
#' @param inital_state_vec initial state vector
#' @param dead_state_num number of dead states
#' @param cycle_length cycle length in years
#' @param ini_t current cycle
#'
#' @return matrix
#'
#' @examples
sd_sim_function<-function(trans_prob_mat,time_horizon,inital_state_vec,dead_state_num,cycle_length,ini_t){
  num_of_states<-dim(trans_prob_mat[,,1])[1]
  #creates the results matrix according to the number of states
  res_matrix<-matrix(data = 0,nrow=0,ncol = num_of_states)
  #adds the inital state vector to the results matrix
  res_matrix<-rbind(res_matrix,inital_state_vec)
  cyc<-ini_t
  flag<-TRUE
  #runs while the sum of dead states isn't 1 or until cycle = time horizon
  while(sum(inital_state_vec[(length(inital_state_vec)-dead_state_num+1):length(inital_state_vec)])<=0.99 & flag){
    if(!is.null(time_horizon)){
      # if time horizon is defined and it was reached change flag to false
      if(cyc>=time_horizon){
        flag<-FALSE
      }
    }
    #multiply the state vector with the transition matrix
    inital_state_vec<-(inital_state_vec %*% trans_prob_mat[,,cyc])
    #add the result to the result matrix
    res_matrix<-rbind(res_matrix,inital_state_vec)
    cyc<-cyc+1
    #calculate the cycle total discount rate

  }
  #calculate final results
  return(res_matrix)
}


