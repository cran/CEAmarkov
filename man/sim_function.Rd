% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_function.R
\name{sim_function}
\alias{sim_function}
\title{Title
This function runs the simulation using the given parameters}
\usage{
sim_function(trans_prob_mat, util, cost, time_horizon, inital_state_vec,
  dead_state_num, discount_rate, cycle_length)
}
\arguments{
\item{trans_prob_mat}{transition matrix}

\item{util}{utility vector}

\item{cost}{cost vector}

\item{time_horizon}{null for lifetime horizon or number of years}

\item{inital_state_vec}{initial state vector}

\item{dead_state_num}{number of dead states}

\item{discount_rate}{discount rate in percentages}

\item{cycle_length}{cycle length in years}
}
\value{
list
}
\description{
Title
This function runs the simulation using the given parameters
}
\examples{
prob_vec<-extract_prob(initial_age=45,selected_country='Ethiopia',time_horizon=NULL,
gender = 'Female',selected_year=2015,cycle_length = (1/12),prob_data = sample_data)

drug_arm_5_years<-matrix(data = c(0,0.029,0.087*0.59,0,0,0.1,0,0.261*0.59,
0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
both_arms_after_20_years<-matrix(data = c(0,0,0,0,0,0.1,0,0.261,0,0,0,0,
0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
after_5_years_drug_arm_no_drug_arm_5<-matrix(data = c(0,0.029,0.087,0,0,
0.1,0,0.261,0,0,0,0,0,0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)

drug_arm_5_years<-change_current_prob_to_cycle_prob(drug_arm_5_years,1,(1/12))
both_arms_after_20_years<-change_current_prob_to_cycle_prob(both_arms_after_20_years,1,(1/12))
after_5_years_drug_arm_no_drug_arm_5<-change_current_prob_to_cycle_prob(
after_5_years_drug_arm_no_drug_arm_5,1,(1/12))

drug_arm<-array(c(replicate(n=5*12,drug_arm_5_years),replicate(n=15*12,
after_5_years_drug_arm_no_drug_arm_5),replicate(n=40*12,
both_arms_after_20_years)),dim = c(5,5,60*12))

up_mat<-prepare_trans_matrix(drug_arm,prob_vec,(1/12))

cost_vec_ethiopia<-c(0,1323/12,1841/12,0,0)
cost_vec_12_months<-c(20000,1323/12,1841/12,0,0)

util_vec<-c(0.94,0.82,0.58,0,0)

total_util_vec_drug<-t(replicate(n=60*12,util_vec))
total_cost_drug_arm<-t(cbind(replicate(n=1,cost_vec_12_months),
replicate(n=59*12+11,cost_vec_ethiopia)))

ini_vec<-c(1,0,0,0,0)

sim_function(up_mat,total_util_vec_drug,total_cost_drug_arm,NULL,ini_vec,2,(1+0.03)^(1/12)-1,(1/12))
}
