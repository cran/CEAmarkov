% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SD_main_func.R
\name{sd_sim_function}
\alias{sd_sim_function}
\title{Title
simulation function used in the sd calculation}
\usage{
sd_sim_function(trans_prob_mat, time_horizon, inital_state_vec,
  dead_state_num, cycle_length, ini_t)
}
\arguments{
\item{trans_prob_mat}{transition matrix}

\item{time_horizon}{null for lifetime horizon or number of years}

\item{inital_state_vec}{initial state vector}

\item{dead_state_num}{number of dead states}

\item{cycle_length}{cycle length in years}

\item{ini_t}{current cycle}
}
\value{
matrix
}
\description{
Title
simulation function used in the sd calculation
}
