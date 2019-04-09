
#' Title
#' This function extract all-cause mortality probability based on the model's parameters
#' @param initial_age initial age
#' @param selected_country country
#' @param time_horizon null for lifetime horizon or number of years
#' @param gender gender
#' @param selected_year the year the data is extracted from
#' @param cycle_length cycle length in years
#' @param prob_data data containing probabilities, country and age
#'
#' @import WHO
#' @import xml2
#'
#' @return five_year_prob
#' @export
#'
#' @examples
#' @ get the relevent probabilities according to:
#' @ inital age, country, time horizon gender, year to extract data from, cycle length, WHO API data
#' prob_vec<-extract_prob(initial_age=45,selected_country='Ethiopia',time_horizon=NULL,
#' gender = 'Female',selected_year=2015,cycle_length = (1/12),prob_data = sample_data)

extract_prob <- function(initial_age, selected_country, time_horizon, gender, selected_year, cycle_length,prob_data){
  filtered_data <- subset(prob_data, prob_data$sex == gender & prob_data$country == selected_country)
  #if a year wasn't chosen get the most recent year
  if (is.null(selected_year)){
    selected_year <- max(filtered_data$year)
  }
  filtered_data <- subset(filtered_data, filtered_data$year == selected_year)
  #exteact by inital age and time horizon
  five_year_prob <- extract_by_age(initial_age, time_horizon, filtered_data)
  #converts the 5 year prob to the cycle prob (in yearly terms)
  five_year_prob$cycle_prob <- change_current_prob_to_cycle_prob(five_year_prob[, "value"], 5, cycle_length)
  return(five_year_prob)
}



#' Title
#' This function extract the all-cause mortality probabilities
#' using the inital age and time-horizon
#' @param initial_age initial age
#' @param time_horizon null for lifetime horizon or number of years
#' @param data_table data frame containing probabilities
#'
#' @return needed_data
#'
#' @import stringr
#'
#' @examples
extract_by_age <- function(initial_age, time_horizon, data_table){
  #trasform the age into two separate columns
  first_split <- gsub("<", "0-", data_table$agegroup)
  second_split <- gsub("[+]", "-120", first_split)
  third_split <- str_split_fixed(second_split, " ", 2)
  forth_split <- str_split_fixed(third_split[, 1], "-", 2)
  #converts to numeric
  final_col <- cbind(as.numeric(forth_split[, 1]), as.numeric(forth_split[, 2]))
  colnames(final_col) <- c("from", "to")
  new_data <- cbind(final_col, data_table)
  #order by age
  ordered_data <- new_data[order(new_data$from), ]
  #get the number of the row we need to start from
  row_to_start <- which(ordered_data$from <= initial_age & ordered_data$to >= initial_age)
  row_to_finish <- which(ordered_data$from <= (initial_age + time_horizon) & ordered_data$to >= (initial_age + time_horizon))
  if (is.null(time_horizon)){
    row_to_finish <- nrow(ordered_data)
  }
  needed_data <- ordered_data[row_to_start:row_to_finish, c("from", "to", "value")]
  if (needed_data[nrow(needed_data), "value"] == 1){
    needed_data[nrow(needed_data), "value"] <- 0.999
  }
  return(needed_data)
}



#' Title
#' This function converts the given probabilities to match the model's cycle length
#' @param prob probability
#' @param year_prob current terms of the given probability in years
#' @param cycle_length cycle length in years
#'
#' @return cycle_prob
#' @export
#'
#' @examples
#' drug_arm_5_years<-matrix(data = c(0,0.029,0.087*0.59,0,0,0.1,0,0.261*0.59,0,0,0,0,0,
#' 0.325,0,0,0,0,1,0,0,0,0,0,1),nrow = 5,ncol = 5,byrow = TRUE)
#' drug_arm_5_years<-change_current_prob_to_cycle_prob(drug_arm_5_years,1,(1/12))

change_current_prob_to_cycle_prob <- function(prob, year_prob, cycle_length){
  #converts the given prob to the cycle prob (in yearly terms)
  cycle_prob <- 1 - (1 - prob) ^ ( (cycle_length) / (year_prob))
  return(cycle_prob)
}
