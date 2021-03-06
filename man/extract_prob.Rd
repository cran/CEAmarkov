% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_prob.R
\name{extract_prob}
\alias{extract_prob}
\title{Title
This function extract all-cause mortality probability based on the model's parameters}
\usage{
extract_prob(initial_age, selected_country, time_horizon, gender,
  selected_year, cycle_length, prob_data)
}
\arguments{
\item{initial_age}{initial age}

\item{selected_country}{country}

\item{time_horizon}{null for lifetime horizon or number of years}

\item{gender}{gender}

\item{selected_year}{the year the data is extracted from}

\item{cycle_length}{cycle length in years}

\item{prob_data}{data containing probabilities, country and age}
}
\value{
five_year_prob
}
\description{
Title
This function extract all-cause mortality probability based on the model's parameters
}
