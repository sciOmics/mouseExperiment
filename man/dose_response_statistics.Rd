\name{dose_response_statistics}
\alias{dose_response_statistics}
\title{Analyze Dose-Response Relationship in Tumor Growth Data}
\usage{
dose_response_statistics(
  df,
  dose_column = "Dose",
  treatment_column = "Treatment",
  volume_column = "Volume",
  day_column = "Day",
  id_column = "ID",
  time_point = NULL,
  control_group_name = "Control"
)
}
\arguments{
\item{df}{A data frame containing tumor growth and dose data}

\item{dose_column}{The name of the column containing dose concentrations}

\item{treatment_column}{The name of the column containing treatment names}

\item{volume_column}{The name of the column storing tumor volume measurements}

\item{day_column}{The name of the column with number of days since the beginning of the experiment}

\item{id_column}{The name of the column with individual mouse identifiers}

\item{time_point}{Optional specific time point (day) to analyze. If NULL, uses the last time point for each mouse.}

\item{control_group_name}{The name of the control group in the treatment column}
}
\value{
A list containing:
  \item{linear_model}{Linear regression model of dose vs. response}
  \item{non_linear_models}{List of fitted non-linear models (EC50, Emax, etc.)}
  \item{anova_model}{ANOVA model comparing dose groups}
  \item{plots}{List of data visualizations including dose-response curves}
  \item{summary_table}{Data frame summarizing results for each dose level}
  \item{statistics}{Key statistical outputs including p-values and model comparisons}
}
\description{
This function analyzes the relationship between drug dose levels and tumor growth 
using linear and non-linear regression models. It can detect both linear and non-linear
dose-response relationships and provides appropriate statistical tests and visualizations.
}
\details{
The function performs several analyses:
\enumerate{
  \item Linear regression to test for a linear dose-response relationship
  \item Non-linear regression to test for sigmoid/Hill equation relationships
  \item ANOVA to test for differences between dose groups
  \item For time-series data, analysis of dose effects on growth rate
}

Multiple models are compared using AIC/BIC to determine the best fit.
}
\examples{
\dontrun{
# Load and prepare data
data <- read.csv("dose_levels_synthetic_data.csv")
data <- calculate_volume(data)
data <- calculate_dates(data, start_date = "24-Mar", 
                       date_format = "\%d-\%b", year = 2023)
                       
# Analyze dose-response at the end of the experiment
results <- dose_response_statistics(data, 
                                  dose_column = "Dose",
                                  treatment_column = "Treatment")
                                  
# Analyze dose-response at a specific time point (day 14)
day14_results <- dose_response_statistics(data, 
                                        dose_column = "Dose",
                                        treatment_column = "Treatment",
                                        time_point = 14)
}
}