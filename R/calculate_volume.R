#' Calculate tumor volume using length and width measurements
#'
#' @param df A data frame containing tumor measurements
#' @param length_column The name of the column storing tumor length measurements
#' @param width_column The name of the column storing tumor width measurements
#'
#' @return Returns the original data frame with a tumor volume column appended
#' @export
#'
#' @examples calculate_volume(synthetic_data)

calculate_volume = function(df, length_column = "Length", width_column = "Width") {
  length_tmp = ifelse(df[,length_column] > df[,width_column], df[,length_column], df[,width_column])
  width_tmp = ifelse(df[,width_column] < df[,length_column], df[,width_column], df[,length_column])
  df[,length_column] = length_tmp
  df[,width_column] = width_tmp
  Volume = df[,length_column] * (df[,width_column]^2)/2
  cbind(df, Volume)
}
