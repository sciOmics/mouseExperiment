#' Calculate tumor volume using length and width measurements
#'
#' @param df A data frame containing tumor measurements
#' @param length_column The name of the column storing tumor length measurements
#' @param width_column The name of the column storing tumor width measurements
#' @param height_column The name of the column storing tumor height measurements (optional, for 3D formulas)
#' @param formula The formula to use for volume calculation. Options:
#'   \itemize{
#'     \item "ellipsoid" (default): V = (length × width² × π) / 6 (assuming height not provided)
#'     \item "modified_ellipsoid": V = (length × width²) / 2
#'     \item "ellipsoid_3axis": V = (length × width × height × π) / 6 (requires height_column)
#'     \item "cylinder": V = (π × width² × length) / 4
#'     \item "sphere": V = (π × width³) / 6 (uses only the width/diameter)
#'     \item "box": V = length × width × height (requires height_column)
#'   }
#' @param in_place Logical, whether to modify the input data frame (TRUE) or return a new data frame (FALSE, default)
#'
#' @return Returns a data frame with a "Volume" column appended, calculated using the specified formula.
#'   If in_place=FALSE (default), returns a new data frame; if in_place=TRUE, modifies the input data frame and returns it.
#' @export
#'
#' @details
#' This function calculates tumor volume using different formulas common in preclinical oncology research.
#' For two-dimensional measurements (length and width), the function ensures the longer dimension
#' is used as length and the shorter as width, regardless of how they are labeled in the input data.
#' 
#' When 3D formulas are used and the height_column is not provided, the height is estimated as equal to
#' the width (a common approximation in models like ellipsoid_3axis).
#'
#' Common volume estimation formulas in the literature include:
#' \itemize{
#'   \item Ellipsoid (default): V = (length × width² × π) / 6 
#'     - Most accurate for ovoid tumors when true height isn't measured
#'     - Widely used in xenograft studies
#'   \item Modified ellipsoid: V = (length × width²) / 2
#'     - Simplified formula omitting π/3
#'     - Used in many older studies
#'   \item Ellipsoid with 3 axes: V = (length × width × height × π) / 6
#'     - Most accurate when all three dimensions can be measured
#'   \item Cylinder: V = (π × width² × length) / 4
#'     - Used for more cylindrical tumors
#'   \item Sphere: V = (π × width³) / 6
#'     - For spherical tumors where only diameter is measured
#'   \item Box: V = length × width × height
#'     - Simple approximation when tumors have more rectangular shape
#' }
#'
#' @examples
#' \dontrun{
#' # Return a new data frame with volume using default ellipsoid formula
#' df_with_volume <- calculate_volume(synthetic_data)
#' 
#' # Use a different formula (modified ellipsoid)
#' df_with_volume <- calculate_volume(synthetic_data, formula = "modified_ellipsoid")
#' 
#' # Modify the original data frame in place
#' synthetic_data <- calculate_volume(synthetic_data, in_place = TRUE)
#' 
#' # Use a 3D formula with height measurements
#' df_with_volume <- calculate_volume(synthetic_data, 
#'                                   height_column = "Height",
#'                                   formula = "ellipsoid_3axis")
#' }
calculate_volume <- function(df, length_column = "Length", width_column = "Width", 
                            height_column = NULL, formula = "ellipsoid", in_place = FALSE) {
  # Input validation
  if(!length_column %in% colnames(df)) {
    stop(paste("Column", length_column, "not found in data frame"))
  }
  if(!width_column %in% colnames(df)) {
    stop(paste("Column", width_column, "not found in data frame"))
  }
  
  # Validate and check height column if needed for 3D formulas
  need_height <- formula %in% c("ellipsoid_3axis", "box")
  if(need_height && is.null(height_column)) {
    message("Note: Height column not provided. Using width as approximation for height.")
  } else if(need_height && !is.null(height_column) && !height_column %in% colnames(df)) {
    stop(paste("Column", height_column, "not found in data frame but required for", formula, "formula."))
  }
  
  # Choose whether to modify in place or create a copy
  if (in_place) {
    df_result <- df  # Reference to the original dataframe
  } else {
    df_result <- df  # Create a copy of the data frame
  }
  
  # Ensure longer dimension is always stored in length_column 
  # and shorter dimension in width_column
  length_tmp <- ifelse(df[,length_column] > df[,width_column], 
                       df[,length_column], df[,width_column])
  width_tmp <- ifelse(df[,width_column] < df[,length_column], 
                      df[,width_column], df[,length_column])
  
  # Get height if needed, otherwise use width as approximation
  if(!is.null(height_column) && height_column %in% colnames(df)) {
    height_tmp <- df[,height_column]
  } else {
    height_tmp <- width_tmp  # Common approximation when height not measured
  }
  
  # Calculate tumor volume using selected formula
  Volume <- switch(formula,
    "modified_ellipsoid" = (length_tmp * width_tmp^2) / 2,
    "ellipsoid" = (length_tmp * width_tmp^2 * pi) / 6,
    "ellipsoid_3axis" = (length_tmp * width_tmp * height_tmp * pi) / 6,
    "cylinder" = (pi * width_tmp^2 * length_tmp) / 4,
    "sphere" = (pi * width_tmp^3) / 6,
    "box" = length_tmp * width_tmp * height_tmp,
    # Default to ellipsoid if formula not recognized
    (length_tmp * width_tmp^2 * pi) / 6
  )
  
  # Add volume to result dataframe
  result <- cbind(df_result, Volume)
  
  # If in_place is TRUE and we're in a function, update the original df in the parent environment
  if (in_place) {
    parent_frame <- parent.frame()
    df_name <- deparse(substitute(df))
    if (exists(df_name, envir = parent_frame)) {
      assign(df_name, result, envir = parent_frame)
    }
  }
  
  # Return the result
  return(result)
}
