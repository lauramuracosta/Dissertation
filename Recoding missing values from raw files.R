## Set your working directory to where the CSV file is located
#setwd("##")

# Read the data
########df <- read.csv('##.csv')

# Function to replace the specified values
replace_values <- function(df) {
  # For each column in the dataframe
  for(col in names(df)) {
    # Replace "MISSING_DATA" with "0" (as character)
    df[[col]][df[[col]] == "MISSING_DATA"] <- "0"
    
    # Replace "UNASKED_QUESTION" with -99 (as numeric)
    df[[col]][df[[col]] == "UNASKED_QUESTION"] <- -99
    
    # Replace NA values with -99
    df[[col]][is.na(df[[col]])] <- -99
  }
  
  return(df)
}

# Apply the function to your data frame
df_cleaned <- replace_values(df)

# Save the results if needed
#######write.csv(df_cleaned, "##.csv", row.names = FALSE)
