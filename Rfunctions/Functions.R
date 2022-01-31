# Data Preprocessing ----
fix_conflict_multi <- function(conflicted_data, conflicted_vars, based_var, pattern_type, replace_value){
  based_var <- sym(based_var)
  # think about map functions purrr
  fixed_data <- conflicted_data
  for (i in seq(length(conflicted_vars))){
    conflicted_var <- conflicted_vars[i]
    conflicted_var <- sym(conflicted_var)
    fixed_data <- fixed_data %>% mutate(across(!!conflicted_var, ~replace(., str_detect(!!based_var, paste(pattern_type, collapse = '|')) & str_detect(!!conflicted_var, "No"), replace_value))) 
  }
  return(fixed_data)
}
#
fix_conflict_unknown <- function(conflicted_data, conflicted_vars, based_var){
  # conflicted_vars <- sym(conflicted_vars)
  based_var <- sym(based_var)
  fixed_data <- conflicted_data
  for (i in seq(length(conflicted_vars))){
    conflicted_var <- conflicted_vars[i]
    conflicted_var <- sym(conflicted_var)
    fixed_data <- fixed_data %>% mutate(across(!!conflicted_var, ~replace(.,str_detect(!!conflicted_var, "Unknown") & str_detect(!!based_var, "No"), "No"))) 
  }
  return(fixed_data)
}