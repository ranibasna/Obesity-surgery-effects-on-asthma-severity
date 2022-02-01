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

# processing the drug vars

get_oh_drug <- function(drug_data){
  drug_data <- drug_data %>% mutate(across(starts_with("drug"), ~replace(., stri_isempty(.), "None")))
  mydata_drug_oh <- drug_data %>% select(id, starts_with("drug")) %>% mutate(across(where(is.character), factor)) %>% 
                    as.data.table %>%  one_hot %>% select(-c(ends_with("None"))) %>% distinct() %>% 
                    rowwise() %>% mutate(R0AK = sum(c_across(ends_with("AK")))) %>% 
                    rowwise() %>% mutate(R0AL = sum(c_across(ends_with("AL")))) %>% 
                    select(- c(drug1_R03AK, drug1_R03AL, drug2_R03AK, drug2_R03AL, drug3_R03AL)) %>% rename(R03BA = drug1_R03BA)
  return(mydata_drug_oh)
}

## processing Smoking vars ----

get_smok_data <- function(smok_data){
  
smok_data <- smok_data %>% select(starts_with(c("smok","id")))
# we need first to replace the Unkown to Unknown
smok_data <- smok_data %>% mutate(across(starts_with("smok"), ~replace(., str_detect(.,"Unkown"), "Unknown")))
# Then we need to replace the empty string to Unknown
smok_data <- smok_data %>% mutate(across(starts_with("smok"), ~replace(., stri_isempty(.), "Unknown")))
# we start with replacing the empty string with NA for all the smoking cols 
# smok_data <- smok_data %>% mutate(across(starts_with("smok"), ~ na_if(.,"")))

# Now we have the following values in the smoking variables: Yes, No, Stopped, Stopped before the surgery (this is only for the smoking var), Unknown and NA
# Here we want to attach the two patterns mentioned. Namely, if Yes or Stopped comes first  then we should not accept No afterwords.
# either we replace with NA or with Stopped
smok_data_fixed <- smok_data %>% 
  fix_conflict_multi(conflicted_vars = c("smoking1", "smoking2", "smoking5"), based_var = "smoking", pattern_type = c("Yes", "Stopped","Temporary stop before surgery"), replace_value = NA) %>%
  fix_conflict_multi(conflicted_vars = c("smoking2","smoking5"), based_var = "smoking1", pattern_type = c("Yes", "Stopped"), replace_value = NA)  %>% 
  fix_conflict_multi(conflicted_vars = c("smoking5"), based_var = "smoking2", pattern_type = c("Yes", "Stopped"), replace_value = NA)  %>% 
  fix_conflict_unknown(conflicted_vars = c("smoking", "smoking1", "smoking2"), based_var = "smoking5" ) %>% 
  fix_conflict_unknown(conflicted_vars = c("smoking", "smoking1"), based_var = "smoking2") %>% 
  fix_conflict_unknown(conflicted_vars = c("smoking"), based_var = "smoking1")

# replace Unknown with NA
smok_data_fixed <- smok_data_fixed %>% mutate(across(starts_with("smok"), ~ na_if(.,"Unknown"))) %>% distinct()

return(smok_data_fixed)
}

# Processing the outcome vars ----

get_OSC <- function(raw_data){
  
  OSC_df <- raw_data %>%  select(id, date_dispensation, packages, prednisolone, date_asthma_main, date_em_asthma, date_surgery, date_first_other)
  
  OSC_df_years <- OSC_df %>% select(-c(date_asthma_main, date_em_asthma, date_first_other)) %>% 
  # mutate(across(c("date_dispensation", "date_surgery"), ~ as.Date(.x,"%d%b%Y"))) %>% 
  # Subtract dates of dispensation from the date oof surgery for each patient 
  mutate(datediff = as.Date(date_surgery) - as.Date(date_dispensation)) %>% 
  # grouping
  group_by(id) %>% 
  arrange(date_dispensation, .by_group = TRUE)  %>% 
  mutate(y1 = case_when(datediff > 0 & datediff <= 364  ~ sum(prednisolone[datediff > 0 & datediff <= 364 ]))) %>% 
  mutate(y2 = case_when(datediff > 364 & datediff <= 728  ~ sum(prednisolone[datediff > 364 & datediff <= 728])))  %>% 
  mutate(y_3 = case_when(datediff > 728 & datediff <= 1092  ~ sum(prednisolone[datediff > 728 & datediff <= 1092 ]))) %>% 
  mutate(y_1 = case_when(datediff < 0 & datediff >= -364  ~ sum(prednisolone[datediff < 0 & datediff >= -364 ]))) %>% 
  mutate(y_2 = case_when(datediff < -364 & datediff >= -728  ~ sum(prednisolone[datediff < -364 & datediff >= -728])))

  cols_to_check = c("y1","y2", "y_1", "y_2")
  OSC_df_years_df <- OSC_df_years %>% select(id,date_surgery, y1, y2, y_1, y_2, datediff) %>% 
                     mutate(across(cols_to_check, ~replace(., is.na(.), 0) )) %>% select(- datediff)
  OSC_df_years_df <- aggregate(. ~ id + date_surgery, data = OSC_df_years_df, FUN=sum) %>% mutate(across(starts_with("y"), ~ if_else( . > 1825, 1,0), .names = "maint_{.col}"))
  
  return(list(OSC_df_years_df, OSC_df_years))
}

get_burst_var <- function(raw_data, OSC_df_years_data){
  
  df_burst <- OSC_df_years_data %>% select(- c( y_3)) %>% dplyr::distinct()  %>% group_by(id) %>% arrange(date_dispensation, .by_group = TRUE)  %>% 
    mutate(diffDate = difftime(date_dispensation, lag(date_dispensation,1), unit = 'day')) %>% 
    mutate( duration = lead(diffDate, order_by = id)) %>% select(- diffDate) %>%
    mutate(across(starts_with("y"), ~ if_else( . < 1825, 1,0), .names = "burstPck_{.col}")) %>% 
    mutate(across(starts_with("burst"), ~ case_when(duration < 14 ~ . -1, TRUE ~ .)))
  
  
  df_burst_adjusted <- df_burst %>% select(- c(packages, date_dispensation, datediff, y1, y2, y_1, y_2, prednisolone, duration)) %>% 
    mutate(across(starts_with("bur"), ~replace(., is.na(.), 0) ))
  
  df_burst_agg <- aggregate(. ~ id + date_surgery, data = df_burst_adjusted, FUN=sum) 
  
  return(df_burst_agg)
}
  ## Emergency outcome var ----
 
get_Em_vars <- function(raw_data){ 
  em_out_var <- raw_data %>% select(id, date_surgery, date_em_asthma) %>% 
    mutate(datediff = date_surgery - date_em_asthma) %>% # here we will get NA values when there are NA in the date_asthma_main
    mutate(has_em = case_when(is.na(date_em_asthma) ~ 0, TRUE ~ 1)) %>%  # if we have NA at. the date_asthma_main make it 0 in has_hos
    # grouping
    group_by(id) %>% 
    arrange(date_em_asthma, .by_group = TRUE)  %>% 
    mutate(em_y1 = case_when(datediff > 0 & datediff <= 364  ~ sum(has_em[datediff > 0 & datediff <= 364 ]), TRUE ~ 0)) %>% # check one year before
    mutate(em_y2 = case_when(datediff > 364 & datediff <= 728  ~ sum(has_em[datediff > 364 & datediff <= 728 ]), TRUE ~ 0)) %>% # check second year before
    mutate(em_y_1 = case_when(datediff < 0 & datediff >= -364  ~ sum(has_em[datediff < 0 & datediff >= -364 ]),  TRUE ~ 0)) %>% # check one year after the surgery
    mutate(em_y_2 = case_when(datediff < -364 & datediff >= -728  ~ sum(has_em[datediff < -364 & datediff >= -728]),  TRUE ~ 0))
  
  
    em_out_var_agg <- em_out_var %>% select(- c(datediff, date_em_asthma, has_em)) %>% dplyr::distinct()
    em_out_var_agg <- aggregate(. ~ id + date_surgery, data = em_out_var_agg,  FUN=sum) 
    return(em_out_var_agg)
}
  ## Hospitalization outcome var ----
get_Hos_vars <- function(raw_data){ 
  hos_out_vars <- raw_data %>% select(id, date_surgery, date_asthma_main, inclusion) %>% 
    mutate(datediff = date_surgery - date_asthma_main) %>% # here we will get NA values when there are NA in the date_asthma_main
    mutate(has_hos = case_when(is.na(date_asthma_main) ~ 0, TRUE ~ 1)) %>%  # if we have NA at. the date_asthma_main make it 0 in has_hos
    # grouping
    group_by(id) %>% 
    arrange(date_asthma_main, .by_group = TRUE)  %>% 
    mutate(hos_y1 = case_when(datediff > 0 & datediff <= 364  ~ sum(has_hos[datediff > 0 & datediff <= 364 ]), TRUE ~ 0)) %>% # check one year before
    mutate(hos_y2 = case_when(datediff > 364 & datediff <= 728  ~ sum(has_hos[datediff > 364 & datediff <= 728 ]), TRUE ~ 0)) %>% # check second year before
    mutate(hos_y_1 = case_when(datediff < 0 & datediff >= -364  ~ sum(has_hos[datediff < 0 & datediff >= -364 ]),  TRUE ~ 0)) %>% # check one year after the surgery
    mutate(hos_y_2 = case_when(datediff < -364 & datediff >= -728  ~ sum(has_hos[datediff < -364 & datediff >= -728]),  TRUE ~ 0))
  
  
  hos_out_vars_agg <- hos_out_vars %>% select(- c(datediff, date_asthma_main, has_hos)) %>% dplyr::distinct()
  hos_out_vars_agg <- aggregate(. ~ id + date_surgery, data = hos_out_vars_agg %>% select(- inclusion), FUN=sum) 
  return(hos_out_vars_agg)

}
