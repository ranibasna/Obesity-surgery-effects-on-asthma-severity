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
    mutate(y1 = case_when(datediff %in% (0:364) ~ sum(prednisolone[datediff %in% (0:364)]), TRUE ~ 0)) %>% # check one year before
    mutate(y2 = case_when(datediff %in% (364:728) ~ sum(prednisolone[datediff %in% (364:728) ]), TRUE ~ 0)) %>% # check second year before
    mutate(y_3 = case_when(datediff %in% (728:1092) ~ sum(prednisolone[datediff %in% (728:1092)]))) %>% 
    mutate(y_1 = case_when(datediff %in% (-364:0) ~ sum(prednisolone[datediff %in% (-364:0) ]),  TRUE ~ 0)) %>% # check one year after the surgery
    mutate(y_2 = case_when(datediff %in% (-364:-728) ~ sum(prednisolone[datediff %in% (-364:-728)]),  TRUE ~ 0))  


  cols_to_check = c("y1","y2", "y_1", "y_2")
  OSC_df_years_df <- OSC_df_years %>% select(id,date_surgery, y1, y2, y_1, y_2, datediff) %>% 
                     mutate(across(cols_to_check, ~replace(., is.na(.), 0) )) %>% select(- datediff) %>% distinct()
  OSC_df_years_df <- aggregate(. ~ id + date_surgery, data = OSC_df_years_df, FUN=sum) %>% mutate(across(starts_with("y"), ~ if_else( . > 1825, 1,0), .names = "maint_{.col}"))
  # we take a random patient and check that the calculation is vald
  sample_rand_id <-  OSC_df_years_df %>% filter(maint_y1 == 1) %>% sample_n(1) %>% pull(id) 
  if(aggregate(. ~ id , raw_data %>% filter(id == sample_rand_id) %>% select(id, date_dispensation, prednisolone, date_surgery) %>% filter(( as.numeric(date_surgery) - as.numeric(date_dispensation)) %in% (1:365)), FUN=sum)$prednisolone > 1825){
    return(list(OSC_df_years_df, OSC_df_years))  
  } 
  else{
    stop( paste0("there are some proble with the following id", sample_rand_id))
  }
}

get_burst_var <- function(raw_data, OSC_df_years_data){
  
  df_burst <- OSC_df_years_data %>% select(- c( y_3)) %>% dplyr::distinct() %>% group_by(id) %>% arrange(date_dispensation, .by_group = TRUE) %>% 
    mutate(diffDate = difftime(date_dispensation, lag(date_dispensation,1), unit = 'day')) %>% 
    mutate( duration = lead(diffDate, order_by = id)) %>% select(- diffDate) %>%
    mutate(burstPck_y1 = case_when(round(y1) %in% (1:1825) ~ 1, TRUE ~ 0)) %>% 
    mutate(burstPck_y2 = case_when(round(y2) %in% (1:1825) ~ 1, TRUE ~ 0)) %>% 
    mutate(burstPck_y_1 = case_when(round(y_1) %in% (1:1825) ~ 1, TRUE ~ 0)) %>%
    mutate(burstPck_y_2 = case_when(round(y_2) %in% (1:1825) ~ 1, TRUE ~ 0)) %>% 
    mutate(across(starts_with("burstPck_y1"), ~ case_when(y1 %in% (1:1825)  & duration < 14 ~ . -1, TRUE ~ .))) %>% # this is instructed by Hannu so that we do not count the burst twice if the differene between two doses is less then 14 days
    mutate(across(starts_with("burstPck_y2"), ~ case_when(y2 %in% (1:1825)  & duration < 14 ~ . -1, TRUE ~ .))) %>% 
    mutate(across(starts_with("burstPck_y_1"), ~ case_when(y_1 %in% (1:1825)  & duration < 14 ~ . -1, TRUE ~ .))) %>% 
    mutate(across(starts_with("burstPck_y_2"), ~ case_when(y_2 %in% (1:1825)  & duration < 14 ~ . -1, TRUE ~ .))) 
  
  df_burst_pre <- df_burst  %>% select(-c(datediff, prednisolone, duration, date_surgery, packages, date_dispensation)) %>% select(- starts_with("y"))
  df_burst_agg <- aggregate(. ~ id , data = df_burst_pre, FUN=sum) 
  # df_burst_adjusted <- df_burst %>% select(- c(packages, date_dispensation, datediff, y1, y2, y_1, y_2, prednisolone, duration)) %>% 
  # mutate(across(starts_with("bur"), ~replace(., is.na(.), 0) ))
  # we take a random patient and check that the calculation is vald
  sample_rand_id <- df_burst_agg %>% filter(burstPck_y1 > 1) %>% sample_n(1) %>% pull(id)
  validating_filtered_df <- raw_data %>% filter(id == sample_rand_id) %>% select(id, date_dispensation, prednisolone, date_surgery) %>% filter(( as.numeric(date_surgery) - as.numeric(date_dispensation)) %in% (1:365))
  if(round(aggregate(. ~ id, validating_filtered_df, FUN=sum)$prednisolone) %in% (1:1825)) {
    return(df_burst_agg)  
  }else{
    stop( paste0("there are some proble with the following id ", sample_rand_id))
  }
  return(df_burst_agg)  
}
  ## Emergency outcome var ----
 
get_Em_vars <- function(raw_data){ 
  em_out_var <- raw_data %>% select(id, date_surgery, date_em_asthma) %>% 
    mutate(datediff = date_surgery - date_em_asthma) %>% # here we will get NA values when there are NA in the date_asthma_main
    mutate(has_em = case_when(is.na(date_em_asthma) ~ 0, TRUE ~ 1)) %>%  # if we have NA at. the date_asthma_main make it 0 in has_em
    # grouping
    group_by(id) %>% 
    arrange(date_em_asthma, .by_group = TRUE)  %>%
    mutate(em_y1 = case_when(datediff %in% (0:364) ~ sum(has_em[datediff %in% (0:364)]), TRUE ~ 0)) %>% # check one year before
    mutate(em_y2 = case_when(datediff %in% (364:728) ~ sum(has_em[datediff %in% (364:728) ]), TRUE ~ 0)) %>% # check second year before
    mutate(em_y_1 = case_when(datediff %in% (-364:0) ~ sum(has_em[datediff %in% (-364:0) ]),  TRUE ~ 0)) %>% # check one year after the surgery
    mutate(em_y_2 = case_when(datediff %in% (-364:-728) ~ sum(has_em[datediff %in% (-364:-728)]),  TRUE ~ 0))
  
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
    mutate(hos_y1 = case_when(datediff %in% (0:364) ~ sum(has_hos[datediff %in% (0:364)]), TRUE ~ 0)) %>% # check one year before
    mutate(hos_y2 = case_when(datediff %in% (364:728) ~ sum(has_hos[datediff %in% (364:728) ]), TRUE ~ 0)) %>% # check second year before
    mutate(hos_y_1 = case_when(datediff %in% (-364:0) ~ sum(has_hos[datediff %in% (-364:0) ]),  TRUE ~ 0)) %>% # check one year after the surgery
    mutate(hos_y_2 = case_when(datediff %in% (-364:-728) ~ sum(has_hos[datediff %in% (-364:-728)]),  TRUE ~ 0))
  
  hos_out_vars_agg <- hos_out_vars %>% select(- c(datediff, date_asthma_main, has_hos)) %>% dplyr::distinct()
  hos_out_vars_agg <- aggregate(. ~ id + date_surgery, data = hos_out_vars_agg %>% select(- inclusion), FUN=sum) 
  return(hos_out_vars_agg)

}

# plot patients ----
patient_follow_vis <- function(first_cut, second_cut, status_var){
  
  # status_var <- enquo(status_var)
  status_var <- sym(status_var)
  
  timeVis_df[first_cut:second_cut,] %>%
    
    # ploting
    ggplot(aes(x = newID, y = observedMonths, fill = dose_cat)) +
    geom_segment(aes(x = newID, xend = newID, y = 0, yend = lastMonth - 0.5), lty = 3) +
    geom_tile(width = 0.5) +
    scale_fill_manual("prednisolone level", values = c("#7D3A2C", "#AA3B2F", "#D36446")) +
    geom_linerange(mapping=aes(x=newID, ymin=monthSurg+0.3, ymax=monthSurg-0.3), width=0.2, size=1, color="blue") +
    geom_point(aes(x = newID, y = lastMonth - 0.5, shape = event)) +
    scale_shape_manual("Event", values = c(15, 5), labels = c("Death", "Alive", "")) +
    # geom_point(aes(x = newID, y = -3, color = UQ(status_var))) +
    geom_point(aes(x = newID, y = -3, color = !!status_var)) +
    scale_color_brewer(palette="BuPu") +
    # scale_color_viridis(discrete=TRUE, option="viridis") +
    # scale_color_manual(UQ(status_var), values = c("steelblue", "lightblue3", "gray30")) +
    guides(fill = guide_legend(override.aes = list(shape = NA), order = 1)) +
    coord_flip() +
    labs(y = "Month", x = "Patient\nNo.") +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(angle = 0))
}

# Define a second function

patient_follow_vis_2 <- function(first_cut, second_cut){
  # status_var <- sym(status_var)
  
  timeVis_df_NewId <- timeVis_df  %>% ungroup() %>% select(-id)
  
  timeVis_df_NewId[first_cut:second_cut,] %>%
    mutate_at(vars(newID, dose_cat), factor) %>% group_by(newID) %>% 
    mutate(month_death = case_when(event == "died" ~ lastMonth, TRUE ~ NA_real_)) %>% ungroup() %>%
    mutate(newID = fct_rev(fct_reorder(newID, lastMonth))) %>%
    ggplot(aes(x = observedMonths, y = newID, group = newID, col = dose_cat)) +
    geom_line() +
    geom_point(shape = 15) +
    geom_point(aes(x = month_death, y = newID), col = "black", shape = 4) +
    theme_bw() + 
    labs(x = "Months of Follow-Up", y = "Patient ID", col = "prednisolone level", title = "Patient Treatment Timeline", subtitle = "x indicates month of patient death") +
    # edit legend box and make patient ids small
    theme(axis.text.y = element_text(size=6), legend.position = c(.6, .9), legend.direction = "horizontal", legend.background = element_rect(linetype="solid", colour ="black")) + 
    # remove extra space around timeline
    scale_x_continuous(expand=c(0.01,0.01)) + # set the color of the lines and points
    scale_color_manual(values=c("dodgerblue","firebrick1","deepskyblue4"))
}

# BN analysis ----
# dis_data <- discretize(data = data_modeling_dis , method = "hartemink", breaks = 3, ordered = FALSE, ibreaks=60, idisc="quantile")
## Preparing the data ----

binner_hos <- function(x){
  # x <- cut(x, breaks =  c(-24, -14, -4, 0, 5, 11, 18), include.lowest = TRUE)
  x <- cut(x, breaks =  c(-24, -5, -0.5, 0, 1, 9, 18), include.lowest = TRUE)
  # now return the group number
  as.numeric(x)
}

binner_em <- function(x){
  x <- cut(x, breaks =  c(-2,-1, 0, 1, 2), include.lowest = TRUE)
  # now return the group number
  as.numeric(x)
}

inc_hos <- c("verybad", "Bad", "light",  "med", "Better", "Good")
inc_em <- c("Bad",  "med", "Better", "Good")

get_BN_processed_data <- function(BN_ready_data, Is_change = TRUE, Is_accChange = FALSE, Is_acc_main = FALSE, Not_burst = TRUE, Not_inclusion = TRUE, discretize_outcome_vars = FALSE, discretize_all_vars_xgb = TRUE, Remove_weight_vars=TRUE, CreatNumericalLevels=TRUE){
  mydata_ready_2 <- BN_ready_data %>% mutate(across(where(is.character), factor)) %>%  mutate(across(starts_with(c("bur","hos","em")), as.double)) %>% 
    mutate(across(starts_with(c("charl","R0","maint")), factor)) %>% mutate(across(starts_with(c("wais","age")), as.double))
  
  # Prepare the vars
  aux_vars <- c("bleeding_amount", "sleep_apnea", "hypertension", "dyslipidemia", "dyspepsia", "diabetes", "waist", "fp_glucose", "b_hba1c") 
  if(Is_change){
    # y1 is before. and y_1 is after
    mydata_ready_processed <- mydata_ready_2 %>% mutate(Diffy = (abs(y1)- abs(y_1)), hos_Diff = hos_y1 - hos_y_1, emDiff = em_y1 - em_y_1,  
                                                 burstDiff = burstPck_y1 - burstPck_y_1) 
    Y_var <- colnames(mydata_ready_processed %>% select(starts_with("Diffy")))
    HosVar <- colnames(mydata_ready_processed %>% select(starts_with("hos_D")))
    EmVar <- colnames(mydata_ready_processed %>% select(starts_with("em_D")))
    binner_hos <- function(x){
      x <- cut(x, breaks =  c(-24, -5, -0.5, 0, 1, 9, 18), include.lowest = TRUE)
      # now return the group number
      as.numeric(x)
    }
    
    binner_em <- function(x){
      x <- cut(x, breaks =  c(-2,-1, 0, 1, 2), include.lowest = TRUE)
      # now return the group number
      as.numeric(x)
    }
  }
  if (Is_accChange){
    # adding the accumulation for -2 +2 years
    mydata_ready_processed <- mydata_ready_2 %>% mutate(y_acc = ((y1+y2)-(y_1+y_2)), hos_acc = ((hos_y1+hos_y2) - (hos_y_1+hos_y_2)), 
                                                 em_acc = ((em_y1+em_y2) - (em_y_1+em_y_2)),  
                                                 burst_acc = ((burstPck_y1+ burstPck_y2) - (burstPck_y_1+burstPck_y_2)) )  
    Y_var <- colnames(mydata_ready_processed %>% select(starts_with("y_acc")))
    HosVar <- colnames(mydata_ready_processed %>% select(starts_with("hos_a")))
    EmVar <- colnames(mydata_ready_processed %>% select(starts_with("em_acc")))
  }
  if(Is_acc_main){
    # adding the accumulation for -2 +2 years
    mydata_ready_processed <- mydata_ready_2 %>% mutate(y_acc_b = y1+y2, y_acc_a = y_1+y_2, hos_acc_b = hos_y1+hos_y2, hos_acc_a = hos_y_1+hos_y_2, 
                                              em_acc_b = em_y1+em_y2, em_acc_a = em_y_1+em_y_2,  
                                              burst_acc_b = burstPck_y1+ burstPck_y2, burst_acc_a =  burstPck_y_1+burstPck_y_2)
    
    Y_var <- colnames(mydata_ready_processed %>% select(starts_with("y_acc")))
    HosVar <- colnames(mydata_ready_processed %>% select(starts_with("hos_a")))
    EmVar <- colnames(mydata_ready_processed %>% select(starts_with("em_acc")))
    #
    binner_hos <- function(x){
      x <- cut(x, breaks =  c(0, 1, 4, 8, 11, 19, 36), include.lowest = TRUE)
      as.numeric(x)
    }
    binner_em <- function(x){
      x <- cut(x, breaks =  c(0, 1, 2, 3, 4), include.lowest = TRUE)
      as.numeric(x)
    }
    inc_hos <- rev(inc_hos)
    inc_em <- rev(inc_em)
  }
  
  # remove the original osc medication vars
  mydata_ready_processed <- mydata_ready_processed %>% select(-c(y1, y_1, y2, y_2, hos_y1, hos_y_1, hos_y2, hos_y_2, em_y1, em_y_1, em_y2, em_y_2, burstPck_y1, burstPck_y_1, burstPck_y2, burstPck_y_2))
  
  # set a formula
  f <- as.formula(paste(HosVar[1], "~ ."))
  # set the recipe
  data_modeling_tran_imp_change_re <- 
    recipes::recipe( formula = f, data = mydata_ready_processed) %>% # replacing Diffy with hos_Diff
    update_role(id, date_surgery, found, missing_death_info, new_role = "ID") %>% 
    update_role(bleeding_amount, sleep_apnea, hypertension, dyslipidemia, dyspepsia, diabetes, waist, fp_glucose, b_hba1c, new_role = "auxiliary_vars") %>% 
    # add_role(starts_with(c("Diff", "hos")), new_role = "numeric outcome feature") %>% 
    # add_role(starts_with(c("Diff", "hos")), new_role = "outcome feature") %>% 
    step_zv(all_predictors()) %>% # remove variables that contain only a single value.
    # step_YeoJohnson(has_role(match = "numeric outcome feature")) %>% 
    step_YeoJohnson(all_of(Y_var)) %>% 
    step_impute_bag( c(all_predictors(), has_role(match = "outcome")), impute_with = imp_vars(c(all_predictors(), has_role(match = "auxiliary_vars")))) %>% 
    step_num2factor( all_of(HosVar), transform = binner_hos, levels = inc_hos, ordered = TRUE ) %>% 
    step_num2factor(all_of(EmVar), transform = binner_em, levels = inc_em, ordered = TRUE ) %>% 
    step_mutate(BMI_6w = bmi - (100*wght_loss6w)*(bmi-25), BMI_1y = bmi - (100*wght_loss1y)*(bmi-25), role = "predictor")
  
  
  # if(Is_change){
  #   data_modeling_tran_imp_change_re <- step_num2factor(recipe = data_modeling_tran_imp_change_re, all_of(HosVar), transform = binner_hos, levels = inc_hos, ordered = TRUE ) %>% 
  #   step_num2factor(all_of(EmVar), transform = binner_em, levels = inc_em, ordered = TRUE ) 
  # }
  
  if(discretize_all_vars_xgb){
    data_modeling_tran_imp_dis_change_re <- step_normalize(recipe = data_modeling_tran_imp_change_re) %>% step_discretize_xgb(all_numeric(), outcome = paste(HosVar[1])) 
    data_modeling_tran_imp_change_juice <- data_modeling_tran_imp_dis_change_re %>% prep() %>% juice()  
  } else{
    data_modeling_tran_imp_change_juice <- data_modeling_tran_imp_change_re %>% prep() %>% juice()
  }
    
  
  # remove auxiliary vars
  data_modeling_tran_imp_change <- data_modeling_tran_imp_change_juice %>% select(- all_of(aux_vars)) %>% select(-c(id, date_surgery, found, missing_death_info))
  # var selection
  data_modeling_tran_imp_change <- data_modeling_tran_imp_change %>% select(-c(maint_y1, maint_y2, maint_y_1, maint_y_2, smoking5, smoking2)) %>% as.data.frame()
  if(Not_burst){
    data_modeling_tran_imp_change <- data_modeling_tran_imp_change %>% select(- starts_with("burs"))
  }
  if(Not_inclusion){
    data_modeling_tran_imp_change <- data_modeling_tran_imp_change %>% select(- inclusion)
  }
  if(Remove_weight_vars){
    data_modeling_tran_imp_change <- data_modeling_tran_imp_change %>% select(- c(wght_loss6w, wght_loss1y))
  }
  if(CreatNumericalLevels){
    oldLevels_b <- sort(unique(data_modeling_tran_imp_change$y_acc_b))
    newLevels_b <- order(sort(unique(data_modeling_tran_imp_change$y_acc_b)))
    oldLevels_a <- sort(unique(data_modeling_tran_imp_change$y_acc_a))
    newLevels_a <- order(sort(unique(data_modeling_tran_imp_change$y_acc_a)))
    
    data_modeling_tran_imp_change <- data_modeling_tran_imp_change %>% mutate(y_acc_b = forcats::fct_recode(y_acc_b, !!! setNames(as.character(oldLevels_b), newLevels_b))) %>% 
      mutate(y_acc_a = forcats::fct_recode(y_acc_a, !!! setNames(as.character(oldLevels_a), newLevels_a)))
  }
  if(discretize_outcome_vars){
    # data_modeling_tran_imp_change <- data_modeling_tran_imp_change %>% mutate(across(starts_with(c("hos","em")), factor))
    dis_data <- discretize(data = data_modeling_dis , method = "hartemink", breaks = 3, ordered = FALSE, ibreaks=60, idisc="quantile")
    data_modeling_tran_imp_change
  }
  return(data_modeling_tran_imp_change)
  # 
} 

get_time_df <- function(ready_data){
    my_clean_time_data <- ready_data %>% mutate(across(starts_with(c("study", "date", "death")) , dmy)) %>% 
     select(id, date_dispensation, date_surgery, prednisolone, deathdate, inclusion, sex, surgery, packages, age) %>%
    # Subtract dates of dispensation from the date oof surgery for each patient
    mutate(numberMonth = interval(date_surgery, date_dispensation) %/% months(1)) %>%
    group_by(id) %>%
    arrange(numberMonth, .by_group = TRUE) %>% ungroup() %>%
    mutate(event = if_else( is.na(deathdate), "died","alive")) %>% group_by(id) %>%
    mutate(observedMonths = numberMonth - first(numberMonth) , monthSurg = first(abs(numberMonth)), lastMonth = max(observedMonths)) %>%
    mutate(dose_cat = case_when(prednisolone < 125 ~  "low", prednisolone >= 125 & prednisolone <= 416 ~ "meduim", prednisolone > 416 ~ "High"))
    # adding a new ordered ID 
    my_clean_time_data$newID <- with(my_clean_time_data, as.numeric(factor(id, levels = unique(id))))
    my_clean_time_data$sex <- as.factor(my_clean_time_data$sex)
    return(my_clean_time_data)
}

get_ITS_data <- function(timeVis_data){
  intrupted_df <- timeVis_data
  intrupted_df$NumberPeople <- 1
  intrupted_df <- intrupted_df %>% ungroup() %>% select(-id) %>% distinct() %>% filter(numberMonth %in% (-48:48)) %>% select(numberMonth, prednisolone, packages, NumberPeople, age, sex) 
  intrupted_df <- intrupted_df %>% mutate(genderPerc = if_else(sex == "Female", 1, 0)) %>% select(-sex)
  intrupted_df_agg <- aggregate(.~ numberMonth, intrupted_df, FUN = sum) 
  intrupted_df_agg  <- intrupted_df_agg %>% mutate(sergury_state = case_when(numberMonth <= 0 ~ "Before surgery", numberMonth > 0 ~ "after surgery"))
  
  timeVis_df48_summary <- timeVis_data %>% filter(numberMonth %in% (-48:48)) %>% select(id, date_surgery, deathdate) %>% distinct() %>% 
    mutate(DeathMonth = interval(date_surgery, deathdate) %/% months(1)) %>%  
    mutate(DeathMonthCount = if_else(DeathMonth < 49, DeathMonth, 0 )) %>% mutate(DeathMonthCount = tidyr::replace_na(DeathMonthCount, 0)) %>%                                                                  group_by(DeathMonthCount) %>% summarise(DeathNumber = n())
  #
  # insert the values of the missing death monthes
  index_set <- setdiff(intrupted_df_agg$numberMonth[-(1:49)], timeVis_df48_summary$DeathMonthCount)
  newDf <- timeVis_df48_summary
  for (i in index_set){
    newrow <- c(i,0)
    newDf <- insertRow(newDf, newrow, i)
  }
  # we will replace the first number by zero as this reflect the number of people who died in that month after the surgery.
  newDf$DeathNumber[1] <- 0
  timeVis_df48_cum <- newDf %>% mutate(cumsum_col = cumsum(DeathNumber))
  deathMonthCount_vec <- rep(0, 48)
  deathMonthCount_vec <- append(deathMonthCount_vec, timeVis_df48_cum$cumsum_col)
  intrupted_df_agg$deathMonthAccCount <- deathMonthCount_vec
  # adding elapsed time
  intrupted_df_agg$elapsed_time <- intrupted_df_agg$numberMonth + 49
  # mean age
  intrupted_df_agg <- intrupted_df_agg %>% mutate(avg_Age = age/NumberPeople) %>% select(-age)
  # adding the elapsed time after surgery (this is similar to the ramp variable)
  ramp <- append(rep(0,49), seq(1,48,1))
  intrupted_df_agg$timeAfterSurg <- ramp 
  #
  intrupted_df_agg <- intrupted_df_agg %>% mutate(step_change = case_when(numberMonth <= 0 ~ 0, numberMonth > 0 ~ 1)) %>% mutate(step_change = as.integer(step_change))
  intrupted_df_agg <- intrupted_df_agg %>% mutate(predAvg = prednisolone/NumberPeople)
  return(intrupted_df_agg)
}


drop_vars_df <- function(data_drop, drop_var){
  assertthat::assert_that(is.vector(drop_var))
  data_drop <- data_drop[ , !(names(data_drop) %in% drop_var)]
  return(data_drop)
}

## Blacklist and whightlist ----

get_bl <- function(final_BN_data, bmi_var){
  # age 
  bl_f_col_age <- colnames(final_BN_data %>% select(-age))
  bl_age = data.frame("from" = bl_f_col_age, "to" = "age")
  # gender
  bl_f_col_gender <- colnames(final_BN_data %>% select(-sex))
  bl_gender = data.frame("from" = bl_f_col_gender, "to" = "sex")
  # black listing backwords in time
  bl_smok <- tiers2blacklist(list("smoking","smoking1"))
  bl_y <- tiers2blacklist(list("y_acc_b", "y_acc_a"))
  bl_hos <- tiers2blacklist(list("hos_acc_b", "hos_acc_a"))
  bl_em <- tiers2blacklist(list("em_acc_b", "em_acc_a"))
  if(bmi_var){
    bl_sur_bmi_wl <- tiers2blacklist(list("bmi","surgery","BMI_6w", "BMI_1y"))  
  }else{
    bl_sur_bmi_wl <- tiers2blacklist(list("bmi","surgery","wght_loss6w", "wght_loss1y"))
  }
  bl  <- rbind(bl_age, bl_gender, bl_sur_bmi_wl, bl_smok)
  return(bl)
}

# whitelist
wl = matrix(c("wght_loss6w", "hos_Diff", "surgery", "wght_loss6w"), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))





## Plot network ----

plot.network <- function(structure, ht = "400px"){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}

## Simulation for validation ----
sim_validation <- function(fitted_BN_obj, BN_ready_data, n_pop){
  # compare the synthetic and original data frames
  df <- BN_ready_data %>% 
    mutate(type = "orig") %>% 
    bind_rows(
      rbn(fitted_BN_obj, n_pop) %>% 
        mutate(type = "sim")
    ) # %>% 
  gg_list <- list()
  grp_var <- "type"
  
  vars <- colnames(df)[colnames(df) != grp_var]
  for(k in 1:length(vars)){
    var_k <- vars[k]
    gg_list[[k]] <- ggplot(df, aes_string(x = var_k, fill = grp_var, col = grp_var))
    if(is.numeric(df[[var_k]])){
      gg_list[[k]] <- gg_list[[k]] + geom_density(alpha = 0.85, size = 0)
    }else{
      gg_list[[k]] <- gg_list[[k]] + geom_bar(position = "dodge")
    }
    gg_list[[k]] <- gg_list[[k]] + 
      theme(
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank()
      ) +
      labs(title = var_k)
  }
  return(gg_list)
}

## conditional probability function ----
cpq_effe_modif <- function(.data, vars, outcome, state, model, repeats = 500000) {
  all.levels <- if (any(length(vars) > 1)) {
    lapply(.data[, (names(.data) %in% vars)], levels)
  } else {
    all.levels <- .data %>%
      select(all_of(vars)) %>%
      sapply(levels) %>%
      as_tibble()
  }
  combos <- do.call("expand.grid", c(all.levels, list(stringsAsFactors = FALSE))) # al combiations
  
  # generate character strings for all combinations
  str1 <- ""
  for (i in seq(nrow(combos))) {
    str1[i] <- paste(combos %>% names(), " == '",
                     combos[i, ] %>% sapply(as.character), "'",
                     sep = "", collapse = "& "
    )
  }
  
  # repeat the string for more than one outcome
  str1 <- rep(str1, times = length(outcome))
  # str1 <- paste("list(", str1, ")", sep = "")
  str1 <- paste("(" ,str1, ")", sep = "")
  # repeat loop for outcome variables (can have more than one outcome)
  all.levels.outcome <- if (any(length(outcome) > 1)) {
    lapply(.data[, (names(.data) %in% outcome)], levels)
  } else {
    all.levels <- .data %>%
      select(all_of(outcome)) %>%
      sapply(levels) %>%
      as_tibble()
  }
  combos.outcome <- do.call("expand.grid", c(all.levels.outcome))
  
  # repeat each outcome for the length of combos
  str3 <- rep(paste("(", outcome, " == '", state, "')", sep = ""), each = length(str1) / length(outcome))
  
  # fit the model
  # fitted <- bn.fit(cextend(model), .data)
  fitted <- model
  
  # join all elements of string together
  cmd <- paste("replicate(200,cpquery(fitted, ", str3, ", ", str1, ", method = 'ls', n = ", repeats, "))", sep = "")
  # cmd <- paste("replicate(200,cpquery(fitted, ", str3, ", ", str1, ", method='lw' n = ", repeats, "))", sep = "")
  prob <- rep(0, length(str1)) # empty vector for probabilities
  q05 <- rep(0,length(str1))
  q975 <- rep(0,length(str1))
  for (i in seq(length(cmd))) {
    prop_vec = eval(parse(text =  cmd[i]))
    q05[i] = quantile(prop_vec, 0.05, na.rm = TRUE)
    q975[i] = quantile(prop_vec,0.975, na.rm = TRUE)
    prob[i] <- mean(prop_vec, na.rm = TRUE)
  } # for each combination of strings, what is the probability of outcome
  res <- cbind(combos, prob, q05, q975)
  
  return(res)
}
# get conditional probability plot for one variable
get_cpq_plot_one_var <- function(res_data, effe_modif_vars, original_raw_data, final_data){
  # get the var cases
  var1Cases <- unique(final_data[,effe_modif_vars[1]])
  #var2Cases <- unique(final_data[,effe_modif_vars[2]])
  # get the var labels
  # Var1Labels <- names(attributes(unique(original_raw_data[,effe_modif_vars[1]]))$labels)
  Var1Labels <- unique(original_raw_data[,effe_modif_vars[1]])
  #
  var_1_sym <- sym(effe_modif_vars[1])
  # define conditions
  conditions_1 <- purrr::map2(var1Cases, Var1Labels, ~quo( !!var_1_sym == !!.x ~ !!.y))
  
  # mutate the new coloumns
  res_data <- res_data %>% mutate("{effe_modif_vars[1]}_Cases" := case_when(!!!conditions_1))
  res_colnames <- colnames(res_data)
  var1 <- sym(res_colnames[grepl("_Cases",res_colnames)])
  #prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var1), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var1), position = position_dodge(0.3)) + scale_color_brewer() + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  #prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var1), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var1), position = position_dodge(0.3)) + scale_color_viridis(discrete = TRUE, option = "D") + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var1), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var1), position = position_dodge(0.3)) + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  res_list <- list()
  res_list[[1]] <- res_data
  res_list[[2]] <- prop_p
  return(res_list)
}
# get the result table and the conditional probability plots
get_cpq_plot <- function(res_data, effe_modif_vars, original_raw_data, final_data){
  # get the var cases
  var1Cases <- unique(final_data[,effe_modif_vars[1]])
  var2Cases <- unique(final_data[,effe_modif_vars[2]])
  # get the var labels
  # Var1Labels <- names(attributes(unique(original_raw_data[,effe_modif_vars[1]]))$labels)
  Var1Labels <- unique(original_raw_data[,effe_modif_vars[1]])
  # Var2Labels <- names(attributes(unique(original_raw_data[,effe_modif_vars[2]]))$labels)
  Var2Labels <- unique(original_raw_data[,effe_modif_vars[2]])
  #
  var_1_sym <- sym(effe_modif_vars[1])
  var_2_sym <- sym(effe_modif_vars[2])
  # define conditions
  conditions_1 <- purrr::map2(var1Cases, Var1Labels, ~quo( !!var_1_sym == !!.x ~ !!.y))
  conditions_2 <- purrr::map2(var2Cases, Var2Labels, ~quo( !!var_2_sym == !!.x ~ !!.y))
  # mutate the new coloumns
  res_data <- res_data %>% mutate("{effe_modif_vars[1]}cases" := case_when(!!!conditions_1)) %>% mutate("{effe_modif_vars[2]}_Cases" := case_when(!!!conditions_2))
  res_colnames <- colnames(res_data)
  var1 <- sym(res_colnames[grepl("cases",res_colnames)])
  var2 <- sym(res_colnames[grepl("_Cases",res_colnames)])
  # browser()
  # prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var2), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var2), position = position_dodge(0.3)) + scale_color_manual(values = c("#00AFBB", "#E7B800",'#999999')) + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var2), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var2), position = position_dodge(0.3)) + scale_color_manual(values = c("#CC0000", "#006600", "#669999", "#00CCCC", "#660099", "#CC0066", "#FF9999", "#FF9900", "black", "black", "black", "black", "black")) + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  res_list <- list()
  res_list[[1]] <- res_data
  res_list[[2]] <- prop_p
  return(res_list)
}

# Processing for Intrupted time sereis ----

insertRow <- function(data, new_row, r) {
  data_new <- rbind(data[1:r, ],            
                    new_row,                
                    data[- (1:r), ])
  data_new <- as.data.frame(data_new)  
  rownames(data_new) <- 1:nrow(data_new)
  data_new <- as_tibble(data_new)
  return(data_new)
}
