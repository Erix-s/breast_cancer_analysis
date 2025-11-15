# Load Libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(caret)
library(pROC)
library(MASS)
library(tidyverse)
library(table1)
library(psych)
library(gridExtra)
library(caret)
library(performance)
library(MASS)
library(car)
library(officer)
library(flextable)
library(ggeffects)
library(gtsummary)
library(glmnet)
library(sjPlot)
library(pROC)



# Final Project - Problem 2
p2_data <- read.csv("data/Project_2_data.csv")
summary(p2_data)
colSums(is.na(p2_data))

# Categorical Variables to Factors
p2_data$Race <- as.factor(p2_data$Race)
p2_data$Marital.Status <- as.factor(p2_data$Marital.Status)
p2_data$T.Stage <- as.factor(p2_data$T.Stage)
p2_data$N.Stage <- as.factor(p2_data$N.Stage)
p2_data$X6th.Stage <- as.factor(p2_data$X6th.Stage)
p2_data$differentiate <- as.factor(p2_data$differentiate) # DM edit: fixed column capitalization, as.factor should run now
p2_data$Grade <- as.factor(p2_data$Grade)
p2_data$A.Stage <- as.factor(p2_data$A.Stage)
p2_data$Estrogen.Status <- as.factor(p2_data$Estrogen.Status)
p2_data$Progesterone.Status <- as.factor(p2_data$Progesterone.Status)
p2_data$Status <- as.factor(p2_data$Status) 

####### DM edit: adjust misspelled col name ####### 
p2_data <- 
  p2_data |> 
  rename(Regional.Node.Positive = Reginol.Node.Positive)
###################################################

summary(p2_data)

# Check for missing values
colSums(is.na(p2_data))

# Pairwise relationship (for numerical variables)
pairs(p2_data[, c("Age", "Tumor.Size", "Regional.Node.Examined", "Regional.Node.Positive")])

# Visualize survival time
ggplot(p2_data, aes(x = Survival.Months, fill = Status)) + 
  geom_histogram(binwidth = 5, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Alive" = "darkblue", "Dead" = "darkgreen"))
labs(title = "Distribution of Survival Months", x = "Months", y = "Frequency") +
  theme_minimal()

# Survival Analysis
# Kaplan-Meier survival curve 
km_fit <- survfit(Surv(Survival.Months, Status == "Dead") ~ 1, data = p2_data)
ggsurvplot(km_fit, data = p2_data, conf.int = TRUE, risk.table = TRUE)
summary(km_fit)

#Kaplan-Meier by Race
km_fit_race <- survfit(Surv(Survival.Months, Status == "Dead") ~ Race, data = p2_data)
ggsurvplot(km_fit_race, data = p2_data, conf.int = FALSE, risk.table = TRUE)
summary(km_fit_race)

# Cox Proportional Hazards Model
# Fit a Cox proportional hazards model
############# DM edit: adjust variable names if necessary #############
cox_model <- coxph(Surv(Survival.Months, Status == "Dead") ~ Age + Race + Marital.Status + 
                     T.Stage + N.Stage + X6th.Stage + differentiate + Grade +
                     A.Stage + Tumor.Size + Estrogen.Status + 
                     Progesterone.Status + Regional.Node.Examined + 
                     Regional.Node.Positive, data = p2_data)
summary(cox_model)
# Logistic Regression
data= read.csv("data/Project_2_data.csv", na.strings = c("NA", "", ".")) |>
  janitor::clean_names() |>
  mutate(
    survival_months = as.integer(as.numeric(survival_months)), 
    reginol_node_positive = as.numeric(reginol_node_positive),
    age = as.numeric(age),
    race = factor(race, levels = c("White", "Black", "Other")),
    grade = factor(grade, levels = c("1", "2", "3", "anaplastic; Grade IV")) |> relevel(ref = "1"),
    marital_status = factor(marital_status, levels = c("Married", "Single", "Divorced", "Separated", "Widowed")) |> relevel(ref = "Married"),
    differentiate = relevel(factor(differentiate), ref = "Well differentiated"),
    a_stage = relevel(factor(a_stage), ref = "Regional"),
    progesterone_status = relevel(factor(progesterone_status), ref = "Positive"), 
    estrogen_status = relevel(factor(estrogen_status), ref = "Positive"),
    t_stage = relevel(factor(t_stage), ref = "T1"),
    n_stage = relevel(factor(n_stage), ref = "N1"), 
    status_bin = ifelse(status=="Dead", 1, 0),
    status_bin = as.numeric(status_bin)
  )

#age is not linear 
m= glm(status_bin ~ age, data=data, family = binomial)

m2= glm(status_bin ~poly(age,2), data=data, family = binomial)
m3= glm(status_bin ~poly(age,3), data=data, family = binomial)
m4= glm(status_bin ~poly(age,4), data=data, family = binomial)

AIC(m,m2,m3,m4)
print(tab_model(m4))
plot_model(m2,type= "eff", terms= "age [all]")

#Summary Statistics Tables 
data= data|> 
  mutate(groupvar = ifelse(status == "Alive", "Alive", "Not Alive"))


descriptive_stat= table1::table1(~age+ race + marital_status+  survival_months  | groupvar, data = data, na.rm = TRUE, digits = 1, format.number = TRUE)

descriptive_stat|>
  as.data.frame()|>
  flextable()|>
  save_as_image( path = "Descriptive_stat_table.png")

descriptive_stat
 # for race
data_race= data|> 
  mutate(groupvar = race)

stats_by_race= table1::table1(~ t_stage+ n_stage + x6th_stage+ differentiate +grade+ a_stage + tumor_size + regional_node_examined + reginol_node_positive + estrogen_status+ progesterone_status + survival_months  | groupvar, data = data_race, na.rm = TRUE, digits = 1, format.number = TRUE)

stats_by_race|>
  as.data.frame()|>
  flextable()|>
  save_as_image( path = "race_stat_table.png")

stats_by_race
 #for clincial variables 
stat= table1::table1(~t_stage+ n_stage + x6th_stage+ differentiate +grade+ a_stage + tumor_size + regional_node_examined + reginol_node_positive + estrogen_status+ progesterone_status + survival_months  | groupvar, data = data, na.rm = TRUE, digits = 1, format.number = TRUE)

stat|>
  as.data.frame()|>
  flextable()|>
  save_as_image( path = "stat_table.png")

stat

#Stepwise variable selection 

full_model = glm(status_bin ~ age + race + marital_status + t_stage + x6th_stage + differentiate + grade + a_stage + tumor_size + estrogen_status + progesterone_status + regional_node_examined + reginol_node_positive, family = binomial, data= data)

stepAIC(full_model)
best_model = stepAIC(full_model)
check_model(best_model)

#Lasso method for variable selection
lambda_seq <- 10^seq(-3, 0, by = .1)
set.seed(2022)


cv_object= cv.glmnet(as.matrix(data[1:14]), data$status_bin, 
                     lambda = lambda_seq, 
                     nfolds = 5)
cv_object 
fit_bestcv= glmnet(as.matrix(data[1:14]), data$status_bin, lambda = cv_object$lambda.min)
coef(fit_bestcv)

# Comparing lasso and stepwise models 
lasso_model = glm(status_bin ~age+  grade +tumor_size + regional_node_examined + reginol_node_positive, family = binomial, data= data)

final_model= glm(formula = status_bin ~ age + race + marital_status + t_stage + 
                   n_stage + differentiate + estrogen_status + progesterone_status + 
                   regional_node_examined + +reginol_node_positive, family = binomial, 
                 data = data)

anova(lasso_model, final_model, test = "Chisq")

#final logistic regression model indcluding age as polynomial variables 
poly = glm(formula = status_bin ~ poly(age,2) + race + marital_status + t_stage  + differentiate + estrogen_status + progesterone_status + 
             regional_node_examined + +reginol_node_positive, family = binomial, 
           data = data)

#model diagnostics 


check_mode(poly)
#model goodness of fit 
performance(poly)
 #AUC
roc(status_bin ~ fitted.values(poly), data= data,
    plot = TRUE, legacy.axes= TRUE,print.auc=TRUE, ci=TRUE)

  #confusion matrix 
d1 = data |> 
  mutate(pred_probs = as.vector(predict(final_model, type = "response"))) |> 
  mutate(
    predicted_outcome = case_when(
      pred_probs > 0.5 ~ "yes",
      TRUE ~ "no"
    ),
    predicted_outcome = factor(predicted_outcome, levels = c("no", "yes")),
    status_bin = factor(status_bin, levels = c(0, 1))
  )|>
  mutate(
    status_bin= recode_factor(status_bin, `0` = "no", `1` = "yes")
  )
 #confusion matrix, specfificty & sensitivity 
cross_table = table(d1$status_bin, d1$predicted_outcome)
caret::confusionMatrix(cross_table)

table = 
  tbl_regression(
    poly, 
    exponentiate = TRUE,  
    add_pairwise_contrasts = TRUE, 
    contrasts_adjust = "bonferroni", 
    pairwise_reverse = FALSE, 
    pvalue_fun = ~style_pvalue(.x, digits = 2)
  ) |>
  add_significance_stars(hide_p = FALSE, hide_se = TRUE, hide_ci = FALSE) |>
  bold_p()
 #odds ratio table 
table|>
  as_flex_table()|>
  save_as_image(path= "OR_table_2.png")
table 

 # age probabilities 
m_age= glm(formula = status_bin ~ age + race + marital_status + t_stage + 
             n_stage + differentiate + estrogen_status + progesterone_status + 
             regional_node_examined + +reginol_node_positive, family = binomial, 
           data = data)

m2_age= glm(formula = status_bin ~ poly(age,2) + race + marital_status + t_stage + 
              n_stage + differentiate + estrogen_status + progesterone_status + 
              regional_node_examined + +reginol_node_positive, family = binomial, 
            data = data)
m3_age= glm(formula = status_bin ~ poly(age,2) + race + marital_status + t_stage + 
              n_stage + differentiate + estrogen_status + progesterone_status + 
              regional_node_examined + +reginol_node_positive, family = binomial, 
            data = data)

AIC(m_age,m2_age,m3_age)

#displays degrees 
library(sjPlot)
print(tab_model(m3_age))

plot_model(m2_age,type= "eff", terms= "age [all]")
b= emmip(m2_age, ~ age, CIs= T, type = "response", at = list(age= c(30,50,70)))+ 
  scale_y_continuous(labels= scales::percent)

b

emmeans(m2_age, pairwise~age, type = "response", at= list(age= c(30,50,70)), infer= T)


 #likelihood ratio test
car::Anova(poly)

############# DM edit: cox results #############
# stat sig variables: age, race (other), race (white), T.StageT4, N.StageN2, N.StageN3, 
# ... differentiate (pooly, undiff, well), Estrogen.StatusPositive, Progesterone.StatusPositive, 
# ... Regional.Node.Examined, Regional.Node.Positive
# positive hazard coefs (higher or true values means higher risk of death): age, T.StageT4, N.StageN2, N.StageN3, 
# ... differentiate (poor, undiff), Regional.Node.Positive
# all the others have negative hazard coefs

# Check proportional hazards assumptions
cox.zph(cox_model)

# View model violation: Visualize Schoenfeld Residuals to better understand the violation

# plot schoenfeld residuals for all predictors 
plot(cox.zph(cox_model))

# plot schoenfeld residuals for residuals, Estrogen Status 
plot(cox.zph(cox_model), var = "Estrogen.Status")

# plot schoenfeld residuals for Progesterone Status
plot(cox.zph(cox_model), var = "Progesterone.Status")

############# DM edit: adjust variable names if necessary #############
# Stratify the Cox Model by Estrogen and Progesterone Status
cox_model_strat <- coxph(Surv(Survival.Months, Status == "Dead") ~ Age + Race + Marital.Status + 
                           T.Stage + N.Stage + X6th.Stage + differentiate + Grade +
                           A.Stage + Tumor.Size + Regional.Node.Examined + 
                           Regional.Node.Positive + strata(Estrogen.Status, Progesterone.Status), data = p2_data)
summary(cox_model_strat)

# check hazards assumption for strat data
cox.zph(cox_model_strat)

# Strat data plots
plot(cox.zph(cox_model_strat))


# Evaluate Model Performance
# Split data into training and testing sets
set.seed(123)
train_index <- createDataPartition(p2_data$Status, p = 0.7, list = FALSE)
train_data <- p2_data[train_index, ]
test_data <- p2_data[-train_index, ]

# Predict risk scores on test data
risk_scores <- predict(cox_model, newdata = test_data, type = "risk")

# ROC curve
roc_obj <- roc(test_data$Status, risk_scores, levels = c("Alive", "Dead"))
plot(roc_obj, main = "ROC Curve")
auc(roc_obj)

# Fairness Evaluation by Race
# Stratify predictions by race
for(r in unique (p2_data$Race)) {
  race_data <- test_data[test_data$Race == r, ]
  race_risk_scores <- predict(cox_model, newdata = race_data, type = "risk")
  roc_race <- roc(race_data$Status, race_risk_scores, levels = c("Alive", "Dead"))
  print(paste("AUC for race:", r))
  print (auc(roc_race))
}

# Example approach: Recalibrate model using balanced sampling or adjusted loss funtions
balanced_train <- train_data %>%
  group_by(Race) %>%
  sample_n(size = min(table(train_data$Race)))

############# DM edit: adjust variable names if necessary #############
balanced_cox_model <- coxph(Surv(Survival.Months, Status == "Dead") ~ Age + Race + Marital.Status + 
                              T.Stage + N.Stage + X6th.Stage + differentiate + Grade +
                              A.Stage + Tumor.Size + Regional.Node.Examined +
                              Regional.Node.Positive, data = balanced_train)

# Evaluate fairness again
for(r in unique(p2_data$Race)) {
  race_data <- test_data[test_data$Race == r, ]
  race_risk_scores <- predict(balanced_cox_model, newdata = race_data, type = "risk")
  roc_race <- roc(race_data$Status, race_risk_scores, levels = c("Alive", "Dead"))
  print(paste("AUC for race:", r))
  print(auc(roc_race))
}

# --------- end

##### DM edit: check cox model results when only using significant variables determined with cox_model #####
cox_model_lim <- coxph(Surv(Survival.Months, Status == "Dead") ~ 
                             Age + Race + T.Stage + N.Stage + differentiate + Estrogen.Status + 
                             Progesterone.Status + Regional.Node.Examined + Regional.Node.Positive, 
                           data = p2_data)
summary(cox_model_lim)



