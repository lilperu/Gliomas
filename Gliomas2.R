if(!require(tidymodels)) install.packages("tidymodels", repos = "http://cran.us.r-project.org")
if(!require(rpart)) install.packages("rpart", repos = "http://cran.us.r-project.org")
if(!require(rsample)) install.packages("rsample", repos = "http://cran.us.r-project.org")
if(!require(tidyr)) install.packages("tidyr", repos = "http://cran.us.r-project.org")
if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
if(!require(magrittr)) install.packages("magrittr", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(baguetts)) install.packages("baguette", repos = "http://cran.us.r-project.org")
if(!require(ranger)) install.packages("ranger", repos = "http://cran.us.r-project.org")
if(!require(vip)) install.packages("vip", repos = "http://cran.us.r-project.org")
if(!require(xgboost)) install.packages("xgboost", repos = "http://cran.us.r-project.org")
if(!require(rpart.plot)) install.packages("rpart.plot", repos = "http://cran.us.r-project.org")
if(!require(knitr)) install.packages("knitr", repos = "http://cran.us.r-project.org")
library(tidymodels)
library(magrittr)
library(dplyr)
library(baguette)
library(vip)
library(rpart.plot)

# open and preview data
genes<- read.csv("TCGA_GBM_LGG_Mutations_all.csv")
head(genes)

# convert non-numeric attributes to factors
genes<- genes %>% mutate(across(.cols= where(is.character), .fns= as.factor))
genes<- genes %>% mutate(across(.cols= where(is.integer), .fns= as.factor))
str(genes)

# remove four patients with empty Gender, Diagnosis, Race
genes<- genes %>% filter(Gender != "--")

# verify zero NAs in each column
colSums(is.na(genes))
str(genes)

# split data in to train and test sets
set.seed(1998)
#strata ensures even distribution of Grade
genes_split<- initial_split(genes, strata=Grade)
genes_train<- training(genes_split)
genes_test<- testing(genes_split)
nrow(genes_train)/nrow(genes)


#############################################################################
# EXPLORATORY DATA ANALYSIS

# table of factor counts
genes_count<- genes[,-2:-7] # exclude demographics
result<- apply(genes_count, 2, function(x) table(x)) # factor counts for each column
result

# convert to data frame
result<- as.data.frame(result)
# convert from wide
result_long <- result %>%
  rownames_to_column(var = "Factor") %>%
  pivot_longer(cols = 3:22, names_to = "Gene", values_to = "Ratio")
# convert Gene column to factor and specify levels
result_long$Gene <- factor(result_long$Gene, levels = colnames(result))
# retain original order
result_long <- result_long[order(result_long$Gene), ]

# Create a bar plot of the mutation frequencies
ggplot(result_long, aes(x = Gene, y = Ratio, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Ratios of Mutations", x = "Gene", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_discrete(labels=c('No Mutation','Mutation'))

# set specifications and fit to create model
model_base<- decision_tree() %>%
  set_mode("classification") %>%
  set_engine("rpart") %>% 
  fit(Grade~.-Project-Case_ID-Primary_Diagnosis, data= genes_train)

# predict on test set
preds_base<- predict(model_base, genes_test, type="prob")
preds_base

# add predictions to test set
preds_base<- data.frame(Grade= genes_test$Grade, preds_base)
# label with new model column 
preds_base<- preds_base %>% mutate(model= "Tree")
preds_base

# visualize Decision Tree
rpart_model_base<- rpart(Grade~. -Project -Case_ID -Primary_Diagnosis, genes, cp=0.02)
rpart.plot(rpart_model_base)

# plot variable importance for predictors
vip(rpart_model_base)

# create filtered Classification Tree
rpart_model<- rpart(Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, genes, cp=0.02)
# customize appearance of visualization (see milbo.org/rpart-plot/prp.pdf)
rpart.plot(rpart_model, type=3, clip.right.labs= FALSE, branch= .3, under= TRUE)


#############################################################################
# CLASSIFICATION TREE #
# model on just genes, removing IDH1 and Age, using the same process as above
model <- decision_tree() %>%
  set_mode("classification") %>%
  set_engine("rpart") %>% 
  fit(Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, data= genes_train)
preds_tprob<- predict(model, genes_test, type="prob")
preds_tprob
preds_tree<- data.frame(Grade= genes_test$Grade, preds_tprob)
preds_tree<- preds_tree %>% mutate(model= "Tree")
preds_tree

# calculate AUC
roc_auc(preds_tree_prob, Grade, .pred_GBM)

# change to class predictions for additional metrics
preds_tclass<- predict(model, genes_test, type="class")
preds_tclass

# combine actual test Grades with class predictions
preds_tree_class<- data.frame(Grade= genes_test$Grade, preds_tclass)
preds_tree_class<- preds_tree_class %>% mutate(model= "Tree")
preds_tree_class

# calculate class metrics
accuracy(preds_tree_class, Grade, .pred_class)
precision(preds_tree_class, Grade, .pred_class)
recall(preds_tree_class, Grade, .pred_class)
f_meas(preds_tree_class, Grade, .pred_class)

# create confusion matrix
conf_mat(preds_tree_class, Grade, .pred_class)


###############################################
# BAGGING #
# create the specification
spec_bagged<- bag_tree() %>% set_mode("classification") %>% set_engine("rpart", times= 100)

# fit model to the training data
model_bagged<- fit(spec_bagged, formula= Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, data= genes_train)
model_bagged

# predicted probabilities of test data
preds_bag<- predict(model_bagged, genes_test, type="prob")
preds_bagging<- data.frame(Grade= genes_test$Grade, preds_bag)
preds_bagging<- preds_bagging %>% mutate(model= "Bagging")
preds_bagging

# change to class predictions for additional metrics
preds_bagclass<- predict(model_bagged, genes_test, type="class")
preds_bagclass

# combine actual test Grades with class predictions
preds_bag_class<- data.frame(Grade= genes_test$Grade, preds_bagclass)
preds_bag_class<- preds_bag_class %>% mutate(model= "Tree")
preds_bag_class

# bagging AUC
roc_auc(preds_bagging, Grade, .pred_GBM)

# bagging class metrics
accuracy(preds_bag_class, Grade, .pred_class)
precision(preds_bag_class, Grade, .pred_class)
recall(preds_bag_class, Grade, .pred_class)
f_meas(preds_bag_class, Grade, .pred_class)


###############################################################################
# RANDOM FOREST #
# specify random forest
spec<- rand_forest() %>% set_mode("classification") %>% 
  # specify algorithm that controls node split
  set_engine("ranger", importance= "impurity")

# train random forest
forest_model<- spec %>% fit(Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, data= genes_train)

# run model on test
preds_f<- forest_model %>% predict(genes_test, type="prob")

# create data frame with actual test Grade and modeled test predictions
preds_forest<- data.frame(Grade= genes_test$Grade, preds_f)
preds_forest<- preds_forest %>% mutate(model= "Random Forest")
preds_forest

# change to class predictions for additional metrics
preds_fclass<- predict(forest_model, genes_test, type="class")
preds_fclass

# combine actual test Grades with class predictions in a data frame
preds_forest_class<- data.frame(Grade= genes_test$Grade, preds_fclass)
preds_forest_class<- preds_forest_class %>% mutate(model= "Tree")
preds_forest_class

# calculate AUC
roc_auc(preds_forest, Grade, .pred_GBM)

# calculate class metrics
accuracy(preds_forest_class, Grade, .pred_class)
precision(preds_forest_class, Grade, .pred_class)
recall(preds_forest_class, Grade, .pred_class)
f_meas(preds_forest_class, Grade, .pred_class)

# create confusion matrix
conf_mat(preds_forest_class, Grade, .pred_class)


###############################################################################
# BOOSTED #

# specify the model class
boost_spec<- boost_tree() %>% set_mode("classification") %>% set_engine("xgboost")

# fit a boosted model
boost_model <- fit(boost_spec, Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, genes_train)

# assign number of folds
set.seed(1998)
folds<- vfold_cv(genes_train, v=10)

# out-of-the-box performance for untuned model
cv_results<- fit_resamples(boost_spec, Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, resamples= folds, metrics= metric_set(roc_auc))
collect_metrics(cv_results, summarize= FALSE)
collect_metrics(cv_results)

# create tuning spec with placeholders
boost_spec<- boost_tree(trees= 100, learn_rate= tune(), tree_depth= tune(), sample_size= tune()) %>% set_mode("classification") %>% set_engine("xgboost")

# create tuning grid
tunegrid_boost<- grid_regular(extract_parameter_set_dials(boost_spec), levels= 4)

# tune along the grid
tune_results<- tune_grid(boost_spec, Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, folds, tunegrid_boost, metric_set(roc_auc))

# visualize the results
autoplot(tune_results)

# select the best parameters for tuning
best_params<- select_best(tune_results)
best_params

#finalize the model
final_spec<- finalize_model(boost_spec, best_params)
final_spec

# train the final model
final_model<- final_spec %>% fit(Grade~.-Project-Case_ID-Primary_Diagnosis-IDH1-Age_at_diagnosis, data= genes_train)
final_model

preds_b<- predict(final_model, genes_test, type="prob")
preds_b

preds_boosting<- data.frame(Grade= genes_test$Grade, preds_b)
preds_boosting<- preds_boosting %>% mutate(model= "Boosting")
preds_boosting

# calculate roc for final model
roc_auc(preds_boosting, Grade, .pred_GBM)

# change to class predictions for additional metrics
preds_bclass<- predict(final_model, genes_test, type="class")
preds_bclass

# combine actual test Grades with class predictions
preds_boost_class<- data.frame(Grade= genes_test$Grade, preds_bclass)
preds_boost_class<- preds_boost_class %>% mutate(model= "Boost")
preds_boost_class

# calculate AUC
roc_auc(preds_boosting, Grade, .pred_GBM)

# calculate additional class metrics
accuracy(preds_boost_class, Grade, .pred_class)
precision(preds_boost_class, Grade, .pred_class)
recall(preds_boost_class, Grade, .pred_class)
f_meas(preds_boost_class, Grade, .pred_class)

# create confidence matrix
conf_mat(preds_boost_class, Grade, .pred_class)


###########################################################
# MODEL COMPARISONS #

# Calculate AUC for each model and add a model type column
auc_preds_tree <- preds_tree %>%
  roc_auc(truth = Grade, .pred_GBM) %>%
  mutate(model = "Tree")

auc_preds_bagging <- preds_bagging %>%
  roc_auc(truth = Grade, .pred_GBM) %>%
  mutate(model = "Bagging")

auc_preds_forest <- preds_forest %>%
  roc_auc(truth = Grade, .pred_GBM) %>%
  mutate(model = "Forest")

auc_preds_boosting <- preds_boosting %>%
  roc_auc(truth = Grade, .pred_GBM) %>%
  mutate(model = "Boosting")

# Combine AUC results with model names using bind_rows()
auc_results <- bind_rows(
  auc_preds_tree,
  auc_preds_bagging,
  auc_preds_forest,
  auc_preds_boosting)

# Print the combined AUC results
print(auc_results)

preds_combined<- bind_rows(preds_tree, preds_bagging, preds_forest, preds_boosting)
preds_combined

cutoffs<- preds_combined %>% group_by(model) %>% roc_curve(Grade, .pred_GBM)
autoplot(cutoffs)