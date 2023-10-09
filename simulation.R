################## SIMULATION STUDY ###########################

rm(list=ls())
#packages
library(MASS)
library(base)
library(ggplot2)
library(caret)
library(dplyr)
library(geigen)
library(CVglasso)
library(penalizedLDA)
library(matrixStats)
library(pcaPP)
library(klaR)
library(robustbase)
library(rrcovHD)
#importing all the functions
setwd('C:/Users/altyn/Desktop/Thesis/code - final version/Chapter 2. Simulations ')
source('functions.R')

### Simulation study ###
#Generating data
set.seed(2023)

#setting parameters
n_train <- 20
k <- 5
n_test <- 20
p <- 200

############ Generating train and test data ###########
#initializing accuracy tables
#cols <- c('GL-LDA', 'rLDA','rGL-LDA', 'PLDA', 'RPLDA')
nrows <- 100
ncols <- 5
table_clean <- matrix(data = NA, nrow=nrows, ncol = ncols)
table_clean_kld <- table_clean
table_clean_time <- table_clean
table_clean_feat <- table_clean

table_5_cont_10 <- table_clean
table_5_cont_10_c <- table_clean
table_5_cont_10_kld <- table_clean
table_5_cont_10_time <- table_clean
table_5_cont_10_feat <- table_clean

table_5_cont_20 <- table_clean
table_5_cont_20_c <- table_clean
table_5_cont_20_kld <- table_clean
table_5_cont_20_time <- table_clean
table_5_cont_20_feat <- table_clean

table_10_cont_10 <- table_clean
table_10_cont_10_c <- table_clean
table_10_cont_10_kld <- table_clean
table_10_cont_10_time <- table_clean
table_10_cont_10_feat <- table_clean

table_10_cont_20 <- table_clean
table_10_cont_20_c <- table_clean
table_10_cont_20_kld <- table_clean
table_10_cont_20_time <- table_clean
table_10_cont_20_feat <- table_clean



for(i in 1:nrows){
  print(paste0('Iteration #', i))
  #generating data and train/test samples
  m <- cov_mat(p)
  train_wo_out <- simulated_data_generator(m, p, k, n_train)
  test_wo_out <- simulated_data_generator(m, p, k, n_test)
  results <- lda_model_results(train_wo_out, test_wo_out,test_wo_out, m)
  table_clean[i, ] <- results$acc
  table_clean_kld[i, ] <- results$kld
  table_clean_time[i, ] <- results$time
  table_clean_feat[i, ] <- results$num_feat

  #contaminated 10% with mean 5
  train_5_cont_10 <- incorp_outs(train_wo_out, k,0.1, 5, 1)
  test_5_cont_10 <- incorp_outs(test_wo_out, k, 0.1, 5, 1)
  results <- lda_model_results(train_5_cont_10, test_wo_out,test_5_cont_10, m)
  table_5_cont_10[i, ] <- results$acc
  table_5_cont_10_c[i, ] <- results$acc_c
  table_5_cont_10_kld[i, ] <- results$kld
  table_5_cont_10_time[i, ] <- results$time
  table_5_cont_10_feat[i, ] <- results$num_feat

  #contaminated 30% with mean 5
  train_5_cont_20 <- incorp_outs(train_wo_out, k, 0.2, 5, 1)
  test_5_cont_20 <- incorp_outs(test_wo_out, k, 0.2, 5, 1)
  results <- lda_model_results(train_5_cont_20, test_wo_out, test_5_cont_20, m)
  table_5_cont_20[i, ] <- results$acc
  table_5_cont_20_c[i, ] <- results$acc_c
  table_5_cont_20_kld[i, ] <- results$kld
  table_5_cont_20_time[i, ] <- results$time
  table_5_cont_20_feat[i, ] <- results$num_feat

  #contaminated 10% with mean 10
  train_10_cont_10 <- incorp_outs(train_wo_out, k, 0.1, 10, 1)
  test_10_cont_10 <- incorp_outs(test_wo_out,k, 0.1, 10, 1)
  results <- lda_model_results(train_10_cont_10, test_wo_out,test_10_cont_10, m)
  table_10_cont_10[i, ] <- results$acc
  table_10_cont_10_c[i, ] <- results$acc_c
  table_10_cont_10_kld[i, ] <- results$kld
  table_10_cont_10_time[i, ] <- results$time
  table_10_cont_10_feat[i, ] <- results$num_feat

  #contaminated 30% with mean 10
  train_10_cont_20 <- incorp_outs(train_wo_out, k, 0.2, 10, 1)
  test_10_cont_20 <- incorp_outs(test_wo_out, k, 0.2, 10, 1)
  results <- lda_model_results(train_10_cont_20, test_wo_out,test_10_cont_20, m)
  table_10_cont_20[i, ] <- results$acc
  table_10_cont_20_c[i, ] <- results$acc_c
  table_10_cont_20_kld[i, ] <- results$kld
  table_10_cont_20_time[i, ] <- results$time
  table_10_cont_20_feat[i, ] <- results$num_feat
  
  
}

#saving tables
write.csv(table_clean, "tables\\clean_data_acc.csv", row.names=FALSE)
write.csv(table_clean_time, "tables\\clean_data_time.csv", row.names=FALSE)
write.csv(table_clean_kld, "tables\\clean_data_kld.csv", row.names=FALSE)
write.csv(table_clean_feat, "tables\\clean_data_feat.csv", row.names=FALSE)

write.csv(table_5_cont_10, "tables\\cln_5_10_acc.csv", row.names=FALSE)
write.csv(table_5_cont_10_time, "tables\\cln_5_10_time.csv", row.names=FALSE)
write.csv(table_5_cont_10_kld, "tables\\cln_5_10_kld.csv", row.names=FALSE)
write.csv(table_5_cont_10_feat, "tables\\cln_5_10_feat.csv", row.names=FALSE)

write.csv(table_5_cont_10_c, "tables\\cont_5_10_acc.csv", row.names=FALSE)

write.csv(table_5_cont_20, "tables\\cln_5_20_acc.csv", row.names=FALSE)
write.csv(table_5_cont_20_time, "tables\\cln_5_20_time.csv", row.names=FALSE)
write.csv(table_5_cont_20_kld, "tables\\cln_5_20_kld.csv", row.names=FALSE)
write.csv(table_5_cont_20_feat, "tables\\cln_5_20_feat.csv", row.names=FALSE)

write.csv(table_5_cont_20_c, "tables\\cont_5_20_acc.csv", row.names=FALSE)

write.csv(table_10_cont_10, "tables\\cln_10_10_acc.csv", row.names=FALSE)
write.csv(table_10_cont_10_time, "tables\\cln_10_10_time.csv", row.names=FALSE)
write.csv(table_10_cont_10_kld, "tables\\cln_10_10_kld.csv", row.names=FALSE)
write.csv(table_10_cont_10_feat, "tables\\cln_10_10_feat.csv", row.names=FALSE)


write.csv(table_10_cont_10_c, "tables\\cont_10_10_acc.csv", row.names=FALSE)


write.csv(table_10_cont_20, "tables\\cln_10_20_acc.csv", row.names=FALSE)
write.csv(table_10_cont_20_time, "tables\\cln_10_20_time.csv", row.names=FALSE)
write.csv(table_10_cont_20_kld, "tables\\cln_10_20_kld.csv", row.names=FALSE)
write.csv(table_10_cont_20_feat, "tables\\cln_10_20_feat.csv", row.names=FALSE)


write.csv(table_10_cont_20_c, "tables\\cont_10_20_acc.csv", row.names=FALSE)


#plots
p <- ggplot(train_wo_out, aes(train_wo_out[,1],train_wo_out[,2], 
                        color = train_wo_out$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")

p <- ggplot(train_wo_out, aes(train_wo_out[,31],train_wo_out[,32], 
                             color = train_wo_out$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")


p <- ggplot(train_5_cont_10, aes(train_5_cont_10[,1],train_5_cont_10[,2], 
                                 color = train_5_cont_10$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(5,1) with eps = 0.1")

p <- ggplot(train_5_cont_10, aes(train_5_cont_10[,31],train_5_cont_10[,32], 
                                 color = train_5_cont_10$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(5,1) with eps = 0.1")

p <- ggplot(train_5_cont_20, aes(train_5_cont_20[,1],train_5_cont_20[,2], 
                                 color = train_5_cont_20$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(5,1) with eps = 0.2")

p <- ggplot(train_5_cont_20, aes(train_5_cont_20[,31],train_5_cont_20[,32], 
                                 color = train_5_cont_20$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(5,1) with eps = 0.2")



p <- ggplot(train_10_cont_10, aes(train_10_cont_10[,1],train_10_cont_10[,2], 
                                 color = train_10_cont_10$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(10,1) with eps = 0.1")

p <- ggplot(train_10_cont_10, aes(train_10_cont_10[,31],train_10_cont_10[,32], 
                                 color = train_10_cont_10$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(10,1) with eps = 0.1")

p <- ggplot(train_10_cont_20, aes(train_10_cont_20[,1],train_10_cont_20[,2], 
                                 color = train_10_cont_20$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(10,1) with eps = 0.2")

p <- ggplot(train_10_cont_20, aes(train_10_cont_20[,31],train_10_cont_20[,32], 
                                 color = train_10_cont_20$target))+geom_point()
p+xlab("First disc feature") + ylab("Second disc feature") + labs(color = "Classes")+
  ggtitle("Outliers coming from N(10,1) with eps = 0.2")
