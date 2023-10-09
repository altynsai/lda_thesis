######### GENERATING DATA #############3
cov_mat <- function(p){
  # Dividing variables into blocks
  p1 <- p*0.3
  p2 <- p*0.3
  p3 <- p*0.4
  
  #setting covariance matrix
  cor1 <- matrix(1,nrow = p1, ncol = p1)
  cor1[upper.tri(cor1)] <- rep(0.8, p1*(p1-1)/2)
  cor1[lower.tri(cor1)] <- rep(0.8, p1*(p1-1)/2)
  
  cor2 <- matrix(1,nrow = p2, ncol = p2)
  cor2[upper.tri(cor2)] <- rep(0.6, p2*(p2-1)/2)
  cor2[lower.tri(cor2)] <- rep(0.6, p2*(p2-1)/2)
  
  
  cor3 <- matrix(1,nrow = p3, ncol = p3)
  cor3[upper.tri(cor3)] <- rep(0, p3*(p3-1)/2)
  cor3[lower.tri(cor3)] <- rep(0, p3*(p3-1)/2)
  
  m <- matrix(0, nrow = p, ncol = p)
  m[1:p1,1:p1] <- cor1
  m[(p1+1):(p1+p2),(p1+1):(p1+p2)] <- cor2
  m[(p1+p2+1):p,(p1+p2+1):p] <- cor3
  return(m)
}


simulated_data_generator <- function(m, p, k, n_obs){
  p1 <- p*0.3
  df <- data.frame()
  # Generating data
  for(i in 1:k){
    mu <- rep(0,p)
    mu[i] <- 1.8
    mu[i+p1] <- -1.8
    df_class <- as.data.frame(mvrnorm(n = n_obs, mu = mu, Sigma = m))
    df_class <- cbind(df_class, rep(i,n_obs))
    df <- rbind(df, df_class)
  }
  
  df[,p+1] <- as.factor(df[,p+1])
  names(df)[p+1] <- "target"
  df <- df[sample(1:nrow(df)), ] #shuffling rows
  return(df)
}


########### INCORPORATING CONTAMINATION INTO DATA ############
incorp_outs <- function(data, nclass, epsilon, k, sdev_contam){
  p1 <- (ncol(data)-1)*0.3
  for (i in 1:nclass){
    numb_out <- round(nrow(data)*epsilon)
    index_out <- sample(1:nrow(data), numb_out, replace=FALSE)
    data[index_out,i] <- rnorm(numb_out, mean=k, sd=sdev_contam)
    index_out <- sample(1:nrow(data), numb_out, replace=FALSE)
    data[index_out,i+p1] <- rnorm(numb_out, mean=k, sd=sdev_contam)
  }
  return(data)
}


##### NORMALIZING DATA ######
med_matrix <- function(data){
  meds <- data %>% group_by(target)%>%summarise_all(median)
  meds_df <- data.frame()
  for(i in data$target){
    meds_df <- rbind(meds_df, filter(meds, target==i))
  }
  return(meds_df[,2:ncol(meds_df)])
}



########### SOLVING LDA PROBLEM FOR GIVEN PRECISION MATRIX ESTIMATOR ##############
lda_by_scatmat <- function(train, test, test_c, Omega){
  p <- ncol(train)-1
  n_train <- nrow(train)
  n_test <- nrow(test)
  k <- nlevels(train$target)
  
  #calculating between-class covariance matrix
  total_med <- apply(train[,1:p], 2, median)
  meds_total <- matrix(rep(total_med,n_train),nrow=n_train, byrow=TRUE)
  meds_by_class <- med_matrix(train)
  Sigma_b <- 1/n_train*t(meds_total - as.matrix(meds_by_class))%*%(meds_total - as.matrix(meds_by_class))
  
  #solving eigenvalue problem
  disc_vect <- Re(eigen(Omega%*%Sigma_b)$vectors[,1:(k-1)])
  train_Y<-as.matrix(train[,1:p])%*% disc_vect
  model <- NaiveBayes(train_Y, train$target)
  
  #testing the model
  test_Y <- as.matrix(test[,1:p])%*% disc_vect
  colnames(test_Y) <- model$varnames
  pred_test <- predict(model, test_Y)
  cm_test <- confusionMatrix(pred_test$class, reference = test$target)
  acc_test <- cm_test$overall['Accuracy']
  
  #testing the model
  test_Y <- as.matrix(test_c[,1:p])%*% disc_vect
  colnames(test_Y) <- model$varnames
  pred_test <- predict(model, test_Y)
  cm_test <- confusionMatrix(pred_test$class, reference = test$target)
  acc_test_c <- cm_test$overall['Accuracy']
  
  return(list(acc = acc_test, dv = disc_vect, acc_c = acc_test_c))
}

KLD <- function(sigma_true, sigma_est, k){
  kl_dist <- (-log(det(sigma_est%*%solve(sigma_true)))+sum(diag(sigma_est%*%solve(sigma_true))))*k - k*ncol(sigma_true)
  return (kl_dist)
}

########## FUNCTION TAKEN FROM ARTICLES #############
kendall.transformed<-function(x){
  ########
  # INPUT
  ########
  # x : Data matrix of dimension N X p (N: number of observations, p: number of variables)
  ########
  # OUTPUT
  ########
  # Cellwise robust covariance matrix estimate based on pairwise Kendall correlation and the Qn scale estimate
  
  x.q=apply(x,2,Qn)
  cor.kendall<- cor.fk(x)
  sigma.kendall=diag(x.q)%*% cor.kendall%*%diag(x.q)
  return(easy.psd(sigma.kendall))
}

easy.psd<-function(sigma)
{
  # returns the nearest semipositive definite matrix from sigma
  eig<-eigen(sigma, symmetric=T)
  d<-pmax(eig$values,0)
  sigma.psd<-eig$vectors%*%diag(d)%*%t(eig$vectors)
  return(sigma.psd)
}


####### COMPARING THE MODELS #########
lda_model_results <- function(train, test, test_c, m){
  acc_results <- c()
  acc_cont_results <- c()
  kld_results <- c()
  time_results <- c()
  num_results <- c()
  p <- ncol(train)-1
  k <- nlevels(train$target)
  nrows <- nrow(train)
  
  
  # Median centered data
  start <- Sys.time()
  train_med <- med_matrix(train)
  train_med_cent <- train[,1:p]-train_med
  tot_med <- matrix(rep(apply(train[,1:p], 2, median),nrows),byrow=T, nrow = nrows)
  end <- Sys.time()
  med_time <- end-start
  
  #Regularized LDA
  start <- Sys.time()
  cov_gl <- CVglasso(train_med_cent, nlam=10, lam.min.ratio = 0.1, start = 'cold', K=3)
  model_gl <- lda_by_scatmat(train, test, test_c, cov_gl$Omega/(nrows-1)*(nrows-k))
  acc <- model_gl$acc
  acc_c <- model_gl$acc_c
  dv <- apply(as.matrix(model_gl$dv),1,sum)
  num_feat <- sum(dv!=0)
  end <- Sys.time()
  
  #saving info
  acc_results <- append(acc_results, acc) #GL-LDA
  acc_cont_results <- append(acc_cont_results, acc_c) #GL-LDA
  kld_results <- append(kld_results, KLD(m, cov_gl$Sigma*(nrows-1)/(nrows-k), k))
  time_results <- append(time_results, as.numeric(end-start+med_time))
  num_results <- append(num_results, num_feat)
  
  #robust LDA
  start <- Sys.time()
  rcov <- kendall.transformed(train_med_cent)*(nrows-1)/(nrows-k) #cellwise robust cov
  model_rlda <- lda_by_scatmat(train, test, test_c, solve(rcov)) #r-LDA
  end <- Sys.time()
  acc <- model_rlda$acc
  acc_c <- model_rlda$acc_c
  dv <- apply(as.matrix(model_rlda$dv),1,sum)
  num_feat <- sum(dv!=0)
  
  #saving info
  acc_results <- append(acc_results, acc)
  acc_cont_results <- append(acc_cont_results, acc_c)
  kld_results <- append(kld_results, KLD(m, rcov, k))
  time_results <- append(time_results, as.numeric(end-start+med_time))
  num_results <- append(num_results, num_feat)
  
  start <- Sys.time()
  cov_rgl <- CVglasso(train_med_cent, S = rcov, nlam=10, K=3, lam.min.ratio = 0.1,
                      start = 'cold')
  model_rgl <- lda_by_scatmat(train, test, test_c,cov_rgl$Omega/(nrows-1)*(nrows-k)) #rGL-LDA
  acc <- model_rgl$acc
  acc_c <- model_rlda$acc_c
  dv <- apply(as.matrix(model_rgl$dv),1,sum)
  num_feat <- sum(dv!=0)
  end <- Sys.time()
  acc_results <- append(acc_results, acc)
  acc_cont_results <- append(acc_cont_results, acc_c)
  kld_results <- append(kld_results, KLD(m, cov_rgl$Sigma*(nrows-1)/(nrows-k), k))
  time_results <- append(time_results, as.numeric(end-start+med_time))
  num_results <- append(num_results, num_feat)
  
  #Penalized LDA
  start <- Sys.time()
  lam_grid <- seq(2,7,length.out = 10)
  class.cent <- train[,1:p] - med_matrix(train)
  wcsd <- apply(class.cent,2,sd)*sqrt((nrows-1)/(nrows-k))
  model_plda <- RPLDA(train, test, test_c, wcsd, lam_grid)
  acc <- model_plda$acc
  acc_c <- model_plda$acc_c
  dv <- apply(as.matrix(model_plda$dv),1,sum)
  num_feat <- sum(dv!=0)
  end <- Sys.time()
  
  acc_results <- append(acc_results, acc)
  acc_cont_results <- append(acc_cont_results, acc_c)
  kld_results <- append(kld_results, KLD(m, diag(wcsd**2), k))
  time_results <- append(time_results, as.numeric(end-start))
  num_results <- append(num_results, num_feat)
  
  #Penalized robust LDA
  start <- Sys.time()
  class.cent <- train[,1:p] - med_matrix(train)
  wcsd <- apply(class.cent,2,Qn)*sqrt((nrows-1)/(nrows-k))
  model_rplda <- RPLDA(train, test, test_c, wcsd, lam_grid)
  acc <- model_rplda$acc
  acc_c <- model_plda$acc_c
  dv <- apply(as.matrix(model_rplda$dv),1,sum)
  num_feat <- sum(dv!=0)
  end <- Sys.time()

  acc_results <- append(acc_results, acc)
  acc_cont_results <- append(acc_cont_results, acc_c)
  kld_results <- append(kld_results, KLD(m, diag(wcsd**2), k))
  time_results <- append(time_results, as.numeric(end-start))
  num_results <- append(num_results, num_feat)
  
  return(list(acc = acc_results, kld = kld_results, time = time_results, 
              num_feat = num_results, acc_c = acc_cont_results))
}



RPLDA <- function(train, test, test_c, wcsd, lam_grid){
  p <- ncol(train)-1
  k <- nlevels(train$target)
  nrows <- nrow(train)
  
  #standardizing the data
  tot.cent <- apply(train[,1:p],2,median)
  train_stand <- cbind(scale(train[,1:p], center=tot.cent, scale = wcsd),
                       target = train$target)
  test_stand <- scale(test[,1:p], center=tot.cent, scale = wcsd)
  test_c_stand <- scale(test[,1:p], center=tot.cent, scale = wcsd)
  sigma.bet <- data.matrix(as.data.frame(train_stand) %>% group_by(target) 
                           %>% summarise_all(median))[,2:(p+1)]
  best_acc <- 0
  best_para <- NULL
  for (para in lam_grid){
    pca <- elasticnet::arrayspc(sigma.bet, K=k-1, para = rep(para,k-1))
    disc_vect <- as.matrix(pca$loadings)
    train_new <- as.matrix(train_stand[,1:p])%*%disc_vect
    colnames(train_new) <- sprintf("DV%s",seq(1:(k-1)))
    model_cv <- train(train_new,train$target,'nb',
                   trControl=trainControl(method='repeatedcv',number=3, repeats = 5))
    if(model_cv$results$Accuracy[2]>= best_acc){
      best_para <- para
      best_acc <- model_cv$results$Accuracy[2]
    }
  }
  print(best_para)
  pca <- elasticnet::arrayspc(sigma.bet, K=k-1, para = rep(best_para,k-1))
  disc_vect <- as.matrix(pca$loadings)
  train_new <- as.matrix(train_stand[,1:p])%*%disc_vect
  colnames(train_new) <- sprintf("DV%s",seq(1:(k-1)))
  model <- NaiveBayes(train_new, train$target, usekernel = T)
  
  #testing the model
  test_new <- as.matrix(test_stand)%*% disc_vect
  colnames(test_new) <- model$varnames
  pred_test <- predict(model, test_new)
  cm_test <- confusionMatrix(pred_test$class, reference = test$target)
  acc_test <- cm_test$overall['Accuracy']
  
  #testing the model
  test_new <- as.matrix(test_c_stand)%*% disc_vect
  colnames(test_new) <- model$varnames
  pred_test <- predict(model, test_new)
  cm_test <- confusionMatrix(pred_test$class, reference = test$target)
  acc_cont_test <- cm_test$overall['Accuracy']
  
  
  
  return(list(acc = acc_test, dv = disc_vect, acc_c = acc_cont_test))
  
}
