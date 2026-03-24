setwd("J:/V/R")
rm(list = ls())
f<-read.csv("ACC cox.csv", header = T)
class(f)
library(tidyverse)
library(survival)
Ftable_rownames <- rownames(f)
# transpose
Ftable <- t(f)
Ftable<-data.frame(Ftable)
# assign names from the stored object
#colnames(Ftable) <- Ftable_rownames
colnames(Ftable)<-Ftable [1,]
Ftable<-Ftable[-1,]
rownames(Ftable) <- NULL

y<-Ftable$os
y<-as.numeric(y)
event<-Ftable$event
event<-as.numeric(event)
x<-Ftable[,-1:-2]
rownames(x) <- NULL

MODEL: 1 # Univarite analysis
dat1<-read.csv("TCGA3GEOCU_match.csv", header = T) # inptfile is matrix of gene normalized expression profiles and survival information {Refer TCGA dataset given the link below} 
dat1<-read.csv("Cancer RSEM.csv", header = T)
dat1<-read.csv("ACC_cox.csv", header = T)
#dat1<-t(dat1)
#dat1<-as.data.frame(dat1)
#colnames(dat1)<-dat1[1,]
#dat1<-dat1[-1,]
#row.names(dat1)<-NULL

dat2<-dat1[,-1]
y<-dat2$os
y<-as.numeric(y)
event<-dat2$event
event<-as.numeric(event)
x<-dat2[,-1:-2]
x<-log2(x+1) # if not normalized then use log2
#row.names(x)<-NULL
#**********************************************************************
#Fit the coxph model 
ans = apply(x,2,function(x,y,event)coxph(Surv(y,event)~x),y=y,event=event)

# Extract data 
univ_results <- lapply(ans,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=4)
                         wald.test<-signif(x$wald["test"], digits=4)
                         beta<-signif(x$coef[1], digits=4);#coeficient beta
                         HR <-signif(x$coef[2], digits=4);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
result<-as.data.frame(res)
write.csv(result,"Cancer RSEM.3.csv")
# results of univariate Cox Regression Analysis 
# Next, survival related genes used as inputfile for multivariate Cox Regression with Penalized Model

MODEL: 2 # Univariate and multivariable regression
#https://www.r-bloggers.com/2019/09/survival-analysis-with-strata-clusters-frailties-and-competing-risks-in-in-finalfit/
#https://www.r-bloggers.com/2018/05/elegant-regression-results-tables-and-plots-in-r-the-finalfit-package/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2394368/
#https://finalfit.org/articles/missing.html 
######Cox proportional hazards regression
#install.packages("finalfit")
dat1<-read.csv("Cancer RSEM.csv", header = T)
y<-dat1$os
y<-as.numeric(y)
event<-dat1$event
event<-as.numeric(event)
#install.packages("finalfit")
#install.packages("dplyr")
library(finalfit)
library(survival)
library(dplyr)
names(dat2)
coxph(Surv(y,event) ~ hpv_test_result +grade+RS.2+RS.1+RS1.1+RS2.2+smoking_history+number_of_pregnancies +age+stage+RS2.2.1+RS1.1.1+clinical_stage+neo.cancer_status+race+grade.1+M.stage+N.stage+T.stage, data = dat1) %>% 
  summary()
dependent_os  <- "Surv(y,event)"
explanatory   <- c("hpv_test_result","grade","RS.2","RS.1","RS1.1","RS2.2", "smoking_history","number_of_pregnancies" ,"age","stage","RS2.2.1","RS1.1.1","clinical_stage" ,"neo.cancer_status","race","grade.1","M.stage","N.stage","T.stage", na.omit(T))
p<-dat1 %>% 
  finalfit(dependent_os, explanatory)
dat1 %>%
  finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
  ff_remove_ref() %>% 
  dependent_label(dat1, dependent_os)-> t
# # Cox Proportional Hazards univariable analysis.
p1<-dat1 %>%
  coxphuni(dependent_os, explanatory) %>%
  fit2df(estimate_suffix=" (univariable)")
# Cox Proportional Hazards multivariable analysis.
p2<-dat1  %>%
  coxphmulti(dependent_os, explanatory) %>%
  fit2df(estimate_suffix=" (multivariable)")
# Dataframe. Missing values must be coded NA, Character vector of any length: name(s) of explanatory variables. to a missing
#data pattern (1=observed, 0=missing).
p3<-dat1 %>%
  finalfit(dependent_os, explanatory, na_to_missing=T)
# coxph
fit = coxph(Surv(time, status) ~ age.factor + sex.factor + obstruct.factor + perfor.factor,
            data = colon_s)
fit %>%
  fit2df(estimate_suffix=" (multivariable)")

write.csv(t,"Uni multi.csv")

dat1 %>% hr_plot(dependent_os, explanatory)

MODEL: 3 #Reduced model: a reduced model can be directly specified and compared
explanatory_multi = c("RS.1","RS.2","stage", na.omit(T))
p2<-dat1 %>% 
  finalfit(dependent_os, explanatory, explanatory_multi, 
           keep_models = TRUE)

MODEL: 4: # Testing for proportional hazards
  explanatory = c("age","cancer_status","RS.1","RS1.1","RS2.2","stage", "os")
p3<-dat1 %>% 
  coxphmulti(dependent_os, explanatory) %>% 
  cox.zph() %>% 
  {zph_result <<- .} %>% 
  plot(var=5)  







# ** Multivariate Cox regression with three penalities including least absolute shrinkage and selection operator(Lasso), Adaptive lasso and Elastic net algorithms for informative prognostic-related genes selection.**
# we need two Package in R 
library(glmnet)
library(survival)
dat<-read.csv("TCGA3GEOCU_match.csv", row.names = 1) # inptfile is matrix of survival related DEGs normalized expression profile and survival information 
dat<-read.csv("Cancer RSEM.csv", row.names = 1)
set.seed(1234)
# Prognostic gene with Elastic net algorithm selection 
y <- Surv(dat$os, dat$event) # os-overall survival, event-overall survival status
which(any(is.na(y)))
x<-model.matrix(y~., dat[,c(-1:-2)])
cv.fit <- cv.glmnet(x,y, family="cox",nfold = 10, alpha= 0.5)
plot(cv.fit)
fit <- glmnet(x,y, family = "cox", alpha =0.5)
plot(fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
pp=coef(cv.fit, s = "lambda.min") # extracting selecting prognostic genes
pp
# extracting selecting prognostic genes
# Prognostic gene with Lasso algorithm selection 
y <- Surv(dat$os, dat$event) 
x<-model.matrix(y~., dat[,c(-1:-2)])
cv.fit1 <- cv.glmnet(x,y, family="cox",nfold = 10, alpha= 1)
plot(cv.fit1)
fit1 <- glmnet(x,y, family = "cox", alpha = 1)
plot(fit1)
cv.fit1$lambda.min

Coefficients <- coef(fit1, s = fit1$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
TT=coef(cv.fit1, s = "lambda.min") # extracting selecting prognostic genes
TT
#adaptive lasso
cv.fit2 <- cv.glmnet(x,y, family="cox", nfold = 10, alpha= 0) # rigid Cox Model
plot(cv.fit2)
fit2 <- glmnet(x,y, family = "cox",
              nfold = 10, alpha= 0)
plot(fit2)
cv.fit2$lambda.min

Coefficients <- coef(fit2, s = fit2$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
aa=coef(cv.fit2, s = "lambda.min")
aa
# addaptive lasso
## Extract coefficients at the error-minimizing lambda
cv.fit2$lambda.min
coef(cv.fit2, s = "lambda.min")
## The intercept estimate should be dropped.
best_ridge_coef <- as.numeric(coef(cv.fit2, s = "lambda.min"))
best_ridge_coef
## Perform adaptive LASSO
alasso1 <- glmnet(x,y, family = "cox", alpha = 1,
                  penalty.factor = 1 / abs(best_ridge_coef))
plot(alasso1, xvar = "lambda")

## Perform adaptive LASSO with 10-fold CV
alasso1_cv <- cv.glmnet(x,y, family = "cox",
                        ## type.measure: loss to use for cross-validation.
                        ## K = 10 is the default.
                        nfold = 10,
                        ## 'alpha = 1' is the lasso penalty, and 'alpha = 0' the ridge penalty.
                        alpha = 1,
                        penalty.factor = 1 / abs(best_ridge_coef),
                        ## prevalidated array is returned
                        keep = TRUE)
## Penalty vs CV MSE plot
plot(alasso1_cv)

## Extract coefficients at the error-minimizing lambda
alasso1_cv$lambda.min
## s: Value(s) of the penalty parameter 'lambda' at which
##    predictions are required. Default is the entire sequence used
##    to create the model.
best_alasso_coef1 <- coef(alasso1_cv, s = alasso1_cv$lambda.min)
best_alasso_coef1 # extracting prognostic gene selected by Adaptive lasso 
# Best subset Cox regression model construction 
# Combine all gene selected by three algorithms and use glumti package to identify which combination of genes would produce optimal prognostic signture 

#best subset regression analysis 
library(survival)
#install.packages("glmulti")
library(glmulti)
dat3<-read.csv("TCGA3GEOCU_match.csv", row.names = 1) # inputfile containing the identifed prognostic genes from penalized models and survival time&event information 
dat <- within(dat1, {
  survival.vector    <- Surv(dat1$os, dat1$event)})
#*************************************************************
#FAM83D +	CDC20 +	TPX2 + 	LECT2 + ANXA10 +	DNASE1L3 +	PON1 +	CD5L +	CYP2C9 +	ADH4 +	CFHR3 +	GHR + LCAT
names(dat)
glmulti.coxph.out <-glmulti(survival.vector ~hpv_test_result +grade+RS.2+RS.1+RS1.1+RS2.2+smoking_history+number_of_pregnancies +age+stage+RS2.2.1+RS1.1.1+clinical_stage+neo.cancer_status+race+grade.1+M.stage+N.stage+T.stage, data = dat,
                            level = 1,               # No interaction considered
                            method = "h",            # Exhaustive approach
                            crit = "aic",            # AIC as criteria
                            confsetsize = 5,         # Keep 5 best models
                            plotty = F, report = F,  # No plot or interim reports
                            fitfunction = "coxph")   # coxph function

## Show result for the best model
summary(glmulti.coxph.out@objects[[1]])
## Show 5 best models (Use @ instead of $ for an S4 object)
glmulti.coxph.out@formulas
res=glmulti.coxph.out
res
plot(res)
plot(res, type="s")
summary(res@objects[[1]])
print(res)
top <- weightable(res)
top

# Seperation of Highrisk and Lowrisk by median risk score
dat1<-read.csv("Cancer RSEM.csv")
z<-dat1$RS.1
str(z)
summary(z)
Hirisk<-dat1[dat1$RS.1>-1.0821,] 
lowrisk<-dat1[dat1$RS.1<= -1.0821,] 
write.csv(lowrisk, "LO.csv")
write.csv(Hirisk, "HI.csv")
# Risk score model construction based Best subset  prognostic genes 
model<- coxph(Surv(os,event) ~dat$Grade+dat$Group, data =dat)
summary(model)
library(tidyverse)
library(survminer)
library(ggplot2)
dat<-read.csv("design matrix.csv")
dat<-read.csv("TCGA3GEOCU_match - 2.csv") # Input file containing Risk score calculated from best subset prognostic genes
dat<-read.csv("Cancer RSEM.csv")
#RS<-(dat$RS)
#class(RS)
#summary(RS) #notedown median of RS and devise HR and LR
#High_risk<-dat[dat$RS>= -1.3894,] #UPREGULATED GENES
#Low_risk<-dat[dat$RS<= -1.3894,] #DOWNREGULATED GENES
#write.csv(High_risk,"HI.csv")
#write.csv(Low_risk,"LO.csv")
# in fit Risk score must be high=1, low=0, RS is either 1 or 0
fit<- survfit(Surv(dat$os,dat$event) ~dat$RS3.3, data = dat)# Rs-Risk score
#fit<- survfit(Surv(dat$os,dat$event) ~dat$RS, data = dat)
fit
# normal KM plot
ggsurvplot(fit, data = dat, risk.table = T, conf.int = T,palette = NULL, ggtheme = theme_minimal())
# advanced KM plot 
p<-ggsurvplot(fit, data = dat, pval = T,conf.int = T,
           legend.labs = c("Low Risk", "High Risk"),
           risk.table = F,
           font.x = c(12, "bold", "black"), font.y = c(12, "bold", "black"),
           font.main = c(15, "bold", "black"),
           xlab = "time(months)", 
           ylab = "survival probablity",
           font.tickslab = c(12, "bold", "black"),
           legend.title = "hsa.mir.137+hsa.mir.1251",
           gtheme = theme_classic2(base_size = 12, base_family = "Times New Roman"))
print(p)
ggpar(p, 
      font.main = c(17, "bold"),
      font.x = c(17, "bold"),
      font.y = c(17, "bold"),
      font.caption = c(17, "bold"), 
      font.legend = c(17, "bold"), 
      font.tickslab = c(17, "bold"))
ggsave("J:\\V\\R\\KMsurvival.jpeg", print(p), height = 4, width = 5, dpi = 300)
#http://www.sthda.com/english/wiki/saving-high-resolution-ggplots-how-to-preserve-semi-transparency 
p1 <- ggsurvplot(fit, data = dat,
                surv.median.line = "hv", # Add medians survival
                pval = TRUE,             # Add p-value and tervals
                
                conf.int = TRUE,        # Add the 95% confidence band
                risk.table = TRUE,      # Add risk table
                tables.height = 0.2,
                tables.theme = theme_cleantable(),
                palette = "jco",
                ggtheme = theme_bw()
)
print(p1)
ggsurvplot(fit, data = dat, pval = T,conf.int = T,
           legend.labs = c("High-risk", "Low-risk"),
           risk.table = T,
           xlab = "time(months)", 
           ylab = "survival probablity",
           legend.title = "Risk",
           risk.table.y.text.col = T,
           risk.table.y.text = T,
           gtheme = theme_minimal())

# Prepare dataset having OS time and status and Risk socre each HCC patients
dat<-read.csv("inputfile", row.names = 1) # example CMUH dataset: Survival_4-gene expression-Profiles.csv

#install.packages("timeROC")
library(timeROC)
library(survival)
#1) Kaplan-Meier estimator(default)
# https://zhuanlan.zhihu.com/p/42296548
ROC.bili.marginal<-timeROC(T=pbc$time,
                           delta=pbc$status,marker=pbc$bili,
                           cause=1,weighting="marginal",
                           times=quantile(pbc$time,probs=seq(0.2,0.8,0.1)),
                           iid=TRUE)
ROC.bili.marginal

# 2) Cox model With competing risks
# evaluate DDST cognitive score as a prognostic tool for dementia onset, accounting for death without dementia competing risk.
ROC.DSST<-timeROC(T=dat$os,delta=dat$event,
                  marker=dat$RS3.3,cause=1,
                  weighting="cox",
                  times=c(12,36,60,84,104),ROC=TRUE) # 1st year, 2nd year and 3rd year risk classification performance 
print(ROC.DSST)
plot(ROC.DSST,time=12,title=FALSE,lwd=2)
plot(ROC.DSST,time=24,col="dodgerblue4",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=36,col="darkgreen",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=48,col="brown",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=60,col="black",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=72,col="pink",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=84,col="dodgerblue4",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=104,col="brown",add=TRUE,title=FALSE,lwd=2)
legend(title ="CESC hsa.mir.526b + hsa.mir.508_RS timeROC",legend = c("1-year(AUC=57.69%)","3-year(AUC=63.95%)","5-year(AUC=72.07%)","7-year(AUC=76.57%)","9-year(AUC=73.44%)"),col = 1:2,
  )
legend("bottomright",c("1-year(AUC=71.29%)","3-year(AUC=64.30%)","5-year(AUC=57.27%)","7-year(AUC=58.08%)"),
       col=c("red","darkgreen", "black","dodgerblue4"),lwd=2, cex = 1)
legend("bottomright",c("1-year(AUC=65.40%)","3-year(AUC=69.11%)","5-year(AUC=74.51%)","7-year(AUC=80.81%)","9-year(AUC=74.56%)"),
       col=c("red","darkgreen", "black","dodgerblue4","brown"),lwd=1.5, cex = 1)
lwd=2 # line width in line
cex = 0.6 # size of the box and content size
inset = 0.01 # it is the how much distance away from boarder lines
lty=1:2 # line is breaks like dotted lines 
legend("bottomright",c("bilirubin","albumin","cholesterol"),
       col=c("red","blue","black"),lty=1:2)
ggsave("ggpubrsave.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
#******************************************************
# Diagnostic performance of Risk model h
library(pROC)
library(parallel)
#install.packages("survivalROC")
library(survivalROC)
dat<-read.csv("Cancer RSEM.1.csv") # example CMUH dataset:4-gene signture expression profiles .csv
class(dat)
#***********************=========================================
gene1_ROC <- roc(dat$os ~ dat$RS.1, plot=TRUE,print.auc=TRUE,col="blue",lwd = 4,print.auc.y=0.4,legacy.axes=TRUE,add = TRUE)

model1<-plot.roc(dat$event, dat$RS2.2,       # data
                 percent = FALSE,  # show all values in percent
                 auc=c(0, 100), 
                 print.auc=TRUE,show.thres=TRUE,                    
                 main = "4-gene signatures",
                 ci = TRUE)

glm.fit=glm(dat$RS ~ dat$FAM83D.2, family=binomial)
lines(dat$FAM83D.2, glm.fit$fitted.values)
roc(dat$FAM83D.2, glm.fit$fitted.values, plot=TRUE)
par(pty = "s")

model2<-plot.roc(dat1$RS, dat1$FAM83D,       # data
                 percent = TRUE,  # show all values in percent
                 auc=c(0, 100), 
                 colorspaces=T,
                 print.auc=TRUE,show.thres=TRUE,                    
                 main = "4-gene signatures",
                 ci = TRUE)

library(pROC) # install with install.packages("pROC")
library(randomForest) # install with install.packages("randomForest")

#######################################
##
## Generate weight and obesity datasets.
##
#######################################
set.seed(420) # this will make my results match yours

num.samples <- 100

## genereate 100 values from a normal distribution with
## mean 172 and standard deviation 29, then sort them
weight <- sort(rnorm(n=num.samples, mean=172, sd=29))
weight
## Now we will decide if a sample is obese or not. 
## NOTE: This method for classifying a sample as obese or not
## was made up just for this example.
## rank(weight) returns 1 for the lightest, 2 for the second lightest, ...
##              ... and it returns 100 for the heaviest.
## So what we do is generate a random number between 0 and 1. Then we see if
## that number is less than rank/100. So, for the lightest sample, rank = 1.
## This sample will be classified "obese" if we get a random number less than
## 1/100. For the second lightest sample, rank = 2, we get another random
## number between 0 and 1 and classify this sample "obese" if that random
## number is < 2/100. We repeat that process for all 100 samples
obese <- ifelse(test=(runif(n=num.samples) < (rank(weight)/num.samples)), 
                yes=1, no=0)
obese ## print out the contents of "obese" to show us which samples were
## classified "obese" with 1, and which samples were classified
## "not obese" with 0.

## plot the data
plot(x=weight, y=obese)

## fit a logistic regression to the data...
glm.fit=glm(obese ~ weight, family=binomial)
lines(weight, glm.fit$fitted.values)


#######################################
##
## draw ROC and AUC using pROC
##
#######################################

## NOTE: By default, the graphs come out looking terrible
## The problem is that ROC graphs should be square, since the x and y axes
## both go from 0 to 1. However, the window in which I draw them isn't square
## so extra whitespace is added to pad the sides.
roc(obese, glm.fit$fitted.values, plot=TRUE)

## Now let's configure R so that it prints the graph as a square.
##
par(pty = "s") ## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region
##                "m" - (the default) creates a maximal plotting region
roc(obese, glm.fit$fitted.values, plot=TRUE)

## NOTE: By default, roc() uses specificity on the x-axis and the values range
## from 1 to 0. This makes the graph look like what we would expect, but the
## x-axis itself might induce a headache. To use 1-specificity (i.e. the 
## False Positive Rate) on the x-axis, set "legacy.axes" to TRUE.
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE)

## If you want to rename the x and y axes...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage")

## We can also change the color of the ROC line, and make it wider...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4)

## If we want to find out the optimal threshold we can store the 
## data used to make the ROC graph in a variable...
roc.info <- roc(obese, glm.fit$fitted.values, legacy.axes=TRUE)
str(roc.info)

## and then extract just the information that we want from that variable.
roc.df <- data.frame(
  tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
  fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
  thresholds=roc.info$thresholds)

head(roc.df) ## head() will show us the values for the upper right-hand corner
## of the ROC graph, when the threshold is so low 
## (negative infinity) that every single sample is called "obese".
## Thus TPP = 100% and FPP = 100%

tail(roc.df) ## tail() will show us the values for the lower left-hand corner
## of the ROC graph, when the threshold is so high (infinity) 
## that every single sample is called "not obese". 
## Thus, TPP = 0% and FPP = 0%

## now let's look at the thresholds between TPP 60% and 80%...
roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]

## We can calculate the area under the curve...
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

## ...and the partial area under the curve.
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, print.auc.x=45, partial.auc=c(100, 90), auc.polygon = TRUE, auc.polygon.col = "#377eb822")


#######################################
##
## Now let's fit the data with a random forest...
##
#######################################
rf.model <- randomForest(factor(obese) ~ weight)

## ROC for random forest
roc(obese, rf.model$votes[,1], plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#4daf4a", lwd=4, print.auc=TRUE)


#######################################
##
## Now layer logistic regression and random forest ROC graphs..
##
#######################################
roc(obese, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)

plot.roc(obese, rf.model$votes[,1], percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Logisitic Regression", "Random Forest"), col=c("#377eb8", "#4daf4a"), lwd=4)


#######################################
##
## Now that we're done with our ROC fun, let's reset the par() variables.
## There are two ways to do it...
##
#######################################
par(pty = "s")

