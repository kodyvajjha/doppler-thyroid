setwd("~/Downloads")
library("caret")
library("stargazer")
library("gmodels")
library("car")
library("tables")
data<-read.csv("MS- 15-9- HALO, NAMES AGE DITED1.csv")


# Removing the extra rows at the end #
datanew <- data[-c(135,136,137,138,139,140),]
View(datanew)


### 1. To assess sensitivity and specificity of doppler ###
###    in diagnosing malignant thyroid nodules.         ###
datanew$DOPPLER.DIAG
preds1=factor(datanew$DOPPLER.DIAG)
truth1=factor(datanew$FINAL.DIAGNOSIS)
sensitivity(preds1,truth1,positive = "MALIGNANT") #prop.of people among all who have disease who test positive
specificity(preds1,truth1,negative = "BENIGN") #prop.of people among all who don't have disease who test negative
posPredValue(preds1,truth1)
negPredValue(preds1,truth1)

### 2. To assess sensitivity and specificity of vascularity ###
###    in diagnosing malignant thyroid nodules.             ###
###    Optimal cutoff has been found to be 4 by trial.      ###

datanew$VASC.DIAG <- ifelse(as.integer(paste(datanew$VASCULARITY)) >= 4 , "MALIGNANT", "BENIGN")
preds2=factor(datanew$VASC.DIAG)
truth2=factor(datanew$FINAL.DIAGNOSIS)
sensitivity(preds2,truth2,positive= "MALIGNANT")
specificity(preds2,truth2,negative= "BENIGN")
datanew$FINAL.DIAGNOSIS

### 3. To assess sensitivity and specificity of resistive index ###
###    in diagnosing malignant thyroid nodules.                 ###
sens=c()
specs=c()
dist=c()
for (i in seq(0,1,0.0001)){
datanew$RI.DIAG <- ifelse(as.numeric(paste(datanew$RI))>i,"MALIGNANT","BENIGN")
preds3=factor(datanew$RI.DIAG)
truth3=factor(datanew$FINAL.DIAGNOSIS)
sens=append(sens,sensitivity(preds3,truth3,positive="MALIGNANT"))
specs=append(specs,specificity(preds3,truth3,negative="BENIGN"))
}

# ROC curve
plot(1-specs,sens,type="b",col="red")
lines(seq(0,1,0.0001),seq(0,1,0.0001),type="l")

# Finding the optimal cutoff. We compute distances of sens,1-specs from (1,0)
# which is the perfect classification point in the ROC space.
d=sqrt((1-sens)^2 + (1-specs)^2)
which.min(d) # equals 6701

optcutoff=seq(0,1,0.0001)[6701] # equals 0.67

# Required maximal sensitivity and specificity #
sens[6701] # equals 77 percent
specs[6701]# equals 84 percent

### 5. To find association between calcification and final diagnosis ###

CrossTable(datanew$CAL,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$CAL,datanew$FINAL.DIAGNOSIS,simulate.p.value = TRUE)

# p-value is very small, so reject the null hypothesis.
# There IS an association between Calcification and the Final Diagnosis.

### 6. Association between variables and final diagnosis. ###

# Gender vs. Final Diagnosis
datanew$SEX.INDICATOR <- ifelse(datanew$SEX == "M",0,1)
CrossTable(datanew$SEX.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$SEX.INDICATOR,datanew$FINAL.DIAGNOSIS, simulate.p.value = TRUE)

## Accept the null hypothesis. No association.


# No. of nodules vs. Final Diagnosis. 
datanew$NODULE.INDICATOR <- ifelse(datanew$NO.OF.NODULES == 1, 1,0)
CrossTable(datanew$NODULE.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$NODULE.INDICATOR,datanew$FINAL.DIAGNOSIS, simulate.p.value = TRUE)

## Accept the null hypothesis. No association.


# Echo genecity vs. Final Diagnosis. 
datanew$ECHO.INDICATOR <- ifelse(datanew$ECHO == "HYPO", 1,0)
CrossTable(datanew$ECHO.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$ECHO.INDICATOR,datanew$FINAL.DIAGNOSIS,simulate.p.value = TRUE)

## Reject the null hypothesis. There is association.


# Calcification vs. Final Diagnosis.
datanew$CAL.INDICATOR <- ifelse(datanew$CAL == "FINE",1,0)
CrossTable(datanew$CAL.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$CAL.INDICATOR,datanew$FINAL.DIAGNOSIS,simulate.p.value = TRUE)

## Reject the null hypothesis. There is association.

# Peripheral Halo vs Final Diagnosis.
datanew$HALO.INDICATOR <- ifelse(datanew$HALO == "-",1,0)
CrossTable(datanew$HALO.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$HALO.INDICATOR,datanew$FINAL.DIAGNOSIS,simulate.p.value = TRUE)

## Accept the null hypothesis. There is no association.

# Vascularity vs Final Diagnosis.
datanew$VASCULARITY.INDICATOR <- ifelse(as.integer(paste(datanew$VASCULARITY)) > 3, 1,0)
CrossTable(datanew$VASCULARITY.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$VASCULARITY.INDICATOR,datanew$FINAL.DIAGNOSIS,simulate.p.value = TRUE)

## Reject the null hypothesis. There is association

# AP Trans vs. Final Diagnosis
datanew$AP.TRANS.INDICATOR <- ifelse(datanew$AP.TRANS == ">1",1,0)
CrossTable(datanew$AP.TRANS.INDICATOR,datanew$FINAL.DIAGNOSIS)
chisq.test(datanew$AP.TRANS.INDICATOR,datanew$FINAL.DIAGNOSIS,simulate.p.value = TRUE)

## Reject the null hypothesis. There is association. 


### 8. Identifying predictors of malignancy. ###
### Dropping variables NODULE.INDICATOR, HALO.INDICATOR, SEX ###
### because they don't show association with the response ###

datanew$FINAL.DIAGNOSIS.INDICATOR <- ifelse(datanew$FINAL.DIAGNOSIS == "MALIGNANT",1,0)
mylogit <- glm(FINAL.DIAGNOSIS.INDICATOR ~ AGE +ECHO.INDICATOR + CAL.INDICATOR + VASCULARITY.INDICATOR + AP.TRANS.INDICATOR + RI, data = datanew, family = "binomial", maxit=100) 
vif(mylogit)
sqrt(vif(mylogit))>2

mylogit2 <- glm(FINAL.DIAGNOSIS.INDICATOR ~ VASCULARITY.INDICATOR + RI, data = datanew, family = "binomial", maxit=100) 
summary(mylogit2)
