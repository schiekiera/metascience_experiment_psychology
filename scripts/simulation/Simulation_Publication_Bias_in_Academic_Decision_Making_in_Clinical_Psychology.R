# Remove data
rm(list = ls())

# Load necessary libraries
library(lme4)
library(mediation)
library(report)
library(faux)
library(effectsize)
library(misty)
# do not load lmerTest, since it is not compatible with mediation package

# Simulate data
set.seed(2023)  # For reproducibility

###########################
## Introductory Comment  ##
###########################

# This code demonstrates the process of simulating and analyzing 
# multilevel and multilevel mediation data using R for an experimental study
# on publication bias in academic decision making of clinical psychology
# researchers.

##########################
##                      ##
##                      ##
## 1. Multilevel Model  ##
##                      ##
##                      ##
##########################




## Multilevel Model
# IR = gamma +  (beta+V_1j+U_i1)*treatment + V_0j + U_i0 + error


# Simulate data for id, stimulus, experimental_factor
n_id <- 75          # Number of unique IDs
n_stimulus <- 16    # Number of unique Stimuli
n<-n_id*n_stimulus  # Number of rows

## factor for variance partions
f_residual<-0.4
f_stimulus_int<-0.2
f_id_int<-0.2
f_stimulus_slope<-0.1
f_id_slope<-0.1


# Intercepts
## γ_00: General intercept
gamma<-0

## V_0j: Weight random intercept id
V_0j<-rnorm(n_id)*f_id_int

## U_i0: Weight random intercept stimulus
U_i0<-rnorm(n_stimulus)*f_stimulus_int

## V_1j: Weight random slopes id
V_1j<-rnorm(n_id)*f_id_slope

## U_i1: Weight random slope stimulus
U_i1<-rnorm(n_stimulus)*f_stimulus_slope


# Create Ordered Identifier Variables
## i: id
## Order: 1,1,1,1, .... 30,30,30,30
id<-sort(rep(1:n_id,n_stimulus))

# j: stimulus
## Order: 1,2, .. 29,30,1,2, ... 29,30
stimulus <- rep(1:n_stimulus,n_id)

# treatment: Experimental factor
## Levels: 0 = control; 1 =  Experimental

### Urn 15 x 0 and 15 x 1
####  Every participant will be exposed to an equal number of control 
####  and experimental stimuli 
urn<-rep(c(0,1),n_stimulus/2)

### Create experimental factor: Sample different combination for each id
treatment<-rep(sample(urn,replace=FALSE),n_id)

# beta_1: effect
## The effect assumed for our power analysis is cohens d = 0.55
## --> transform to standardized beta = pearson's r
beta<-d_to_r(-0.55)

# Error: R
## Normally distributed random variable
error=rnorm(n)*f_residual

# IR = gamma +  (beta+V_1j+U_i1)*treatment + V_0j + U_i0 + error
IR<-c()
for(i in 1:n_id){
  for (j in 1:n_stimulus){
    IR[j+n_stimulus*(i-1)]<-gamma+(beta+V_1j[i]+U_i1[j])*
      treatment[j+n_stimulus*(i-1)]+V_0j[i]+U_i0[j]+error[j+n_stimulus*(i-1)]
  }
}

## combining the dataframe
df<-data.frame(id=id,
               stimulus=stimulus,
               treatment=treatment,
               IR=IR)

##############
# Null Model #
##############
model_X.1.1 <- lmer(IR ~ 1 + (1|id) + (1|stimulus), data = df)
report(model_X.1.1)

##########################
# Random Intercept Model #
##########################
model_X.1.2 <- lmer(IR ~ treatment + (1|id) + (1|stimulus), data = df)
report(model_X.1.2)

##################################
# Main Model: Random Slope Model #
##################################
model_X.1.3 <- lmer(IR ~ treatment + (treatment|id) + (1|stimulus), data = df)
report(model_X.1.3)

#########################
# Exploratory Variables #
#########################

## recode treatment for interaction effects
df$treatment<-recode(df$treatment, `0` = -1L, `1` = 1L, )
df$treatment[1:20]

## Familiarity with Publication Bias
df$familiarity_PB<-rnorm(n)++df$IR*0.4*runif(n)

## Pressure to Publish
df$pub_pressure<-rnorm(n)+df$IR*0.1*runif(n)

# Professional Status
df$professor<-sample(0:1,n,replace = TRUE,prob=c(0.9,0.1))
df$postdoc<-rep(0,n)

for (i in 1:n){
  sample<-sample(0:1,1,prob=c(0.55,0.45))
  if(df$professor[i]==0){
    df$postdoc[i]<-sample
  }
}

#position
df$position<-rep(sample(1:n_stimulus,n_stimulus),n_id)

######################
# Exploratory Models #
######################

model_X.1.4 <- lmer(IR ~ treatment * familiarity_PB +
                      treatment * pub_pressure +
                      postdoc + professor + position + 
                      (treatment|id) + (1|stimulus), data = df)
report(model_X.1.4)






####################################
##                                ##
##                                ##
## 2. Multilevel Mediation Model  ##
##                                ##
##                                ##
####################################

# Intercepts
## γ_00: General intercept
gamma_a<-0
gamma_b<-0
gamma_c<-0


## V_0j: Weight random intercept id
random_intercepts_id <- rnorm_multi(
  n = n_id, 
  mu = c(V_0j_a = 0, V_0j_b = 0, V_0j_c = 0),
  r = c( 1,   0.1, 0.1, 
         0.1, 1,    0.05, 
         0.1, 0.05,  1)
)*f_id_int

## U_i0: Weight random intercept stimulus
random_intercepts_stimulus <- rnorm_multi(
  n = n_stimulus, 
  mu = c(U_i0_a = 0, U_i0_b = 0, U_i0_c = 0),
  r = c( 1,   0.05, 0.05, 
         0.05, 1,    0.05, 
         0.05, 0.05,  1)
)*f_stimulus_int

## V_1j: Weight random slopes id
V_1j_a<-rnorm(n_id)*f_id_slope
V_1j_b<-rnorm(n_id)*f_id_slope
V_1j_c<-rnorm(n_id)*f_id_slope

## U_i1: Weight random slope stimulus
U_i1_a<-rnorm(n_stimulus)*f_stimulus_slope
U_i1_b<-rnorm(n_stimulus)*f_stimulus_slope
U_i1_c<-rnorm(n_stimulus)*f_stimulus_slope


# Create Ordered Identifier Variables
## i: id
## Order: 1,1,1,1, .... 30,30,30,30
id<-sort(rep(1:n_id,n_stimulus))


# j: stimulus
## Order: 1,2, .. 29,30,1,2, ... 29,30
stimulus <- rep(1:n_stimulus,n_id)

# treatment: Experimental factor
## Levels: 0 = control; 1 =  Experimental

### Urn 15 x 0 and 15 x 1
### Every participant will be exposed to an equal number of control and experimental stimuli 
urn<-rep(c(0,1),n_stimulus/2)

### Create experimental factor: Sample different combination for each id
treatment<-rep(sample(urn,replace=FALSE),n_id)

# beta: effect
# assuming effects of d=0.4 between all variables
beta_a<-d_to_r(-0.4)
beta_b<-d_to_r(-0.4)
beta_c<-d_to_r(0.4)


# Error: R
## Normally distributed random variable
error_a=rnorm(n)*f_residual
error_b=rnorm(n)*f_residual
error_c=rnorm(n)*f_residual

# FOR = gamma_a +  (beta_a+V_1j_a+U_i1_a)*treatment + V_0j_a + Um_i0_a + error_a      ## a part
FOR<-c()
for(i in 1:n_id){
  for (j in 1:n_stimulus){
    FOR[j+n_stimulus*(i-1)]<-gamma_a+(beta_a+V_1j_a[i]+U_i1_a[j])*treatment[j+n_stimulus*(i-1)]+random_intercepts_id$V_0j_a[i]+random_intercepts_stimulus$U_i0_a[j]+error_a[j+n_stimulus*(i-1)]
  }
}

#     CR_IR[j+n_stimulus*(i-1)]<-gamma_b+(beta_b+V_1j_b+U_i1_b)*FOR+V_0j_b+U_i0_b+error_b+
#                               gamma_c+(beta_c+V_1j_c+U_i1_c)*treatment+V_0j_c+error_c
CR_IR<-c()
for(i in 1:n_id){
  for (j in 1:n_stimulus){
    CR_IR[j+n_stimulus*(i-1)]<-gamma_b+(beta_b+V_1j_b[i]+U_i1_b[j])*FOR[j+n_stimulus*(i-1)]+random_intercepts_id$V_0j_b[i]+random_intercepts_stimulus$U_i0_b[j]+error_b[j+n_stimulus*(i-1)]+
      gamma_c+(beta_c+V_1j_c[i]+U_i1_c[j])*treatment[j+n_stimulus*(i-1)]+random_intercepts_id$V_0j_c[i]+random_intercepts_stimulus$U_i0_c[j]+error_c[j+n_stimulus*(i-1)]
  }
}

## combining the dataframe
df<-data.frame(id=id,
               stimulus=stimulus,
               FOR=FOR,
               treatment=treatment,
               CR_IR=CR_IR)

# Centering: No Centering is conducted, since FOR is also the dependent variable
# in model_X.2.1_med, which can result in singular fits. "treatment" is binary, so
# also no centering is conducted.

###############
# Main Models #
###############
model_X.2.1_med <- lmer(FOR ~ treatment + (treatment|id), data = df)
model_X.2.1_out <- lmer(CR_IR ~ FOR + treatment + (treatment|id), data = df)
model_X.2.1 <- mediation::mediate(model_X.2.1_med, model_X.2.1_out, treat = "treatment", mediator = "FOR", sims = 1000)
summary(model_X.2.1)


#########################
# Exploratory Variables #
#########################
## recode treatment for interaction effects
df$treatment<-recode(df$treatment, `0` = -1L, `1` = 1L, )
df$treatment[1:20]

## Familiarity with Publication Bias
df$familiarity_PB<-rnorm(n)+df$FOR*0.3*runif(n)+df$CR_IR*0.4*runif(n)

## Pressure to Publish
df$pub_pressure<-rnorm(n)+df$FOR*0.2*runif(n)+df$CR_IR*0.1*runif(n)

# Professional Status
df$professor<-sample(0:1,n,replace = TRUE,prob=c(0.9,0.1))
df$postdoc<-rep(0,n)

for (i in 1:n){
  sample<-sample(0:1,1,prob=c(0.55,0.45))
  if(df$professor[i]==0){
    df$postdoc[i]<-sample
  }
}

#position
df$position<-rep(sample(1:n_stimulus,n_stimulus),n_id)

#Conscientiousness
conscientiousness<-rnorm(n)+df$CR_IR*-0.3*runif(n)

######################
# Exploratory Models #
######################
model_X.2.2_med <- lmer(FOR ~ treatment * familiarity_PB + 
                          treatment * pub_pressure + (treatment|id), data = df)
model_X.2.2_out <- lmer(CR_IR ~ treatment * familiarity_PB +
                          treatment * pub_pressure + FOR + conscientiousness +
                          position + (treatment|id), data = df)
model_X.2.2 <- mediation::mediate(model_X.2.2_med, model_X.2.2_out,
                                  treat = "treatment", mediator = "FOR", 
                                  sims = 1000)
summary(model_X.2.2)