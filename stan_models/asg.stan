
// In this STAN model, ILI data is modeled according to the assymmetric
// Gaussian function where there is a hierarchy of parameters over flu
// seasons. The ASG function is taken as the expected value of a beta
// distribution and the scale parameter of the beta distribution is 
// considered independent over seasons.
// A discrepancy reverse random walk is included to the mean to capture
// any pattern that is shared by all or most flu seasons. 
// Flu hospitalizations are then modeled as a linear function of the ILI
// model, the squared ILI model, and one lagged hospitalization. The 
// parameters in the linear model between the two seasons, 2022 and 2023,
// are independent.


functions {

   real asg(row_vector theta, real t) {
    
    real beta;
    real eta;
    real mu; 
    real sig1;
    real sig2;
    real ASG;
    beta = theta[1];
    eta = exp(theta[2]);
    mu = theta[3];
    sig1 = exp(theta[4]);
    sig2 = exp(theta[5]);
    
    ASG = ((beta + (eta)*exp(-((t-mu)^2)/(2*sig1^2)))*(t < mu) +
       (beta + (eta)*exp(-((t-mu)^2)/(2*sig2^2)))*(mu <= t));
    return ASG;
  }
}

data {
int<lower=0> n_params;
int<lower=0> n_weeks;                                 // weeks in a season
array[n_weeks] real<lower=0, upper=1> ili;
array[n_weeks] int<lower=0> weeks;                         // number of parameters ASG
vector[n_params] m0;                            // prior means
matrix[n_params,n_params] C0;                   // prior sds

real<lower=0> sigma_kappa; //1000

}
parameters {
  
row_vector[n_params] theta;
real<lower=0> kappa;

}

model {

theta ~ multi_normal(m0,C0);
kappa ~ normal(0,sigma_kappa);
for (i in 1:n_weeks) ili[i] ~ beta_proportion(
                                    inv_logit(asg(theta,
                                    i)),
                                    kappa); 
}

generated quantities {
    array[n_weeks + 4] real<lower=0,upper=1> pred_ili;

    for (i in 1:(n_weeks + 4)) {
        pred_ili[i] = beta_proportion_rng(
                              inv_logit(asg(theta,
                              i)),
                              kappa);
    }
}
