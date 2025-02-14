
// This Stan program defines a dual model for ILI data and for 
// hospitalizations. The ILI data is modeled according to an SIR
// compartmental model where ILI is treated as the observed 
// infected. The infected portion is the expected value of a beta
// distribution. The mean and scale parameters are independent 
// across seasons.
// Flu hospitalizations are then modeled as a linear function of 
// the ILI model, the squared ILI model, and one lagged 
// hospitalization. The parameters in the linear model between 
// the two seasons, 2022 and 2023, are independent.

functions {
  
  
  vector sir(real t, 
             vector y, 
             real beta,
             real gamma) {
    
    real S = y[1];
    real I = y[2];
    real R = y[3];
    
    vector[3] dydt;
    
    dydt[1] = -beta * I * S;
    dydt[2] =  beta * I * S - gamma * I;
    dydt[3] =  gamma * I;
    
    return dydt;
  }
}
data {
  int<lower=1> n_weeks;
  array[n_weeks] int<lower=1> weeks;
  real S0;
  real t0;
  array[n_weeks] real<lower=0> ts;
  array[n_weeks] real<lower=0> ili;
  
  //prior hyperparameters
  real<lower=0> sigma_kappa;
  real rho_mu; // .68
  real<lower=0> rho_sigma; //.08
  real beta_mu; //.8
  real<lower=0> beta_sigma; //.3
  real I0_mu; //.005
  real I0_sigma;
  
  
  
}

parameters {
  real<lower=0,upper=S0> rho;
  real<lower=0> beta;
  real<lower=0,upper=1 - S0> I0;
  real<lower=0> kappa;
}

transformed parameters{
  vector[3] y0;
  array[n_weeks] vector[3] y;
  array[n_weeks] real IS;
  real<lower=0> R0;
  real<lower=0> gamma;
  //array[HM] real ili_ps;
  
  R0 = 1 - S0 - I0;
  gamma = rho*beta;
  
    {
      y0[1] = S0;
      y0[2] = I0;
      y0[3] = R0;
      
      
      y[1:n_weeks,] =
        ode_rk45(sir, y0, t0,
                 segment(ts, 1, n_weeks),
                 beta, gamma);
      
      IS = y[1:n_weeks, 2];
    }
  
}

model {
  rho ~ normal(rho_mu, rho_sigma);
  beta ~ normal(beta_mu, beta_sigma);
  I0 ~ normal(I0_mu, I0_sigma);
  kappa ~ normal(0, sigma_kappa);
  
  for (i in 1:n_weeks) {
        
    ili[i] ~ beta_proportion(
      (IS[i]+.00001), 
      kappa);
      
  }
  
}
generated quantities {

  array[n_weeks + 4] real pred_ili;
  array[n_weeks + 4] real tp;
  array[n_weeks + 4] vector[3] ypm;
  array[n_weeks + 4] real IP;
  vector[3] yp;
  yp[1] = S0;
  yp[2] = I0;
  yp[3] = 1 - S0 - I0;
  for (i in 1:(n_weeks + 4)) tp[i] = i;
  ypm = ode_rk45(sir, yp, t0, tp, beta, gamma);
  IP = ypm[,2];


    for (i in 1:(n_weeks + 4)) {
      pred_ili[i] = beta_proportion_rng(
        (IP[i]+.00001),
        kappa);
    }
}



