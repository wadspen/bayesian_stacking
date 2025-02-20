//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.

functions {
  
 
  
    real empCRPSs(int N, int num_comp, vector omega, data matrix mae, 
                       data array[] matrix absdiff) {
    
      // vector[T] crpss;
      real crps = dot_product(omega, mae[,N]) + quad_form(absdiff[N], omega);
      
    
    return crps;
  }
  
  real riskCRPSs(int T, int num_comp, array[] vector omegas, data matrix mae, 
                       data array[] matrix absdiff) {
    // real rcrpss = mean(empCRPSs(T, num_comp, omega, mae, absdiff));
    real rcrpss;
    vector[T] crpss;
    
    for (t in 1:T) crpss[t] = empCRPSs(t, num_comp, omegas[t], mae, absdiff);
    rcrpss = mean(crpss);
    
    return rcrpss;
  }
  
  
  
  
}

  
data {
  int<lower=0> T;
  int<lower=0> num_comp;
  real<lower=0> eta;
  vector<lower=0>[num_comp] alpha;
  matrix[num_comp,T] mae;
  array[T] matrix[num_comp, num_comp] absdiff;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // simplex[num_comp] omega;
  // vector[num_comp] omega;
  real<lower=0> sigma_z;
  array[T] vector[num_comp] zetas;
  real<lower=0,upper=1> beta;
}

transformed parameters {
  // vector[num_comp] omegat = exp(omega)/sum(exp(omega));
  array[T] simplex[num_comp] omegas;
  for (t in 1:T) omegas[t] = softmax(zetas[t]);
  real<lower=0> risk_crps = riskCRPSs(T, num_comp, omegas, mae, absdiff);
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // omega ~ dirichlet(alpha);
  sigma_z ~ normal(0, 10);
  beta ~ normal(0, 1);
  zetas[1] ~ normal(0, sigma_z);
  for (t in 2:T) {zetas[t] ~ normal(beta*zetas[t-1], sigma_z);}
  
  target += exponential_lupdf(risk_crps | eta*T);
}

generated quantities {
  vector[num_comp] zetaT1;
  for (c in 1:num_comp) zetaT1[c] = normal_rng(beta*zetas[T][c], sigma_z);
  simplex[num_comp] omegaT1 = softmax(zetaT1);
}

