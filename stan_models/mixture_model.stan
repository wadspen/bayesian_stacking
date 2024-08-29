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



  
data {
  int<lower=0> T;
  vector[T] y;
  int<lower=0> num_comp;
  vector[num_comp - 1] sigma_s;
  vector[num_comp] mus;
  vector<lower=0>[num_comp] sigmas;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // simplex[num_comp] omega;
  // vector[num_comp] omega;
  vector<lower=0>[num_comp - 1] sigma_z;
  // real<lower=0> sigma_z;
  array[T] vector[num_comp - 1] zetas;
  // real<lower=0,upper=1> beta;
  real beta;
}

transformed parameters {
  // vector[num_comp] omegat = exp(omega)/sum(exp(omega)); 
  
  array[T] simplex[num_comp] omegas;
  for (t in 1:T) {
    omegas[t][1:(num_comp - 1)] = exp(zetas[t])/(1 + sum(exp(zetas[t])));
    omegas[t][num_comp] = 1 - sum(omegas[t][1:(num_comp - 1)]);
  }
  
  
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // omega ~ dirichlet(alpha);
  // omega ~ normal(1, 1);
  
  sigma_z ~ normal(0, 1);
  beta ~ normal(0, 1);
  zetas[1] ~ multi_normal(log(sigma_s) + 7,
  diag_matrix(sigma_s));
  
  target += log_sum_exp({log(omegas[1][1]) + normal_lpdf(y[1] | mus[1], sigmas[1]),
                          log(omegas[1][2]) + normal_lpdf(y[1] | mus[2], sigmas[2]),
                          log(omegas[1][3]) + normal_lpdf(y[1] | mus[3], sigmas[3]),
                          log(omegas[1][4]) + normal_lpdf(y[1] | mus[4], sigmas[4]),
                          log(omegas[1][5]) + normal_lpdf(y[1] | mus[5], sigmas[5]),
                          log(omegas[1][6]) + normal_lpdf(y[1] | mus[6], sigmas[6])});
                          
  for (t in 2:T) {
    zetas[t] ~ multi_normal(beta + zetas[t-1], diag_matrix(sigma_z));
    
    target += log_sum_exp({log(omegas[t][1]) + normal_lpdf(y[t] | mus[1], sigmas[1]),
                          log(omegas[t][2]) + normal_lpdf(y[t] | mus[2], sigmas[2]),
                          log(omegas[t][3]) + normal_lpdf(y[t] | mus[3], sigmas[3]),
                          log(omegas[t][4]) + normal_lpdf(y[t] | mus[4], sigmas[4]),
                          log(omegas[t][5]) + normal_lpdf(y[t] | mus[5], sigmas[5]),
                          log(omegas[t][6]) + normal_lpdf(y[t] | mus[6], sigmas[6])});
  }
  
  
  
}


generated quantities {
  vector[num_comp - 1] zetaT1;
  
  zetaT1 = multi_normal_rng(beta + zetas[T], diag_matrix(sigma_z));
  simplex[num_comp] omegaT1;
  omegaT1[1:(num_comp -1)] = exp(zetaT1)/(1 + sum(exp(zetaT1)));
  omegaT1[num_comp] = 1 - sum(omegaT1[1:(num_comp - 1)]);
}

