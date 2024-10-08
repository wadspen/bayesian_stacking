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
  
 
  
    real mixnormCRPS(int num_comp, int T, row_vector betas, matrix alphas,
                       vector omega) {
    
      real crps;
      crps = dot_product(omega, betas) + quad_form(alphas, omega);
    
    return crps;
  }
  
  real riskCRPS(int num_comp, int T, matrix betas, array[] matrix alphas,
                       array[] vector omegas, vector wts, real lambda) {
    vector[T] crpss;
    vector[T] reg;
    real rcrps;
    for (t in 1:T) {
      crpss[t] = mixnormCRPS(num_comp, T, betas[t,], alphas[t], omegas[t]);
      reg[t] = distance(omegas[t], wts);
    }
    rcrps = mean(crpss) + lambda*mean(reg);
    return rcrps;
  }
  
  
  
  
}

  
data {
  int<lower=0> T;
  vector[T] y;
  int<lower=0> num_comp;
  real<lower=0> eta;
  array[T] matrix[num_comp, num_comp] alphas;
  matrix[T, num_comp] betas;
  vector[num_comp - 1] sigma_s;
  vector[num_comp] wts;
  real<lower=0> lambda;
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
  // real<lower=0> lambda;
}

transformed parameters {
  // vector[num_comp] omegat = exp(omega)/sum(exp(omega)); 
  
  array[T] simplex[num_comp] omegas;
  for (t in 1:T) {
    omegas[t][1:(num_comp - 1)] = exp(zetas[t])/(1 + sum(exp(zetas[t])));
    omegas[t][num_comp] = 1 - sum(omegas[t][1:(num_comp - 1)]);
  }
  
  
  
  // array[T] simplex[num_comp] omegas;
  // for (t in 1:T) omegas[t] = softmax(zetas[t]);
  
  real<lower=0> risk_crps = riskCRPS(num_comp, T, betas, alphas,
                                     omegas, wts, lambda);
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // omega ~ dirichlet(alpha);
  // omega ~ normal(1, 1);
  // lambda ~ normal(0,3);
  // sigma_z ~ normal(0, 1);
  beta ~ normal(0, 1);
  zetas[1] ~ multi_normal(log(sigma_s) + 7,
  diag_matrix(sigma_s));
  for (t in 2:T) {zetas[t] ~ multi_normal(beta + zetas[t-1], diag_matrix(sigma_z));}
  
  sigma_z ~ normal(0, 10);
  // beta ~ normal(0, 1);
  // zetas[1] ~ normal(0, sigma_z);
  // for (t in 2:T) {zetas[t] ~ normal(beta*zetas[t-1], sigma_z);}
  
  // target += gibbsLik(risk_crps, eta, N);
  // target += normal_lpdf(risk_crps | 0, 1);
  target += exponential_lupdf(risk_crps | eta*T);
}


generated quantities {
  vector[num_comp - 1] zetaT1;
  // for (c in 1:num_comp) zetaT1[c] = normal_rng(zetas[T][c], sigma_z);
  
  
  
  // vector[num_comp] zetaT1;
  // for (c in 1:num_comp) zetaT1[c] = normal_rng(beta*zetas[T][c], sigma_z);
  // simplex[num_comp] omegaT1 = softmax(zetaT1);
  
  zetaT1 = multi_normal_rng(beta + zetas[T], diag_matrix(sigma_z));
  simplex[num_comp] omegaT1;
  omegaT1[1:(num_comp -1)] = exp(zetaT1)/(1 + sum(exp(zetaT1)));
  omegaT1[num_comp] = 1 - sum(omegaT1[1:(num_comp - 1)]);
}

