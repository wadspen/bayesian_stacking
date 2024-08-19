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
                       vector omega, real tweight) {
    vector[T] crpss;
    real rcrps;
    for (t in 1:T) {
      crpss[t] = tweight^((T-t))*mixnormCRPS(num_comp, 
                                           T, betas[t,], alphas[t], omega);
    }
    rcrps = mean(crpss);
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
  vector<lower=0>[num_comp] alpha;
  real<lower=0,upper=1> tweight;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  simplex[num_comp] omega;
}

transformed parameters {
  // vector[num_comp] omegat = exp(omega)/sum(exp(omega)); 
  real<lower=0> risk_crps = riskCRPS(num_comp, T, betas, alphas,
                                     omega, tweight);
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // omega ~ dirichlet(alpha);
  // omega ~ normal(1, 1);
  omega ~ dirichlet(alpha);
  
  // target += gibbsLik(risk_crps, eta, N);
  // target += normal_lpdf(risk_crps | 0, 1);
  target += exponential_lupdf(risk_crps | eta*T);
}


generated quantities {
  simplex[num_comp] omegaT1 = dirichlet_rng(alpha);
}

