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
  
 
  
    vector mixnormCRPS(int num_comp, int N, matrix betas, matrix alphas,
                       vector omega) {
    
      vector[N] crps;
      
      for (n in 1:N) {
        crps[n] = dot_product(omega, betas[n,]) + quad_form(alphas, omega);
      }
    
    return crps;
  }
  
  real riskCRPS(int num_comp, int N, matrix betas, matrix alphas,
                       vector omega) {
    real rcrps = mean(mixnormCRPS(num_comp, N, betas, alphas,
                                  omega));
    return rcrps;
  }
  
  
  
  
}

  
data {
  int<lower=0> N;
  vector[N] y;
  int<lower=0> num_comp;
  vector[num_comp] mu;
  vector<lower=0>[num_comp] sigma;
  real<lower=0> eta;
  vector<lower=0>[num_comp] alpha;
  matrix[num_comp, num_comp] alphas;
  matrix[N, num_comp] betas;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  simplex[num_comp] omega;
  // vector[num_comp] omega;
}

transformed parameters {
  // vector[num_comp] omegat = exp(omega)/sum(exp(omega)); 
  real<lower=0> risk_crps = riskCRPS(num_comp, N, betas, alphas,
                                     omega);
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  omega ~ dirichlet(alpha);
  // omega ~ normal(1, 1);
  
  // target += gibbsLik(risk_crps, eta, N);
  // target += normal_lpdf(risk_crps | 0, 1);
  target += exponential_lupdf(risk_crps | eta*N);
}

