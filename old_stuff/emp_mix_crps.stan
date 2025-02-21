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
  
 
  
    vector empCRPSs(int N, int num_comp, vector omega, data matrix mae, 
                       data array[] matrix absdiff) {
    
      vector[N] crpss;
      
      for (n in 1:N) {
        crpss[n] = dot_product(omega, mae[,n]) 
                   + quad_form(absdiff[n], omega);
      }
    
    return crpss;
  }
  
  real riskCRPSs(int N, int num_comp, vector omega, data matrix mae, 
                       data array[] matrix absdiff) {
    real rcrpss = mean(empCRPSs(N, num_comp, omega, mae, absdiff));
    return rcrpss;
  }
  
  
  
  
}

  
data {
  int<lower=0> N;
  int<lower=0> num_comp;
  real<lower=0> eta;
  vector<lower=0>[num_comp] alpha;
  matrix[num_comp,N] mae;
  array[N] matrix[num_comp, num_comp] absdiff;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  simplex[num_comp] omega;
  // vector[num_comp] omega;
}

transformed parameters {
  // vector[num_comp] omegat = exp(omega)/sum(exp(omega)); 
  real<lower=0> risk_crps = riskCRPSs(N, num_comp, omega, mae, absdiff);
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  omega ~ dirichlet(alpha);
  target += exponential_lupdf(risk_crps | eta*N);
}

