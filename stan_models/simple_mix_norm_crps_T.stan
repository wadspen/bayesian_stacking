
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

parameters {
  simplex[num_comp] omega;
}

transformed parameters {
  
  real<lower=0> risk_crps = riskCRPS(num_comp, T, betas, alphas,
                                     omega, tweight);
  
}


model {
  omega ~ dirichlet(alpha);
  target += exponential_lupdf(risk_crps | eta*T);
}



