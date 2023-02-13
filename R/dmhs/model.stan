functions {
  // for likelihood estimation
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + lgamma(sum(y)+1) + sum(lgamma(alpha + to_vector(y)))
      - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha))-sum(lgamma(to_vector(y)+1));
    }
  }

data {
  int<lower=1> N; // total number of observations
  int<lower=1> M; // total number of observations in the prediction set
  int<lower=2> J; // number of categories
  int<lower=2> P; // number of predictor levels
  matrix[N,P] X; // predictor design matrix
  matrix[M,P] Xp; // predictor design matrix
  int <lower=0> Y[N,J]; // data // response variable
  //real<lower=0> sd_prior;
  real<lower=0> sd_prior;
  real<lower=0> psi;
}

parameters {
  matrix[P, J] beta_raw; // coefficients (raw)
  vector[J] beta0; // intercept
  matrix<lower=0>[P,J] lambda_tilde; // truncated local shrinkage
  vector<lower=0>[J] tau; // global shrinkage
}

transformed parameters{
  matrix[P,J] beta; // coefficients
  matrix<lower=0>[P,J] lambda; // local shrinkage
  lambda = diag_post_multiply(lambda_tilde, tau);
  beta = beta_raw .* lambda;
}

model {
  // prior:
  //for(k in 1:J){
  //  beta0[k] ~ normal(0, 10);
  //}
  beta0 ~ normal(0, 10);

  for (k in 1:P) {
    //for (l in 1:J) {
      //tau[l] ~ cauchy(0.1, 1); // flexible
      lambda_tilde[k,] ~ cauchy(0, 1);
      beta_raw[k,] ~ normal(0,sd_prior);
    //}
  }
  tau ~ cauchy(0.1, 1);

  for (i in 1:N) {
    vector[J] logits;
    for (j in 1:J){
      logits[j] = beta0[j]+X[i,] * beta[,j];
    }
    Y[i,] ~ dirichlet_multinomial(softmax(logits)*(1-psi)/psi);
  }
}

generated quantities {
  matrix[M, J] pipred;
  for (i in 1:M) {
    vector[J] logits;
    for (j in 1:J){
      logits[j] = beta0[j]+Xp[i,] * beta[,j];
    }
    pipred[i,] = transpose(softmax(logits));
  }
}
