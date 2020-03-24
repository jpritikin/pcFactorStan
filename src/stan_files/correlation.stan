#include /pre/license.stan
functions {
#include /functions/pairwise.stan
}
data {
  // dimensions
  int<lower=1> NPA;             // worths or players or objects or things
  int<lower=1> NCMP;            // unique comparisons
  int<lower=1> N;               // observations
  int<lower=1> numRefresh;      // when change in item/pa1/pa2
  int<lower=1> NITEMS;
  int<lower=1> NTHRESH[NITEMS];         // number of thresholds
  int<lower=1> TOFFSET[NITEMS];
  vector[NITEMS] scale;
  // response data
  int<lower=1, upper=NPA> pa1[numRefresh];
  int<lower=1, upper=NPA> pa2[numRefresh];
  int weight[NCMP];
  int pick[NCMP];
  int refresh[numRefresh];
  int numOutcome[numRefresh];
  int item[numRefresh];
}
transformed data {
  int totalThresholds = sum(NTHRESH);
  int rcat[NCMP];
  {
    int cmpStart = 0;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      for (cmp in 1:refresh[rx]) {
        rcat[cmpStart + cmp] = pick[cmpStart + cmp] + NTHRESH[ix] + 1;
      }
      cmpStart += refresh[rx];
    }
  }
}
parameters {
  vector[totalThresholds] threshold;
  vector<lower=0>[NITEMS] alpha;
  matrix[NPA,NITEMS]      rawTheta;
  cholesky_factor_corr[NITEMS] rawThetaCorChol;
}
transformed parameters {
  vector[totalThresholds] cumTh;
  matrix[NPA,NITEMS]      theta;

  // non-centered parameterization due to thin data
  for (pa in 1:NPA) {
    theta[pa,] = (rawThetaCorChol * rawTheta[pa,]')';
  }
  for (ix in 1:NITEMS) {
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  rawThetaCorChol ~ lkj_corr_cholesky(2);
  for (pa in 1:NPA) {
    rawTheta[pa,] ~ std_normal();
  }
  threshold ~ normal(0, 2.0);
  alpha ~ exponential(0.1);
  {
    int cmpStart = 1;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      target += pairwise_logprob(rcat, weight, cmpStart, refresh[rx],
                                 scale[ix], alpha[ix], theta[pa1[rx], ix],
                                 theta[pa2[rx], ix], cumTh[from:to]);
      cmpStart += refresh[rx];
    }
  }
}
generated quantities {
  corr_matrix[NITEMS] thetaCor;
  thetaCor = multiply_lower_tri_self_transpose(rawThetaCorChol);
}
