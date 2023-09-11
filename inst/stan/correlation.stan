#include /pre/license.stan
functions {
#include /functions/pairwise.stan
}
data {
  real alphaScalePrior;
  // dimensions
  int<lower=1> NPA;             // worths or players or objects or things
  int<lower=1> NCMP;            // unique comparisons
  int<lower=1> N;               // observations
  int<lower=1> numRefresh;      // when change in item/pa1/pa2
  int<lower=1> NITEMS;
  array[NITEMS] int<lower=1> NTHRESH;         // number of thresholds
  array[NITEMS] int<lower=1> TOFFSET;
  vector[NITEMS] scale;
  real corLKJPrior;
  // response data
  array[numRefresh] int<lower=1, upper=NPA> pa1;
  array[numRefresh] int<lower=1, upper=NPA> pa2;
  array[NCMP] int weight;
  array[NCMP] int pick;
  array[numRefresh] int refresh;
  array[numRefresh] int numOutcome;
  array[numRefresh] int item;
}
transformed data {
  int totalThresholds = sum(NTHRESH);
  array[NCMP] int rcat;
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
  vector<lower=0,upper=1>[totalThresholds] rawThreshold;
  vector<lower=0>[NITEMS] alpha;
  matrix[NPA,NITEMS]      rawTheta;
  cholesky_factor_corr[NITEMS] rawThetaCorChol;
}
transformed parameters {
  vector[totalThresholds] threshold;
  vector[totalThresholds] rawCumTh;
  matrix[NPA,NITEMS]      theta;

  // non-centered parameterization due to thin data
  for (pa in 1:NPA) {
    theta[pa,] = (rawThetaCorChol * rawTheta[pa,]')';
  }
  for (ix in 1:NITEMS) {
    real maxSpan = max(theta[,ix]) - min(theta[,ix]);
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    threshold[from:to] = maxSpan * rawThreshold[from:to];
    rawCumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  rawThetaCorChol ~ lkj_corr_cholesky(corLKJPrior);
  for (pa in 1:NPA) {
    rawTheta[pa,] ~ std_normal();
  }
  rawThreshold ~ beta(1.1, 2);
  for (ix in 1:NITEMS) alpha[ix] ~ normal(1.749, alphaScalePrior) T[0,];
  {
    int cmpStart = 1;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      target += pairwise_logprob(rcat, weight, cmpStart, refresh[rx],
                                 scale[ix], alpha[ix], theta[pa1[rx], ix],
                                 theta[pa2[rx], ix], rawCumTh[from:to]);
      cmpStart += refresh[rx];
    }
  }
}
generated quantities {
  corr_matrix[NITEMS] thetaCor;
  thetaCor = multiply_lower_tri_self_transpose(rawThetaCorChol);
}
