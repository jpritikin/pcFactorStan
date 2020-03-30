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
  real alphaShape;
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
  vector<lower=0,upper=1>[totalThresholds] rawThreshold;
  vector<lower=0>[NITEMS] alpha;
  matrix[NPA,NITEMS]      rawTheta;
  cholesky_factor_corr[NITEMS] rawThetaCorChol;
}
transformed parameters {
  vector[totalThresholds] threshold;
  vector[totalThresholds] cumTh;
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
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  rawThetaCorChol ~ lkj_corr_cholesky(2);
  for (pa in 1:NPA) {
    rawTheta[pa,] ~ std_normal();
  }
  rawThreshold ~ beta(1.1, 1.1);
  alpha ~ inv_gamma(alphaShape, 1.749*(1+alphaShape));
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
  real log_lik[N];

  corr_matrix[NITEMS] thetaCor;
  thetaCor = multiply_lower_tri_self_transpose(rawThetaCorChol);

  {
    int cmpStart = 1;
    int cur = 1;
    for (rx in 1:numRefresh) {
      int ix = item[rx];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      int nout = numOutcome[rx];
      int last = cur + nout - 1;
      log_lik[cur:last] =
        pairwise_loo(rcat, weight, nout, cmpStart, refresh[rx],
                     scale[ix], alpha[ix],
                     theta[pa1[rx],ix], theta[pa2[rx],ix], cumTh[from:to]);
      cmpStart += refresh[rx];
      cur += numOutcome[rx];
    }
  }
}
