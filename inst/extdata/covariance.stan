functions {
  vector cmp_probs(real scale, real pa1, real pa2, vector thr) {
    int nth = num_elements(thr);
    vector[1+nth*2] unsummed;
    real paDiff = scale * (pa1 - pa2);
    unsummed[1] = 0;
    for (tx in 1:nth) {
      real t1 = thr[nth+1-tx];
      unsummed[1+tx] = paDiff - t1;
      unsummed[1+2*nth+1-tx] = paDiff + t1;
    }
    return softmax(cumulative_sum(unsummed));
  }
}
data {
  // dimensions
  int<lower=1> NPA;             // number of players or objects or things
  int<lower=1> NCMP;            // number of unique comparisons
  int<lower=1> N;               // number of observations
  int<lower=1> NITEMS;
  int<lower=1> NTHRESH[NITEMS];         // number of thresholds
  int<lower=1> TOFFSET[NITEMS];
  real scale;
  // response data
  int<lower=1, upper=NPA> pa1[NCMP];        // PA1 for observation N
  int<lower=1, upper=NPA> pa2[NCMP];        // PA2 for observation N
  int weight[NCMP];
  int pick[NCMP];
  int refresh[NCMP];
  int item[NCMP];
}
transformed data {
  int totalThresholds = sum(NTHRESH);
  int rcat[NCMP];
  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH[item[cmp]] + 1;
  }
}
parameters {
  vector[totalThresholds] threshold;
  matrix[NPA,NITEMS]      rawTheta;
  vector<lower=0>[NITEMS] sigma;
  cholesky_factor_corr[NITEMS] rawThetaCorChol;
}
transformed parameters {
  vector[totalThresholds] cumTh;
  matrix[NPA,NITEMS]      theta;

  // non-centered parameterization due to thin data
  for (pa in 1:NPA) {
    theta[pa,] = (sigma .* (rawThetaCorChol * rawTheta[pa,]'))';
  }
  for (ix in 1:NITEMS) {
    int from = TOFFSET[ix];
    int to = TOFFSET[ix] + NTHRESH[ix] - 1;
    cumTh[from:to] = cumulative_sum(threshold[from:to]);
  }
}
model {
  vector[max(NTHRESH)*2 + 1] prob;
  int probSize;

  rawThetaCorChol ~ lkj_corr_cholesky(2);
  for (pa in 1:NPA) {
    rawTheta[pa,] ~ normal(0,1);
  }
  threshold ~ normal(0,2);
  sigma ~ lognormal(1, 1);
  for (cmp in 1:NCMP) {
    if (refresh[cmp]) {
      int ix = item[cmp];
      int from = TOFFSET[ix];
      int to = TOFFSET[ix] + NTHRESH[ix] - 1;
      probSize = (2*NTHRESH[ix]+1);
      prob[:probSize] = cmp_probs(scale,
               theta[pa1[cmp], ix],
               theta[pa2[cmp], ix], cumTh[from:to]);
    }
    if (weight[cmp] == 1) {
      target += categorical_lpmf(rcat[cmp] | prob[:probSize]);
    } else {
      target += weight[cmp] * categorical_lpmf(rcat[cmp] | prob[:probSize]);
    }
  }
}
generated quantities {
  corr_matrix[NITEMS] thetaCor;
  thetaCor = multiply_lower_tri_self_transpose(rawThetaCorChol);
}
