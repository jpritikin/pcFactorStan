functions {
  vector cmp_probs(real scale, real pa1, real pa2, vector thr) {
    int nth = num_elements(thr);
    vector[1+nth*2] unsummed;
    real paDiff = scale * (pa2 - pa1);
    unsummed[1] = 0;
    for (tx in 1:nth) {
      real t1 = thr[nth+1-tx];
      unsummed[1+tx] = paDiff + t1;
      unsummed[1+2*nth+1-tx] = paDiff - t1;
    }
    return softmax(cumulative_sum(unsummed));
  }
}
data {
  // dimensions
  int<lower=1> NPA;             // number of players or objects or things
  int<lower=1> NCMP;            // number of unique comparisons
  int<lower=1> N;               // number of observations
  int<lower=1> NTHRESH;         // number of thresholds
  real varCorrection;
  // response data
  int<lower=1, upper=NPA> pa1[NCMP];        // PA1 for observation N
  int<lower=1, upper=NPA> pa2[NCMP];        // PA2 for observation N
  int weight[NCMP];
  int pick[NCMP];
  int refresh[NCMP];
}
transformed data {
  int rcat[NCMP];

  for (cmp in 1:NCMP) {
    rcat[cmp] = pick[cmp] + NTHRESH + 1;
  }
}
parameters {
  vector[NPA] theta;
  vector[NTHRESH] threshold;
  real<lower=0> sigma;
}
transformed parameters {
  vector[NTHRESH] cumTh = cumulative_sum(threshold);
  real scale = (sigma * sigma) ^ varCorrection;

}
model {
  vector[NTHRESH*2 + 1] prob;
  sigma ~ lognormal(1, 1);
  theta ~ normal(0, sigma);
  threshold ~ normal(0, 2.0);
  for (cmp in 1:NCMP) {
    if (refresh[cmp]) {
      prob = cmp_probs(scale, theta[pa1[cmp]], theta[pa2[cmp]], cumTh);
    }
    if (weight[cmp] == 1) {
      target += categorical_lpmf(rcat[cmp] | prob);
    } else {
      target += weight[cmp] * categorical_lpmf(rcat[cmp] | prob);
    }
  }
}
