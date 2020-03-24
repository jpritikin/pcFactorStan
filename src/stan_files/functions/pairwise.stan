vector cmp_probs(real scale, real alpha, real pa1, real pa2, vector thr) {
  int nth = num_elements(thr);
  real at[nth*2];
  real pr[2+nth*2];
  vector[1 + nth*2] out;
  real paDiff = scale * (pa1 - pa2);
  for (tx in 1:num_elements(at)) {
    if (tx <= nth) {
      at[tx] = paDiff - thr[nth+1-tx];
    } else {
      at[tx] = paDiff + thr[tx-nth];
    }
  }
  pr[1] = 0;
  pr[2+nth*2] = 1;
  for (tx in 1:num_elements(at)) {
    // TODO hoist mult by alpha
    pr[1+tx] = 1.0/(1.0+exp(-at[tx]*alpha));
  }
  for (tx in 1:num_elements(out)) {
    out[tx] = pr[tx+1] - pr[tx];
  }
  return out;
}

real pairwise_logprob(int[] rcat, int[] weight, int cmpStart, int len,
                      real scale, real alpha, real pa1, real pa2, vector cumTh)
{
  real lp = 0;
  int nth = num_elements(cumTh);
  vector[1+nth*2] prob;
  prob = cmp_probs(scale, alpha, pa1, pa2, cumTh);
  for (ox in cmpStart:(cmpStart + len - 1)) {
    if (weight[ox] == 1) {
      lp += log(prob[rcat[ox]]);
    } else {
      lp += weight[ox] * log(prob[rcat[ox]]);
    }
  }
  return lp;
}

real[] pairwise_loo(int[] rcat, int[] weight, int numOutcome, int cmpStart, int len,
                  real scale, real alpha, real pa1, real pa2, vector cumTh)
{
  real lp[numOutcome];
  int cur = 1;
  int nth = num_elements(cumTh);
  vector[1+nth*2] prob;
  prob = cmp_probs(scale, alpha, pa1, pa2, cumTh);
  for (ox in cmpStart:(cmpStart + len - 1)) {
    real lp1 = log(prob[rcat[ox]]);
    for (wx in 1:weight[ox]) {
      lp[cur] = lp1;
      cur += 1;
    }
  }
  return lp;
}
