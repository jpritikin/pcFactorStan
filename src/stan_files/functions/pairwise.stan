vector cmp_probs(real scale, real alpha, real pa1, real pa2, vector thr, int[] want) {
  int nth = num_elements(thr);
  real at[nth*2];
  real pr[1+nth*2];
  vector[1 + nth*2] out;
  real paDiff = alpha * scale * (pa1 - pa2);
  vector[nth] thrAlpha = thr * alpha;
  for (tx in 1:num_elements(at)) {
    if (tx <= nth) {
      at[tx] = paDiff - thrAlpha[nth+1-tx];
    } else {
      at[tx] = paDiff + thrAlpha[tx-nth];
    }
  }
  pr[1+nth*2] = 1;
  for (tx in 1:num_elements(at)) {
    if (want[tx] || want[tx+1]) {
      pr[tx] = 1.0/(1.0+exp(-at[tx]));
    } else {
      pr[tx] = 0; // not used
    }
  }
  out[1] = pr[1];
  for (tx in 2:num_elements(out)) {
    out[tx] = pr[tx] - pr[tx-1];
  }
  return out;
}

real pairwise_logprob(int[] rcat, int[] weight, int cmpStart, int len,
                      real scale, real alpha, real pa1, real pa2, vector cumTh)
{
  real lp = 0;
  int nth = num_elements(cumTh);
  vector[1+nth*2] prob;
  int want[1+nth*2];
  for (ox in 1:num_elements(want)) want[ox] = 0;
  for (ox in cmpStart:(cmpStart + len - 1)) {
    want[ rcat[ox] ] = 1;
  }
  prob = cmp_probs(scale, alpha, pa1, pa2, cumTh, want);
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
  int want[1+nth*2];
  for (ox in 1:num_elements(want)) want[ox] = 0;
  for (ox in cmpStart:(cmpStart + len - 1)) {
    want[ rcat[ox] ] = 1;
  }
  prob = cmp_probs(scale, alpha, pa1, pa2, cumTh, want);
  for (ox in cmpStart:(cmpStart + len - 1)) {
    real lp1 = log(prob[rcat[ox]]);
    for (wx in 1:weight[ox]) {
      lp[cur] = lp1;
      cur += 1;
    }
  }
  return lp;
}
