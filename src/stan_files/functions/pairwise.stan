vector cmp_probs(real scale, real alpha, real pa1, real pa2, vector thr) {
  int nth = num_elements(thr);
  vector[1+nth*2] unsummed;
  real paDiff = scale * (pa2 - pa1);
  unsummed[1] = 0;
  for (tx in 1:nth) {
    real t1 = thr[nth+1-tx];
    unsummed[1+tx] = paDiff + t1;
    unsummed[1+2*nth+1-tx] = paDiff - t1;
  }
  return softmax(cumulative_sum(alpha * unsummed));
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
