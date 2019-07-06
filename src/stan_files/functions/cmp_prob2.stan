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
