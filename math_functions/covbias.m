function bias = covbias(alpha2,e_alpha2,J1)

bias = e_alpha2 - J1 * alpha2;