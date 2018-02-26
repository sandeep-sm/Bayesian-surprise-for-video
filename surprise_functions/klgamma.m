function d = klgamma(a1,a2,b1,b2)

d = -a2 + a2.*log(b1./b2) + log(gamma(a2)./gamma(a1)) + b2.*(a1./b1) + (a1 - a2).*digamma(a1,1000);