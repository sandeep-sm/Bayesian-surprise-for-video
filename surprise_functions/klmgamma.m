function d = klmgamma(a1,a2,b1,b2)

% Notice that the form of the gamma is very similar to the
% standard KL gamma expressed as:
% d = a2 * log(b1/b2) + log(gamma(a2)/gamma(a1)) + b2*(a1/b1) + (a1 - a2) * digamma(a1,1000);
N     = size(a1,2);
blog  = log(b1 ./ b2) / N;
broot = ((b2 ./ b1)^(1/N)) / N;

d = (1/N) * sum(a2 .* blog)      + ... 
    prod(gamma(a2) ./ gamma(a1)) + ...
    sum(a1 .* broot)             + ...
    sum(-1 * a1 + (a1 - a2) .* digamma(a1));






  
