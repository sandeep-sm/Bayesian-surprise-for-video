

function [b1,b2] = betavalues(d,u)

% Check bounds on decay, it must be less than 1 and greater than 0
if d <= 0 || d >= 1
    error('D must be bounded such that 0 < D < 1\n');
end

% compute beta and beta' at their asymptotic values
b2 = -u ./ (d - 1);
b1 = b2;
