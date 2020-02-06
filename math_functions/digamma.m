function d = digamma(a,k)

if k < 1 || mod(k,1) ~= 0
    error('The value t in digamma must be a whole number greater than or equal to 1\n');
elseif a == 0
    error('The digamma is only defined for numbers where `a` is not 0\n');
end

emc = eulermasch(k);

if size(a,2) > 1 || size(a,1) > 1
    Z   = zeros(size(a,1),size(a,2));
    d   = Z - 1 * emc;
    ONE = ones(size(a,1),size(a,2));
    for i=1:k
        Ii = ONE * (1/i);
        I  = ONE * i;
        d  = d + Ii - (ONE./(a + I - ONE));
    end
else
    d = -1 * emc;
    for i=1:k
        d = d + (1/i) - (1/(a + i - 1));
    end
end

