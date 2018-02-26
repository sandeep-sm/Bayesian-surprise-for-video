

function d = trigamma(a,k)

if k < 1 || mod(k,1) ~= 0
    error('The value t in trigamma must be a whole number greater than or equal to 1\n');
elseif a == 0
    error('The trigamma is only defined for numbers where `a` is not 0\n');
end

if size(a,2) > 1 | size(a,1) > 1
    d   = zeros(size(a,1),size(a,2));
else
    d = 0;
end

for i=1:k
    d = d + 1/((a + i)^2);
end