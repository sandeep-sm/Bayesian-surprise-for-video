

function gamma = eulermasch(t)

if t < 1 || mod(t,1) ~= 0
    error('The value t in Euler-Mascheroni must be a whole number greater than or equal to 1\n');
end

gamma = 0;
for i=1:t
    gamma = gamma + (1/i) - log(1 + (1/i));
end