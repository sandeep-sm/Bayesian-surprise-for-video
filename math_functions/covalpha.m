function J = covalpha(alpha,e_alpha,mode)

if mode == 0
    J = e_alpha * inv(alpha);
elseif mode == 1
    [Q,R] = qr(alpha,0);    
    J = R\(Q' * e_alpha); % Least Sqaures Approximation. 
    %NOTE: We can also use J = inv(R)*Q' * e_alpha which is slower but more
    %accurate
else
    error('Unknown mode given for computing J passed as an option. Only 0 and 1 supported');
end


