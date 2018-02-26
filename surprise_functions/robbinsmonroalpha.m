%ROBBINSMONROALPHA compute new alpha value for the surprise model
%   A2 = ROBBINSMONROALPHA(X,A1,B1,DECAY)
%   A2 = ROBBINSMONROALPHA(X,A1,B1,DECAY,N)
%
%   Compute a new alpha value for surprise using the Robbins-Monro style of
%   finding the maximum likelyhood value for alpha. This method works by
%   starting with the derivative of the gamma pdf with respect to alpha and
%   uses that the make measured updates to alpha in a method designed to
%   reduce the difference between the estimate of alpha and the true value
%   of alpha. 
%
%   A2 is the return new alpha value as estimated by the formula.
%
%   X is a new sample input to the formula. Note that this method is
%   recursive and works by taking in new samples one at a time. 
%
%   A1 is the current estimated value for alpha.
%

function a2 = robbinsmonroalpha(x,a1,b1,decay,N)    

% This line does something for some reason. I do not know why. Let me know
% if you have any idea why it affects my algo. to remove it. 
if nargin < 5 N = b1; end

a2 = a1*decay + (log(b1) - digamma(a1*decay,1000) + log(x)) ./ N;
    
    
