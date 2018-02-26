
function smod = newsm(decay,inita,initb,x,initlamb,updatefac,dim)

% depending on the number of input arguments, set some values to their
% default values. Otherwise accept the input values
if nargin < 7 dim       = 1; end
if nargin < 5 updatefac = 1; end
if nargin < 5 initlamb  = 0; end
if nargin < 4 x         = linspace(-10,10,101); end
if nargin < 3 initb     = 1; end
if nargin < 2 inita     = 1; end
if nargin < 1 
    error('DECAY value must be specified in calling newsm');
end

% check bounds on decay rate
if decay <= 0 || decay >= 1
    error('Error in DECAY value. It must be set 0 < DECAY < 1');
end

% Set up structures and containers for surprise parameters and values. 
smod = struct('Description','Surprise Model Container');
smod.decay     = trimf(x,[decay decay decay]);
smod.updatefac = trimf(x,[updatefac updatefac updatefac]); 
smod.smodisset = 1;
smod.options   = struct('debug',1,'graph','surprise','setbetamax','no','factordecay','yes','robbins_monro','no','robbins_monro_2','no',...
                            'robbins_monro_3','no');
smod.options.jointmodel = 'none'; % Other options include 'linear' and 'overdet'
% Use standard univariate model
if dim == 1 
    A = trimf(x,[1 1 1]);
    [d1, ix] = min(abs(x-inita));
    [d2, ix] = min(abs(x-initb));
    smod.alpha0    = A;
    smod.alpha1    = A;
    smod.alpha2    = trimf(x,[inita inita inita]);
    smod.beta1     = A;
    smod.beta2     = trimf(x,[initb initb initb]);
    smod.xbar1     = A;
    smod.xbar2     = A;
    smod.surprise  = 0;
    smod.epoch     = 0;
    smod.data0     = 1;
    smod.dim       = 1;
    smod.temp_prod = A;
    smod.max       = struct('Description','Maximum and upper bound limits on model');
    % Obtain asymptotic maximum values for beta and beta'
    % [smod.max.beta1,smod.max.beta2] = betavalues(decay,updatefac);
    
% Use experimental multi-variate model
else
    mat            = ones(1,dim);
    mat2           = ones(dim,dim);
    smod.alpha0    = mat;
    smod.alpha1    = mat;
    smod.alpha2    = mat;
    smod.e_alpha1  = mat;
    smod.e_alpha2  = mat;
    smod.m_alpha1  = mat2;
    smod.m_alpha2  = mat2;
    smod.me_alpha1 = mat2;
    smod.me_alpha2 = mat2;
    smod.cov       = mat2;
    smod.beta0     = mat;
    smod.beta1     = mat;
    smod.beta2     = mat;
    smod.xbar1     = mat;
    smod.xbar2     = mat; 
    smod.J0        = mat2;
    smod.J1        = mat2;
    smod.surprise  = mat * 0;
    smod.epoch     = mat * 0;
    smod.data0     = 1;
    smod.dim       = dim;
    smod.max       = struct('Description','Maximum and upper bound limits on model');
    % Obtain asymptotic maximum values for beta and beta'
    [smod.max.beta1,smod.max.beta2] = betavalues(decay,updatefac);
    smod.joint          = struct('Description','Contains values for joint surprise if computed');
    smod.joint.beta1    = 1;
    smod.joint.beta2    = 1;
    smod.joint.alpha1   = 1;
    smod.joint.alpha2   = 1;
    smod.joint.surprise = 0;
end
    
    