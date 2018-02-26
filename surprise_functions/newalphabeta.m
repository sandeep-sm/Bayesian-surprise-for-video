

function smod = newalphabeta(data,smod,x)

%%%%%%%%%%%%%
% UPDATE ALPHA
% Here we can try a few different alpha updates for experimentation
% purposes
smod.alpha0 = smod.alpha1;
smod.alpha1 = smod.alpha2;

% Robbins-Monro Style Update
if strcmp(smod.options.robbins_monro,'yes') 
    smod.alpha2 = robbinsmonroalpha(data,smod.alpha1,smod.beta1,smod.decay);   
elseif strcmp(smod.options.robbins_monro_2,'yes') 
    smod.alpha2 = robbinsmonroalpha(data,smod.alpha1,smod.max.beta1,smod.decay); 
elseif strcmp(smod.options.robbins_monro_3,'yes') 
    smod.alpha2 = robbinsmonroalpha(data,smod.alpha1,smod.beta1,smod.decay,smod.max.beta1); 
else
    smod.temp_prod = fuzarith(x,smod.alpha1,smod.decay,'prod')';
    smod.alpha2 = fuzarith(x,smod.temp_prod,data,'sum');
end

% Update Beta based on whether this is a univariate model or a multivariate
% model

if smod.dim == 1
    smod = uni_beta(data,smod,x);
elseif smod.dim == 2
    %smod = mult_ln_beta(data,smod);
    smod = biv_beta(data,smod);
elseif smod.dim > 2
    smod = mult_parts_beta(data,smod);
else
    error('Error in dimensions specified as `smod.dim`');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Univariate Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call this to run surprise on univariate data

function smod = uni_beta(data,smod,x) 
%%%%%%%%%%%%%%
% UPDATE BETA
% Here we can try a few different beta updates for experimentation
% purposes
 smod.beta1  = smod.beta2;

% Use Kalman Style Decay Factorization
if strcmp(smod.options.factordecay,'yes') 
    if smod.options.debug > 1
        fprintf('Using option `factordecay` in surprise model\n');
    end
    
    % The expected value of alpha is different depending on how we compute
    % it
    if strcmp(smod.options.robbins_monro,'yes') 
        expectedAlpha2 = robbinsmonroalpha(smod.xbar2,smod.alpha1,smod.beta1,smod.decay);   
    elseif strcmp(smod.options.robbins_monro_2,'yes') 
        expectedAlpha2 = robbinsmonroalpha(smod.xbar2,smod.alpha1,smod.max.beta1,smod.decay); 
    elseif strcmp(smod.options.robbins_monro_3,'yes') 
        expectedAlpha2 = robbinsmonroalpha(smod.xbar2,smod.alpha1,smod.beta1,smod.decay,smod.max.beta1); 
    else
        smod.temp_prod = fuzarith(x,smod.alpha1,smod.decay,'prod');
        expectedAlpha2 = fuzarith(x,smod.temp_prod,smod.xbar2,'sum');      
    end
    %fuzzy maths
    temp_prod = fuzarith(x,smod.beta1,smod.decay,'prod');
    fuz_1 = trimf(x,[1 1 1]);
    temp_sum = fuzarith(x,temp_prod,fuz_1,'sum');
    temp_prod = fuzarith(x,temp_sum,expectedAlpha2,'prod');
    
    %final result for beta2
    %smod.beta2  = expectedAlpha2*(smod.beta1*smod.decay + 1)./smod.alpha2
    smod.beta2  = fuzarith(x,temp_prod,smod.alpha2,'div'); 
    
    if smod.options.debug > 1
        fprintf('Values [%f,%f,%f,%f]\n',smod.xbar1,smod.xbar2,smod.beta1,smod.beta2);
    end
    
    
    %smod.xbar2  = (data + smod.xbar1*smod.decay) ./ (1 + smod.decay);  
    %smod.xbar1  = smod.xbar2;
    
% Base Surprise Beta Update    
elseif strcmp(smod.options.setbetamax,'no') 
    if smod.options.debug > 1
        fprintf('Updating beta based on standard model\n');
    end
    temp_prod = fuzarith(x,smod.beta1,smod.decay,'prod');
    smod.beta2 = fuzarith(x,temp_prod,smod.updatefac,'sum');
else
    fprintf('Using option `setbetamax` in surprise model\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bivariate Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call this to run surprise on bivariate data

% I am putting in the bivariate model to make it stable, then I will put in
% the multivaraite model using it as a template
function smod = biv_beta(data,smod) 

smod = mult_parts_beta(data,smod);

if smod.options.debug > 1
    fprintf('PRIOR BETA2(1) = %f\n',smod.beta1(:,1));
    fprintf('PRIOR BETA2(2) = %f\n',smod.beta1(:,2));
    fprintf('BETA VALUES (1): [E(A2):A2:E(A1):A1:COV] = [ %f : %f : %f : %f : %f]\n',...
            smod.e_alpha2(:,1),smod.alpha2(:,1),smod.e_alpha1(:,1),smod.alpha1(:,1),smod.cov(1,2));
    fprintf('BETA VALUES (2): [E(A2):A2:E(A1):A1:COV] = [ %f : %f : %f : %f : %f]\n',...
            smod.e_alpha2(:,2),smod.alpha2(:,2),smod.e_alpha1(:,2),smod.alpha1(:,2),smod.cov(2,1));
    fprintf('BETA2(1) = %f\n',smod.beta2(:,1));
    fprintf('BETA2(2) = %f\n',smod.beta2(:,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multivariant Model - Computes surprise as parts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call this to run surprise on multivariant data

function smod = mult_parts_beta(data,smod)

%%%%%%%%%%%%%%
% UPDATE BETA
% Here we can try a few different beta updates for experimentation
% purposes
smod.beta0 = smod.beta1;
smod.beta1 = smod.beta2;

% The expected value of alpha is different depending on how we compute
% it
% This is not yet optimized. It's in a more visible form which may be
% better so that it is understandable. Opinions?
if strcmp(smod.options.robbins_monro,'yes') 
    smod.e_alpha1 = robbinsmonroalpha(smod.xbar1,smod.alpha0,smod.beta0,smod.decay); 
    smod.e_alpha2 = robbinsmonroalpha(smod.xbar2,smod.alpha1,smod.beta1,smod.decay); 
elseif strcmp(smod.options.robbins_monro_2,'yes') 
    smod.e_alpha1 = robbinsmonroalpha(smod.xbar1,smod.alpha0,smod.max.beta1,smod.decay);
    smod.e_alpha2 = robbinsmonroalpha(smod.xbar2,smod.alpha1,smod.max.beta1,smod.decay);
elseif strcmp(smod.options.robbins_monro_3,'yes') 
    smod.e_alpha1 = robbinsmonroalpha(smod.xbar1,smod.alpha0,smod.beta0,smod.decay,smod.max.beta1);
    smod.e_alpha2 = robbinsmonroalpha(smod.xbar2,smod.alpha1,smod.beta1,smod.decay,smod.max.beta1);
else
    smod.e_alpha1 = smod.alpha0*smod.decay + smod.xbar1;     
    smod.e_alpha2 = smod.alpha1*smod.decay + smod.xbar2;
end

% Find the general covariance matrix
for i = 1:smod.dim
    for j = 1:smod.dim
        if i ~= j
            smod.cov(i,j) = (smod.alpha1(:,j)*(smod.e_alpha2(:,i)/smod.e_alpha2(:,j) - smod.e_alpha1(:,i)/smod.e_alpha1(:,j)))/(smod.dim - 1);
        else
            smod.cov(i,j) = smod.e_alpha2(:,i);
        end
    end
end

% Compute the expected value of alpha' based on the covariance matrix
expectedAlpha2 = sum(smod.cov,2)';

% Compute the final beta update based on the ratio of E(alpha2)/alpha2 
smod.beta2 = (expectedAlpha2./smod.alpha2).*(smod.beta1*smod.decay + 1);

% Recompute the basic prediction for the expected value of alpha as xbar
smod.xbar2  = (data      + smod.xbar1*smod.decay) / (1 + smod.decay);  
smod.xbar1  = smod.xbar2;

if strcmp(smod.options.jointmodel,'blind')
    smod.joint.beta1   = smod.joint.beta2; 
    smod.joint.beta2   = sum(sum(smod.beta2))/smod.dim;    
    smod.joint.xbar    = sum(sum(data))/smod.dim;
    smod.joint.alpha2  = robbinsmonroalpha(smod.joint.xbar,smod.joint.alpha1,smod.joint.beta2,smod.decay);
    smod.joint.alpha1  = smod.joint.alpha2;
elseif strcmp(smod.options.jointmodel,'mux')
    smod.joint.beta1   = smod.joint.beta2;
    smod.joint.beta2   = sum(sum(smod.beta2));
    smod.joint.xbar    = sum(sum(data))/smod.dim;
    smod.joint.alpha2  = robbinsmonroalpha(smod.joint.xbar,smod.joint.alpha1,smod.joint.beta2,smod.decay);
    smod.joint.alpha1  = smod.joint.alpha2;
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multivariant Model - Computes surprise as a joint model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call this to run surprise on multivariant data

function smod = mult_joint_beta(data,smod)