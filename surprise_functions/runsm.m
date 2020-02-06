function smod = runsm(data,smod,x,iterator)

% We require that smod be reset each run. This keeps values clean
if isfield(smod,'smodisset')
    if smod.smodisset == 0
        fprintf('NOTICE: SMOD has been set once in newsm, but is no longer usable\n');
        fprintf('Most likely you have tried to call runsm in batch mode twice on\n');
        fprintf('the same smod. You must make a new smod each time you call a\n');
        fprintf('batch mode run on your data\n');
        error('SMOD is not valid');
    end
else
    error('You must create SMOD using the function newsm before you call this function');
end

% Use the asymptotic max value for beta
if strcmp(smod.options.setbetamax,'yes')
    smod.beta1    = smod.max.beta1;
    smod.beta2    = smod.max.beta2;
end

% Figure out if we are running sample by sample or in batch mode. If
% running in batch mode, we may need to transpose the vector matrix
if size(data,1) > 1
    if smod.options.debug == 2
        fprintf('Running in batch mode. Input is a column vector\n')
    end
    smod = runsmbatch(data,smod,x,iterator);
elseif size(data,2) > 1  
    if smod.dim == 1
        if smod.options.debug == 1
            fprintf('Running in batch mode (Univariate). Input is a row vector\n')
        end
        data = data';
        smod = runsmbatch(data,smod,x,iterator);
    else
        if smod.options.debug == 2
            fprintf('Running in single step mode (Multivariate). Input is a row vector\n')
        end
        smod = runsmsingle(data,smod,x);
    end  
else
    if smod.options.debug == 2
        fprintf('Running in single step mode (Univaraite)\n')
    end
    smod = runsmsingle(data,smod,x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch Run Function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call this to run surprise on a sample vector

function smod = runsmbatch(data,smod,x,iterator)

SX = size(data,1);
if smod.dim == 1
    SY = 1;
else
    SY = size(data,2);
end

% Create surprise values output matrix
smod.surprise = zeros(SX,SY);
smod.surprise = [];
smod.iterator = iterator;
% create some debug information if requested
if smod.options.debug > 0
    smod.debugdata        = struct('Description','Debug values from runsm');
    smod.debugdata.alpha1 = zeros(SX,SY);
    smod.debugdata.alpha2 = zeros(SX,SY);
    smod.debugdata.beta1  = zeros(SX,SY);
    smod.debugdata.beta2  = zeros(SX,SY);
end

% For each data item run the surprise model on it. This is the core of the
% batch mode surprise model. We compute new beta and alpha values fromt the
% sample value, then we compute the KL distance between the two Gamma PDF's
% using klgamma. We take the absolute value to support negative data
% values. However, using negative values as inputs may not in fact make
% sense. 
for n = 1:SX
    Sy = size(data,2);
    for var = 1:Sy
        smod = newalphabeta(data(n,Sy),smod,x);
        %smod.surprise = [smod.surprise abs(klgamma(smod.alpha1,smod.alpha2,smod.beta1,smod.beta2,x))];
        smod.surprise = [smod.surprise abs(bhattacharya_distance(smod.alpha1,smod.alpha2,smod.beta1,smod.beta2,x))];
        smod.iterator = smod.iterator + 1;
        fprintf('iteration number : %d\n',smod.iterator);
    end
    if ~strcmp(smod.options.jointmodel,'none')
        smod.joint.surprise(n,:) = abs(klgamma(smod.joint.alpha1,smod.joint.alpha2,smod.joint.beta1,smod.joint.beta2,x));
    end
    
    if smod.options.debug > 0
        smod.debugdata.alpha1(n,:) = smod.alpha1;
        smod.debugdata.alpha2(n,:) = smod.alpha2;
        smod.debugdata.beta1(n,:)  = smod.beta1;
        smod.debugdata.beta2(n,:)  = smod.beta2;
        if smod.options.debug > 1
            fprintf('RUNNING INPUT %f LOOP %d ALPHA [%f,%f] BETA [%f,%f] SURPRISE VALUE %f\n',data(n,1),n,smod.alpha1,smod.alpha2,smod.beta1,smod.beta2,smod.surprise(n,1));
            fprintf('\n');    
        end
    end
    
    smod.epoch  = smod.epoch + 1;
    smod.smodisset = 0;
end

% Graph the surprise results if requested.
if strcmp(smod.options.graph,'surprise')
    if strcmp(smod.options.setbetamax,'yes')
        Title = 'setbetamax';
    elseif strcmp(smod.options.factordecay,'yes')
        Title = 'factordecay';
    else
        Title = 'basic';
    end
    
    res = basicsurprisegraph(smod.surprise,data,'Surprise Value','Input Value   ',Title,smod);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single Run Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call this to run surprise on one sample at a time. 

function smod = runsmsingle(data,smod,x)
% For each data item run the surprise model on it. This is the core of the
% batch mode surprise model. We compute new beta and alpha values fromt the
% sample value, then we compute the KL distance between the two Gamma PDF's
% using klgamma. We take the absolute value to support negative data
% values. However, using negative values as inputs may not in fact make
% sense. 
smod = newalphabeta(data,smod,x);
smod.surprise = abs(klgamma(smod.alpha1,smod.alpha2,smod.beta1,smod.beta2));
if ~strcmp(smod.options.jointmodel,'none')
    smod.joint.surprise = abs(klgamma(smod.joint.alpha1,smod.joint.alpha2,smod.joint.beta1,smod.joint.beta2));
end
%smod.alpha1   = smod.alpha2; 
%if strcmp(smod.options.setbetamax,'no')
%    smod.beta1    = smod.beta2;
%end
smod.epoch    = smod.epoch + 1;

