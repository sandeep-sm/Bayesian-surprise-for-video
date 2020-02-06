function d = bhattacharya_distance(a1,a2,b1,b2,x)
a1_a = (find(a1,1)-1)*200/4000000;
a1_b = find(a1 == max(a1),1)*200/4000000;
a1_c = (find(a1,1,'last')+1)*200/4000000;

a2_a = (find(a2,1)-1)*200/4000000;
a2_b = find(a2 == max(a2),1)*200/4000000;
a2_c = (find(a2,1,'last')+1)*200/4000000;

b1_a = (find(b1,1)-1)*200/4000000;
b1_b = (find(b1 == max(b1),1))*200/4000000;
b1_c = (find(b1,1,'last')+1)*200/4000000;

b2_a = (find(b2,1)-1)*200/4000000;
b2_b = (find(b2 == max(b2),1))*200/4000000;
b2_c = (find(b2,1,'last')+1)*200/4000000;

%{
temp = trimf(x,[2 2 2]);

% this is alpha'
alpha0_temp = fuzarith(x,a1,a2,'sum');
alpha0 = fuzarith(x,alpha0_temp,temp,'div');

alpha0_a = (find(alpha0,1)-1)*200/4000000;
alpha0_b = (find(alpha0 == max(alpha0),1))*200/4000000;
alpha0_c = (find(alpha0,1,'last')+1)*200/4000000;

% this is beta'
beta0_temp = fuzarith(x,b1,b2,'sum');
beta0 = fuzarith(x,beta0_temp,temp,'div');

beta0_a = (find(beta0,1)-1)*200/4000000;
beta0_b = (find(beta0 == max(beta0),1))*200/4000000;
beta0_c = (find(beta0,1,'last')+1)*200/4000000;

% this is log(beta')
log_beta0_temp = sort([log(beta0_a) log(beta0_b) log(beta0_c)]);
log_beta0 = trimf(x,log_beta0_temp);

% this is log(gammma(alpha'))
log_gamma_alpha0_temp = sort([gamma(log(alpha0_a)) gamma(log(alpha0_b)) gamma(log(alpha0_c))]);
log_gamma_alpha0 = trimf(x,log_gamma_alpha0_temp);

%}

%{
fprintf('%f -- %f -- %f',a2_c,a2_b,a2_a);
gamma_a1 = trimf(x,[ gamma(a1_c) gamma(a1_b) gamma(a1_a)]);
gamma_a2 = trimf(x,[a2_c a2_b a2_a]);
%}

%d = (1/2)*(a1+a2)log((b1+b2)/2) - (a1/2)*log(b1) - (a2/2)*log(b2) + 0.5*(log(gamma(a1)) + log(gamma(a2))) - log(gamma((a1+a2)/2));

d_a = (1/2)*(a1_a+a2_a)*log((b1_a+b2_a)/2) - (a1_a/2)*log(b1_a) - (a2_a/2)*log(b2_a) + 0.5*(log(gamma(a1_a)) + log(gamma(a2_a))) - log(gamma((a1_a+a2_a)/2));
d_b = (1/2)*(a1_b+a2_b)*log((b1_b+b2_b)/2) - (a1_b/2)*log(b1_b) - (a2_b/2)*log(b2_b) + 0.5*(log(gamma(a1_b)) + log(gamma(a2_b))) - log(gamma((a1_b+a2_b)/2));
d_c = (1/2)*(a1_c+a2_c)*log((b1_c+b2_c)/2) - (a1_c/2)*log(b1_c) - (a2_c/2)*log(b2_c) + 0.5*(log(gamma(a1_c)) + log(gamma(a2_c))) - log(gamma((a1_c+a2_c)/2));

d1 =[round(d_a,4),round(d_b,4),round(d_c,4)];
d1 = sort(d1);
%fprintf('%f -- %f -- %f\n\n',round(d1(1),4),round(d1(2),4),round(d1(3),4));
%pause;
d = d1(2);
%plot(x,d);