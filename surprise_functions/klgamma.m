function d = klgamma(a1,a2,b1,b2,x)
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
a2_a = gamma(a2_a);
a2_b = gamma(a2_b);
a2_c = gamma(a2_c);

fprintf('%f -- %f -- %f',a2_c,a2_b,a2_a);
gamma_a1 = trimf(x,[ gamma(a1_c) gamma(a1_b) gamma(a1_a)]);
gamma_a2 = trimf(x,[a2_c a2_b a2_a]);

%{
figure;
plot(x,gamma_a1);
figure;
plot(x,gamma_a1);

pause;
%}

div1 = fuzarith(x,gamma_a2,gamma_a1,'div');
div1_a = find(div1,1)-1;
div1_b = find(div1 == 1,1);
div1_c = find(div1,1,'last')+1;
term3 = trimf(x,[log(div1_a) log(div1_b) log(div1_c)]);

div1 = fuzarith(x,b1,b2,'div');
div1_a = find(div1,1)-1;
div1_b = find(div1 == 1,1);
div1_c = find(div1,1,'last')+1;
temp = trimf(x,[log(div1_a) log(div1_b) log(div1_c)]);
term2 = fuzarith(x,a2,temp,'prod');

temp = fuzarith(x,a1,b1,'div');
term4 = fuzarith(x,b2,temp,'prod');

d = -a2 + term2 + term3 + term4 + (a1 - a2).*digamma(a1,1000);

%}

if a1_a == 0
    a1_a = a1_b;
end
if a2_a == 0
    a2_a = a2_b;
end

fprintf('%f -- %f -- %f\n',a1_a,a1_b,a1_c);
fprintf('%f -- %f -- %f\n',a2_a,a2_b,a2_c);
fprintf('%f -- %f -- %f\n',b1_a,b1_b,b1_c);
fprintf('%f -- %f -- %f\n',b2_a,b2_b,b2_c);
%d = 1;
%pause;


d_a = -a2_a + a2_a.*log(b1_a./b2_a) + log(gamma(a2_a)./gamma(a1_a)) + b2_a.*(a1_a./b1_a) + (a1_a - a2_a).*digamma(a1_a,1000);
d_b = -a2_b + a2_b.*log(b1_b./b2_b) + log(gamma(a2_b)./gamma(a1_b)) + b2_b.*(a1_b./b1_b) + (a1_b - a2_b).*digamma(a1_b,1000);
d_c = -a2_c + a2_c.*log(b1_c./b2_c) + log(gamma(a2_c)./gamma(a1_c)) + b2_c.*(a1_c./b1_c) + (a1_c - a2_c).*digamma(a1_c,1000);

if ~d_a
    if ~d_b
        if ~d_c
            d_c = trimf(x,[1 1 1]);
        end
        d_b = d_c;
    elseif ~d_c
        d_c = d_b;
    end
    d_a = round((d_b+d_c)/2,4);
elseif ~d_b
    if ~d_a
        if ~d_c
            d_c = trimf(x,[1 1 1]);
        end
        d_a = d_c;
    elseif ~d_c
        d_c = d_a;
    end
    d_b = round((d_a+d_c)/2,4);
elseif ~d_c
    if ~d_a
        if ~d_b
            d_b = trimf(x,[1 1 1]);
        end
        d_a = d_b;
    elseif ~d_b
        d_b = d_a;
    end
    d_c = round((d_a+d_b)/2,4);
end

d1 =[round(d_a,3),round(d_b,3),round(d_c,3)];
d1 = sort(d1);
fprintf('\n%f -- %f -- %f\n\n',round(d1(1),4),round(d1(2),4),round(d1(3),4));
d = d1(2);
%plot(x,d);