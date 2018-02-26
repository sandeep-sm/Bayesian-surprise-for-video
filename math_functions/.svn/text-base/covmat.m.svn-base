function Cov = covmat(Mat)

SX = size(Mat,2);

Cov = ones(SX,SX);

for i = 1:SX
    Cov(i,:) = Cov(i,:) * Mat(1,i);
    Cov(:,i) = Cov(:,i) * Mat(1,i);
end
 
Cov = sqrt(Cov);