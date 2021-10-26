%Import and load monthly returns data
load('returnsdata.mat')
[rows,cols] = size(returnsdata);

%Q1 - do for specified time periods
a=zeros(cols,1);
for i=1:cols
    a(i) = mean(returnsdata(:,i));
end

%Q1 - do for specified time periods
b1=zeros(6,1);
b2=zeros(6,1);
b3=zeros(6,1);
b4=zeros(6,1);
b5=zeros(6,1);
b6=zeros(6,1);
for i=1:6
    b1(i)=covariance(1,i);
    b2(i)=covariance(2,i);
    b3(i)=covariance(3,i);
    b4(i)=covariance(4,i);
    b5(i)=covariance(5,i);
    b6(i)=covariance(6,i);
end
% Covariance matrix
A = [b1,b2,b3,b4,b5,b6];

% Diagonal symmetry check
sym = issymmetric(A);
%~1(true) ~0(false)

%Positive semi-definite check
[d,lambda] = eig(A);
count=0;
count1=0;
for i=1:length(lambda)
    if lambda(i,i) < 0
        count=count+1;
    elseif lambda(i,i) >=0
        count1=count1+1;
    end
end 
if count>0
    disp('Covariance matrix is NOT Positive Semi-Definite')
elseif count1==6
    disp('Covariance matrix is Positive Semi-Definite')
end

%Visual representation of covariance matrix
figure(1)
imagesc(A)
colormap flag;
title('Covariance Matrix');xlabel('Asset Class');ylabel('Asset Class');