opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "B4:H248";

% Specify column names and types
opts.VariableNames = ["DATE", "AustralianEquitiesG", "DevelopedEquitiesG", "EmergingMarketEquitiesG", "AustralianFixedInterestD", "GlobalGovernmentBondsD", "CashD"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "DATE", "InputFormat", "");

% Import the data
data = readtable("rawdata.xlsx", opts, "UseExcel", false);

%Prepare Data

keepcols = {'DATE','AustralianEquitiesG', 'DevelopedEquitiesG', 'EmergingMarketEquitiesG', 'AustralianFixedInterestD', 'GlobalGovernmentBondsD', 'CashD'};
data = data(:,keepcols);

colnames = {'Date','AU_Equities','Dev_Equities','Em_Equities','AU_Fixed','Dev_Gov_Bonds','Cash'};
data.Properties.VariableNames = colnames;

%Remove First Row
data(1,:) = [];

%PARAMETER ESTIMATION

%Question 1

%Isolate the Relevant Time Intervals

%Interval A
date_A_1 = find(data.Date == '29-Dec-2006');
date_A_2 = find(data.Date == '31-Dec-2010');

interval_A = data(date_A_1:date_A_2,:);

%Interval B
date_B_1 = find(data.Date == '31-Dec-2010');
date_B_2 = find(data.Date == '31-Dec-2014');

interval_B = data(date_B_1:date_B_2,:);

%Second Time Interval 
[rows,cols] = size(interval_B);

matrix2 = table2array(interval_B(:, 2:cols));

for i=1:(cols - 1)
    matrix2(:, i) = log(matrix2(:, i) + 1); 
end 

%First Time Interval 
[rows,cols] = size(interval_A);

matrix1 = table2array(interval_A(:, 2:cols));

for i=1:(cols - 1)
    matrix1(:, i) = log(matrix1(:, i) + 1); 
end 


%Q3
%Finding 2-year returns for the first time interval 
interval_A_2year = zeros(2, 6); 
for i=1:(cols-1) 
    for j=1:24
        interval_A_2year(1, i) = sum(matrix1(j, i)); 
    end
    for j=25:49 
        interval_A_2year(2, i) = sum(matrix1(j, i));
    end 
end

%Find 2-year returns for the second time interval 
interval_B_2year = zeros(2, 6); 
for i=1:(cols-1) 
    for j=1:24 
        interval_B_2year(1, i) = sum(matrix2(j, i)); 
    end 
    for j=25:49 
        interval_B_2year(2, i) = sum(matrix2(j, i)); 
    end 
end


%Finding mu for each time interval for R2 
interval_A_mu2 = zeros(1, 6);
interval_B_mu2 = zeros(1, 6);

for i=1:6 
    interval_A_mu2(1, i) = (interval_A_2year(1, i) + interval_A_2year(2, i)) / 2; 
end 

for i=1:6 
    interval_B_mu2(1, i) = (interval_B_2year(1, i) + interval_B_2year(2, i)) / 2; 
end 


%Creating covariance matrices
covariance_A2 = cov(interval_A_2year); 
covariance_B2 = cov(interval_B_2year); 


%Creating correlation coefficient matrices (returns a matrix of 1's and -1's, something seems wrong) 
corr_A2 = zeros(6, 6); 
corr_B2 = zeros(6, 6); 

for j=1:6 
    for i=1:6 
        corr_A2(i, j) = (covariance_A2(i, j)) ./ ((sqrt(covariance_A2(i, i))) * (sqrt(covariance_A2(j, j)))); 
    end 
end 

for j=1:6 
    for i=1:6 
        corr_B2(i, j) = (covariance_B2(i, j)) ./ ((sqrt(covariance_B2(i, i))) * (sqrt(covariance_B2(j, j)))); 
    end 
end 


%Question 4 (Results in empty sym, there's probably a much better way to do this) 
syms w1 w2 w3 w4 w5 w6 lambda 
A = exp(-(interval_A_mu2(:, 1)).*w1 -(interval_A_mu2(:,2)).*w2 -(interval_A_mu2(:,3)).*w3 -(interval_A_mu2(:,4)).*w4 -(interval_A_mu2(:,5)).*w5 -(interval_A_mu2(:,6)).*w6); %Objective
V = w1+w2+w3+w4+w5+w6 == 1; %Constraint 
L = A + lambda * lhs(V); %Lagrange 
dL_dw1 = diff(L, w1) == 0; % derivative of L with respect to wl
dL_dw2 = diff(L, w2) == 0; % derivative of L with respect to w2
dL_dw3 = diff(L, w3) == 0; % derivative of L with respect to w3
dL_dw4 = diff(L, w4) == 0; % derivative of L with respect to w4
dL_dw5 = diff(L, w5) == 0; % derivative of L with respect to w5
dL_dw6 = diff(L, w6) == 0; % derivative of L with respect to w6
dL_dlambda = diff(L,lambda) == 0; % derivative of L with respect to lambda
system = [dL_dw1;dL_dw2;dL_dw3;dL_dw4;dL_dw5;dL_dw6;dL_dlambda;]; % build the system of equations
[w1_val, w2_val, w3_val, w4_val, w5_val, w6_val, lambda_val] = solve(system, [w1 w2 w3 w4 w5 w6 lambda ], 'Real', true);
results_numeric = [w1_val, w2_val, w3_val, w4_val, w5_val, w6_val, lambda_val];
results_numeric 

%Question 5
%Efficeint Frontier 
portopt(interval_A_mu2, covariance_A2, 2000)    %This is fucked, not too sure why 

%5b) (this doesn't work, not too sure why)
syms w1 w2 w3 w4 w5 w6 lambda 
A = [w1 w2 w3 w4 w5 w6] * covariance_A2 * transpose([w1 w2 w3 w4 w5 w6]); %Objective
V = w1+w2+w3+w4+w5+w6 == 1; %Constraint (Missing other constraint, no point putting it in rn if base case doesn't work) 
L = A + lambda * lhs(V); %Lagrange 
dL_dw1 = diff(L, w1) == 0; % derivative of L with respect to wl
dL_dw2 = diff(L, w2) == 0; % derivative of L with respect to w2
dL_dw3 = diff(L, w3) == 0; % derivative of L with respect to w3
dL_dw4 = diff(L, w4) == 0; % derivative of L with respect to w4
dL_dw5 = diff(L, w5) == 0; % derivative of L with respect to w5
dL_dw6 = diff(L, w6) == 0; % derivative of L with respect to w6
dL_dlambda = diff(L,lambda) == 0; % derivative of L with respect to lambda
system = [dL_dw1;dL_dw2;dL_dw3;dL_dw4;dL_dw5;dL_dw6;dL_dlambda;]; % build the system of equations
[w1_val, w2_val, w3_val, w4_val, w5_val, w6_val, lambda_val] = solve(system, [w1 w2 w3 w4 w5 w6 lambda ], 'Real', true);
results_numeric1 = [w1_val, w2_val, w3_val, w4_val, w5_val, w6_val, lambda_val];
results_numeric1
