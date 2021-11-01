% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/thomaspapasavvas/Desktop/data.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 26-Oct-2021 03:34:51
% Set up the Import Options and import the data
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
data = readtable("/Users/Dan/Downloads/FMAT3888/Data.xlsx", opts, "UseExcel", false);

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
date_A_1 = find(data.Date == '31-Jan-2007');
date_A_2 = find(data.Date == '31-Dec-2010');

interval_A = data(date_A_1:date_A_2,:);

%Interval B
date_B_1 = find(data.Date == '31-Jan-2011');
date_B_2 = find(data.Date == '31-Dec-2014');

interval_B = data(date_B_1:date_B_2,:);

%Isolate Raw Data
raw_interval_A = interval_A(:,2:7);
raw_interval_A = table2array(raw_interval_A);
raw_interval_B = interval_B(:,2:7);
raw_interval_B = table2array(raw_interval_B);

%Estimate a_i's
%Interval A
A_normal = log(raw_interval_A+1); %Convert to X's X = ln(alpha + 1)
a_A = mean(A_normal); %Calculate a's

%Interval B
B_normal = log(raw_interval_B+1); %Convert to X's X = ln(alpha + 1)
a_B = mean(B_normal); %Calculate a's

%Estimate b_ij's

%Interval A
C_A = cov(A_normal); %Find covariance matrix

%Interval B
C_B = cov(B_normal); %Find covariance matrix

%
%Question 3

%Find Annual Returns

%Interval A
annual_A = zeros(3,6); %Create an empty 3by6 matrix to store annual returns
for i=1:3 %Cycle through rows of annual matrix
    raw = 1+raw_interval_A(12*(i-1)+1,:); 
    for j = ((i-1)*12+1)+1:i*12 %Cycle through each set of 12 rows (raws)
        raw = raw.*(1+raw_interval_A(j,:));  
    end
    annual_A(i,1:6) = raw;
end
annual_A = annual_A-1;

%Interval B
annual_B = zeros(3,6); %Create an empty 3by6 matrix to store annual returns
for i=1:3 %Cycle through rows of annual matrix
    raw = 1+raw_interval_B(12*(i-1)+1,:); 
    for j = ((i-1)*12+1)+1:i*12 %Cycle through each set of 12 rows (raws)
        raw = raw.*(1+raw_interval_B(j,:));
    end
    annual_B(i,1:6) = raw;
end
annual_B = annual_B-1;

%Calculate 2 Year Returns
%Interval A
two_year_A = zeros(2,6);
for i = 1:2
    row = 1+annual_A(i,:);
    for j = i+1:i+1
        row = row.*(1+annual_A(j,:));
    end
    two_year_A(i,1:6) = row;
end
two_year_A = two_year_A-1;

%Interval B
two_year_B = zeros(2,6);
for i = 1:2
    row = 1+annual_B(i,:);
    for j = i+1:i+1
        row = row.*(1+annual_B(j,:));
    end
    two_year_B(i,1:6) = row;
end
two_year_B = two_year_B-1;

%Estimate Parameters for R1 and R2

%Interval A
mu_R1_A = mean(annual_A);
mu_R2_A = mean(two_year_A);

C_R1_A = cov(annual_A);
C_R2_A = cov(two_year_A);

corr_R1_A = corrcoef(annual_A);
corr_R2_A = corrcoef(two_year_A);

%Interval B
mu_R1_B = mean(annual_B);
mu_R2_B = mean(two_year_B);

C_R1_B = cov(annual_B);
C_R2_B = cov(two_year_B);

corr_R1_B = corrcoef(annual_B);
corr_R2_B = corrcoef(two_year_B);

%
%Question 4

%Using interval_A data

clear w1 w2 w3 w4 w5 w6 lambda

syms f(w1,w2,w3,w4,w5,w6) g(w1,w2,w3,w4,w5,w6)
f(w1,w2,w3,w4,w5,w6) = w1*mu_R2_A(1)+w2*mu_R2_A(2)+w3*mu_R2_A(3)+w4*mu_R2_A(4)+w5*mu_R2_A(5)+w6*mu_R2_A(6);
g(w1,w2,w3,w4,w5,w6) = w1^2*C_R2_A(1,1)+w2^2*C_R2_A(2,2)+w3^2*C_R2_A(3,3)+w4^2*C_R2_A(4,4)+w5^2*C_R2_A(5,5)+w6^2*C_R2_A(6,6);

mu_A = f(w1,w2,w3,w4,w5,w6);
sigma_A = g(w1,w2,w3,w4,w5,w6);

syms w1 w2 w3 w4 w5 w6 lambda
S = solve(w1+w2+w3+w4+w5+w6 == 1,lambda == (-2*mu_R2_A(1)-2*w1*C_R2_A(1,1))*exp(-2*mu_A-sigma_A), lambda == (-2*mu_R2_A(2)-2*w2*C_R2_A(2,2))*exp(-2*mu_A-sigma_A),lambda == (-2*mu_R2_A(3)-2*w3*C_R2_A(3,3))*exp(-2*mu_A-sigma_A),lambda == (-2*mu_R2_A(4)-2*w4*C_R2_A(4,4))*exp(-2*mu_A-sigma_A),lambda == (-2*mu_R2_A(5)-2*w5*C_R2_A(5,5))*exp(-2*mu_A-sigma_A),lambda == (-2*mu_R2_A(6)-2*w6*C_R2_A(6,6))*exp(-2*mu_A-sigma_A));

%Calculate Expected Utility from Solution

utility_A = double(exp(-2*f(S.w1,S.w2,S.w3,S.w4,S.w5,S.w6)-g(S.w1,S.w2,S.w3,S.w4,S.w5,S.w6)));

%Using Interval B data

clear w1 w2 w3 w4 w5 w6 lambda

syms f(w1,w2,w3,w4,w5,w6) g(w1,w2,w3,w4,w5,w6)
f(w1,w2,w3,w4,w5,w6) = w1*mu_R2_B(1)+w2*mu_R2_B(2)+w3*mu_R2_B(3)+w4*mu_R2_B(4)+w5*mu_R2_B(5)+w6*mu_R2_B(6);
g(w1,w2,w3,w4,w5,w6) = w1^2*C_R2_B(1,1)+w2^2*C_R2_B(2,2)+w3^2*C_R2_B(3,3)+w4^2*C_R2_B(4,4)+w5^2*C_R2_B(5,5)+w6^2*C_R2_B(6,6);

mu_B = f(w1,w2,w3,w4,w5,w6);
sigma_B = g(w1,w2,w3,w4,w5,w6);

syms w1 w2 w3 w4 w5 w6 lambda
SS = solve(w1+w2+w3+w4+w5+w6 == 1,lambda == (-2*mu_R2_B(1)-2*w1*C_R2_B(1,1))*exp(-2*mu_B-sigma_B), lambda == (-2*mu_R2_B(2)-2*w2*C_R2_B(2,2))*exp(-2*mu_B-sigma_B),lambda == (-2*mu_R2_B(3)-2*w3*C_R2_B(3,3))*exp(-2*mu_B-sigma_B),lambda == (-2*mu_R2_B(4)-2*w4*C_R2_B(4,4))*exp(-2*mu_B-sigma_B),lambda == (-2*mu_R2_B(5)-2*w5*C_R2_B(5,5))*exp(-2*mu_B-sigma_B),lambda == (-2*mu_R2_B(6)-2*w6*C_R2_B(6,6))*exp(-2*mu_B-sigma_B));

%Calculate Expected Utility from Solution

utility_B = double(exp(-2*f(SS.w1,SS.w2,SS.w3,SS.w4,SS.w5,SS.w6)-g(SS.w1,SS.w2,SS.w3,SS.w4,SS.w5,SS.w6)));


%Efficeint Frontier 
stdev_R2_A = zeros(1, 6); 
for i=1:6
    stdev_R2_A(:, i) = sqrt(C_R2_A(i, i)); 
end 
stdev_R2_A; 
%portopt(mu_R2_B, C_R2_B, 10);

%Base case without expected return constraint, this worked fine
syms w1 w2 w3 w4 w5 w6 lambda 
A = [w1 w2 w3 w4 w5 w6] * C_R2_A * transpose([w1 w2 w3 w4 w5 w6]); %Objective
V = w1+w2+w3+w4+w5+w6 -1 == 0; %Constraint 
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



