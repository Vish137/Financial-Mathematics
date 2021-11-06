%% Set up the Import Options and import the data
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



%%
%Question 3

%Find the mean annual return Interval A & B
One_A=12*a_A;
One_B=12*a_B;
%Find Annual Returns
%Find the Covariance Matrix for annual Return A & B
Cov_One_A=12*C_A;
Cov_One_B=12*C_B;
%%Since the e^Y-1=Y, Y is close to zero, we use naive approach.
%Find the Ceoffience corelation for annual return A & B
Corr_One_A=corrcov(Cov_One_A);
Corr_One_B=corrcov(Cov_One_B);

%Find the mean two year return Interval A & B
Two_A=24*a_A;
Two_B=24*a_B;
%Interval B
%Find the Covariance Matrix for two year Return A & B
Cov_Two_A=24*C_A;
Cov_Two_B=24*C_B;
%Find the  Ceoffience corelation for two year return A & B
Corr_Two_A=corrcov(Cov_Two_A);
Corr_Two_B=corrcov(Cov_Two_B);

% Q6
% a

R_1 = exp(mvnrnd(12*a_B,12*C_B,12));
R_2 = exp(mvnrnd(24*a_B,24*C_B,12));
fun = @(x) -1 * mean(-1 * exp(-1 * reshape(((1 + R_1 * x(1:6)) * (1 + R_2 * x(7:12))' - 1),[],1)));
x0 = [0.1,0.2,0.3,0.2,0.1,0.1,0.1,0.2,0.3,0.2,0.1,0.1]';
% x1 = [1,0,0,0,0,0,1,0,0,0,0,0]';
Aeq = [1,1,1,1,1,1,0,0,0,0,0,0;0,0,0,0,0,0,1,1,1,1,1,1];
Beq = [1;1];
LB = 0.0 * ones(12);
UB = 1.0 * ones(12);
[x_6_a,fval_6_a] = fmincon(fun,x0,[],[],Aeq,Beq,LB,UB)









