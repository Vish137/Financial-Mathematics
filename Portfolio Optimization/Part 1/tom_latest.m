%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/thomaspapasavvas/Desktop/data.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 26-Oct-2021 03:34:51

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
data = readtable("/Users/thomaspapasavvas/Desktop/data.xlsx", opts, "UseExcel", false);

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

%Find Annual Returns

%Interval A
annual_A = zeros(3,6); %Create an empty 3by6 matrix to store annual returns
for i=1:3 %Cycle through rows of annual matrix
    raw = raw_interval_A(12*(i-1)+1,:); 
    for j = (i-1)*12+1:i*12 %Cycle through each set of 12 rows (raws)
        raw = raw.*(1+raw_interval_A(j,:)); 
    end
    annual_A(i,1:6) = raw;
end

%Interval B
annual_B = zeros(3,6); %Create an empty 3by6 matrix to store annual returns
for i=1:3 %Cycle through rows of annual matrix
    raw = raw_interval_B(12*(i-1)+1,:); 
    for j = (i-1)*12+1:i*12 %Cycle through each set of 12 rows (raws)
        raw = raw.*(1+raw_interval_B(j,:)); 
    end
    annual_B(i,1:6) = raw;
end

disp(annual_A);
disp(annual_B);

    


