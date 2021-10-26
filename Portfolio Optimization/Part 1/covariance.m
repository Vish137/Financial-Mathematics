%Title: covariance.m
%Date: 26/10/21 16:10
%Author: Vish137
%Description: calculates covariance (sample and population) of assets (i,j)
%
function [covariances,covariancep] = covariance(i,j)
arguments
    i (1,1) double {mustBeNonnegative}
    j (1,1) double {mustBeNonnegative}
end
   
%TODO
%Create KWARG for specified time period - or do it manually
%
%Import and load monthly returns data
load('returnsdata.mat')
[rows,~] = size(returnsdata); 
isum=0;
nsum=0;
for index=1:rows
    isum = isum + (returnsdata(index,i) - mean(returnsdata(:,i)))*(returnsdata(index,j) - mean(returnsdata(:,j)))/(rows-1);
    nsum = nsum + (returnsdata(index,i) - mean(returnsdata(:,i)))*(returnsdata(index,j) - mean(returnsdata(:,j)))/(rows);
end
covariances=isum; %sample covariance
covariancep=nsum; %population covariance
end

