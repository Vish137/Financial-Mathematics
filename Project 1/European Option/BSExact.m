function [put, call] = BSExact(S, K, r, sigma, T)
%FuncName: BSExact.m
%Description: this computes the exact value of European options in the
%Black Scholes model
%
%---------
%ARGUMENTS
%
% - S: the initial price of the stock
% - K: the strike price in the option contract
% - r: the risk-free interest rate on the option contract
% - sigma: volatilty of the option price
% - T: the time to expiry (maturity) of the option contract
%
%---------
%OUTPUTS
%
% put: the price of a European put option in the Black Scholes Model
% call: the price of a European call option in the Black Scholes Model
%---------

d1 = (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
N1 = 0.5*(1+erf(-1*d1/sqrt(2)));
N2 = 0.5*(1+erf(-1*d2/sqrt(2)));

%Put option exact price
put = K.*exp(-r*T).*N2 - S.*N1;

%Call option exact price
call = -K.*exp(-r*T).*N2 + S.*N1;
end
