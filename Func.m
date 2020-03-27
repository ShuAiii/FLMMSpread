function [sigma11,sigma12] = Func(dtau,S1,S2,sigma1,sigma2,rho,pS1)
%% Greeks %%
[gamma11,gamma12]=Greek(S1,S2,sigma1,sigma2,rho,dtau); 
%% Price Impact Function %%
lambda=Lambda(dtau,S1,pS1);
%% Drift Functions %%
% Ft=lambda*charm1;
% F2=lambda*gamma12*drift2*S2;
% F11=0.5*(speed111*(sigma1^2*S1^2+lambda^2*gamma11^2*sigma2^2*S2^2+2*lambda*gamma12*rho*sigma1*sigma2*S1*S2))/(1-lambda*gamma11)^2;
% F12=speed112*(rho*sigma1*sigma2*S1*S2+lambda*gamma12*sigma2^2*S2^2)/(1-lambda*gamma11);
% F22=0.5*speed122*sigma2^2*S2^2;
% mu1=1/(1-lambda*gamma11) * (drift1*S1+Ft+F2+F11+F12+F22);
% mu2=drift2;
%% Volatility Function %%
sigma11=sigma1 / (1-lambda*gamma11);
sigma12=sigma2*lambda*gamma12*S2 / (S1*(1-lambda*gamma11));
end

