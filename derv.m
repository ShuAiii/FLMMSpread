function [sigma11,sigma12,d1sigma11,d2sigma11,d1sigma12,d2sigma12] = derv(dtau,S1,S2,sigma1,sigma2,rho,pS1)
%% Greeks %%
[gamma11,gamma12,speed111,speed112,speed122]=Greek(S1,S2,sigma1,sigma2,rho,dtau); 
%% Price Impact Function %%
[lambda]=Lambda(dtau,S1,pS1);
%% Volatility Function %%
sigma11=sigma1 / (1-lambda*gamma11);
sigma12=sigma2*lambda*gamma12*S2 / (S1*(1-lambda*gamma11));
%% Derivative Function %%
d1sigma11=sigma1*(1/(1-lambda*gamma11)+S1*lambda*speed111/(1-lambda*gamma11)^2);
d2sigma11=sigma1*S1*lambda*speed112/(1-lambda*gamma11)^2;
d1sigma12=sigma2*S2*(lambda*speed112/(1-lambda*gamma11)+lambda^2*gamma12*speed111/(1-lambda*gamma11)^2);
d2sigma12=sigma2*lambda*((S2*speed122+gamma12)/(1-lambda*gamma11)+lambda*S2*gamma12*speed112/(1-lambda*gamma11)^2);
end