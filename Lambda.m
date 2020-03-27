function lambda = Lambda(tau,S1,pS1)
%% Function Parameters %%
e=0.04;
b=100;
floor=0.6*pS1;
cap=1.4*pS1;
%% Price Impact Function %%
if (S1>floor) && (S1<cap)
    lambda=e*(1-exp(-b*tau^(3/2)));
else 
    lambda=0;
end
