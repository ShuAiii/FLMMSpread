function [value,Var,CI95,CI99]=MC(N,m,int_S1,int_S2,K,T,t,sigma1,sigma2,rho,r)
%% T partation %%
tau=T-t;
GBm=zeros(N,2);
%% Loading Trained Networks %%
pathGamma11 = '\Gamma11_Net.h5';
pathGamma12 = '\Gamma12_Net.h5';
GammaNet11 = importKerasNetwork(pathGamma11);
GammaNet12 = importKerasNetwork(pathGamma12);
%% Discretization %%
dt=tau/m;
S1=zeros(1,m+1);
S1(1)=int_S1;
S2=zeros(1,m+1);
S2(1)=int_S2;
%% Sampling %%
tic
for j=1:N
    %% Sampling %%
    meanmatrix=zeros(2,1);
    covmatrix=[dt,0;0,dt];
    bm=mvnrnd(meanmatrix,covmatrix,m)';
    %% Markov Chain %%
    pS1=int_S1;
    for i = 1:m
        dtau=tau-(i-1)*dt;
        gamma11 = predict(GammaNet11,[S1(i),S2(i),K,r,sigma1,sigma2,rho,tau]);
        gamma12 = predict(GammaNet12,[S1(i),S2(i),K,r,sigma1,sigma2,rho,tau]);
        lambda=Lambda(dtau,S1(i),pS1);
        sigma11=sigma1 / (1-lambda*gamma11);
        sigma12=sigma2*lambda*gamma12*S2 / (S1*(1-lambda*gamma11));
        S1(i+1)=S1(i)+r*S1(i)*dt+sigma11*S1(i)*bm(1,i)+sigma12*S1(i)*bm(2,i);
        S2(i+1)=S2(i)+r*S2(i)*dt+sigma2*rho*S2(i)*bm(1,i)+sigma2*sqrt(1-rho^2)*S2(i)*bm(2,i);
        pS1=S1(i);
    end
    a=S1(m+1);
    b=S2(m+1);
    GBm(j,:)=[a,b];
end
%% Monte Carlo %%
    Payout = GBm(:,1)-GBm(:,2)-K;
    Ind = Payout > 0;
    V = Payout .* Ind;
    MCMu = 1/N * sum(V);
    value = exp(-r*tau) * MCMu;
%% Variance %%
Var = exp(-2*r*tau) * sum((V - MCMu).^2) / N;
std=sqrt(Var);
CI95=[value-1.96/sqrt(N)*std,value+1.96/sqrt(N)*std];
CI99=[value-2.576/sqrt(N)*std,value+2.576/sqrt(N)*std];
toc
end
