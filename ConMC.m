function [ConValue,ConVar,CI95,CI99]=ConMC(N,m,int_S1,int_S2,K,T,t,sigma1,sigma2,rho,r)
%% T partation %%
tau=T-t;
ConGBm=zeros(N,4);
%% Loading Trained Networks %%
pathGamma11 = 'C:\Users\Jkzhang\Desktop\MsC\Gamma11_Net.h5';
pathGamma12 = 'C:\Users\Jkzhang\Desktop\MsC\Gamma12_Net.h5';
GammaNet11 = importKerasNetwork(pathGamma11);
GammaNet12 = importKerasNetwork(pathGamma12);
%% Discretization %%
dt=tau/m;
S1=zeros(1,m+1);
S1(1)=int_S1;
S2=zeros(1,m+1);
S2(1)=int_S2;
cS1=zeros(1,m+1);
cS1(1)=int_S1;
cS2=zeros(1,m+1);
cS2(1)=int_S2;
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
            cS1(i+1)=cS1(i)+r*cS1(i)*dt+sigma1*cS1(i)*bm(1,i);
            cS2(i+1)=cS2(i)+r*cS2(i)*dt+sigma2*rho*cS2(i)*bm(1,i)+sigma2*sqrt(1-rho^2)*cS2(i)*bm(2,i);
        end
        a = S1(m+1);
        b = S2(m+1);
        c = cS1(m+1);
        d = cS2(m+1);
        ConGBm(j,:)=[a,b,c,d];
    end
%% Control Monte Carlo %%
    tmean=Margrabe(int_S1,int_S2,tau,sigma1,sigma2,rho);
    Payout = ConGBm(:,1)-ConGBm(:,2)-K;
    Ind = Payout > 0;
    V = Payout .* Ind;
    Mu = 1/N * sum(V);
    ConPayout = ConGBm(:,3)-ConGBm(:,4);
    ConInd = ConPayout > 0;
    ConV = ConPayout .* ConInd;
    ConMu = 1/N * sum(ConV);
    
    Var=1/N*sum((V-Mu).^2);
    ConVar=1/N*sum((ConV-ConMu).^2);
    Cov=1/N*sum((V-Mu).*(ConV-ConMu));
    c=-Cov/ConVar;
    
    ConValue=exp(-r*T)/N*sum(V+c*ConV)-c*tmean;
toc
%% Variance %%
ConVar=exp(-2*r*tau)*Var*(1-(Cov/(sqrt(Var)*sqrt(ConVar)))^2);

std=sqrt(ConVar);
CI95=[ConValue-1.96/sqrt(N)*std,ConValue+1.96/sqrt(N)*std];
CI99=[ConValue-2.576/sqrt(N)*std,ConValue+2.576/sqrt(N)*std];
end