function  value = Margrabe(int_S1,int_S2,tau,sigma1,sigma2,rho)


sigma=sqrt(sigma1^2+sigma2^2-2*rho*sigma1*sigma2);

d1=(log(int_S1/int_S2)+1/2*sigma^2*tau)/(sigma*sqrt(tau));
d2=d1-sigma*sqrt(tau);

value=int_S1*normcdf(d1)-int_S2*normcdf(d2);

end

