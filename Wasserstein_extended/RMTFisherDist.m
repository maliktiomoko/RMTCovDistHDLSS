function [int_test] = RMTFisherDist(X,Y,a,a0)
%Function that compute the Wasserstein distance between Gaussian centered
%distribution based on the article Random  Matrix-Improved  Estimation  of  the  Wasserstein  Distance
%between  two  Centered  Gaussian  Distribution (Malik TIOMOKO & Romain Couillet)
%Input Need the samples from the first class X of dimension p*n and the
%samples from the second class Y of size p*n
%Return the estimate est proposed in the article and the classical esthat
%Define the dimensions
p=size(X,1);
n1=size(X,2);
n2=size(Y,2);
c1=p/n1;
c2=p/n2;
%Sample covariance estimate
hatC1=X*X'/n1;hatC2=Y*Y'/n2;
 lambda=sort(eig(hatC1\hatC2));
 m=@(z) mean(1./(lambda-z));
mp=@(z) mean(1./(lambda-z).^2);
phi=@(z) z+c1*z.^2.*m(z);
phip=@(z) 1+2*c1*z.*m(z)+c1*z.^2.*mp(z);
psi=@(z) 1-c2-c2*z.*m(z);
psip=@(z) -c2*m(z)-c2*z.*mp(z);
eta=sort(real(eig(diag(lambda)-(1/(p-n1))*sqrt(lambda)*sqrt(lambda)')));
zeta=sort(real(eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
%g=@(x) a(6)*x.^6;
 maxV = max(lambda)*1.5;
lambdac=lambda(lambda>0.001);
etac=eta(eta>0.001);
zetac=zeta(zeta>0.001);

 phi_test=@(z) (1-c1)*z.*(z-etac(1))./(z-lambdac(1));
 phi_psi=@(z) (1-c1)*(z-etac(1));
 psi_test=@(z) z./(z-lambdac(1));
 rest=@(z) 1./(z-etac(1));
 g=@(x) a0;
 for i=1:length(a)
 g=@(x) g(x)+a(i)*x.^i;
 end
 tro1=@(z) log((-c1/c2.*z)+(-c1/c2).*z.*(1/c1+1/c2-1)./(-1/c2+1+z.*m(z)));
 tro2=@(z) log(abs(-c1/c2.*z))+log(1+abs((1/c1+1/c2-1)./(-1/c2+1+z.*m(z))));
  %tro=@(z) log((-c1/c2.*z))+log(1+(1/c1+1/c2-1)./(-1/c2+1+z.*m(z)));
 %tro=@(z) (-c1/c2)*(1+(1/c1+1/c2-1)/(-1/c2)+1+z.*m(z));
 for i=1:length(zetac)
     phi_psi=@(z) phi_psi(z).*(z-etac(i+1))./(z-zetac(i));
     phi_test=@(z) phi_test(z).*(z-etac(i+1))./(z-lambdac(i+1));
     psi_test=@(z) psi_test(z).*(z-zetac(i))./(z-lambdac(i+1));
     rest=@(z) rest(z)+1./(z-etac(i+1))-1./(z-zetac(i));
 end
integrand=@(z) psi(z).*((phi(z)./psi(z)).^n).*rest(z);


minV = -1;altitude=1;step=1e-4;
contour = fliplr([altitude*1i+minV+(0:step:maxV-minV),altitude*1i+maxV+(-step:-step:-2)*1i*altitude,-1i*altitude+maxV+(-step:-step:-maxV+minV),minV-1i*altitude+(0:step:2)*1i*altitude]); 
 integrand = @(z) g(phi(z)./psi(z)).*(psi(z)/c2).*(phip(z)./phi(z)-psip(z)./psi(z));
 int2=real(trapz(contour,(1/(2*pi*1j))*integrand(contour)))-a0*(1-c2)/c2;
 int_test=0;
   for i=1:length(a)
       int_test=int_test+a(i)*residu_calculus_fisher(lambdac,zetac,etac,i,c1,c2,p,psip,phi)/c2;
   end
int_test=int_test+a0
end

