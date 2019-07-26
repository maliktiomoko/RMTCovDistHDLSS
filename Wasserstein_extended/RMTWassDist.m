function [int,esthat] = RMTWassDist(X,Y,a,a0,n1,n2)
%Function that compute the Wasserstein distance between Gaussian centered
%distribution based on the article Random  Matrix-Improved  Estimation  of  the  Wasserstein  Distance
%between  two  Centered  Gaussian  Distribution (Malik TIOMOKO & Romain Couillet)
%Input Need the samples from the first class X of dimension p*n and the
%samples from the second class Y of size p*n
%Return the estimate est proposed in the article and the classical esthat
%Define the dimensions
p=size(X,1);
%n1=size(X,2);
%n2=size(Y,2);
c1=p/n1;c2=p/n2;
%Sample covariance estimate
hatC1=X*X'/n1;hatC2=Y*Y'/n2;
lambda=sort(eig(hatC1*hatC2));
m=@(z) mean(1./(lambda-z));
mp=@(z) mean(1./(lambda-z).^2);
phi=@(z) z./(1-c1-c1.*z.*m(z));
phip=@(z) (1-c1+c1*z.^2.*mp(z))./((1-c1-c1*z.*m(z)).^2);
psi=@(z) 1-c2-c2*z.*m(z);
psip=@(z) -c2*m(z)-c2*z.*mp(z);
%f=@(z) sqrt(z);
eta=sort(real(eig(diag(lambda)-(1/n1)*sqrt(lambda)*sqrt(lambda)')));
zeta=sort(real(eig(diag(lambda)-(1/n2)*sqrt(lambda)*sqrt(lambda)')));
% phi_test=@(z) z;
% psi_test=@(z) 1;
% phipsi=@(z) sqrt(z)/(c2);
% for i=1:length(lambda)
%     phi_test=@(z) phi_test(z).*((z-lambda(i))./(z-eta(i)));
%     psi_test=@(z) psi_test(z).*(z-zeta(i))./(z-lambda(i));
%     phipsi=@(z) phipsi(z).*sqrt((z-zeta(i))./(z-eta(i)));
% end
% Distinguish the case where n1<n2 to the case where n1>n2
% if eta(1)<zeta(1)
%     my_eta=zeta;
%     my_zeta=eta;
% else
%     my_zeta=zeta;
%     my_eta=eta;
% end
% other=@(z) 2*sum(1./(z-zeta))-2*sum(1./(z-lambda));
% integrand_real=@(z) (1/(2*pi))*2*f(-(phi(z)./psi(z))).*other(z).*(psi(z)/c2);
% %Computing the second term (real_integral)
% real_integral=0;
% % for i=1:length(my_zeta)
% %     real_integral=real_integral+integral(integrand_real,my_zeta(i),my_eta(i));
% % end
% %Computing the first term (pole in lambda)
% pole=2*(sqrt(c2/c1))*sum(sqrt(lambda))/c2;
% esty=pole+real_integral;
% est=(1/p)*trace(hatC1+hatC2)-2*esty;
% %Distinguish the case n1=n2
% if n1==n2
%     est=(1/p)*trace(hatC1+hatC2)-2*(sum(sqrt(lambda))-sum(sqrt(zeta)))*(2*n1/p);
% end
% %Classical estimate
esthat=(1/p)*trace(hatC1)+(1/p)*trace(hatC2)-2*(1/p)*(trace((hatC1^(1/2)*hatC2*hatC1^(1/2))^(1/2)));
%a{1}=@(z)1;
%b{1}=@(z)1;
% for i=2:10
%     a{i}=@(z)1/2*(a{i-1}(z)+z./a{i-1}(z));
%     %b{i}=@(z) b{i-1}(z)-1+z.*exp(-b{i-1}(z));
% end
maxV = max(lambda)*1.2;
lambdac=lambda(lambda>0.001);
etac=eta(eta>0.001);
zetac=zeta(zeta>0.001);
% choose=1e-5:1e-4:etac(1);
%                      [maximum,maxi]=max(phi(choose)./psi(choose));
%                      ind_max=choose(maxi);
minV = -1;altitude=0.5;step=1e-3;
contour = fliplr([altitude*1i+minV+(0:step:maxV-minV),altitude*1i+maxV+(-step:-step:-2)*1i*altitude,-1i*altitude+maxV+(-step:-step:-maxV+minV),minV-1i*altitude+(0:step:2)*1i*altitude]);        
%g=@(z) 0.00051*z.^3-0.022*z.^2+0.43*z+0.56;
%g=@(z) -5.3*1e-5*z.^4+0.0026*z.^3-0.049*z.^2+0.55*z+0.44;
%g=@(z) a(4)*z.^4+a(3)*z.^3+a(2)*z.^2+a(1)*z;
%g=@(z) a0;
%g=@(x) 6.4*10^(-6)*x.^5-0.00037*x.^4+0.0083*x.^3-0.092*x.^2+0.67*x+0.36;
%g=@(x) a0;
g=@(x) a0;
 for i=1:length(a)
     g=@(x) g(x)+a(i).*x.^(i);
 end
%g=@(x) a(6)*x.^6+a(5)*x.^5+a(4)*x.^4+a(3)*x.^3+a(2)*x.^2+a(1)*x+a0;
   integrand = @(z) g(phi(z)./psi(z)).*(psi(z)/c2).*(phip(z)./phi(z)-psip(z)./psi(z));
   int2=real(trapz(contour,(1/(2*pi*1j))*integrand(contour)))-g(0)*(1-c2)/c2;
%  int_test=0;
%  phip_inv=@(z) -phi(z).^2./phip(z);
%   for i=10
%       i
%       int_test=int_test+a(i)*residu_calculus(lambdac,zetac,etac,i,phi,psi,phip_inv,psip)/c2;
%   end
% int_test=int_test+a0;
%int=(1/p)*trace(hatC1)+(1/p)*trace(hatC2)-2*int_test;
int=(1/p)*trace(hatC1)+(1/p)*trace(hatC2)-2*int2;
 
end

