clear all
clc
close all
p=32;c1=1.5;c2=2;
seed=167;
rng(167)
for j=1:6
    p=2*p;
    j
    n1=floor(p/c1);n2=floor(p/c2);
    C1=toeplitz(0.0.^(0:p-1));
    
    %C2=toeplitz(0.5.^(0:p-1));
    N2=2*p;
    W=zeros(p,N2);
    for i=1:N2
        W(:,i)=mvnrnd(zeros(p,1),eye(p));
    end
    C2=(W*W')/((1/p)*trace(W*W'));
    lambdavr=eig(C1*C2);
    n = 20;
    lower = min(lambdavr);
    upper = 100;
    %x = exp(linspace(log(lower),log(upper),100));
    x=logspace(log(min(lambdavr))/log(10),log(50)/log(10),100);
    %x=linspace(lower,upper,100);
    y=sqrt(x);
    a=fliplr(polyfit(x,y,n));
    a0=a(1);a(1)=[];
    g=@(x) a0;
    for i=1:length(a)
        g=@(x) g(x)+a(i)*x.^i;
    end
    n_simulation=10;
for i=1:n_simulation
    X=zeros(p,n1);
    Y=zeros(p,n2);
for k=1:n1
    %k
    X(:,k) = mvnrnd(zeros(1,p),C1);
end
for k=1:n2
    %k
    Y(:,k) = mvnrnd(zeros(1,p),C2);
end
%Proposed Wasserstein distance

[est(i),esthat(i)] = RMTWassDist(X,Y,a,a0);
%[est(i),esthat(i),estvr(i),est_vrai(i)] = RMTFisherDist(X,Y,C1,C2);
%Real Wasserstein distance
est_vrai(i)=(1/p)*trace(C1)+(1/p)*trace(C2)-2*trace((C1^(1/2)*C2*C1^(1/2))^(1/2))/p;
estvr(i)=(1/p)*trace(C1)+(1/p)*trace(C2)-2*mean(g(lambdavr));
end
est_vrai_mean(j)=mean(est_vrai)
estvr_mean(j)=mean(estvr)
esthat_mean(j)=mean(esthat)
est_mean(j)=mean(est)
end
% plot(abs(int_mean-est_vrai_mean))
% hold on
% plot(abs(int_mean-esthat_mean))