clear all
clc
close all
p=64;c1=0.75;c2=8;
for j=1
    %p=2*p;
    %j
    n1=floor(p/c1);n2=floor(p/c2);
    C1=toeplitz(0.2.^(0:p-1));
    C2=toeplitz(0.4.^(0:p-1));
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

[int] = RMTFisherDist(X,Y);
%[est(i),esthat(i),estvr(i),est_vrai(i)] = RMTFisherDist(X,Y,C1,C2);
%Real Wasserstein distance
est_vrai(i)=log((C1^(1/2)*C2*C1^(1/2))^(1/2))/p;
end
est_vrai_mean(j)=mean(est_vrai)
%estvr_mean=mean(estvr)
esthat_mean(j)=mean(esthat)
est_mean=mean(est)
end
plot(abs(int_mean-est_vrai_mean))
hold on
plot(abs(int_mean-esthat_mean))