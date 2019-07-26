function [int] = residu_calculus_fisher(lambda,zeta,eta,n,c1,c2,p,psip,phi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rest=@(z) 1./(z-eta(1));
for i=1:length(zeta)
    rest=@(z) rest(z)+1./(z-eta(i+1))-1./(z-zeta(i));
end
if n==0
    int=0;
    for i=1:length(lambda)
        lambdac=lambda;lambdac(i)=[];
        A1=@(x) x;
        C1=@(x) 1/(x-eta(end));
        for j=1:length(lambdac)
            A1=@(x) A1(x).*(x-zeta(j))./(x-lambdac(j));
            A2=@(x) A1(x).*(x-zeta(j))./(x-lambdac(j));
            C1=@(x) C1(x)+1./(x-eta(j))-1./(x-zeta(j));
        end
        int=int+A1(lambda(i)).*C1(lambda(i))+...
          A1(eta(i))./(eta(i)-lambda(i));
    end
else
int=0;int2=0;
for i=1:length(zeta)
    zetar=zeta;etar=eta;
    zetar(i)=[];etar(i)=[];
    my_cell3=cell(2,length(lambda));
    my_cell4=cell(2,length(lambda)+1);
    my_cell1=cell(2,2*length(lambda)-1);
    my_cell1{1,1}=[0 1];
    my_cell1{1,length(eta)}=[eta(end-1) n];
    my_cell1{1,length(eta)+1}=[eta(end) n];
    my_cell1{2,length(lambda)-1}=[lambda(end-1) 1];
    my_cell1{2,length(lambda)}=[lambda(end) 1];
    my_cell1{2,end}=[0 0];
    my_cell3{2,length(eta)-1}=[0,0];
    my_cell3{2,length(eta)}=[0,0];
    my_cell3{1,length(eta)-1}=[eta(end-1),1];
    my_cell3{1,length(eta)}=[eta(end),1];
    my_cell4{1,length(lambda)+1}=[lambda(end) 1/(n-1)];
    my_cell4{1,length(lambda)}=[lambda(end-1) 1/(n-1)];
    my_cell4{2,length(eta)}=[eta(end) 1/(n-1)];
    my_cell4{2,length(eta)-1}=[eta(end-1) 1/(n-1)];
    my_cell4{2,end}=[0 0];
    my_cell4{1,1}=[0 -1/(n-1)];
    for j=1:length(zetar)
        my_cell1{1,j+1}=[eta(j) n];
        my_cell1{1,length(eta)+j+1}=[zetar(j) 1];
        my_cell1{2,j}=[lambda(j) 1];
        my_cell1{2,length(lambda)+j}=[zetar(j) n];
        my_cell3{1,j}=[eta(j) 1];
        my_cell3{2,j}=[zetar(j) 1];
        my_cell4{1,j+1}=[lambda(j) 1/(n-1)];
        my_cell4{2,j}=[eta(j) 1/(n-1)];
    end
    fun=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell3);
    %fun3=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell4,psip,phi,c1);
    fun_av=@(x) derivee(my_cell1,max([0 n-1]),x);
    %fun_av2=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell1);
    if n>=2
    int=int+(1/(factorial(max([0 n-2]))))*fun(zeta(i))*(1-c1)^n;
    end
    int=int-(1/(factorial(max([0 n-1]))))*fun_av(zeta(i))*(1-c1)^n;
    %test1=(1/(factorial(max([0 n-1]))))*fun_av(zeta(i));
    %test2=(1/(n-1))*(1/(factorial(max([0 n-2]))))*fun_av2(zeta(i));
    %int2=int2+(1/(factorial(max([0 n-2]))))*fun3(zeta(i))*(1-c1)^n
end
for i=1:length(lambda)
    lambdac=lambda;lambdac(i)=[];
    expre=@(z) z.*(z-eta(1)).^n;
    expre2=@(z) (1-c1)^n*(z-eta(1)).^n;
    expre3=@(z) z;
    for j=1:length(lambdac)
        expre=@(z) expre(z).*(z-eta(j+1)).^n./((z-zeta(j)).^(n-1)*(z-lambdac(j)));
        expre2=@(z) expre2(z).*(z-eta(j+1)).^n./((z-zeta(j)).^(n));
        expre3=@(z) expre3(z).*(z-zeta(j))./((z-lambdac(j)));
    end
    int=int+rest(lambda(i))*expre(lambda(i))*(1-c1)^n;
    %int2=int2+rest(lambda(i))*expre(lambda(i))*(1-c1)^n;
    ret4(i)=rest(lambda(i))*expre(lambda(i))*(1-c1)^n;
    %int2=int2+rest(lambda(i))*((-c1/c2)*lambda(i))^n*(c2/p)*lambda(i)
    %int2=int2+(1-(1/c1+1/c2-1)*(p))*((-c1/c2)*lambda(i))^n*(c2/p)
end
end
end

