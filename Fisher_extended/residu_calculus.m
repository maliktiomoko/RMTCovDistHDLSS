function [int] = residu_calculus(lambda,zeta,eta,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if n==0
    int=0;
    for i=1:length(lambda)
        lambdac=lambda;lambdac(i)=[];
        A1=@(x) x;
        B1=@(x) 1./x;
        C1=@(x) -1./(x-eta(end));
        for j=1:length(lambdac)
            A1=@(x) A1(x).*(x-zeta(j))./(x-lambdac(j));
            B1=@(x) B1(x)+1./(x-zeta(j))-1./(x-lambdac(j));
            C1=@(x) C1(x)+2./(x-lambdac(j))-1./(x-zeta(j))-1./(x-eta(j));
        end
        int=int+A1(lambda(i)).*C1(lambda(i))+2*A1(lambda(i))*B1(lambda(i))-...
            A1(eta(i))./(eta(i)-lambda(i));
    end
else
int=0;
for i=1:length(zeta)
    zetar=zeta;etar=eta;
    zetar(i)=[];etar(i)=[];
    my_cell3=cell(2,2*length(lambda));
    my_cell1=cell(2,3*length(lambda)-1);
    my_cell1{1,length(lambda)}=[lambda(end-1) n];
    my_cell1{1,length(lambda)+1}=[lambda(end) n];
    my_cell1{1,2*length(lambda)}=[lambda(end-1) n];
    my_cell1{1,2*length(lambda)+1}=[lambda(end) n];
    %my_cell1{1,3*length(lambda)}=[zeta(end) 1];
    my_cell1{2,length(lambda)}=[eta(end) n];
    my_cell1{2,length(lambda)-1}=[eta(end-1) n];
    my_cell1{2,2*length(lambda)}=[lambda(end) 1];
    my_cell1{2,2*length(lambda)-1}=[lambda(end-1) 1];
    
    my_cell3{1,length(lambda)-1}=[lambda(end-1) 1];
    my_cell3{1,length(lambda)}=[lambda(end) 1];
    my_cell3{1,2*length(lambda)-1}=[lambda(end-1) 1];
    my_cell3{1,2*length(lambda)}=[lambda(end) 1];
    my_cell3{2,length(eta)-1}=[eta(end-1) 1];
    my_cell3{2,length(eta)}=[eta(end) 1];
    my_cell3{2,end-1}=[0 0];
    my_cell3{2,end}=[0 0];
    my_cell1{1,1}=[0 1];
    my_cell1{2,end}=[0 0];
    for j=1:length(zetar)
        my_cell3{1,j}=[lambda(j) 1];
        my_cell3{1,length(lambda)+j}=[lambda(j) 1];
        my_cell3{2,j}=[eta(j) 1];
        my_cell3{2,length(eta)+j}=[zetar(j) 1];
        my_cell1{1,j+1}=[lambda(j) n];
        my_cell1{1,length(lambda)+j+1}=[lambda(j) n];
        my_cell1{1,2*length(lambda)+j+1}=[zetar(j) 1];
        
        my_cell1{2,j}=[eta(j) n];
        my_cell1{2,length(eta)+j}=[lambda(j) 1];
        my_cell1{2,2*length(eta)+j}=[zetar(j) n];
        %my_cell1{2,end-1}=[0 0];
    end
    my_cell2=cell(2,3*length(lambda));
    my_cell4=cell(2,2*length(lambda));
    my_cell4{1,length(lambda)}=[lambda(end) 1];
    my_cell4{1,2*length(lambda)}=[lambda(end) 1];
    my_cell2{1,length(lambda)}=[lambda(end-1) n];
    my_cell2{1,length(lambda)+1}=[lambda(end) n];
    my_cell2{1,2*length(lambda)}=[lambda(end-1) n];
    my_cell2{1,2*length(lambda)+1}=[lambda(end) n];
    my_cell2{1,3*length(lambda)}=[zeta(end) 1];
    my_cell2{2,length(lambda)}=[lambda(end) 1];
    my_cell4{2,end-1}=[0 0];
    my_cell4{2,end}=[0 0];
    for j=1:length(etar)
        my_cell4{1,j}=[lambda(j) 1];
        my_cell4{1,length(lambda)+j}=[lambda(j) 1];
        my_cell4{2,j}=[etar(j) 1];
        my_cell4{2,length(etar)+j}=[zeta(j) 1];
        my_cell2{1,1}=[0 1];
        my_cell2{1,j+1}=[lambda(j) n];
        my_cell2{1,length(lambda)+j+1}=[lambda(j) n];
        my_cell2{1,2*length(lambda)+j+1}=[zeta(j) 1];
        
        my_cell2{2,j}=[lambda(j) 1];
        my_cell2{2,length(lambda)+j}=[etar(j) n];
        my_cell2{2,length(etar)+length(lambda)+j}=[zeta(j) n];
        my_cell2{2,end}=[0 0];
        my_cell2{2,end-1}=[0 0];   
    end
    fun=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell3);
    fun2=@(x) derivee_totale(my_cell2,max([0 n-1]),x,my_cell4);
    master_zeta=@(z) 2./(z-lambda(1))-1./(z-eta(1));
    master_zeta=@(z) master_zeta(z)+2./(z-lambda(end))-1./(z-eta(end));
    for k=1:length(zetar)
    master_zeta=@(z) master_zeta(z)+2./(z-lambda(k+1))-1./(z-eta(k+1))-1./(z-zetar(k));
    end
    master_eta=@(z) 2./(z-lambda(1));
    for k=1:length(etar)
    master_eta=@(z) master_eta(z) +2./(z-lambda(k+1))-1./(z-etar(k))-1./(z-zeta(k));
    end
    fun_av=@(x) derivee(my_cell1,max([0 n-1]),x);
    fun2_av=@(x) derivee(my_cell2,max([0 n]),x);
    if n>=2
        int=int+(1/(factorial(max([0 n-2]))))*fun(zeta(i));
    end
    int=int+...
    1/(factorial(n-1))*fun2(eta(i));
    int=int-(1/(factorial(max([0 n-1]))))*fun_av(zeta(i))-...
    1/(factorial(n))*fun2_av(eta(i));
    %int=int+1/(factorial(n-1))*fun2(eta(i))

    fun2(eta(i))
end
my_cell4=cell(2,2*length(lambda));
my_cell4{1,length(lambda)}=[lambda(end) 1];
my_cell4{1,2*length(lambda)}=[lambda(end) 1];
my_cell2=cell(2,3*length(lambda));
my_cell2{1,length(lambda)}=[lambda(end-1) n];
my_cell2{1,length(lambda)+1}=[lambda(end) n];
my_cell2{1,2*length(lambda)}=[lambda(end-1) n];
my_cell2{1,2*length(lambda)+1}=[lambda(end) n];
my_cell2{1,3*length(lambda)}=[zeta(end) 1];
my_cell2{2,length(lambda)}=[lambda(end) 1];
my_cell4{2,end-1}=[0 0];
my_cell4{2,end}=[0 0];
etar=eta;etar(end)=[];
for j=1:length(etar)
    my_cell4{1,j}=[lambda(j) 1];
    my_cell4{1,length(lambda)+j}=[lambda(j) 1];
    my_cell4{2,j}=[etar(j) 1];
    my_cell4{2,length(etar)+j}=[zeta(j) 1];
    my_cell2{1,1}=[0 1];
    my_cell2{1,j+1}=[lambda(j) n];
    my_cell2{1,length(lambda)+j+1}=[lambda(j) n];
    my_cell2{1,2*length(lambda)+j+1}=[zeta(j) 1];

    my_cell2{2,j}=[lambda(j) 1];
    my_cell2{2,length(lambda)+j}=[etar(j) n];
    my_cell2{2,length(etar)+length(lambda)+j}=[zeta(j) n];
    my_cell2{2,end}=[0 0];
    my_cell2{2,end-1}=[0 0];   
end
fun2=@(x) derivee_totale(my_cell2,n-1,x,my_cell4);
fun2_av=@(x) derivee(my_cell2,n,x);
master_eta=@(z) 2./(z-lambda(1));
    for k=1:length(etar)
    master_eta=@(z) master_eta(z)+2./(z-lambda(k+1))-1./(z-etar(k))-1./(z-zeta(k));
    end
int=int+1/(factorial(n-1))*fun2(eta(end));
int=int-1/(factorial(n))*fun2_av(eta(end));
end
end

