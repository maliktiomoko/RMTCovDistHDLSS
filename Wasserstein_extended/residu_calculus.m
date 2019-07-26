function [int2] = residu_calculus(lambda,zeta,eta,n,phi,psi,phip,psip)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
int2=0;
for i=1:length(zeta)
    zetar=zeta;etar=eta;
    zetar(i)=[];etar(i)=[];
    my_cell5=cell(2,length(lambda)+1);
    my_cell5{1,length(lambda)-1}=[eta(end-1) 1/(n-1)];
    my_cell5{1,length(lambda)}=[eta(end) 1/(n-1)];
    my_cell5{2,length(lambda)}=[lambda(end-1) 1/(n-1)];
    my_cell5{2,length(lambda)+1}=[lambda(end) 1/(n-1)];
    my_cell5{2,1}=[0 1/(n-1)];my_cell5{1,end}=[0 0];
%     my_cell3=cell(2,2*length(lambda));
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
    
%     my_cell3{1,length(lambda)-1}=[lambda(end-1) 1];
%     my_cell3{1,length(lambda)}=[lambda(end) 1];
%     my_cell3{1,2*length(lambda)-1}=[lambda(end-1) 1];
%     my_cell3{1,2*length(lambda)}=[lambda(end) 1];
%     my_cell3{2,length(eta)-1}=[eta(end-1) 1];
%     my_cell3{2,length(eta)}=[eta(end) 1];
%     my_cell3{2,end-1}=[0 0];
%     my_cell3{2,end}=[0 0];
    my_cell1{1,1}=[0 1];
    my_cell1{2,end}=[0 0];
    for j=1:length(zetar)
        my_cell5{1,j}=[eta(j) 1/(n-1)];
        my_cell5{2,j+1}=[lambda(j) 1/(n-1)];
%         my_cell3{1,j}=[lambda(j) 1];
%         my_cell3{1,length(lambda)+j}=[lambda(j) 1];
%         my_cell3{2,j}=[eta(j) 1];
%         my_cell3{2,length(eta)+j}=[zetar(j) 1];
        my_cell1{1,j+1}=[lambda(j) n];
        my_cell1{1,length(lambda)+j+1}=[lambda(j) n];
        my_cell1{1,2*length(lambda)+j+1}=[zetar(j) 1];
        
        my_cell1{2,j}=[eta(j) n];
        my_cell1{2,length(eta)+j}=[lambda(j) 1];
        my_cell1{2,2*length(eta)+j}=[zetar(j) n];
        %my_cell1{2,end-1}=[0 0];
    end
    my_cell2=cell(2,3*length(lambda));
%     my_cell4=cell(2,2*length(lambda));
    my_cell6=cell(2,length(lambda));
    my_cell6{1,length(lambda)}=[lambda(end) 1/n];
    my_cell6{2,end}=[0 1/n];
%     my_cell4{1,length(lambda)}=[lambda(end) 1];
%     my_cell4{1,2*length(lambda)}=[lambda(end) 1];
    my_cell2{1,length(lambda)}=[lambda(end-1) n];
    my_cell2{1,length(lambda)+1}=[lambda(end) n];
    my_cell2{1,2*length(lambda)}=[lambda(end-1) n];
    my_cell2{1,2*length(lambda)+1}=[lambda(end) n];
    my_cell2{1,3*length(lambda)}=[zeta(end) 1];
    my_cell2{2,length(lambda)}=[lambda(end) 1];
%     my_cell4{2,end-1}=[0 0];
%     my_cell4{2,end}=[0 0];
    for j=1:length(etar)
        my_cell6{1,j}=[lambda(j) 1/n];
        my_cell6{2,j}=[zeta(j) 1/n];
%         my_cell4{1,j}=[lambda(j) 1];
%         my_cell4{1,length(lambda)+j}=[lambda(j) 1];
%         my_cell4{2,j}=[etar(j) 1];
%         my_cell4{2,length(etar)+j}=[zeta(j) 1];
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
%     fun=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell3);
%     fun2=@(x) derivee_totale(my_cell2,max([0 n-1]),x,my_cell4);
    fun3=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell5,phi,psi,phip,psip,'zeta');
    fun4=@(x) derivee_totale(my_cell2,max([0 n-1]),x,my_cell6,phi,psi,phip,psip,'eta');
%     master_zeta=@(z) 2./(z-lambda(1))-1./(z-eta(1));
%     master_zeta=@(z) master_zeta(z)+2./(z-lambda(end))-1./(z-eta(end));
%     for k=1:length(zetar)
%     master_zeta=@(z) master_zeta(z)+2./(z-lambda(k+1))-1./(z-eta(k+1))-1./(z-zetar(k));
%     end
%     master_eta=@(z) 2./(z-lambda(1));
%     for k=1:length(etar)
%     master_eta=@(z) master_eta(z) +2./(z-lambda(k+1))-1./(z-etar(k))-1./(z-zeta(k));
%     end
%     fun_av=@(x) derivee(my_cell1,max([0 n-1]),x);
%     fun2_av=@(x) derivee(my_cell2,max([0 n]),x);
%     fun_av_bis=@(x) derivee_totale(my_cell1,max([0 n-2]),x,my_cell1);
%     fun2_av_bis=@(x) derivee_totale(my_cell2,max([0 n-1]),x,my_cell2);
%     if n>=2
%         int=int+(1/(factorial(max([0 n-2]))))*fun(zeta(i));
%     end
%     int=int+...
%     1/(factorial(n-1))*fun2(eta(i));
%     int=int-(1/(factorial(max([0 n-1]))))*fun_av(zeta(i))-...
%     1/(factorial(n))*fun2_av_bis(eta(i));
    if n>=2
        int2=int2+(1/(factorial(max([0 n-2]))))*fun3(zeta(i)); 
        test1(i)=(1/(factorial(max([0 n-2]))))*fun3(zeta(i));
    else
        int2=int2-derivee(my_cell1,0,zeta(i),phi,psi,phip,psip,'zeta');
    end
    test2(i)=(1/(factorial(max([0 n-1]))))*fun4(eta(i));
    int2=int2+(1/(factorial(max([0 n-1]))))*fun4(eta(i));
end
my_cell6=cell(2,length(lambda));
my_cell6{1,length(lambda)}=[lambda(end) 1/n];
my_cell6{2,end}=[0 1/n];
% my_cell4=cell(2,2*length(lambda));
% my_cell4{1,length(lambda)}=[lambda(end) 1];
% my_cell4{1,2*length(lambda)}=[lambda(end) 1];
my_cell2=cell(2,3*length(lambda));
my_cell2{1,length(lambda)}=[lambda(end-1) n];
my_cell2{1,length(lambda)+1}=[lambda(end) n];
my_cell2{1,2*length(lambda)}=[lambda(end-1) n];
my_cell2{1,2*length(lambda)+1}=[lambda(end) n];
my_cell2{1,3*length(lambda)}=[zeta(end) 1];
my_cell2{2,length(lambda)}=[lambda(end) 1];
% my_cell4{2,end-1}=[0 0];
% my_cell4{2,end}=[0 0];
etar=eta;etar(end)=[];
for j=1:length(etar)
    my_cell6{1,j}=[lambda(j) 1/n];
    my_cell6{2,j}=[zeta(j) 1/n];
%     my_cell4{1,j}=[lambda(j) 1];
%     my_cell4{1,length(lambda)+j}=[lambda(j) 1];
%     my_cell4{2,j}=[etar(j) 1];
%     my_cell4{2,length(etar)+j}=[zeta(j) 1];
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
% fun2=@(x) derivee_totale(my_cell2,n-1,x,my_cell4);
% fun2_av=@(x) derivee(my_cell2,n,x);
fun5=@(x) derivee_totale(my_cell2,max([0 n-1]),x,my_cell6,phi,psi,phip,psip,'eta');
% master_eta=@(z) 2./(z-lambda(1));
%     for k=1:length(etar)
%     master_eta=@(z) master_eta(z)+2./(z-lambda(k+1))-1./(z-etar(k))-1./(z-zeta(k));
%     end
% int=int+1/(factorial(n-1))*fun2(eta(end));
% int=int-1/(factorial(n))*fun2_av(eta(end));
int2=int2+1/(factorial(n-1))*fun5(eta(end));
test_last=1/(factorial(n-1))*fun5(eta(end));
r=1;


end

