function [poly] = derivee(my_cell,n,x,phi,psi,phip,psip,in)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%binome de newton;
% [m,o]=size(my_cell);
% A= 1;
% for i=1:o
%     A= A.*(x-my_cell{1,i}(1)).^my_cell{1,i}(2)./((x-my_cell{2,i}(1)).^my_cell{2,i}(2));
% end
if strcmp(in,'zeta')
    A=((phi(x)./psip(x))^my_cell{1,2}(2)*psip(x));
else
    A=((phip(x)./psi(x))^my_cell{1,2}(2)*psi(x));
end
if n==0
    poly=A;
else
    k_vec=0:n-1;
    poly=0;
for i=k_vec
    poly= poly+derivee(my_cell,i,x,phi,psi,phip,psip,in).*nchoosek(n-1,i)*derivee_poly(my_cell,n-i-1,x);
end
end
end

