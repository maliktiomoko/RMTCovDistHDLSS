function [poly] = derivee_totale(my_cell,n,x,my_cell2,phi,psi,phip,psip,in)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%binome de newton;
k_vec=0:n;
poly=0;
for i=k_vec
    poly= poly+derivee(my_cell,i,x,phi,psi,phip,psip,in).*nchoosek(n,i)*derivee_poly(my_cell2,n-i,x);
end
end

