function [out] = derivee_poly(my_cell,k,x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[m,o]=size(my_cell);
out=0;
for i=1:o
    out=out+(-1)^(k).*(factorial(k)).*my_cell{1,i}(2)./((x-my_cell{1,i}(1)).^(k+1));
    out=out-(-1)^(k).*(factorial(k)).*my_cell{2,i}(2)./((x-my_cell{2,i}(1)).^(k+1));
end
end

