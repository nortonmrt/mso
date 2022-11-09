function [shape,nderiv] = shapeFunct_Truss(xi)
%shape é uma matriz com N1 e N2

shape = [1-xi;1+xi]*1/2;

nderiv = [-1;1]*1/2;
end

