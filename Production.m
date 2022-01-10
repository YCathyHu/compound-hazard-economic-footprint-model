function [x] = Production(StockMat,E_CZ,VA_t,E_VA,OrderMat,RD,xx,S,mat_trivial)
% Calculate the actual production considering both capacity and orders
% StockMat - amount of available intermediate input
% E_CZ - input coefficients to produce one unit of x
% VA_t - amount of available production factor input at time step t
% E_VA - input coefficients to produce one unit of x
% OrderMat -  order matrix
% RD - reconstruction demand
% xx - pre-disaster level of output
% S - number of sectors
% mat_trivial - identification of trivial intermediate input
ZX = StockMat./E_CZ;
ZX(isnan(ZX)) = 1e20;
temp = 1.25 .* repmat(xx,S,1);
ZX(mat_trivial) = temp(mat_trivial);
VX = VA_t./E_VA;
OX = sum([OrderMat,RD],2)';
x = min([ZX;VX;OX],[],1);
end