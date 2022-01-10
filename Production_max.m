function [x_max] = Production_max(StockMat,E_CZ,VA_t,E_VA,xx,S,mat_trivial)
% Calculate the maximum production capacity with available intermediate stock and production factor input
% StockMat - amount of available intermediate input
% E_CZ - input coefficients to produce one unit of x
% VA_t - amount of available production factor input at time step t
% E_VA - input coefficients to produce one unit of x
% xx - pre-disaster level of output
% S - number of sectors
% mat_trivial - identification of trivial intermediate input
ZX = StockMat./E_CZ;
temp = 1.25 .* repmat(xx,S,1);
ZX(mat_trivial) = temp(mat_trivial);
VX = VA_t./E_VA;
x_max = min([ZX;VX],[],1);
end