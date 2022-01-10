function [OverProdSign] = OverProdSignFun(StockMat,E_CZ,VA_t,E_VA,OrderMat,RD,xx,F,S,R)
% Determine the sign or direction of temporal production expansion
% StockMat - amount of available intermediate input
% E_CZ - input coefficients to produce one unit of x
% VA_t - amount of available production factor input at time step t
% E_VA - input coefficients to produce one unit of x
% OrderMat -  order matrix
% RD - reconstruction demand
% xx - pre-disaster level of output
% F - number of industtial sectors
% S - number of sectors
% R - number of regions
Result = zeros(F,S*R);
ZX = StockMat./E_CZ; % S*(S*R)
ZX(isnan(ZX)) = 1e20;
VX = VA_t./E_VA;  % F*(S*R)
OX = sum([OrderMat,RD],2)';  % 1*(S*R)
ZOX = [ZX;OX];
VCU = repmat(VX(1,:),size(ZOX,1),1)<ZOX;
Result(1,sum(VCU)==size(ZOX,1)) = 1;
VCD = repmat(VX(1,:),size(ZOX,1),1)>ZOX;
Result(1,sum(VCD)>0) = -1;
VLU = repmat(VX(2,:),size(ZOX,1),1)<ZOX;
Result(2,sum(VLU) == size(ZOX,1)) = 1;
VLD = repmat(VX(2,:),size(ZOX,1),1)>ZOX;
Result(2,sum(VLD)>0) = -1;
temp = abs(repmat(OX,F,1)-VX)./xx;
OverProdSign = Result .* temp;  % F*(S*R)
OverProdSign(isnan(OverProdSign))=0;
end