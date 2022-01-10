function [z_dis] = substi(zz,zz_c,I_sum)
% Calculate the regional contribution of unique products
% I_sum - indicate the number of unique products and each product is provided by which regions 
z_dis = zeros(size(zz));
for i = 1:size(I_sum,1)
    z_dis(I_sum(i,:)==1,:) = zz(I_sum(i,:)==1,:)./zz_c(i,:);
end
z_dis(isnan(z_dis)) = 0;
end