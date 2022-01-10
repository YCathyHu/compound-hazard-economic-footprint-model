function [zz] = substi_inv(z_dis,order_c,I_sum)
% Allocate orders of unique products to regions involved
% I_sum - indicate the number of unique products and each product is provided by which regions 
zz = zeros(size(z_dis));
for i = 1:size(I_sum,1)
    zz(I_sum(i,:)==1,:) = z_dis(I_sum(i,:)==1,:).*order_c(i,:);
end
end