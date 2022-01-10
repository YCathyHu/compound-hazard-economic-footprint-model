function [x_sorted] = resultsort3(x,Scenarios,Flood_set_int)
x_sorted = zeros(size(Scenarios,1)*size(Flood_set_int,1),1);
for i = 1:size(Flood_set_int,1)
    for j = 1:size(Scenarios,1)
        x_sorted(j+(i-1)*size(Scenarios,1)) = x(i+(j-1)*size(Flood_set_int,1));
    end
end
end