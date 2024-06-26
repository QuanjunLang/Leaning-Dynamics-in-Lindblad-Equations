function result = combine_Lindbladian(all_C, rho)


result = zeros(size(rho));
N = length(all_C);
for i = 1:N
    C = all_C{i};
    temp = 2*C*rho*C' - rho*(C'*C) - (C'*C)*rho;
    result = result + temp;
end
result = result/2;


end