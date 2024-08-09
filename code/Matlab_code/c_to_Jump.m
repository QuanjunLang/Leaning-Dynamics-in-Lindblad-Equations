function Jump = c_to_Jump(F, c)
    

[n, ~] = size(F{1});
Jump = zeros(n, n);
for i = 1:n^2-1
    Jump = Jump + c(i) * F{i};
end

end