function PAULI = allpaulis(n)
% make the vec version of all the Paulis
% by chosing the m rows at random, 
% we can sample Paulis without replacement
%
% I changed the code slightly. I believe there was an error that results in
% computing the rows of PAULI incorrectly needed a sign flip. This should
% work appropriately now, but the changes are highlighted so we can revert
% to the previous version if necessary

d=2^n; m = d^2;

% X=1, Y=2, Z=3, I=4;

PAULI = sparse([],[],[],d^2,m,m*d );
fprintf('Taking measurements...      ');
for i=1:m
    fprintf('\b\b\b\b\b\b%5.1f%%', 100*i/m );
    list = int2base(i,4,n); list(find(list == 0)) = 4;
    E_i = explicitPauliTensor( list );    
    %E_i = E_i.'; %%?? could be used instead of the congugate transpose
    PAULI(:,i) = E_i(:);  % access via column is MUCH faster
end
fprintf('\n');
%PAULI = PAULI.' originally 
PAULI = PAULI';
% transpose it, since we were implicitly dealing with transpose 
% earlier (since updating the columns of a sparse matrix is 
% much more efficient than updating the rows )   

end