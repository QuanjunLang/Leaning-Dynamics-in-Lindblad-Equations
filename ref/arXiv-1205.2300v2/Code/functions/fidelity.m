function F = fidelity(U,S,V,R,depol)
% F = fidelity(A,B)
% calculate the fidelity between two
% density operators A and B
% 
% F = fidelity(U,S,V,R)
% calculate the fidelity between the factorized density operators 
% A = U*S*U'; B = V*R*V';
% 
% F = fidelity(U,S,V,R,depol)
% calculate the fidelity between the factorized density operators 
% A = (1-depol)*U*S*U' + depol*I/d;   B = V*R*V';
% where
% I/d = eye(length(U))/length(U);
% is the maximally mixed state.
% Note that in this mode, the fidelity is NOT SYMMETRIC
% as calculated by this function.
% 
% We use the convention that we square the fidelity here!!!


if nargin == 2
    F = fullfidelity(U,S);
    return
elseif (nargin == 4 || nargin == 5)
    FACTORIZED = 1;
else
    F=0;
    fprintf('fidelity only takes 2, 4 or 5 arguments.');
    return
end

% 4 arguments is the same as zero depolarizing
if nargin == 4
    depol = 0;
end

% Factorized
if FACTORIZED
    Z = V'*U;
    % define the new "density operators", 
    % where B is subnormalized and lives on the support of A
    A = R; B = Z*S*Z'; 
    % add depolarizing noise, if necessary
    B = (1-depol)*B + depol*eye(length(R))/length(U);
    F = fullfidelity(A,B);
end

end

function F = fullfidelity(A,B)
% do the full fidelity calculation
    x = sqrtm(A);
    F = x*B*x';
    F = trace(sqrtm(F));
    F = real(F)^2; 
end
