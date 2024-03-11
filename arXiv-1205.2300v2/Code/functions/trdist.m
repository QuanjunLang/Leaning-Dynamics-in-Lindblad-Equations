function F = trdist(U,S,V,R,depol)
% calculate the trace distance between
% two operators A and B

% F = trdist(A,B)
% calculate the trace distance between two
% operators A and B
% 
% F = trdist(U,S,V,R)
% calculate the trace distance between the factorized operators 
% A = U*S*U'; B = V*R*V';
% 
% F = trdist(U,S,V,R,depol)
% calculate the trace distance between the factorized operators 
% A = (1-depol)*U*S*U' + depol*I/d;   B = V*R*V';
% where
% I/d = eye(length(U))/length(U);
% is the maximally mixed state.
% Note that in this mode, the trace distance is NOT SYMMETRIC
% as calculated by this function.

if nargin == 2
    F = fulltrdist(U,S); % full
	return
elseif (nargin == 4 || nargin == 5)
    FACTORIZED = 1;
else
    F = 0;
    fprintf('trdist only takes 2, 4 or 5 arguments.');
    return
end

% 4 arguments is the same as zero depolarizing
if nargin == 4
    depol = 0;
end

% Factorized
if FACTORIZED
    W = horzcat(U,V);
    [Q,X] = qr(W,0); % Q is an orthonormal basis for the column space of W
    QU = (Q'*U);
    QV = (Q'*V);
    % isometry into subspace spanned by the support.
    A = QU*S*QU'; B = QV*R*QV';
    k = min(size(Q)); % dimension of joint support ( == rank(Q) )
    d = length(U);
    % add depolarizing noise, if necessary
    A = (1-depol)*A+depol*eye(k)/d;
    F = fulltrdist(A,B) + depol*(d-k)/d;
end

end

function F = fulltrdist(A,B)
% do the full trace distance calculation
F = sum(svd(A-B))/2;
end
