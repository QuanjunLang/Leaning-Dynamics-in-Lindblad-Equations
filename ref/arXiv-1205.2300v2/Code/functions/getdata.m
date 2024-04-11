function [U,S] = getdata(n,r)
% 
% [U,S] = getdata(n,r)
% 
% Sample a density operator from the following distribution.  
% Take a Haar-random pure state on a composite n*r system then trace out the 
% r-dimensional environment.  This state is guaranteed to have rank at most r.
% See Bengtsson and Zyczkowski "Geometry of Quantum States" for details.
% 
% The output is given in a factorized form, so that 
%       rho = U*S*U'
% where U is the first r columns of a unitary matrix and S is a diagonal
% positive matrix with unit trace.
% 
% To enable (disable) complex number support, set the global 
% variable COMPLEX to 1 (0).
% 
% Steve Flammia, June 27, 2009

% choose unit normal compex numbers
A = randn(n,r)+1i*randn(n,r); 
[Q,R] = qr(A,0); % thin QR decomposition
S = R*R'; S = S/trace(S);
% now diagonalize S
[U,S] = eig(S);
U=Q*U; % and redefine U

% now we have the factorized eigendecomposition
%   rho = U*S*U'

end
