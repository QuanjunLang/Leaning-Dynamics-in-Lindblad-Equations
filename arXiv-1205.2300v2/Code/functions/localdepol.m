function rho = localdepol(M,p)
% take an initial state M of n qubits and apply the local depolarizing
% channel with local noise strength p, i.e., the qubit channel
%       E(M) = (1-p) M + p I/2
% where I/2 is the maximally mixed state.

s = size(M); n = log2(s(1));
errorcheck(s,n);

rho = reshape(M,[s(1)*s(2),1]); % vec rho
for j=1:n,
    rho = depol(j,n,p)*rho; % depolarize the jth qubit
end
rho = reshape(rho,[s(1),s(2)]); % unvec rho

end

function errorcheck(s,n)

if s(1)~=s(2)
    error('Input should be square.');
end

if mod(n,1)~=0
    error('Size of input should be a power of 2.');
end

end

function chan = depol(j,n,p)
% local depolarizing noise on the jth qubit out of n, with strength p
q = p/4;

% define the Paulis
X = sparse([0,1;1,0]); Y = sparse([0,-i;i,0]); Z = sparse([1,0;0,-1]);

% put each Pauli in the jth slot
X = kron(speye(2^(j-1)),X); X = kron(X,speye(2^(n-j)));
Y = kron(speye(2^(j-1)),Y); Y = kron(Y,speye(2^(n-j)));
Z = kron(speye(2^(j-1)),Z); Z = kron(Z,speye(2^(n-j)));

% the depolarizing channel on the jth qubit
chan = (1-3*q)*speye(4^n) + q*( kron(X,X) + kron(Y,conj(Y)) + kron(Z,Z) );

end
