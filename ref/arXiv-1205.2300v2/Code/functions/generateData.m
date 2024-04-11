function [PAULI,data,A,D] = generateData(d,r,m,depol)
% Generate a random state with depolarizing noise, together with a random
% set of Pauli expectation values.

global allP


n = log(d)/log(2);

% Compute all the Paulis.
% This takes a while, but you only have to do it once.
if (~exist('allP','var')) || (d^2~=length(allP))
    fprintf('Computing all the Pauli matrices on %d qubits.\n',n);
    fprintf('This is only done each time you redefine n.\n');
    allP = allpaulis(n);
end

%randn('state',2012); rand('state',2012); % Fix for debug
[A,D] = getdata(d,r);
rho = A*D*A';
% add depolarizing noise
rho = localdepol(rho,depol);
[A,D] = svd(rho);

% select m Paulis at random (without replacement!)
rp = randperm(d^2);
PAULI = allP(rp(1:m),:);

data = PAULI*vec(rho);
% even if measurements are complex, data should still be real
data = real(data);

end

