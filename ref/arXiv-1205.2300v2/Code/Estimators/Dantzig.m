function [ x ] = Dantzig(n,~,m,T,c,PAULI,data,~,~)
% FUNCTION Dantzig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       This takes as input everything possible. This should be done for
%       any user-created estimator, as it preserves flexibility in creating
%       multiple estimators that may use different parameters.
%
%
%   OUTPUTS:
%       X, the reconstructed state
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = floor(T-m*c);
eta = floor(N/m);
t = eta*m;
d = 2^n;

% estimator parameter
lambda = 3*d/sqrt(t);

% make noisy data by simulating measurements
s = Simulate_Measurement_Noise(data,eta);

ss = reshape(PAULI'*s,d,d); % data+noise in matrix form

% cvx sdp quiet
% cvx sdp
cvx_begin sdp
    variable x(d,d) hermitian;
    x == hermitian_semidefinite(d);
    xx = reshape(PAULI'*PAULI*vec(x),d,d);
    
    minimize(trace(x));
    subject to
    (d/m)*norm(xx-ss) <= lambda;
    x >= 0;
cvx_end


end


function s = Simulate_Measurement_Noise(data,eta)

% make noisy data by simulating measurements
p = (1+data)/2; % Probability of getting +1 on a measurement
s = rand(numel(data),eta); % Make some measurements
P = repmat(p,1,eta);
s = P > s;
s = sum(s,2); % Determine how many +1s we got
s = s/eta;
s = 2*s-1; % This is now data+noise

end
