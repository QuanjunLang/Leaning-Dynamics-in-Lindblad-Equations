function [ x ] = Lasso(n,~,m,T,c,PAULI,data,~,~)
% FUNCTION Lasso %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
d = 2^n;

% estimator parameter
mu = 4*sqrt(m/eta);

% make noisy data by simulating measurements
s = Simulate_Measurement_Noise(data,eta); % This is now data+noise

cvx_begin sdp quiet
    variable x(d,d) hermitian;
    x == hermitian_semidefinite(d);
    minimize( 2*mu*trace(x) + pow_pos(norm(PAULI*vec(x)-s),2) );
    subject to
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
