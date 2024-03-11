function [ X ] = MaxLikelihood(n,~,m,T,c,PAULI,data,~,~)
% FUNCTION MaxLikelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

p = (1+data)/2; % Probability of getting +1 on a measurement
s = rand(numel(data),eta); % Make some measurements
P = repmat(p,1,eta);
s = P > s;
s = sum(s,2); % Determine how many +1s we got
aPlus = s;
aMinus = eta - aPlus;

% Using the method and notation of arXiv:quant-ph/0611244

% make each measured Pauli into a projector, and stack them in rows (vec'd)
veye = repmat(vec(eye(d))',m,1);
Pip = (veye+PAULI)/2;
Pim = (veye-PAULI)/2;
Pi = [Pip;Pim]; % all the positive operators for our measurements

w = [aPlus;aMinus]; % weights (a function of the data; "f" in the paper)

TOL = 1e-7;
X = MLE(Pi,w,d,TOL);

end
