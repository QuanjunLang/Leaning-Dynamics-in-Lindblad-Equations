function [X_est, outputInfo] = ALS_total_N2(A, b, r, varargin)
%%
% This function learns the reshaped quantum operator from random initial
% states and random observables without assuming any structure. 

% Suppose X is a quantum operator, reshaped to a matrix with (N*N)x(N*N).
% It is the Choi matrix of the channel operator in the quantum process
% tomography. We decompose X to N x N submatrices of size N*N.
%
% We assume that X has rank r << N^2. 
%
%%%%%%%%%%%%%%% Learning Process: %%%%%%%%%%%%%
% 1. Reshape the initial states and observables by kronecker product and
% apply the matrix sensing framework.
% 
% 2. ALS 
%
%
%
% Inputs:
%   - A: Tensor of size (N, N, M), representing the M initial staes, of size N x N.
%   - b: Measurement result of X, from the observable E_11, (E_12+E_21, ..., E_1N+E_N1), (iE_12-iE_21, ..., iE_1N-iE_N1).
%   - r: Target rank for reconstruction
%   - varargin: Optional arguments (debugON, plotON, X_true, etc, see below.)
%
% Outputs:
%   - X_est: Reconstructed reshaped operator X in matrix form
%   - outputInfo: Structure containing reconstruction details (errors, time, etc.)
%
%
%
% copyright - Quanjun Lang, 2024









%%
M = sysInfo.M;
N = sysInfo.n;
all_O = observableInfo.O;
all_rho0 = observableInfo.rho0;

RE = trueInfo.RE_true;
RL = trueInfo.RL_true;

A = cell(M, 1);
for m = 1:M
    O = all_O(:, :, m);
    rho0 = all_rho0(:, :, m);
    A{m} = kron(conj(O), rho0);
end

A_mat = zeros(N^2, N^2, M);
for m = 1:M
    A_mat(:, :, m) = A{m};
end



% fprintf('Estimating E: ')
% b_E = squeeze(all_rho(sysInfo.channel_dt_rate+1, :))';
% r_E = trueInfo.rank_RE_true;
% [RE_est, outputInfo_E] = ALS(A_mat, b_E, r_E, 'X_true', RE, 'debugON', 1, 'operator_name', 'E');
% fprintf('Time: %.3f\n', outputInfo_E.time)


fprintf('Estimating L: ')
b_L = squeeze(all_rho(2, :) - all_rho(1, :))'/sysInfo.dt;
r_L = trueInfo.rank_RL_true;
[RL_est, outputInfo_L] = ALS(A_mat, b_L, r_L, 'X_true', RL, 'debugON', 1, 'operator_name', 'L');
fprintf('Time: %.3f\n', outputInfo_L.time)