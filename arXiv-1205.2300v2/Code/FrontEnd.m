%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This is a front end for calling the Compressed Sensing Tomography
%   estimators code. The following parameters should always be user
%   specified:
%       n = The number of qubits
%       r = The rank of the test matrix
%       iterations = The number of iterations the test is run
%      [c_vals] = A vector specifying the values of the time required to
%           change measurement settings (c). It should be of the following
%           form:
%
%           c1  c2 ... cn
%
%           To use only one time just enter the number.
%      [T_vals] = A matrix specifying the values of the total experiment time
%           length (T). It should be of the following form
%
%           T11 T21... Tn1
%           .           .
%           .           .  
%           .           .
%           T1m .......Tmn
%
%           This is intended to allow for flexibility in use. To use only
%           one time just enter the value.
%       [m] = To use a single value of m, just set this to a number. To use
%               multiple values, let [m] = [m_start m_finish m_step] where
%               this dictates start finish and step size for m
%       estimator = A function handle or an array of function handles to
%               estimator functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all


global COMPLEX
global allP
addPaths
COMPLEX = true;
% cvx_quiet(true);

n = 5;
d = 2^n;
r = 1;
depol = .01;
iterations = 120;
c_vals = [20];
T_vals = [41000; 80000; 270000];
m = [2*r*d d^2 2*r*d];

estimator =  {@Dantzig,@Lasso,@MaxLikelihood};

% This gives us a function to give to the test code to generate fake states 
% for testing. This could be replaced with any function that gets data from 
% a real experiment
dataFunction = @generateData;

% To run multiple estimators add in each at the end of this call. For example:
% runCVX(n,r,k,c,T,m_start,m_finish,m_step,depol,'LeastSquares','MaxLikelihood'); 
tic; !date
testEstimatorParallel(n,r,iterations,c_vals,T_vals,m,depol,dataFunction,estimator);
toc; !date