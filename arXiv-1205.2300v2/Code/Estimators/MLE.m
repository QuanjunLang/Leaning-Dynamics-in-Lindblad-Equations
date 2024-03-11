function rho = MLE(Pi,w,d,TOL)
% FUNCTION MLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       Pi = The projectors from the measurements
%       w = The weights for each outcome (depends on the data)
%       d = the Hilbert space dimension
%       TOL = a tolerance for exiting the optimization
%
%   OUTPUTS
%       rho = The maximum likelihood estimate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MIGHT NEED TO CONJ THIS!!!
Gp = pinv(mat(sum(conj(Pi)))); % pseudoinverse of G

rho = eye(d)/d; % initialize rho as the maximally mixed state
pr = conj(Pi)*vec(rho); % initialize the probabilites
llf(1) = L(pr,w); % compute the initial log-likelihood function

j=2;
rho = step(Pi,pr,w,rho,Gp); % next rho
pr = conj(Pi)*vec(rho); % compute the new probabilites
llf(j) = L(pr,w); % compute the log-likelihood function (llf)

while llf(j-1)-llf(j)>TOL
    j=j+1;
    rho = step(Pi,pr,w,rho,Gp); % next rho
    pr = conj(Pi)*vec(rho); % compute the new probabilites
    llf(j) = L(pr,w); % compute the log-likelihood function (llf)
end

%plot(llf);
rho=conj(rho);
%display(rho);


end

function rhonew = step(Pi,pr,w,rho,Gp)
% FUNCTION step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       Pi = The projectors from the measurements
%       pr = The probabilities from the measurements (depends on the state)
%       w = The weights for each outcome (depends on the data)
%       rho = the current state
%       Gp = a normalizing factor (matrix)
%
%   OUTPUTS
%       rhonew = The next iteration of the algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter so that we don't get zero probabilities
c=w./pr; c(isnan(c))=0; c(isinf(c))=0;% replace NaNs and Infs with 0.
R = Gp*mat(c'*Pi)/sum(w);
rhonew = R*rho*R'; % iterate rho and normalize
rhonew = rhonew/trace(rhonew);

% delta=norm(rho-rhonew,'fro');

end

function llf = L(pr,w)
% FUNCTION L %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%       pr = The probabilities from the measurements(depends on the state)
%       w = The weights for each outcome (depends on the data)
%
%   OUTPUTS
%       llf = The (negative) log-likelihood function of pr(X)|w
%
%   NOTES
%       Since we work with the negative version, we want to *minimize*.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pr(pr==0)=1; % eliminate impossible outcomes
llf  = -w'*log(pr)/sum(w);
llf = real(llf);

end