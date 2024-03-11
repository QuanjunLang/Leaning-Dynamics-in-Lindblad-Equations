function [] = testEstimatorParallel(n,r,iterations,c_vals,T_vals,m,depol,dataFunction,estimators)
% FUNCTION testEstimatorParallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ANNOYING GLOBAL INPUTS:
%   I hard-coded my directory and the number of cores on my computer.
%   Sorry about that.
% 
% INPUTS:
%   n = number of qubits
%   r = rank of generated test matrix
%   iterations = number of iterations to run CVX
%   c_vals = a vector (or single number) of the values of c to use
%   T_vals = a matrix (or single number) of the values of T to use
%   m = information about which measurement settings to use
%   dataFunction = function to be used to either get or generate data
%   estimator1,2,3 = which estimators to use
%
% OUTPUTS
%   None. Saves to file(s) in form of intersave_n_r_T_c_estimator.mat.
%   
%   This function requires the parallelization package
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 2^n;

% Choose a directory to save a file
%[fname,pname] = uiputfile('*.mat','Save file name');
%savepath = strcat(pname,fname);
pname = '/Users/sflammia/Dropbox/Tomo Code (Paulis)/Data/';

% if ( not(matlabpool('size') > 0) )
%     matlabpool(12); % 12 is the number of cores I have on my machine.
% end

% Used to determine if m is a single number or a series of them.
multi = false;
if numel(m) == 3
    m_start = m(1);
    m_finish = m(2);
    m_step = m(3);
    m_iter = ceil(1+(m_finish-m_start)/m_step);
    fprintf('Running from m = %d to %d in increments of %d for a total of %d iterations\n',...
        m_start,m_finish,m_step,m_iter);
    multi = true;
elseif numel(m) == 1
    m_start = m(1);
    m_finish = m(1);
    m_step = 1;
    m_iter = 1;
else
   fprintf('Error in inputs for m'); 
end

%Determine type of structure we are given for estimators
multiEstimators = true;
if isa(estimators,'function_handle')
   multiEstimators = false; 
end

for i = 1:size(T_vals,2)
    c = c_vals(i);
    for j = 1:size(T_vals,1)
        T = T_vals(j,i);
        k=1;
        while k <= numel(estimators)
            out.m = zeros(m_iter-1,1);
            out.Fid = zeros(m_iter-1,iterations);
            out.Tr = zeros(m_iter-1,iterations);
            out.FroDist = zeros(m_iter-1,iterations);
            
            if(multiEstimators)
                E = estimators{k};
            else
                E = estimators; 
            end
            m = m_start;
            while m <= m_finish
                % Create save location string
                fprintf('Running m = %d, c = %d, T = %d, type = %s\n',m,c,T,char(E));
                if(multi)
                    interSave = strcat(int2str(n),'_',int2str(r),'_',...
                     int2str(T),'_',int2str(c),'_',char(E),'_',int2str(m_start),...
                     '_',int2str(m_finish),'_',int2str(m_step),'.mat');
                else
                    interSave = strcat(int2str(n),'_',int2str(r),'_',...
                     int2str(T),'_',int2str(c),'_',char(E),'_',int2str(m_start),'.mat');
                end
                interSave = strcat(pname,interSave);
                
                % Call the estimator and save to out
                for loop = 1:iterations
                    [PAULI,data,A,D] = feval(dataFunction,d,r,m,depol);
                    x = feval(E,n,r,m,T,c,PAULI,data,A,D);
                    [U,S] = svd(x);
                    S = S/trace(S); % renormalize the trace
                    t = trdist(A,D,U,S);
                    f = fidelity(A,D,U,S);
                    fro_dist = norm(A*D*A'-U*S*U','fro');
                    
                    % The reason for this contorted way of defining things
                    % is so that the parfor loop doesn't complain.
                    Fidpara(loop) = f;
                    Trpara(loop) = t;
                    FroDistpara(loop) = fro_dist;
                    
                end 
                mval = ceil(1+(m-m_start)/m_step);
                out.m(mval,1) = m;
                out.Fid(mval,:) = Fidpara;
                out.Tr(mval,:) = Trpara;
                out.FroDist(mval,:) = FroDistpara;
                
                m = m + m_step;
            end
            % Check to see if we have run this test before. If we have we
            %  can just append the original data with the new data we have
            %  collected.
            try
                prevRun = load(interSave);
                out.Fid = [prevRun.out.Fid out.Fid];
                out.Tr = [prevRun.out.Tr out.Tr];
                out.FroDist = [prevRun.out.FroDist out.FroDist];
                save( interSave , 'out');
            catch
                save( interSave , 'out');
            end
            k = k+1;
        end
    end
end


end