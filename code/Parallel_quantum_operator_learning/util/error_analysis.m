function outputInfo = error_analysis(outputInfo, X_est, X_true, varargin)


%% Parse inputs using inputParser
p = inputParser;

% Required inputs
addRequired(p, 'outputInfo');
addRequired(p, 'X_est');
addRequired(p, 'X_true');


% Optional arguments with default values
addOptional(p, 'X_est_blocks', 0);          % Estimated blocks
addOptional(p, 'X_true_sub_blocks', 0)      % True blocks
addOptional(p, 'compute_dnorm', 0)          % Diamond norm is expensive
addOptional(p, 'displayON', 0)          % Diamond norm is expensive

% Parse arguments
parse(p, outputInfo, X_est, X_true, varargin{:});

% Extract parsed inputs
X_est_blocks            = p.Results.X_est_blocks;
X_true_sub_blocks       = p.Results.X_true_sub_blocks;
compute_dnorm           = p.Results.compute_dnorm;
displayON               = p.Results.displayON;

if ismatrix(X_true)

    % Frobenius errors for the entire matrix
    error_fro_X = norm(X_est - X_true, 'fro');
    rel_error_fro_X = error_fro_X / norm(X_true, 'fro');

    outputInfo.error_fro_X = error_fro_X;
    outputInfo.rel_error_fro_X = rel_error_fro_X;


    % Nuclear norm error of reshaped matrix X
    error_nuc_X     = sum(svd(X_est - X_true));
    rel_error_nuc_X = error_nuc_X./sum(svd(X_true));

    outputInfo.error_nuc_X = error_nuc_X;
    outputInfo.rel_error_nuc_X = rel_error_nuc_X;


    % Nuclear norm error of operator
    error_nuc_op     = sum(svd(rearrangement_R(X_est) - rearrangement_R(X_true)));
    rel_error_nuc_op = error_nuc_op./sum(svd(rearrangement_R(X_true)));

    outputInfo.error_nuc_op = error_nuc_op;
    outputInfo.rel_error_nuc_op = rel_error_nuc_op;

    if compute_dnorm
        dnorm_tic = tic;
        if displayON
            fprintf('computing Diamond norm ...')
        end
        % Diamond norm error of Choi matrix
        Operator_est = rearrangement_R(X_est);
        Choi_est = RR(Operator_est);

        Operator_true = rearrangement_R(X_true);
        Choi_true = RR(Operator_true);

        error_diamond = dnorm(Choi_est - Choi_true);
        norm_diamond = dnorm(Choi_true);
        
        rel_error_diamond = error_diamond/norm_diamond;
        outputInfo.rel_error_diamond = rel_error_diamond;

        outputInfo.error_diamond = error_diamond;
       

        dnorm_time = toc(dnorm_tic);
        fprintf('finished in %.2f s\n', dnorm_time)
        outputInfo.dnorm_time = dnorm_time;
    end


    if (~isscalar(X_est_blocks)) && (~isscalar(X_true_sub_blocks))
        [NN, ~] = size(X_est);
        N = floor(sqrt(NN));
        % Frobenius norm for each submatrix
        error_fro_blocks = zeros(N, N);
        for i = 1:N
            for j = 1:N
                error_fro_blocks(i, j) = norm(X_est_blocks{i, j} - X_true_sub_blocks{i, j}, 'fro');
            end
        end

        % Frobenius norms of true blocks
        norm_blocks = zeros(N, N);
        for i = 1:N
            for j = 1:N
                norm_blocks(i, j) = norm(X_true_sub_blocks{i, j}, 'fro');
            end
        end

        % Relative Frobenius norm of each submatrix
        rel_error_fro_blocks = error_fro_blocks ./ norm_blocks;


        outputInfo.error_fro_blocks = error_fro_blocks;
        outputInfo.rel_error_fro_blocks = rel_error_fro_blocks;
        outputInfo.norm_blocks = norm_blocks;
    end



    if displayON
        fprintf('Relative Frobenius error \t= %.8f, \tRelative Nuclear error on Choi \t= %.8f, \tRelative Nuclear error on operator \t= %.8f \n',rel_error_fro_X, rel_error_nuc_X, rel_error_nuc_op);
        fprintf('Frobenius error \t\t= %.8f, \tNuclear error on Choi \t\t= %.8f, \tNuclear error on operator \t\t= %.8f \n',error_fro_X, error_nuc_X, error_nuc_op);
        if compute_dnorm
            fprintf('Relative Diamond norm error \t= %.8f, \nDiamond norm error \t\t= %.8f, ', rel_error_diamond, error_diamond);
            % fprintf('Diamond norm error \t\t= %.8f\n', error_diamond);
        end
    end

end

