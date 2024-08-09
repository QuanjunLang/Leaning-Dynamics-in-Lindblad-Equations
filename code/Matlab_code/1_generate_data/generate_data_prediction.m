function pred_traj = generate_data_prediction(tgrid, Hamiltonian, C, varargin)
% Use Python package to generate trajectories

%% Input parser
p = inputParser;
addRequired(p, 'tgrid');
addRequired(p, 'Hamiltonian');
addRequired(p, 'Jump');
addOptional(p, 'plotON', 0);

parse(p, tgrid, Hamiltonian, C, varargin{:});

%% Load parameters

M           = 10;
p           = length(C);
n           = length(Hamiltonian);
%% Generate data using python packages
Jump        = zeros(n, n, p);
for i = 1:p
    Jump(:, :, i) = C{i};
end

a = pyrunfile("Lindblad_validate.py", 'a', Hamiltonian = Hamiltonian, Jump = Jump, tgrid = tgrid, M = M);



%% Load data from Python code
result_temp     = cell(a);

all_rho_temp        = double(result_temp{1});
pred_traj             = permute(all_rho_temp, [3, 4, 2, 1]);

end


