%% Add code path to the environment
folder = fileparts(which(mfilename));
addpath(genpath(folder));
clear folder

% if ~exist('chebfun')
%     fprintf('Please download chebfun for improved performance\n');
% end

%% Please set user specified data saving folder here
% USERNAME = char(java.lang.System.getProperty('user.name'));
% switch USERNAME
%     case 'langquanjun'
%         sysInfo.data_saving_folder = '/Users/langquanjun/QuanjunFei_data/Partial_Observation/data';
%     otherwise
%         SAVE_DIR = [getenv('HOME'),'/IPS_Graph_data/'];
%         if ~exist(SAVE_DIR,'dir'); mkdir(SAVE_DIR); end
%         sysInfo.data_saving_folder = SAVE_DIR;
% end
% 
% assert(isfield(sysInfo,'data_saving_folder'), 'Please specify the data saving folder in system_settings.m')


