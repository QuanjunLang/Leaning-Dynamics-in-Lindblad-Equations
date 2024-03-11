function MakeGraphs(key)

% FUNCTION MakeGraphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS
%       key = a num with the common value of the plot that you want 
%             alternatively, it should be a string 
%             of the form 'Da', 'La', or 'Ma'
%             (enough to uniquely identify an estimator).
%   OUTPUTS
%       none. Produces a graph which is good enough for testing purposes.
%
%   Written: Steve Flammia 4/12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set y-axis limits for both plots.
set_Fid_y_lim = [0.48 1.0];
set_Tr_y_lim = [0 .61];

files = Get_Filenames(key);

% how many plots will we merge, and what are their parameters?
L = length(files);
n=cell(1,L); T=n; c=n; est=n;
[n,T,c,est] = Parse_Filename(files);

% Note: either T or est is constant, but not both!
% Also: n and c should be constant.

% make the title
titlestr = Make_Title(key,T{1},c{1});

% load this first so we can use it as axes lables
mat = load(files{1});
M = mat.out.m;

% set limits on the x axis.
if strcmp(n{1},'4'); setxlim = [24 262]; else setxlim = [20 1060]; end

% determine properties of the main figure (Fidelity)
close all;
fig1 = figure(1);
set(fig1,'Position',[1 1 850 470]);

axes1 = axes('ylim',set_Fid_y_lim,'xlim',setxlim,...
    'XTick',M,'XTickLabel',M,'FontSize',12);
xlabel(axes1,'Number of sampled Pauli operators');
ylabel(axes1,'Fidelity');
title(titlestr);
% tighten bounding box.
tight = get(gca,'TightInset')*[-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1];
display(get(gca,'Position'));
display(get(gca,'OuterPosition'));
tight = get(gca,'Position')+tight;
set(gca,'OuterPosition',tight);
display(get(gca,'Position'));
display(get(gca,'OuterPosition'));


hold('all');
box on;

% determine properties of the inset figure (Tr Distance)
axes2 = axes('pos',[0.32 0.17 0.48 0.31],'XTick',M,...
  'xlim',setxlim,'Color','none','ylim',set_Tr_y_lim);
ylabel(axes2,'Trace Distance');
box('on');
hold('all');

s = length(M);
F = zeros(L,s); dF=F; Tr=F; dTr=F; h = zeros(1,2*L);

for j=1:L
    mat = load(files{j});
    F(j,:) = mean(mat.out.Fid,2);
    dF(j,:) = std(mat.out.Fid,1,2);
    Tr(j,:) = mean(mat.out.Tr,2);
    dTr(j,:) = std(mat.out.Tr,1,2);
    
    h(j) = errorbar(M,F(j,:),dF(j,:),'Parent',axes1);
    h(j+L) = errorbar(M,Tr(j,:),dTr(j,:),'Parent',axes2);
    hold('all');
end


Make_Legend(key,L,h(1:L),T,est);

hold('off');

end


function files = Get_Filenames(key)

    % Convert key to a string, if necessary
    if ~ischar(key); key=num2str(key); end

    % Find all .mat files in the "Data" folder
    %filenames = what('Data/4_qubit_data');
    filenames = what('Data');
    filenames = filenames.mat;

    % find those filenames which match key
    ind = strfind(filenames,key);
    ind = cellfun(@length,ind);
    ind = find(ind);
    
    files = filenames(ind);
end

function [n,T,c,est] = Parse_Filename(files)

L = length(files);

n=cell(1,L); T=n; c=n; est=n;
for j=1:L
    data = textscan(files{j},'%s %s %s %s %s %s %s %s','Delimiter','_');
    [n{j},T{j},c{j},est{j}] = data{[1 3 4 5]};
end

end

function titlestr = Make_Title(key,T,c)
% if the key is numeric, then plot estimators at fixed T;
% if it's a string, plot a fixed estimator vs. T
titlestr = 'Fidelity and Trace Distance for ';
if isnumeric(key) 
    titlestr = strcat(titlestr,'T= ',T,' and c=',c);
elseif ischar(key)
    if ~isempty(strfind('Dantzig',key))
        x=' the Matrix Dantzig Selector';
    elseif ~isempty(strfind('MaxLik',key))
        x=' Maximum Likelihood';
    elseif ~isempty(strfind('Lasso',key))
        x=' the Matrix Lasso';
    else
        error('If key is a string, it must be Dantzig, MaxLik, or Lasso.');
    end    
    titlestr = strcat(titlestr,x,' vs. T at c=',c);
else
    error('key must be a num or a str.');
end
    
end

function Make_Legend(key,L,h,T,est)
% This still needs some work! 

% define the quantity we are varying (T or est)
if ischar(key);
	x = T;
elseif isnumeric(key);
    x = est;
else
    error('key should be either a number or a string.');
end

for j=1:L; x(j) = x{j}; end

leg = legend(h,x);
set(leg,'pos',[.18 .165 .01 .04],'FontSize',9);

end