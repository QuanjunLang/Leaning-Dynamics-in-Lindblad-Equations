clc
close all
clear all

addpath(genpath(fileparts(which(mfilename))));
%% generate data
I = system_settings();

%%
h_grid = I.h_grid;
dx = I.dx;
xgrid = I.xgrid;

tau = 0.6;
n = 15;



s = round(tau/dx);
vgrid = zeros(2*n, 1);
y = zeros(2*n, 1);
for v = 1:2*n
    vgrid(v) = xgrid((v-1)*s+1);
    y(v) = h_grid((v-1)*s+1);
end

figure;hold on;
plot(xgrid, h_grid);
plot(vgrid, y, 'o');

%% Prony's method
m = n;
F = zeros(2*n-m, n-m);
for i = 1:m
    F(:, i) = y(m-i+1:2*n-i);
end

Y = y(m+1:2*n);

P = F\Y;
plot(P)



