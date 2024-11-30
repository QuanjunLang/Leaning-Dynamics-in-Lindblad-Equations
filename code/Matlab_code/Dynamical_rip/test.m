clc
close all
clear all


%%
delta_old = 0.02;
m = 1;
n = 1;
A_norm = 1;
e = 10.^linspace(-10, 1, 5000);

top = 4*delta_old + 4*m*n^2*e*A_norm + m*n^2*e.^2;
bot = 4 + m*n^2*e.^2;

delta_new = top ./ bot;
%%
figure;
plot(log10(e),  log10(delta_new))




%%
x = e;
alpha = 0.2;
y = -x.^2 - 2*alpha*(delta_old - 1)*x + alpha;

plot(log10(x), y)
ylim([-3, 3])