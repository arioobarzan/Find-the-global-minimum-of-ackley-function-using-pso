clc;
clear;
close all;

%% set constant value
np = 200; % particle number
lb = -1; % lower boundary value
ub = 1; % upper boundary value

%% find the global minimum of the cost function
[xBest,yMin] = pso(@(x)ackley(x),np,lb,ub);

%% plot the result
x = -1:0.01:1;
y = ackley(x')';
figure(1);
hold on
plot(x,y);
h = plot(xBest,yMin,'gO','MarkerFaceColor','g','MarkerSize',14);
legend(h,'Global Minimum');
xlabel('x');
ylabel('y');
title('Find the global minimum of ackley function using pso');