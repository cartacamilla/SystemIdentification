clear all
close all
clc

flexibleData = load('CE.mat');

%% Initialisation
u = flexibleData.u;
y = flexibleData.y;

[N,M] = size(y);

Te = 0.015;

%% Divide data
N = N/2;

data = iddata(y(1:N),u(1:N),Te);
valid = iddata(y(N+1:end),u(N+1:end),Te);

