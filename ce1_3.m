close all;
clc;

a = -0.2;
b = 0.2;
Te = 0.2;
N = 50/Te;

signal = a + (b-a).*rand(N,1);

simin = struct();
simin.signals = struct('values', signal);
simin.time = linspace(0,N*Te, N);

sim('ce1_1_sim')

