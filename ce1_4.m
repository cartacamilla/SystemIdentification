close all;
hold off;

% input signal
Uprbs = prbs(4,7);

% sample time
Te = 0.2;

sim_time = 100; % seconds
N = length(Uprbs);

simin = struct();
simin.signals = struct('values', Uprbs);
simin.time = linspace(0,N*Te, N);

sim('ce1_1_sim')
plot(simout)

