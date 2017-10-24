% input signal
saturation = 0.5;
noiseVariance = 0.1;

u = 0.5*prbs(7,16);
Te = 0.2; % sample time
N = length(u);
sim_time = N*Te;

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')

plot(simout)
%%

