close all;
%clc;

a = -0.4;
b = -a;
Te = 0.2;
sim_time = 100; % seconds
N = sim_time/Te;

u = a + (b-a).*rand(N,1);

simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);

sim('ce1_1_sim')
% plot(simout)

U = toeplitz(u, zeros(N,1));
Y = simout.Data;

k = 300;
Uk = U(:,1:k);

theta = pinv(Uk)*Y;

%theta = U\Y;
plot(simin.time(1:k), theta);


G = tf([4],[1 1 4]);
Z = c2d(Te*G, Te, 'zoh');

hold on;
impulse(Z);
