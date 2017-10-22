close all;
%clc;

Te = 0.2;
sim_time = 70;
N = sim_time/Te;

% random input signal
a = -0.4; % scaling to avoid saturation
b = -a;
u = a + (b-a).*rand(N,1);

% simulation
simin = struct();
simin.signals = struct('values', u);
simin.time = linspace(0,N*Te, N);
sim('ce1_1_sim')

% deconvolution
U = toeplitz(u, [u(1);zeros(N-1,1)]);
k = 300;
Uk = U(:,1:k); % truncate U

theta = pinv(Uk)*simout;

%theta = U\Y;
plot(simin.time(1:k), theta);


G = tf([4],[1 1 4]);
Z = c2d(Te*G, Te, 'zoh');

hold on;
impulse(Z);
