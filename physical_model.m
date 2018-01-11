%% double integrator model with two resonance frequencies
w1 = 14.1;
w2 = 26.6;
s = tf('s');
k = 10;
G = k * 1/s^2 * 1/(((s/w1)^2+1)*((s/w2)^2+1));
%c2d(G, Te)
freq = (0:1450).*(omega_s/N);
figure
subplot(1,2,1)
bode(spectral_analysis_model); hold on;
bode(G, freq); hold off;
%% ... model with lowpass instead of doule integrator
wc = 1.71;
k = 3.16;
G = k * 1/(s/wc + 1)^2 * 1/(((s/w1)^2+1)*((s/w2)^2+1));
%c2d(G, Te)
%figure
subplot(1,2,2)
bode(spectral_analysis_model)
hold on;
bode(G, freq)
hold off;
%% Analysis of physical system
syms x1 x2 x3 J1 J2 J3 k1 k2 k3 T s
A=[
    (J1*s^2+k1), -k1, 0;
    -k1, (J2*s^2+k1+k2), -k2;
    0, -k2, (J3*s^2+k2)
];
B=[T;0;0];
x=[x1;x2;x3];
res = solve(A*x-B, x);
G = res.x3/T;

% G has the structure K * 1/s^2 * 1/((s/w1)^2+1) * 1/((s/w2)^2+1);
% J1 = 1; J2 = 1; J3 = 1; k1 = 1; k2 = 1; K3 = 1;
% [n,d] = numden(G);
% n = sym2poly(subs(n));
% d = sym2poly(subs(d));
% G = tf(n,d);
%bode(G)
%% state space model
syms x1 x2 x3 dx1 dx2 dx3 J1 J2 J3 k1 k2 k3 km u
A=[zeros(3), eye(3);
    [
        -k1/J1, k1/J1, 0;
        k1/J2, -(k1+k2)/J2, k2/J2;
        0, k2/J3, -k2/J3
    ], zeros(3)
];
B = [0 0 0 km 0 0]';
C = [0 0 1 0 0 0];
D = zeros(1,1);
x = [x1; x2; x3; dx1; dx2; dx3];
