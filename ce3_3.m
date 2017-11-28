laserbeamdata = load('laserbeamdataN.mat');
y = laserbeamdata.y;
u = laserbeamdata.u;
Te = 1e-3; % sampling time

N = length(y);

% ARX model: yhat(k) = -a1*y(k-1) -a2*y(k-2) +b1*u(k-2) +b1*u(k-2)
% A(z) = 1 + a1*z^-1 + a2*z^-2
% B(z) = b1*z^-1 + b2*z^-2
% G(q^-1) = (q^-d * B(q^-1) / A(1^-1))

% n=2, m=2, d=0
% phi(k) = [-y(k-1), ..., -y(k-n), u(k-d-1), ..., u(k-d-m)]'
% Phi = [phi(1); phi(2); ..., phi(N)]
% theta = [a1 a2 b1 b2]'
n=2;
m=2;
d=0;

% least squares
Phi = toeplitz([0;y(1:N-1)], [0; 0]);
Phi = [Phi, toeplitz([0;u(1:N-1)], [0; 0])];
theta = pinv(Phi'*Phi)*Phi'*y;

% reconstruct
yhat = Phi*theta;

% loss function
J = sum((y - yhat).^2);

% compare
plot(y);
hold on
plot(yhat);
hold off

%% lsim
% G = z^-d*B(z^-1)/A(z^-1)
A = [1; theta(1:n)];
B = [theta(n+1: end)];
model = tf(B, A, Te);

%% Instrumental Variable method
%todo
