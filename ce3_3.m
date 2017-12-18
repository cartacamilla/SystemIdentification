clear all
close all
clc

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
Phi = toeplitz([0;-y(1:N-1)], [0; 0]);
Phi = [Phi, toeplitz([0;u(1:N-1)], [0; 0])];
% theta_hat = inv(Phi'*Phi)*Phi'*y;
theta_hat = Phi\y


%%

% reconstruct
yhat = Phi*theta_hat;

% compare
x_ax = Te*[1:length(yhat)];

figure
plot(x_ax,y);
hold on
plot(x_ax,yhat);
hold off

xlabel('Time[s]')
ylabel('Amplitude')
legend('Real output','ARX output')

% loss function
J_pred = sum((y - yhat).^2) %Predicition error



%% lsim
% G = z^-d*B(z^-1)/A(z^-1)
A = [1 theta_hat(1:n)'];
B = [theta_hat(n+1:end)'];
ARX_model = tf(B, A, Te);

yARX = lsim(ARX_model,u);

% compare
plot(x_ax,y,'b');
hold on
% plot(yhat,'r');
hold on
plot(x_ax,yARX,'r')
hold on
%% Instrumental Variable method
yM = yARX;

% least squares
Phi_iv = toeplitz([0;-yM(1:N-1)], [0; 0]);
Phi_iv = [Phi_iv, toeplitz([0;u(1:N-1)], [0; 0])];
%theta_hat_iv = Phi_iv\y;

theta_hat_iv = pinv(Phi_iv'*Phi)*Phi_iv'*y;

% reconstruct
yhat_iv = Phi_iv*theta_hat_iv;

J_pred_iv = sum((y-yhat_iv).^2) % error

%% lsim
% G = z^-d*B(z^-1)/A(z^-1)
A_iv = [1 theta_hat_iv(1:n)'];
B_iv = [theta_hat_iv(n+1:end)'];
IV_model = tf(B_iv, A_iv, Te);

yIV = lsim(IV_model,u);

%% compare
plot(yIV,'g'); 
hold off

legend('Real output','ARX output','IV output')













