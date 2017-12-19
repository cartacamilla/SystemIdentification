clear all
close all
clc

laserbeamdata = load('laserbeamdataN.mat');

y = laserbeamdata.y;
u = laserbeamdata.u;
Te = 1e-3; % sampling time

%%
N = length(y);
for r = 1:9; %optimizes error

    Y = [];
    U = [];

    for k=1:N-r
        Y = [Y, y(k:k+r)];
        U = [U, u(k:k+r)];
    end


    %% U orthogonal
    U_orth = eye(N-r,N-r) - U'*pinv(U*U')*U;

    Q = Y*U_orth;

    thres = 0.01;
    n = sum(svd(Q) > thres)

    figure 
    plot(svd(Q))
    xlabel('Order')
    ylabel('Value')

    Or = Q(1:r,1:n);

    C = Or(1,:);

    A = pinv(Or(1:end-1,:))*Or(2:end,:);

    %%

    z = tf('z',Te);

    F = C*inv((z*eye(size(A)))-A);

    U_f = [];

    for k = 1:size(F,2)
        U_f = [U_f lsim(F(k),u)];
    end
    
    %%
    B = pinv(U_f'*U_f)*U_f'*y;

    yhat = U_f*B;

    %loss function
    J(r) = sum((y-yhat).^2)/N

    %%
    figure
    plot(y), hold on
    plot(yhat)

    legend('real','approx')

end

figure
stem(J)

xlabel('r')
ylabel('Loss function')

%% Model
t = 0:Te:Te*(length(y)-1);
[tf_a, tf_b] = ss2tf(A, B, C, 0);
sys_ss = tf(tf_a, tf_b, Te, 'variable', 'z');
y_hat = lsim(sys_ss,u,t);

norm2 = norm(y_hat-y);
fprintf('Norm-2 of the error: %E\n', norm2);

% compare
figure
plot(t,y,'b');
hold on
% plot(yhat,'r');
hold on
plot(t,y_hat,'r')

axis([0 0.5 -0.004 0.005])

xlabel('Time[s]')
ylabel('Amplitude')
legend('Real output','ym for SSID model')

