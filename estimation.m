u = FPdata.u;
y = FPdata.y;
r = FPdata.r;
N = length(y);

data = detrend(iddata(y,u,Te));
estim_data = detrend(iddata(y(1:N/2),u(1:N/2),Te));
valid_data = detrend(iddata(y(N/2+1:end),u(N/2+1:end),Te));

% ARX models for different orders n
%% na estimate
order = 2:10;
loss_func = [];
for na = order
    SYS = arx(data, [na na 0]);
    loss_func = [loss_func, SYS.EstimationInfo.LossFcn];
end

scaling_factor = 0.001;
penalty_BIC = scaling_factor*order*log(N)/N;
global_crit = loss_func + penalty_BIC;
[~,i] = min(global_crit);
na_est = order(i);

figure; hold on;
plot(order, loss_func, 'b');
plot(order, penalty_BIC, '--k')
plot(order, global_crit, '-*r')
hold off;
title('ARX loss function');
xlabel('Model order na');
ylabel('Loss function');
legend('Loss Function', 'Penalty term BIC', 'Global criterion');
na = na_est

%% nk estimate
nb = na; nk = 0;
SYS = arx(data, [na nb nk]);
figure;
errorbar(SYS.b, SYS.db*2) % plot B with a confidence interval of 2 sigma

lower = SYS.b + 2 * SYS.db;
higher = SYS.b - 2 * SYS.db;
nk_est = 0;
for i = 1:length(SYS.b)
    if(lower(i).*higher(i) <= 0)
        nk_est = nk_est + 1;
    else
        break; % if 0 in 2sigma confidence interval
    end
end
nk = nk_est

%% nb estimate
order = 1:10;
loss_func = [];
for nb = order
    orders = [na nb nk];
    SYS = arx(data, orders);
    loss_func = [loss_func, SYS.EstimationInfo.LossFcn];
end

scaling_factor = 1e-4;
penalty_BIC = scaling_factor*order*log(N)/N;
global_crit = loss_func + penalty_BIC;
[~,i] = min(global_crit);
nb_est = order(i);

figure; hold on;
plot(order, loss_func, 'b');
plot(order, penalty_BIC, '--k')
plot(order, global_crit, '-*r');
hold off;
title('ARX loss function');
xlabel('Model order nb'); ylabel('Loss function');
legend('Loss Function', 'Penalty term BIC', 'Global criterion');
nb = nb_est
n = max(na, nb+nk)

%% Validate order with ARMAX model by observing pole-zero-cancellation
figure
order=5:8;
for n=order
    SYS = armax(data, [n n n 0]);
    subplot(2,2,n-min(order)+1);
    h = iopzplot(SYS);
    showConfidence(h, 2);
    title('ARMAX n = '+string(n));
    axis([-1.2 3 -2 2])
end
printpdf(gcf, 'final-report/images/3_zero_pole_cancel.pdf', 1.5, 1.8)

%% compare with selstruc
%NN = struc(5:7,5:7,0:4);
%V = arxstruc(estim_data, valid_data, NN);
%order = selstruc(V)

