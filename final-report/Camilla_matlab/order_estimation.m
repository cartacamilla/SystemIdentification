function order_estimation(y,u,Te,freq_response)

N = length(u);

data = iddata(y(1:N/2),u(1:N/2),Te);
valid = iddata(y(N/2+1:end),u(N/2+1:end),Te);

data = detrend(data);
valid = detrend(valid);

%% Loss function method
nb_order = 25;
loss_function = zeros(nb_order,1);
scaling_factor = 0.001;
penalty_BIC = scaling_factor*(1:nb_order)'*log(N)/N;


for n = 1:nb_order
    nk = 1; na = n; nb = n;
    orders = [na nb nk];
    SYS = arx(data, orders);
    
    loss_function(n) = SYS.Report.Fit.LossFcn;
    
end

global_crit = loss_function + penalty_BIC;
[~,na_est] = min(global_crit)

fig1 = figure(1);
plot(1:nb_order, loss_function(), 'b'); hold on;
plot(1:nb_order, penalty_BIC, '--k'); hold on;
plot(1:nb_order, global_crit, '-*r'); hold off;
title('Loss function vs model order');
ylabel('Loss function value'); xlabel('Model order');
legend('Loss Function', 'Penalty term BIC', 'Global criterion');
axis([1 nb_order 0 2e-4]);

saveas(fig1, 'final-report/images/1_loss_funct_order', 'png');


%% Zero/pôle cancellation method using ARMAX
nb_order = na_est+5;

figure
for n = 1:nb_order
    na = n; nb = n; nc = n; nk = 1;
    orders = [na, nb, nc, nk];
    
    SYS = armax(data, orders);
    
    fig = figure(n)
    h = iopzplot(SYS);
    title('n = '+string(n))
    showConfidence(h,2)
    saveas(fig, sprintf('final-report/images/1_zpc_%d',n), 'png');

end
hold off;

na_est = 6




%% Estimation delay

n = na_est;
na = 1; nb = n; nc = n; 

SYS = arx(data,[na nb  1])

lower = SYS.b + 2 * SYS.db;
higher = SYS.b - 2 * SYS.db;

fig3 = figure(3);
errorbar(SYS.b, SYS.db*2)

title('Estimation delay ');
ylabel('B coefficient value'); xlabel('Coefficient order');

saveas(fig3, 'final-report/images/1_est_delay_order', 'png');

nk_est = 0;
for i = 1:length(SYS.b)
    if(lower(i).*higher(i) <= 0)
        nk_est = nk_est + 1;
    else
        i
        break;
    end
end

nk_est


%% nb estimation
nk = nk_est; na = na_est;
scaling_factor = 1e-4;
penalty_BIC = scaling_factor*(1:nb_order)'*log(N)/N;
loss_function = [];

for n = 1:nb_order
    nb = n;
    orders = [na nb nk];
    SYS = arx(data, orders);
    
    loss_function(n) = SYS.Report.Fit.LossFcn;
    
end

global_crit = loss_function' + penalty_BIC;
[~,nb_est] = min(global_crit)

fig4 = figure(4);
plot(1:nb_order, loss_function(), 'b'); hold on;
plot(1:nb_order, penalty_BIC, '--k'); hold on;
plot(1:nb_order, global_crit, '-*r'); hold off;
title('Loss function vs model order');
ylabel('Loss function value'); xlabel('Model order');
legend('Loss Function', 'Penalty term BIC', 'Global criterion');
axis([1 nb_order 0 1e-5]);

saveas(fig4, 'final-report/images/2_loss_funct_order', 'png');

%% Zero/pôle cancellation method for nb
nb_order = nb_est+5;
na = na_est; nk = nk_est;
figure
for n = 1:nb_order
    nb = n; nc = n; 
    orders = [na, nb, nc, nk];
    
    SYS = armax(data, orders);
    
    fig = figure(n)
    h = iopzplot(SYS);
    title('n = '+string(n))
    showConfidence(h,2)
    %saveas(fig, sprintf('final-report/images/2_zpc_%d',n), 'png');

end
hold off;

%% Matlab structure
na = 5:6; nb = 1:4; nk = 1:3;

M_SYS = struc(na, nb, nk);

NN = selstruc(arxstruc(data, data, M_SYS));


%% Frequency response
na = na_est; nb = nb_est; nc = na_est; nk = nk_est

fig5 = figure(5);

SYS = armax(data,[na nb nc 1]);
bode(SYS); hold on;
%bode(freq_response); 
hold off;

legend('Identified model using ARMAX', 'Location','southwest');

saveas(fig5, 'final-report/images/1_freq_resp_order', 'png');

end