clear all
close all
clc

flexibleData = load('CE.mat');

%% Initialisation
u = flexibleData.u;
y = flexibleData.y;

[N,M] = size(y);

Te = 0.015;

%% Create Data Object

data = iddata(y,u,Te);

[data_d,T] = detrend(data);

plot(data,data_d)


%% Order estimation with ARX
nc = 0; nd = 0; nf = 0; nk = 1;
thres = 0.002;
order_arx = 0;

figure
for n =1:10
    orders = [n n nk];
    SYS = arx(data_d, orders);
    stem(n,SYS.Report.Fit.LossFcn), hold on
    
    if(SYS.Report.Fit.LossFcn > thres)
        order_arx = order_arx + 1;
    end
end

order_arx

%% Validation with ARMAX

order_armax = 0;
figure
for n =1:8
    orders = [n n n nk];
    SYS = armax(data_d, orders);
    stem(n,SYS.Report.Fit.LossFcn), hold on
    
    if(SYS.Report.Fit.LossFcn > thres)
        order_armax = order_armax + 1;
    end
    
    h = figure
    iopzplot(SYS)
    saveas(h, 'CE3-report/images/'+int2str(n), 'png');

end
order_armax










