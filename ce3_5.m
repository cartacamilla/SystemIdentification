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

figure
na = order_arx; nb = order_arx;
SYS = arx(data_d, [na nb 1]);
errorbar(SYS.b, SYS.db*2)

%% Validation with ARMAX
nc = 0; nd = 0; nf = 0; nk = 1;
order_armax = 0;
%figure
for n =1:10
    orders = [n n n nk];
    SYS = armax(data_d, orders);
    %stem(n,SYS.Report.Fit.LossFcn), hold on
    
    if(SYS.Report.Fit.LossFcn > thres)
        order_armax = order_armax + 1;
    end
end
order_armax

%% estimate delay nk
na = order_armax; nb = order_armax; nc = order_armax;
SYS = armax(data_d, [na nb nc 1]);
lower = SYS.b - 2*SYS.db;
upper = SYS.b + 2*SYS.db;
test = lower.*upper <= 0; % test if 0 within 2 sigma
nk = 0;
for i = test
    if i == 0
        break
    else
        nk = nk+1;
    end
end
nk
errorbar(SYS.b, SYS.db*2)

%% Plot Zero/Pole and their confidence interval
figure
h = iopzplot(SYS)
showConfidence(h,2)

%% Divide data
N1 = N/2;

data = iddata(y(1:N1),u(1:N1),Te);
valid = iddata(y(N1+1:end),u(N1+1:end),Te);

%% Compare

NN = struc(1:10,1:10,0:10);
V=arxstruc(data,valid,NN);
[NN, Vm] = selstruc(V, 'PLOT');
NN



