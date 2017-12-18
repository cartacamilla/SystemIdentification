clear all
close all
clc

flexibleData = load('CE.mat');

%% Initialisation
u = flexibleData.u;
y = flexibleData.y;

[N,M] = size(y);

Te = 0.015;

%% Divide data
N = N/2;

data = iddata(y(1:N),u(1:N),Te);
valid = iddata(y(N+1:end),u(N+1:end),Te);

data = detrend(data);
valid = detrend(valid);

%%
na = 6;
nb = 6;
nk = 3;
nc = na;
nd = na;
nf = na;
nx = max(nb + nk, na);

figure
S_arx = arx(data, [na nb nk]);
subplot(3,2,1)
compare(valid, S_arx)
title('ARX')

S_iv4 = iv4(data, [na nb nk]);
subplot(3,2,2)
compare(valid, S_iv4)
title('Instrumental Variables')

S_armax = armax(data, [na nb nc nk]);
subplot(3,2,3)
compare(valid, S_armax)
title('ARMAX')

S_oe = oe(data, [nb nf nk]);
subplot(3,2,4)
compare(valid, S_oe)
title('Output-Error')

S_bj = bj(data, [nb nc nd nf nk]);
subplot(3,2,5)
compare(valid, S_bj)
title('Box-Jenkins')

S_n4sid = n4sid(data, nx);
subplot(3,2,6)
compare(valid, S_n4sid)
title('State-space')


printpdf(gcf, 'ce3_5_2_system_compare.pdf', 1, 1.5)

