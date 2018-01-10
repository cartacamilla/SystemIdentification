function parametric_models(y,u,Te,na,nb,nk)

[N,M] = size(y);

train = 1:N/2;
test = N/2+1:N;
train_data = iddata(y(train),u(train),Te);
test_data = iddata(y(test),u(test),Te);

train_data = detrend(train_data);
test_data = detrend(test_data);

nc = na; nd = na; nf = na; nx = max(nb + nk, na);

figure
S_arx = arx(train_data, [na nb nk]);
subplot(3,2,1)
compare(test_data, S_arx)
title('ARX')

S_iv4 = iv4(train_data, [na nb nk]);
subplot(3,2,2)
compare(test_data, S_iv4)
title('Instrumental Variables')

S_armax = armax(train_data, [na nb nc nk]);
subplot(3,2,3)
compare(test_data, S_armax)
title('ARMAX')

S_oe = oe(train_data, [nb nf nk]);
subplot(3,2,4)
compare(test_data, S_oe)
title('Output-Error')

S_bj = bj(train_data, [nb nc nd nf nk]);
subplot(3,2,5)
compare(test_data, S_bj)
title('Box-Jenkins')

S_n4sid = n4sid(train_data, nx);
subplot(3,2,6)
compare(test_data, S_n4sid)
title('State-space')


printpdf(gcf, 'final-report/images/system_compare.pdf', 1, 1.8)

%%

figure
S_arx = arx(train_data, [na nb nk]);
subplot(3,1,1)
resid(test_data, S_arx)
title('ARX')

S_iv4 = iv4(train_data, [na nb nk]);
subplot(3,1,2)
resid(test_data, S_iv4)
title('Instrumental Variables')

S_armax = armax(train_data, [na nb nc nk]);
subplot(3,1,3)
resid(test_data, S_armax)
title('ARMAX')

printpdf(gcf, 'final-report/images/system_resid1.pdf', 1, 1.8)


figure
S_oe = oe(train_data, [nb nf nk]);
subplot(3,1,1)
resid(test_data, S_oe)
title('Output-Error')

S_bj = bj(train_data, [nb nc nd nf nk]);
subplot(3,1,2)
resid(test_data, S_bj)
title('Box-Jenkins')

S_n4sid = n4sid(train_data, nx);
subplot(3,1,3)
resid(test_data, S_n4sid)
title('State-space')

printpdf(gcf, 'final-report/images/system_resid2.pdf', 1, 1.8)






















end