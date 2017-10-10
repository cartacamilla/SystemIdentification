%x = ones(5,1);
%[R,h] = intcor(x,x);
%comment

close all;
clc;

x = prbs(4,3);
[R,h] = intcor(x,x);

figure 
grid on
stem(h,R, 'o')

title('')
xlabel('')
ylabel('')

