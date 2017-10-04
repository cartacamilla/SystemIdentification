%x = ones(5,1);
%[R,h] = intcor(x,x);

x = prbs(4,1);
[R,h] = intcor(x,x);
plot(h,R)
