function [R,h] = intcor(u,y)

N = length(u);
if length(y) <= N
    y = [y;zeros(N-length(y),1)];
else
    N = length(y);
    u = [u;zeros(N-length(u),1)];
end

R = [];
h = [];

for h0 = -N+1:N-1
    buff = 0;
    for k = 1:N
        buff = buff + 1/N * u(k)*y(mod(k-h0+N-1,N)+1);
    end
    R = [R;buff];
    h = [h;h0];
end
