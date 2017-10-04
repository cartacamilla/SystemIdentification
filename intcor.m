function [R,h] = intcor(u,y)
%if size(u)[1] > size(y)

N = size(u,1);

R = [];
h = [];

for d = -N+1:-1
    x = 0;
    for n = 1:N+d
        x = x + 1/N * u(n)*y(n - d);
    end
    R = [R; x];
    h = [h; d];
end

for d = 0:N-1
    x = 0;
    for n = 1:N-d
        x = x + 1/N * u(n)*y(n + d);
    end
    R = [R; x];
    h = [h; d];
end
