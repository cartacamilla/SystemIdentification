function [R,h] = intcor(u,y)
% intercorrelation function
% supposes u and y to be of same length

N = length(u);
% R = [];
% h = [];
% 
% for h0 = -N+1:N-1
%     buff = 0;
%     for k = 1:N
%         buff = buff + 1/N * u(k)*y(mod(k-h0+N-1,N)+1);
%     end
%     R = [R;buff];
%     h = [h;h0];
% end


y=[y;y;y];
for h = 0:2*N-1
    R(h+1) = sum(u .* y(h+1:N+h));
end

R = R'./N;
