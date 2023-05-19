
%%
% attractor random map
% --------------------------------------------------------------------------

N=1e7;
y=zeros(1,N+1);
y(1) = rand;
delta = 2;
for n=1:N,
    y(n+1) = W(n).^(delta) + (1-W(n).^(delta)).*y(n);
end
