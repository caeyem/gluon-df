%%
% magnetism
% --------------------------------------------------------------------------

k = 3;
N = 1e3;
s=zeros(1,N+1);
% stable core s(t) frac of spins tat do not deped on the initial condition
for t=1:N
for p = 0:k,
    s(t+1) = s(t+1) + (factorial(k)/(factorial(k-p)*factorial(p)))*((s(t))^(k-p))*(((1-s(t))^(p))*(2^(1-2^p)));
end
end
figure;
plot(s(1:200));
xlabel('t'); ylabel('s'); title(['Evolution of s for k = ', num2str(k)]);
figure;
histogram(s)
