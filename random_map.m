%%
% random map
% ------------------------------------------------------------------------------------------------
% new Ws and distrib here
N = 1e7;
bins = 100;
u = rand(1, N);
Finvu = 1 - (1 - u).^2;
w = Finvu;

% checking w distrib
figure;
subplot(1,2,1);
histogram(w,bins);
xlabel('W'); ylabel('count'); title('Histogram of constructed samples');
subplot(1,2,2);
xd = [0.0:0.01:1];
yd = 0.5./sqrt(1-xd);
plot(xd, yd);
xlabel('W'); title('Theoretical distribution');

% probabilities
W = w;
Wmax = zeros(1, 1e7+1);
Wcmax = zeros(1, 1e7+1);
Y = zeros(1, 1e7+1);
for n = 1:1e7,
    W(n) = w(n); % optional as defined in preallocation
    Wmax(n+1) = max([(1-W(n))*Wmax(n), W(n)]);
    Wcmax(n+1) = max([min([(1-W(n))*Wmax(n), W(n)]),(1-W(n))*Wcmax(n)]);
    Y(n+1) = (W(n))^2 + ((1 - W(n))^2)*Y(n);
end

figure;
subplot(1,3,1);
[cnt1, ctr1] = hist(Wmax, bins);
plot(ctr1, bins*cnt1/N); xlim([0 1]); ylim([0 4]);
xlabel('W: W_{max}'); ylabel('P_1(W)'); title('Max Valley Weight Distribution');
subplot(1,3,2);
[cnt2, ctr2] = hist(Wcmax, bins);
plot(ctr2, bins*cnt2/N); xlim([0 1]); ylim([0 4]);
xlabel('W: W_{max}^c'); ylabel('P_2(W)'); title('Second Max Valley Weight Distribution');
subplot(1,3,3);
[cnt3, ctr3] = hist(Y, bins);
plot(ctr3, bins*cnt3/N); xlim([0 1]); ylim([0 4]);
xlabel('Y'); ylabel('\Pi(Y)'); title('Overlap Distribution');
sgtitle(['Distributions for the RM model ', ' N:', num2str(N), ' bins:', num2str(bins)]);
