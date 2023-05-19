%%
% mean field theory of spin glasses
%---------------------------------------------------------------------------------
% verification of intermediate distributions

% alt calc
yList = [0.7 0.7 0.9 0.3 ];
NList = [1e7 3.25e6 1e7 1e7 ];
binsList = [1e2 1e2 1e2 1e3 ];
idx = 4;
y = yList(idx);
N = NList(idx);
bins = binsList(idx);
xn = [0.0:1/N:1.0];
alpha = y;
% beta = 1-y;
% probabilities
Wmax = zeros(1, N+1);
Wcmax = zeros(1, N+1);
rhonx3 = zeros(1, N);
% rhonx3 = betarnd(alpha, beta, [1 N]);
W = zeros(1, N+1);
Y2 = zeros(1, N+1);
Y3 = zeros(1, N+1);
Y4 = zeros(1, N+1);
R = zeros(1, N+1);
Wfac = 1;
S = 1;
% try to gen using W distrib
for n = 1:N,
   beta = n*(1-y);
%     rhonx3(n) = betarnd(alpha, beta);
    rhonx3(n) = betarnd(alpha, (n*(1-y)));
%     W(n) = rhonx3(n);
    W(n) = Wfac*rhonx3(n);
    Wfac = Wfac - W(n);
    R(n) = (gamma(n-n*y+y)/(gamma(y)*gamma(n-n*y)))*((W(n)/(1-S))^(y-1))*((1 - W(n)/(1-S))^(n-n*y-1))*(1/(1-S));
    S = S - W(n);

%     Wmax(n+1) = max([Wmax(n), W(n)]);
%     Wcmax(n+1) = max([min([(1-W(n))*Wmax(n), W(n)]),(1-W(n))*Wcmax(n)]);

    [Wmax(n+1), max_i] = max([Wmax(n), W(n)]);
    if max_i == 1
        Wcmax(n+1) = max([Wcmax(n), W(n)]);
    else
        Wcmax(n+1) = Wmax(n);
    end

    Y2(n+1) = (W(n))^2 + ((1 - W(n))^2)*Y2(n);
    Y3(n+1) = (W(n))^20 + ((1 - W(n))^20)*Y3(n);
    Y4(n+1) = sum(W(1:n).^2);
end
Y = Y4;
figure;
subplot(1,3,1);
[cnt1, ctr1] = hist(Wmax, bins);
plot(ctr1, bins*cnt1/N); xlim([0 1]); ylim([0 4]);
subplot(1,3,2);
[cnt2, ctr2] = hist(Wcmax, bins);
plot(ctr2, bins*cnt2/N); xlim([0 1]); ylim([0 4]);
subplot(1,3,3);
[cnt3, ctr3] = hist(Y, bins);
plot(ctr3, bins*cnt3/N); xlim([0 1]); ylim([0 4]);
sgtitle(['idx:', num2str(idx), ' y:', num2str(y), ' N:', num2str(N), ' bins:', num2str(bins)]);
figure;
scatter(W,R);
ylim([0 1000])
figure;
plot(W);
xlim([0 100]);

bins = 50;
max_idx = max(find(W > 0.01));
Y = Y2;
figure;
subplot(1,3,1);
[cnt1, ctr1] = hist(Wmax(1:max_idx), bins);
plot(ctr1, bins*cnt1/max_idx); xlim([0 1]); ylim([0 4]);
subplot(1,3,2);
[cnt2, ctr2] = hist(Wcmax(1:max_idx), bins);
plot(ctr2, bins*cnt2/max_idx); xlim([0 1]); ylim([0 4]);
subplot(1,3,3);
[cnt3, ctr3] = hist(Y(1:max_idx), bins);
plot(ctr3, bins*cnt3/max_idx); xlim([0 1]); ylim([0 4]);
sgtitle(['idx:', num2str(idx), ' y:', num2str(y), ' N:', num2str(N), ' bins:', num2str(bins)]);
figure;
[cnt4, ctr4] = hist(W(1:max_idx), bins);
plot(ctr4, bins*cnt4/max_idx); xlim([0 1]); ylim([0 4]);
