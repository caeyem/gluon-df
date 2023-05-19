%%
% random energy model
% --------------------------------------------------------------------------

y = 0.7;
w = [0:0.01:1.0];
N = 1e7;
bins = 100;
% P(W)
PW = ((w.^(-1 + y)).*(1-w).^(-y))./(gamma(y)*gamma(1-y));
figure;
plot(w, PW);
xlabel('W'); ylabel('P(W)'); title(['Theoretical Distribution of Weights P(W) for y = ', num2str(y)]);

% P(W1,W2)
y = 0.7;
ny = 101; nx = 101;
w1 = repmat([0:0.01:1.0], ny, 1);
w2 = repmat(flipud([0:0.01:1.0]'), 1, nx);
PW1 = repmat(PW, ny, 1);
delta = 2;
theta = 1;
% assuming theta function is indicator for above 0 - theta is step function
% PW1W2 = w1.*PW1.*delta.*(w1-w2) + ...
%     (1-y)*((w1.*w2).^(-1+y)).*((1-w1-w2).^(1-2*y)).*theta.*(1-w1-w2)./(((gamma(y))^2)*gamma(2-2*y));
PW1W2 = w1.*PW1.*delta.*(w1-w2) + ...
    (1-y)*((w1.*w2).^(-1+y)).*(real((1-w1-w2).^(2-2*y))).*theta./(((gamma(y))^2)*gamma(2-2*y));
figure;
imagesc(PW1W2); colorbar;
PW1W2alt = w1.*PW1.*delta.*(w1-w2) + ...
    (1-y)*((w1.*w2).^(-1+y)).*(((1-w1-w2).^(1-2*y)).*((1-w1-w2) >= 0)).*theta./(((gamma(y))^2)*gamma(2-2*y));
figure;
imagesc(PW1W2alt); colorbar;
colormap bone;
colorbar;
title(['Theoretical Distribution of Weight Pairs P(W_1,W_2) for y = ', num2str(y)]);

rhonx = PW;
rhonxc = rhonx;
rhonxc(isinf(rhonx)) = 0;
%cdf
% \Delta x for the Riemann sum
dx = 1/length(w);
Frho = cumsum(rhonxc).*dx; %check how to remove the 2 correction
% further adjustments
maxFrho = nanmax(Frho);
Frhoc = Frho;
Frhoc(isnan(Frho)) = maxFrho;
Frhoc = Frhoc/maxFrho;
F_dist = makedist('PiecewiseLinear', 'x', w, 'Fx', Frhoc);
v = random(F_dist, 1, N);
% probabilities
W = v;
Wmax = zeros(1, N+1);
Wcmax = zeros(1, N+1);
Y = zeros(1, N+1);
for n = 1:N,
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
sgtitle(['Distributions for the REM-SG/GluonTM model ', ' y:', num2str(y), ' N:', num2str(N), ' bins:', num2str(bins)]);

