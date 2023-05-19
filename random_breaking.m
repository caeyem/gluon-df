%%
% random breaking
%---------------------------------------------------------------------------------
alpha = [0.1:0.1:5.0];
% check overload here
W = [0:0.05:1.0];
rho = zeros(size(alpha, 2), size(W, 2));

figure;
for a = 1:size(alpha, 2),
    for r = 1:size(W, 2),
        rho(a, r) = (alpha(a) + 1)*(1 - W(r))^(alpha(a));
    end

%     plot(W, rho);
%     title(['alpha = ' num2str(alpha(a))]);
%     pause;
end

plot(W, rho); xlabel('W'); ylabel('rho'); title('Distributions for BI corresponding to alpha = 0.1:0.1:5.0');

% new Ws and distrib here
N = 1e7;
u = rand(1, N);

% probabilities
W = u; % check alternate approach
Wmax = zeros(1, 1e7+1);
Wcmax = zeros(1, 1e7+1);
Y = zeros(1, 1e7+1);
for n = 1:1e7,
    W(n) = u(n); % optional as defined in preallocation
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
sgtitle(['Distributions for the BI model ', ' alpha:', num2str(0), ' N:', num2str(N), ' bins:', num2str(bins)]);
