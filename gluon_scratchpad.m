% generate log pdf plots - log scale  going from 10^-4 to 1 on x

%% Preliminaries - Scratchpad
%-------------------------------------------------------------------------------

s = 100;
fs = 1; % free energy
T = 100;
Wl = 100;
delta = 2; % TODO: redo all calcs with delta = 2

% ii
Ws = exp(-fs/T);

% iii
Ws = rand(1, s);

Ws = Ws/sum(Ws);

Y = 0;
for s = 1:size(Ws, 2),
    Y = Y + Ws*Ws;
end

% norm Ws
HW = zeros(1, Wl);
for W = 1/Wl:1/Wl:1,
    HW(W*Wl) = W*sum(floor(Ws*Wl) == floor(W*Wl));
end
% index of HW in reality is 1/W1 times the index above

HWc = HW;
HWcc = HW;
W = 1/Wl:1/Wl:1;
Wc = W;
Wcc = W;

% unclear if this an expectation or a function of W
gW = mean(HW, 2);
gWWc = mean(HW.*HWc, 2);
gWWcWcc = mean(HW.*HWc.*HWcc, 2);

% check these logic
fW = gW./W;
fWWc = gWWc./(W.*Wc) - gWc./(Wc);

% random map
fRMW = (0.5)./W./((1-W).^2);
fRMWWc = (0.25)./(W.*Wc)./((1-W-Wc).^2);

% spin glass
fSGW = ((W.^(y-2)).*((1-W).^(-y)))./(gamma(y)*gamma(1-y));
fSGWWc = ((1-y)*((W.*Wc).^(y-2)).*((1-W-Wc).^(1-2*y)))./(((gamma(y)).^(2))*gamma(2-2*y));
K=100;
% check changes for unequal Wk
fSGWk = ((1-y).^(K-1))*gamma(K)./((gamma(y)^K)*(gamma(K-K*y)));
for k=1:K,
    fSGWk = (fSGWk.*(W.*(y-2)));
end
fSGWk = fSGWk*((1-K*W).^(K-K*y-1));

YP = 0;
P = 2;
for s = 1:size(Ws, 2),
    YP = YP + Ws.^P;
end

% Examples
% Random breaking of an interval
% W changes here
distrib = @uniform;
N = 100;
W = zeros(1, N);
prefac = 1;
for n = 1:N,
    x = distrib(1);
    W(n) = prefac * x;
    prefac = prefac * (1-x);
end

alpha = [0.1:0.1:5.0];
% check overload here
W = [0:0.05:1.0];
rho = zeros(size(alpha, 2), size(W, 2));

figure;
hold on;
for a = 1:size(alpha, 2),
    for r = 1:size(W, 2),
        rho(a, r) = (alpha(a) + 1)*(1 - W(r))^(alpha(a));
    end

%     plot(W, rho);
%     title(['alpha = ' num2str(alpha(a))]);
%     pause;
end

plot(W, rho); xlabel('W'); ylabel('rho'); title('Distributions for BI corresponding to alpha = 0.1:0.1:5.0');

alphaval = 0.2;
gB1W = rho(squeeze(find(alpha == alphaval)));
% looks like this needs to be made into a matrix -tbu
gB1WWc = [(alphaval + 1)*W.*(1 - W).^(alphaval)].*[floor(Ws*Wl) == floor(W*Wl)];
