
%%
% theoretical overlap distributions
% --------------------------------------------------------------------------

% singularity
% random map
w1 = [0.51:0.01:1.00];
P1RMW = 1./(2.*w1) * 1./(1-w1).^(0.5);
w2 = [0.34:0.01:0.49];
P2RMW = (1./((4.*w2).*((1-w2).^(0.5)))) .* (log(1+(1-w2./(1-w2)).^(0.5))-log(1-(1-w2./(1-w2)).^(0.5)));
w = [w2 0.5 w1];
PRMW = [P2RMW nan P1RMW];
figure;
plot(w, PRMW); xlim([0.5 1]); ylim([0 4]);
xlabel('Y'); ylabel('\Pi(Y)'); title(['Theoretical Overlap Distribution for RM']);
% spin glass
y = 0.7;
w1 = [0.51:0.01:1.00];
P1SGW = ((w1.^(y-2)).*(1-w1).^(-y))/(gamma(y)*gamma(1-y));
w2 = [0.34:0.01:0.49];
% complex integral
% P2SGW = 1./(2*w) * 1./(1-w)^(0.5);
P2SGW = nan(size(w2));
w = [w2 0.5 w1];
PSGW = [P2SGW nan P1SGW];
figure;
plot(w, PSGW); xlim([0.5 1]); ylim([0 4]);
xlabel('Y'); ylabel('\Pi(Y)'); title(['Theoretical Overlap Distribution for y= ', num2str(y)]);
% break
alpha = 0.2;
w1 = [0.51:0.01:1.00];
P1BIW = (alpha+1).*((1-w1).^(alpha))./w1;
w2 = [0.34:0.01:0.49];
P2BIW = (1./(w2)) .* log((1-w2)./w2);
w = [w2 0.5 w1];
PBIW = [P2BIW nan P1BIW];
figure;
plot(w, PBIW); xlim([0.5 1]); ylim([0 4]);
xlabel('Y'); ylabel('\Pi(Y)'); title(['Theoretical Overlap Distribution for BI for alpha = ', num2str(alpha)]);
