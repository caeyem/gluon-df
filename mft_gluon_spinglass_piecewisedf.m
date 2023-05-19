
%%
% mean field theory of spin glasses
%---------------------------------------------------------------------------------
% approach 2: using direct piecewise linear distribution construction
% distribution check
xn = [0.0:0.001:1.0];
rhonx = xn;
y = 0.7;
for n=1:size(xn,2),
    rhonx(n) = (gamma(n-n*y+y)/(gamma(n-n*y)*gamma(y)))*((xn(n))^(y-1))*((1-xn(n))^(n*(1-y)-1));
end

% check only 1/2 the distribution
% correction for singularity
rhonxc = rhonx;
rhonxc(isinf(rhonx)) = 0;
%cdf
% \Delta x for the Riemann sum
dx = 1/length(xn);
Frho = 2*cumsum(rhonxc).*dx; %check how to remove the 2 correction
% further adjustments
maxFrho = nanmax(Frho);
Frhoc = Frho;
Frhoc(isnan(Frho)) = maxFrho;
Frhoc = Frhoc/maxFrho;

figure;
subplot(1,2,1);
plot(xn, rhonxc);
subplot(1,2,2);
plot(xn, Frhoc);

F_dist = makedist('PiecewiseLinear', 'x', xn, 'Fx', Frhoc);
rand_nums = random(F_dist, 1, 1e7);
figure;
histogram(rand_nums, 1e2);
