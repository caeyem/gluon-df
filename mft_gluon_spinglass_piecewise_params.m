%%
% mean field theory of spin glasses
%---------------------------------------------------------------------------------
% approach 2: Work on parameter sets using piecewise linear distribution construction

% try using distribution construction
for idx = 1:4,
    xn = [0.0:0.001:1.0];
    rhonx = xn;
    y = yList(idx);
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
    Frho = cumsum(rhonxc).*dx; %check how to remove the 2 correction
    % further adjustments
    maxFrho = nanmax(Frho);
    Frhoc = Frho;
    Frhoc(isnan(Frho)) = maxFrho;
    Frhoc = Frhoc/maxFrho;
    F_dist = makedist('PiecewiseLinear', 'x', xn, 'Fx', Frhoc);

    N = NList(idx);
    v = random(F_dist, 1, N);
    bins = binsList(idx);

    % probabilities
    W = v;
    Wmax = zeros(1, N+1);
    Wcmax = zeros(1, N+1);
    Y = zeros(1, N+1);
    Wfac = 1;
    for n = 1:N,
        W(n) = Wfac*v(n); % optional as defined in preallocation
        Wfac = (Wfac - W(n));
        Wmax(n+1) = max([(1-W(n))*Wmax(n), W(n)]);
        Wcmax(n+1) = max([min([(1-W(n))*Wmax(n), W(n)]),(1-W(n))*Wcmax(n)]);
        Y(n+1) = (W(n))^2 + ((1 - W(n))^2)*Y(n);
    end

    figure; title(['idx: ', num2str(idx)]);
    subplot(1,3,1);
    histogram(Wmax,bins);
    subplot(1,3,2);
    histogram(Wcmax,bins);
    subplot(1,3,3);
    histogram(Y,100);
    pause;
end
