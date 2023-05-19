%%
% mean field theory of spin glasses
%---------------------------------------------------------------------------------
% approach 1: using beta distribution
% distribution check

% try using beta distribution
xx = [0:0.01:1];
for yCtr = 1:size(yList,2),
    figure;
    y = yList(yCtr);
    sgtitle(['Sample Beta distributions: y: ', num2str(y)]);
    for nn=1:1:10
        yy = betapdf(xx,y,(1-y)*nn);
        subplot(5, 2, nn);
        plot(xx, yy);
        title(['n = ', num2str(nn)]);
    end
end

xx = [0:0.01:1];
for yCtr = 1:size(yList,2),
    figure;
    y = yList(yCtr);
    sgtitle(['Sample Beta distributions: y: ', num2str(y)]);
    for nn=10:10:100
        yy = betapdf(xx,y,(1-y)*nn);
        subplot(5, 2, nn/10);
        plot(xx, yy);
        title(['n = ', num2str(nn)]);
    end
end
nn = 1;


for idx = 1:size(yList, 2),
    y = yList(idx);
    N = NList(idx);
    bins = binsList(idx);

    xn = [0.0:1/N:1.0];
    n = 1:size(xn,2);

    alpha = y;
    beta1 = n*(1-y);
    beta2 = nn*(1-y);
    rhonx1 = betarnd(alpha, beta1);
    rhonx2 = betarnd(alpha, beta2, [1 N]);

    verbose_plot = false;
    if(verbose_plot)
        figure;
        subplot(121);
        histogram(rhonx2, 1e2);
        xlabel('x'); ylabel('count'); title('Histogram of constructed samples');
        subplot(122);
        xx = [0:0.01:1];
        yy = betapdf(xx,y,(1-y)*nn);
        subplot(122);
        plot(xx, yy);
        xlabel('x'); ylabel('count'); title('Theoretical Distribution');
    end

    % probabilities
    W = rhonx2;
    Wmax = zeros(1, N+1);
    Wcmax = zeros(1, N+1);
    Y = zeros(1, N+1);
    Wfac = 1;
    for n = 1:N,
%         W(n) = rhonx1(n);
%         W(n) = Wfac*rhonx2(n);
%         Wfac = (Wfac - W(n));
        Wmax(n+1) = max([(1-W(n))*Wmax(n), W(n)]);
        Wcmax(n+1) = max([min([(1-W(n))*Wmax(n), W(n)]),(1-W(n))*Wcmax(n)]);
        Y(n+1) = (W(n))^2 + ((1 - W(n))^2)*Y(n);
    end

    if(verbose_plot)
        figure;
        subplot(1,3,1);
        histogram(Wmax, bins, 'Normalization', 'pdf');
        xlabel('W: W_{max}'); ylabel('H_1(W)'); title('Max Valley Weight Histogram');
        subplot(1,3,2);
        histogram(Wcmax, bins, 'Normalization', 'pdf');
        xlabel('W: W_{max}^c'); ylabel('H_2(W)'); title('Second Max Valley Weight Histogram');
        subplot(1,3,3);
        histogram(Y, bins, 'Normalization', 'pdf');
        xlabel('Y'); ylabel('\Pi(Y)'); title('Overlap Histogram');
        sgtitle(['Histograms for the SG/GluonTM model ', ' y:', num2str(y), ' N:', num2str(N), ' bins:', num2str(bins)]);
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
    sgtitle(['Distributions for the SG/GluonTM model ', ' y:', num2str(y), ' N:', num2str(N), ' bins:', num2str(bins)]);

    log_plot = false;
    if(log_plot)
        figure;
        subplot(1,3,1);
        [cnt1, ctr1] = hist(Wmax, bins);
        semilogx(ctr1, bins*cnt1/N); xlim([0 1]); ylim([0 4]);
        xlabel('log W: W_{max}'); ylabel('P_1(W)'); title('Max Valley Weight Distribution');
        subplot(1,3,2);
        [cnt2, ctr2] = hist(Wcmax, bins);
        semilogx(ctr2, bins*cnt2/N); xlim([0 1]); ylim([0 4]);
        xlabel('log W: W_{max}^c'); ylabel('P_2(W)'); title('Second Max Valley Weight Distribution');
        subplot(1,3,3);
        [cnt3, ctr3] = hist(Y, bins);
        semilogx(ctr3, bins*cnt3/N); xlim([0 1]); ylim([0 4]);
        xlabel('log Y'); ylabel('\Pi(Y)'); title('Overlap Distribution');
        sgtitle(['Distributions for the SG/GluonTM model ', ' y:', num2str(y), ' N:', num2str(N), ' bins:', num2str(bins)]);
    end

end
