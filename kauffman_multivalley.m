%%
% multivalley structure in kauffman model
% --------------------------------------------------------------------------

Nmax = 23;
ybar = zeros(1,Nmax);
y2bar = zeros(1,Nmax);
y3bar = zeros(1,Nmax);
y4bar = zeros(1,Nmax);
ybarsk = zeros(1,Nmax);
y2barsk = zeros(1,Nmax);
y3barsk = zeros(1,Nmax);
y4barsk = zeros(1,Nmax);

for N=1:Nmax,
    w = rand(1,2^N);
    ybar(N) = mean(w);
    y2bar(N) = mean(w.^2);
    y3bar(N) = mean(w.^3);
    y4bar(N) = mean(w.^4);
    % using random free energy picture for the Sherrington-Kirkpatrick model
    y2barsk(N) = 0.3*(ybar(N)+2*ybar(N)^2);
    y3barsk(N) = 0.5*ybar(N)*(ybar(N)+1);
    y4barsk(N) = 0.167*ybar(N)*(ybar(N)+1)*(ybar(N)+2);
end

figure;
hold on;
scatter([1:Nmax], ybar, 'or');
scatter([1:Nmax], y2bar, '^r');
scatter([1:Nmax], y3bar, 'vr');
scatter([1:Nmax], y4bar, '*r');
xlabel('N'); ylabel('Moment value'); title(['MC simulated moments of Kauffman model 1.o 2.^ 3.v 4.*']);

figure;
hold on;
scatter([1:Nmax], ybar, 'or');
scatter([1:Nmax], y2barsk, '^r');
scatter([1:Nmax], y3barsk, 'vr');
scatter([1:Nmax], y4barsk, '*r');
xlabel('N'); ylabel('Moment value'); title(['MC simulated moments of SK model 1.o 2.^ 3.v 4.*']);

figure;
hold on;
scatter([1:Nmax], abs(y2barsk-y2bar), '^r');
scatter([1:Nmax], abs(y3barsk-y3bar), 'vr');
scatter([1:Nmax], abs(y4barsk-y4bar), '*r');
xlabel('N'); ylabel('Moment value'); title(['Error in MC simulated moments 1.o 2.^ 3.v 4.*']);
