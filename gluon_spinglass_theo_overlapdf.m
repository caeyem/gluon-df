%%
% spin glass overlap functions
% --------------------------------------------------------------------------
y=[0.5:0.01:1.00];
alphaval = 0.2;

Py = (alphaval./((1-y).*((2*y-1).^0.5))) .* ((1-(2*y-1).^(0.5))./(1+(2*y-1).^(0.5))).^(alphaval);
y=[0.0:0.01:0.49, y];
Py=[zeros(1,50) Py];
figure;
plot(y,Py);
xlabel('Y'); ylabel('\Pi(Y)'); title(['Overlap Function Distribution for \alpha: ', num2str(alphaval)]);
gammaval = 0.577215664901532860606512090082; % euler constant
Py2 = ((exp(-gammaval*alphaval))/(2^(alphaval) * gamma(alphaval))*(1-y).^(alphaval-1));
figure;
plot(y,Py2);
xlabel('Y'); ylabel('\Pi(Y)'); title(['Overlap Function Distribution for \alpha: ', num2str(alphaval)]);
Fy = 0.5*exp(-alphaval)*((1+alphaval./((alphaval^2+y.^2).^(0.5))).*exp((alphaval^2+y.^2).^0.5) + ...
    (1-alphaval./((alphaval^2+y.^2).^(0.5))).*exp(-(alphaval^2+y.^2).^0.5));
figure;
plot(y,Fy);
xlabel('Y'); ylabel('F(y)'); title(['Function F(y) for \alpha: ', num2str(alphaval)]);
