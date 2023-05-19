%%
% nature of spin glass
% --------------------------------------------------------------------------
y=[0:0.01:1.00];
Keps = (1-y).*gamma((1-y).^2)./gamma(2-2.*y);
figure;
plot(y, Keps);
xlabel('y'); ylabel('P_k'); title(['Distribution of P_k']);

y = 0.7;
y = 0.9;
% The sequence of double factorials for even n = 0, 2, 4, 6, 8,... starts as
% 1, 2, 8, 48, 384, 3840, 46080, 645120,
% The sequence of double factorials for odd n = 1, 3, 5, 7, 9,... starts as
% 1, 3, 15, 105, 945, 10395, 135135
y1 = y;
y2 = (1/3)*(y+2*y^2);
y3 = (1/15) *(3*y+7*y^2+5*y^3);
y4 = (1/105) *(15*y+39*y^2+37*y^3+14*y^4);
y5 = (1/945) *(105*y+296*y^2+326*y^3+176*y^4+42*y^5);
y6 = (1/10395) *(945*y+2838*y^2+3458*y^3+2228*y^4+794*y^5+132*y^6);
y7 = (1/135135) *(10395*y+32859*y^2+43191*y^3+31235*y^4+13553*y^5+3473*y^6+429*y^7);
%mgf, % characteristic func
syms w;
F = 1 + (i*w*y1) + (((i*w)^2)*y2)/factorial(2) + (((i*w)^3)*y3)/factorial(3) + ...
     + (((i*w)^4)*y4)/factorial(4) +  (((i*w)^5)*y5)/factorial(5)  + (((i*w)^6)*y6)/factorial(6) + ...
      + (((i*w)^7)*y7)/factorial(7);
pt = ifourier(F);
% pt = (2*pi*dirac(x) + (7*pi*dirac(1, x))/5 + (14*pi*dirac(2, x))/25 + (161*pi*dirac(3, x))/1000 + (453*pi*dirac(4, x))/12500 + (482512061196773*pi*dirac(5, x))/72057594037927936 + (1699251701497829*pi*dirac(6, x))/1621295865853378560 + (6453806293987741*pi*dirac(7, x))/45396284243894599680)/(2*pi)
x=[0:0.01:1];
ptval=subs(pt);
figure;
plot(x,ptval);
% alternate
wf = 0:0.01:2*pi;
Ff = 1 + (i*wf*y1) + (((i*wf).^2)*y2)/factorial(2) + (((i*wf).^3)*y3)/factorial(3) + ...
     + (((i*wf).^4)*y4)/factorial(4) +  (((i*wf).^5)*y5)/factorial(5)  + (((i*wf).^6)*y6)/factorial(6) + ...
      + (((i*wf).^7)*y7)/factorial(7);
ptt = ifft(Ff);
% pt = (2*pi*dirac(x) + (7*pi*dirac(1, x))/5 + (14*pi*dirac(2, x))/25 + (161*pi*dirac(3, x))/1000 + (453*pi*dirac(4, x))/12500 + (482512061196773*pi*dirac(5, x))/72057594037927936 + (1699251701497829*pi*dirac(6, x))/1621295865853378560 + (6453806293987741*pi*dirac(7, x))/45396284243894599680)/(2*pi)
figure;
plot(abs(ptt));
