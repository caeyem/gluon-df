%%
% mean field theory of spin glasses
%---------------------------------------------------------------------------------
% Parameter Sets

% charts
yList = [0.6275 0.7 0.7 0.707 0.9 0.3 ];
NList = [1e7 1e7 3.25e6 1e7 1e7 5e7 ];
binsList = [1e2 1e2 1e2 1e2 1e2 1e3 ];

% alternate charts
yList = [0.1:0.1:1.0];
NList = 1e7 * ones(1, size(yList, 2));
binsList = 1e2 * ones(1, size(yList, 2));
