function [Esp_Coef, Esp_Sigma Esp_Fit] = ExpDecay_2Cmp_fit(data, k0s, f01, range);

% Calculate best two component exponential fit with a the lsqnonlin non 
%linear least square fitting routine.


% Read input data

% tlist =data(:,1);
% Y = data(:,2);

tlist =data(range(1):min(range(end),size(data,1)),1);
Y = data(range(1):min(range(end),size(data,1)),2);


k01 = k0s(1);
k02 = k0s(2);
f01 = 0.5;
A0 = max(Y);


% define initial values for the parameters
COEF0 = [k01, k02, f01, A0];


% define lower boundaries for the fitted parameters
COEF_LB = [0, 0, 0 , 0];

% define upper boundaries for the fitted parameters
COEF_UB = [Inf, Inf, 1,Inf];

% Define anonymous function for the fitting
fitfun = @(COEF) (ExpDecay_2Cmp_fun(COEF, tlist)) - (Y);

% Select Options for the fitting
options = optimset('FunValCheck','off', 'MaxIter', 1000, 'MaxFunEvals', 1000,'display','off');

% run fitting routine
[Esp_Coef, resNorm, residuals,exitflag,output,lambda,jacobian] = lsqnonlin(fitfun,COEF0,COEF_LB,COEF_UB,options);

ci = nlparci(Esp_Coef,residuals,'jacobian',jacobian);
Esp_Sigma = (ci(:,2) - ci(:,1))/2;
Esp_Sigma = Esp_Sigma';



% Compute output
Esp_Fit(:,1) = data(:,1);
Esp_Fit(:,2) = ExpDecay_2Cmp_fun(Esp_Coef, data(:,1));

%Compute the individual components --DB
Esp_Coef_1 = [Esp_Coef(1) Esp_Coef(4)*Esp_Coef(3)];
Esp_Fit(:,3) = ExpDecay_fun(Esp_Coef_1,data(:,1));

Esp_Coef_2 = [Esp_Coef(2) Esp_Coef(4)*(1 - Esp_Coef(3))];
Esp_Fit(:,4) = ExpDecay_fun(Esp_Coef_2,data(:,1));

