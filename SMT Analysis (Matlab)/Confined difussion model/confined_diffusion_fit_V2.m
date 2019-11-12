%Difussion model fitting for a confined movement in circular shape

function [X_fitted, residual,jacobian,resnorm] = confined_diffusion_fit_V2(mmsd,level,num_points,D0,R0);

%Check that there are enough points on the mmsd:
if size(mmsd,1) < num_points;
num_points = size(mmsd,1);
end

Xdata = mmsd(1:num_points,1);
Ydata = mmsd(1:num_points,2);
Ydata(1,1) = level;
Xdata0 = [R0,D0,level];

msd_confined_circle_diffusion = @(X,Xdata)(X(1)^2)*(1-exp(-4*X(2)*Xdata/(X(1)^2))) + X(3);

% WEIGHTS = mmsd(1:num_points,3) ./ sqrt(mmsd(1:num_points,4));
% ft = fittype('msd_conf_circ_diff(x,R,D,offset)');
% [f , Goodness] = fit( x, y, ft,'Weight', WEIGHTS);

options = optimoptions('lsqcurvefit','Display','off');
lb = [];
ub = [];
[X_fitted,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(msd_confined_circle_diffusion,Xdata0,Xdata,Ydata,lb,ub,options);

end