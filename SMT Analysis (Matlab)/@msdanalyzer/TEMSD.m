function  varargout = TEMSD(obj, TEMSD_points)
%%FITMEANMSD Fit the weighted averaged MSD by a linear function.
%
% obj.fitMeanMSD computes and fits the weighted mean MSD by a
% straight line y = a * x. The fit is therefore valid only for
% purely diffusive behavior. Fit results are displayed in the
% command window.
%
% obj.fitMeanMSD(clip_factor) does the fit, taking into account
% only the first potion of the average MSD curve specified by
% 'clip_factor' (a double between 0 and 1). If the value
% exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.
%
% [fo, gof] = obj.fitMeanMSD(...) returns the fit object and the
% goodness of fit.


if ~obj.msd_valid
    obj = obj.computeMSD;
end

ft = fittype('poly1');
mmsd = obj.getMeanMSD;

%If the size of Tracks is smaller than the number of points, change that
%number.
if (TEMSD_points+1) > size(mmsd,1);
    TEMSD_points = size(mmsd,1)-1;
end

% Clip data, never take the first one dt = 0
t = mmsd(2:TEMSD_points+1,1);
y = mmsd(2:TEMSD_points+1,2);
w = 1./mmsd(2:TEMSD_points+1,3);

     % Thrash bad data (When the track is shorter than TMSD_points-1)
    nonnan = ~isnan(y);  
    t = t(nonnan);
    y = y(nonnan);
    w = w(nonnan);
 
%Normal fit    
[fo, gof] = fit(t, y, ft, 'Weights', w);
ci = confint(fo);

%Least Squares Fit
Initial_values = [0.05,0.01];
Xdata = t;
Ydata = y;
linear_poly = @(X,Xdata)( (X(1)*Xdata) + X(2));
options = optimset('Display','off');
[X_fitted,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(linear_poly,Initial_values,Xdata,Ydata,[],[],options);
conf_int = nlparci(X_fitted,residual,'jacobian',jacobian);



%---------------------------------------
str = sprintf(['Estimating D through Least Squares linear fit of the mean MSD curve.\n', ...
    'D = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n'], X_fitted(1)/2/obj.n_dim, conf_int(1,1)/2/obj.n_dim, conf_int(1,2)/2/obj.n_dim);
disp(str)

fprintf('\n');

str = sprintf([
    'Estimating D through linear weighted fit of the mean MSD curve.\n', ...
    'D = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n', ...
    'Goodness of fit: R² = %.3f.' ], ...
    fo.p1/2/obj.n_dim, ci(1)/2/obj.n_dim, ci(2)/2/obj.n_dim, gof.adjrsquare);
disp(str)

if nargout > 0
    varargout{1} = fo;
    if nargout > 1
        varargout{2} = gof;
        if nargout > 2
            varargout{3} = X_fitted;
            if nargout > 3
                varargout{4} = conf_int;
            end
        end
    end
end

end