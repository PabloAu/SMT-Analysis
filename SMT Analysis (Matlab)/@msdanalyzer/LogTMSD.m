function obj = LogTMSD(obj, LogTMSD_points)
%%FITLOGLOGMSD Fit the log-log MSD to determine behavior.
%
% obj = obj.fitLogLogMSD fits each MSD curve stored in this object
% in a log-log fashion. If x = log(delays) and y = log(msd) where
% 'delays' are the delays at which the msd is calculated, then this
% method fits y = f(x) by a straight line y = alpha * x + gamma, so
% that we approximate the MSD curves by MSD = gamma * delay^alpha.
% By default, only the first 25% of each MSD curve is considered
% for the fit,
%
% Results are stored in the 'loglogfit' field of the returned
% object. It is a structure with 3 fields:
% - alpha: all the values for the slope of the log-log fit.
% - gamma: all the values for the value at origin of the log-log fit.
% - r2fit: the adjusted R2 value as a indicator of the goodness of
% the fit.
%
% obj = obj.fitLogLogMSD(clip_factor) does the fit, taking into
% account only the first potion of each MSD curve specified by
% 'clip_factor' (a double between 0 and 1). If the value
% exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.


if ~obj.msd_valid
    obj = obj.computeMSD;
end
n_spots = numel(obj.msd);


%If the size of Tracks is smaller than the number of points, change that
%number.
if (LogTMSD_points+1) > size(obj.msd{1},1);
    LogTMSD_points = size(obj.msd{1},1)-1;
end

    fprintf('Fitting %d curves of log(MSD) = f(log(t)), taking only the first %d points of each curve... ',...
        n_spots, LogTMSD_points)

alpha = NaN(n_spots, 1);
gamma = NaN(n_spots, 1);
r2fit = NaN(n_spots, 1);
ft = fittype('poly1');

% fprintf('%d/%d', 0, n_spots);



for i_spot = 1 : n_spots
    
%     fprintf('\b\b\b\b\b\b\b\b\b%d/%d', i_spot, n_spots);
    
    msd_spot = obj.msd{i_spot};
    
    % Clip data
    t = msd_spot(2:LogTMSD_points+1,1);
    y = msd_spot(2:LogTMSD_points+1,2);
    w = msd_spot(2:LogTMSD_points+1,4);
    
     % Thrash bad data (When the track is shorter than TMSD_points-1)
    nonnan = ~isnan(y);  
    t = t(nonnan);
    y = y(nonnan);
    w = w(nonnan);    
    if numel(y) < 2
        continue
    end
    
    xl = log(t);
    yl = log(y);
    
    bad_log =  isinf(xl) | isinf(yl);
    xl(bad_log) = [];
    yl(bad_log) = [];
    w(bad_log) = []; 
    if numel(xl) < 2
        continue
    end
    
    [fo, gof] = fit(xl, yl, ft, 'Weights', w);
    
    alpha(i_spot) = fo.p1;
    gamma(i_spot) = exp(fo.p2);
    r2fit(i_spot) = gof.adjrsquare;
    
end
% fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
fprintf('\n')

obj.loglogfit = struct(...
    'alpha', alpha, ...
    'gamma', gamma, ...
    'r2fit', r2fit);

end