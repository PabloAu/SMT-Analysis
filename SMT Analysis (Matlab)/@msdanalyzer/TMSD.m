function obj = TMSD(obj, TMSD_points)
%%FITMSD Fit all MSD curves by a linear function.
%
% obj = obj.fitMSD fits all MSD curves by a straight line
%                      y = a * x + b.
% The fit is therefore rigorously valid only for purely
% diffusive behavior.
%
% Results are stored in the 'fit' field of the returned
% object. It is a structure with 2 fields:
% - a: all the values of the slope of the linear fit.
% - b: all the values for the intersect of the linear fit.
% - r2fit: the adjusted R2 value as a indicator of the goodness
% of the fit.
%
% obj = obj.fitMSD(clip_factor) does the fit, taking into
% account only the first potion of the average MSD curve
% specified by 'clip_factor' (a double between 0 and 1). If the
% value exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.


if ~obj.msd_valid
    obj = obj.computeMSD;
end
n_spots = numel(obj.msd);

%If the size of Tracks is smaller than the number of points, change that
%number.
if (TMSD_points+1) > size(obj.msd{1},1);
    TMSD_points = size(obj.msd{1},1)-1;
end

    fprintf('Fitting %d curves of MSD = f(t), taking only the first %d points of each curve... ',...
        n_spots, TMSD_points )

a = NaN(n_spots, 1);
b = NaN(n_spots, 1);
r2fit = NaN(n_spots, 1);
ft = fittype('poly1');

% fprintf('%d/%d', 0, n_spots);




%Iterate on trough the Tracks
for i_spot = 1 : n_spots
    
%     fprintf('\b\b\b\b\b\b\b\b\b%d/%d', i_spot, n_spots);
    
    msd_spot = obj.msd{i_spot};
    
% Clip data, never take the first one dt = 0
    t = msd_spot(2:TMSD_points+1,1);
    y = msd_spot(2:TMSD_points+1,2);
    w = msd_spot(2:TMSD_points+1,4);
    
     % Thrash bad data (When the track is shorter than TMSD_points-1)
    nonnan = ~isnan(y);  
    t = t(nonnan);
    y = y(nonnan);
    w = w(nonnan);    
        
    if numel(y) < 2
        continue
    end
    
    [fo, gof] = fit(t, y, ft, 'Weights', w);
    
    a(i_spot) = fo.p1;
    b(i_spot) = fo.p2;
    r2fit(i_spot) = gof.adjrsquare;
    
end
% fprintf('\b\b\b\b\b\b\b\b\bDone.\n')
fprintf('\n')

obj.lfit = struct(...
    'a', a, ...
    'b', b, ...
    'r2fit', r2fit);

end