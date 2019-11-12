%functions that models the MSD of a confined circle diffusion

function y = msd_conf_circ_diff_V2(x,R,D,offset);


% With offset
y = R^2*(1-exp(-4*D*x/R^2)) + offset;

% %Without offset
% y = R^2*(1-exp(-4*D*x/R^2));



end