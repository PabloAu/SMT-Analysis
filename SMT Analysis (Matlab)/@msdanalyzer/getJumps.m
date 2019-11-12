function Jumps = getJumps(obj, indices)
%%GETVELOCITIES Generate and return the instantaneous velocities.
%
% Jumps = obj.getJumps returns in v the instantaneous jumps from frame to
% frame calculated over all the particles trajectories stored in this
% object.

% This method returns a cell array, one cell per particle. Arrays
% are N x (Ndim+1) double arrays, with Ndim the dimensionality set
% at object creation. Data is organized as follow:  [ Ti JumpXi JumpYi ... ].
%
% v = obj.getJumps(indices) restrict the calculation over only
% the particles with specified indices. Use an empty array to use
% take all.

if nargin < 2 || isempty(indices)
    indices = 1 : numel(obj.tracks);
end

n_tracks = numel(indices);
Jumps = cell(n_tracks, 1);

for i = 1 : n_tracks
    
    index = indices(i);
    
    t = obj.tracks{index}(:, 1);
    X = obj.tracks{index}(:, 2:end);
    
    
    dX = diff(X, 1);
    Jumps{i} = [ t(1:end-1) dX];
end

end