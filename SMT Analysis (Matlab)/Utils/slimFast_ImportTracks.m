%slimFast track import.

function [trajectories, tracklength] = slimFast_ImportTracks(tracks_file);

table = csvread(tracks_file,1,1);
trajectories = cell(max(table(:,3)), 1);

for i=1:max(table(:,3));
idx = find(table(:,3) == i);
track_temp = table(idx,[2 4 5]);

trajectories{i} = track_temp;
tracklength(i) = size(track_temp,1);

end

end