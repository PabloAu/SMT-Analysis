function [max_dist, min_dist, avg_dist] = distance_scatter(points);

%% Find the maximum, minimum and avegare distances in a group of scatter points

%Points if a 2D matrix with x and y coordinates
x = points(:,1);
y = points(:,2);
numberOfCoords = length(x);

maxDistance = zeros(1, numberOfCoords);
minDistance = zeros(1,numberOfCoords);

indexOfMax = zeros(1, numberOfCoords, 'int32');
indexOfMin = zeros(1, numberOfCoords, 'int32');

%-----------------------------------------------
% Find the furthest away points.
for k = 1 : numberOfCoords;
    
  distances = sqrt((x-x(k)).^2 + (y-y(k)).^2);
  
  maxDistance(k) = max(distances);
  abl = min(nonzeros(distances));
  if isempty(abl);
      abl = 0;
  end
  minDistance(k) = abl;
  avgDistance(k) = mean(distances);


end

max_dist = max(maxDistance);
min_dist = min(minDistance);
avg_dist = mean(avgDistance);



end