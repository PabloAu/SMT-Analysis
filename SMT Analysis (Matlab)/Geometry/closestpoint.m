%Calculates the closes point in a plane to a certain point (3D).

function x = closestpoint(n, d, p);


% n is the vector [A,B,C] that defines the plane (Normalized)
%Ax + By + Cz = D (equation of the plane)
% d is the distance of the plane from the origin
% p is the point  [P,Q,R]
v = (d - sum(p.*n)) / sum(n.*n);
x = p + v * n;

end