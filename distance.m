function distZ1Star = distance(z1)
global z1Star continuumRad

continuumRad = 0;

% find the distance to the minimizer
if (abs(z1 - z1Star) <= continuumRad)
    distZ1Star = 0;
else
    distZ1Star = abs(abs(z1) - (z1Star + continuumRad));
end

end