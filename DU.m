function inside = DU(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: DU.m
%--------------------------------------------------------------------------
% Description: Jump set
% Return 0 if outside of D, and 1 if inside D
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2020 1:07:00

global gamma M cTilde_10 d_10 alpha c_0 V0

% state
z1 = x(1);
z2 = x(2);
q = x(3);

absGradL = abs(GradientL(z1));

if(q == 0)
    V0 = gamma*(alpha/(M^2))*absGradL^2 + (1/2)*z2^2;
elseif(q == 1)
    z2Squared = z2^2;
end

if (q == 0 && V0 >= c_0)||(q == 1 && absGradL <= cTilde_10 && z2Squared <= d_10)  
    inside = 1;
else
    inside = 0;
end

end