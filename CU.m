function [value] = CU(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: CU.m
%--------------------------------------------------------------------------
% Description: Flow set
% Return 0 if outside of C, and 1 if inside C
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2020 1:07:00

global cTilde_0 cTilde_10 d_0 d_10

% state
z1 = x(1);
z2 = x(2);
q = x(3);

absGradL = abs(GradientL(z1));

if(q == 0)
    halfz2Squared = (1/2)*z2^2;
elseif(q == 1)
    z2Squared = z2^2;
end

if (q == 1 || absGradL <= cTilde_0 || halfz2Squared <= d_0) && (q == 0 || absGradL >= cTilde_10 || z2Squared >= d_10) 
    value = 1;
else 
    value = 0;
end

end