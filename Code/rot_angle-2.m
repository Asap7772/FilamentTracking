function [theta] = rot_angle(slope,x)
%Input: Differential of the spline evaluated at x = slope at every x
%Output: Angle of rotation for the image in a region around x

for i=1:length(x)
         theta(i)= atan(slope(i));   % Plug in the x-value and solve for tangent angle
end
