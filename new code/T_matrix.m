function [ T ] = T_matrix( Angle , Ax_rot )
% Function returns a transformation matrix for a 3D beam element given the
% an angle a vector and, the vector is the axis of rotation.
%
% Arthur Schout
% 07/21/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: rotation in radians, axis of rotation
% output: transformation matrix [3x3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT: rotation in radians, what is positive what is negative"?
%       vector which is the axis of rotation, if no rotation vector is
%       given then the axis of rotation is assumed to be [0 0 1] (z-axis)
%
% OUTPUT: the transformation matrix [3x3]
%

if nargin == 1
    % axis of rotation normalized vector
    u = [0 0 1]; 
    % transformation matrix
    T = zeros(3,3);
    T(1,1) = cos(Angle)+u(1)^2*(1-cos(Angle));
    T(2,2) = cos(Angle)+u(2)^2*(1-cos(Angle));
    T(3,3) = cos(Angle)+u(3)^2*(1-cos(Angle));
    T(1,2) = u(1)*u(2)*(1-cos(Angle))-u(3)*sin(Angle);
    T(2,1) = u(1)*u(2)*(1-cos(Angle))+u(3)*sin(Angle);
    T(1,3) = u(1)*u(3)*(1-cos(Angle))+u(2)*sin(Angle);
    T(3,1) = u(1)*u(3)*(1-cos(Angle))-u(2)*sin(Angle);
    T(2,3) = u(2)*u(3)*(1-cos(Angle))-u(1)*sin(Angle);
    T(3,2) = u(2)*u(3)*(1-cos(Angle))+u(1)*sin(Angle);
else
    u = Ax_rot/norm(Ax_rot);
    T = zeros(3,3);
    T(1,1) = cos(Angle)+u(1)^2*(1-cos(Angle));
    T(2,2) = cos(Angle)+u(2)^2*(1-cos(Angle));
    T(3,3) = cos(Angle)+u(3)^2*(1-cos(Angle));
    T(1,2) = u(1)*u(2)*(1-cos(Angle))-u(3)*sin(Angle);
    T(2,1) = u(1)*u(2)*(1-cos(Angle))+u(3)*sin(Angle);
    T(1,3) = u(1)*u(3)*(1-cos(Angle))+u(2)*sin(Angle);
    T(3,1) = u(1)*u(3)*(1-cos(Angle))-u(2)*sin(Angle);
    T(2,3) = u(2)*u(3)*(1-cos(Angle))-u(1)*sin(Angle);
    T(3,2) = u(2)*u(3)*(1-cos(Angle))+u(1)*sin(Angle);
end


end

