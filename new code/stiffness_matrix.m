function [K_el ] = stiffness_matrix(M_prop,Sec_prop,Length,dir,typ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simple element definiton function
%%% Arthur Schout
%%% 27/11/2017
%%% V0.1 27/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function defines an element matrix in the global coordinate system
% the output is a matrix that can be used in the assembly of a the global
% stiffness matrix.
% Ref 1 Friedmand Z. and Kosmatka J.B. : AN IMPROVED TWO-NODE TIMOSHENKO BEAM FINITE ELEMENT 


% input:    M_prop = [ m , E , G ];
%           Sec_prop = [ A , Ix , Iy , PHI ];
%           length = L
%           dir = angle for now we remain in plane so only one angle is
%           given
%           type  = 'string' 

% L = the length of the element
% E = youngs modulus
% G = shear modulus
% A = the area of the beam
% I = moment of inertia of the beam
% PHI = the ratio between shear stiffness and beanding stiffness (only for
% TIMOSHENKO beams
% dir = the angle of the element and the x-axis
% type = choose between euler bernoulli or timoshenko beam.    

L = Length;
A = Sec_prop(1);
I = Sec_prop(2);
PHI = Sec_prop(4);
E = M_prop(2);



% create element matrices
if strcmp(typ,'EB')
    % definition of the euler bernoulli beam element stiffness matrix
    K_basic = E*I/(L^3) * ...
        [   12      6*L     -12     6*L     ; ...
            6*L     4*L^2   -6*L    2*L^2   ; ...
            -12     -6*L    12      -6*L    ; ...
            6*L     2*L^2   -6*L    4*L^2   ];
    % mass matrix
    M_basic = rho*A*L/420* ...
        [   156     22*L    54      -13*L   ;...
            22*L    4*L^2   13*L    -3*L^2  ;...
            54      13*L    156     -22*L   ;...
            -13*L   -3*L^2  -22*L   4*L^2   ];
elseif typ == 'TIM'
    % definition of the upper triangle element stiffnes matrix from Ref 1
    K_basic = E*I/((1+PHI)*(L^3)) * ...
        [   12      6*L             -12     6*L             ; ...
            0       (4+PHI)*L^2     -6*L    (2-PHI)*L^2     ; ...
            0       0               12      -6*L            ; ...
            0       0               0       (4+PHI)*L^2     ];
    % construct full matrix
    K_basic = triu(K_basic) + triu(K_basic,1)';
    % mass matrix
    M_basic = 
else
    msg = 'Wrong type defined, choose between TIM or EB as "type"';
    error(msg)
end

% transformation fo the element to the global coordinates
R = [ cos(dir) sin(dir) ; -sin(dir) cos(dir) ];

T = zeros(4,4);
T(1:2,1:2) = R;
T(3:4,3:4) = R;

K_el  = T*K_basic*T';

% implement symmetry check on elements and PSD check


    
    
    
    