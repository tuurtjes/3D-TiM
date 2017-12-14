function [K_el , M_el] = element_matrix_1D(M_prop,Sec_prop,Length,dir,typ)
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


% input:    M_prop = [ m , E , G , rho ];
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
rho = M_prop(4);
Ir = 1;


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
            54      13*L    156     22*L    ;...
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
    % definition of the upper triangle element mass matrix for translational
    % inertia
    M_basic_trans = rho*A*L/(210*(1+PHI)^2)* ...
        [   (70*PHI^2+147*PHI+78)   (35*PHI^2+77*PHI+44)*L/4    (35*PHI^2+63*PHI+27)        -(35*PHI^2+63*PHI+26)*L/4   ;...
            0                       (7*PHI^2+14*PHI+8)*L^2/4    (35*PHI^2+63*PHI+26)*L/4    -(7*PHI^2+14*PHI+6)*L^2/4   ;...
            0                       0                           (70*PHI^2+147*PHI+78)       -(35*PHI^2+77*PHI+44)*L/4   ;...
            0                       0                           0                           (7*PHI^2+14*PHI+8)*L^2/4    ];

    % construct full matrix
    M_basic_trans = triu(M_basic_trans) + triu(M_basic_trans,1)';

    % definition of the upper triangle element mass matrix for rotational
    % inertia
    M_basic_rot   = rho*Ir / (30*(1+PHI)^2*L)* ...
        [   36      -(15*PHI-3)*L           -36             -(15*PHI-3)*L           ;...
            0       (10*PHI^2+5*PHI+4)*L^2  (15*PHI-3)*L    (5*PHI^2-5*PHI-1)*L^2   ;...
            0       0                       36              (15*PHI-3)*L            ;...
            0       0                       0               (10*PHI^2+5*PHI+4)*L^2  ];
    % construct full matrix
    M_basic_rot = triu(M_basic_rot) + triu(M_basic_rot,1)';
    M_basic   = M_basic_rot+M_basic_trans;
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
M_el  = T*M_basic*T';
% implement symmetry check on elements and PSD check
%%% some sanity checks on the matrices
% check if element matrices are symmetric and print warning if not UNNNECCESARY DUE TO ALTERNATIVE WAY OF CONSTRUCTING ELEMENT MATRICES
sym_K_el = isequal(K_el,K_el');
sym_M_el = isequal(M_el,M_el');

if sym_K_el ~= 1
    fprintf('!!!\n!!!\tElement stiffness matrix is not symmetric\n!!!\n');
end
if sym_M_el ~= 1
    fprintf('!!!\n!!!\tElement mass matrix is not symmetric\n!!!\n');
end

% check if element matrices are positive definite
[ ~ , PSD_K] = chol(K_el);
[ ~ , PSD_M] = chol(M_el);

if PSD_K ~= 0
    fprintf('!!!\n!!!\tElement Stiffness matrix is not PSD\n!!!\n')
end

if PSD_M ~= 0
    fprintf('!!!\n!!!\tElement Mass matrix is not PSD\n!!!\n')
end

    
    
    
    