function [K_el , M_el] = element_matrix_2D(M_prop,Sec_prop,Length,dir,typ)
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
%           dir = angle in radians between the element and the y-axis
%           for now we remain in plane so only one angle is given
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

L   = Length;
A   = Sec_prop(1);
Ix  = Sec_prop(2);
Iy  = Sec_prop(3);
E   = M_prop(2);
rho = M_prop(4);
Ir  = Ix;
nu  = 0.33;
Ks  = 5*(1+nu)/(6+5*nu); 
As  = Ks*A;
G   = E/(2*(1+nu)); 
PHI = 12*E*Ix/(G*As*L^2); % PHI should be different for bending in x and y


% create element matrices in XZ PLANE
if strcmp(typ,'EB')
    % definition of the euler bernoulli beam element stiffness matrix nodal
    % movement in in xz plane
    K_basic_x = E*Ix/(L^3) * ...
        [   12      6*L     -12     6*L     ; ...
            6*L     4*L^2   -6*L    2*L^2   ; ...
            -12     -6*L    12      -6*L    ; ...
            6*L     2*L^2   -6*L    4*L^2   ];
    % mass matrix
    M_basic_x = rho*A*L/420* ...
        [   156     22*L    54      -13*L   ;...
            22*L    4*L^2   13*L    -3*L^2  ;...
            54      13*L    156     -22*L   ;...
            -13*L   -3*L^2  -22*L   4*L^2   ];
elseif strcmp(typ,'TIM')
    % definition of the upper triangle element stiffness matrix from Ref 1
    % in xz plane
    K_basic_x = E*Ix/((1+PHI)*(L^3)) * ...
        [   12      6*L             -12     6*L             ; ...
            0       (4+PHI)*L^2     -6*L    (2-PHI)*L^2     ; ...
            0       0               12      -6*L            ; ...
            0       0               0       (4+PHI)*L^2     ];
    % construct full matrix
    K_basic_x = triu(K_basic_x) + triu(K_basic_x,1)';
    % definition of the upper triangle element mass matrix for translational
    % inertia
    M_basic_trans_x = rho*A*L/(210*(1+PHI)^2)* ...
        [   (70*PHI^2+147*PHI+78)   (35*PHI^2+77*PHI+44)*L/4    (35*PHI^2+63*PHI+27)        -(35*PHI^2+63*PHI+26)*L/4   ;...
            0                       (7*PHI^2+14*PHI+8)*L^2/4    (35*PHI^2+63*PHI+26)*L/4    -(7*PHI^2+14*PHI+6)*L^2/4   ;...
            0                       0                           (70*PHI^2+147*PHI+78)       -(35*PHI^2+77*PHI+44)*L/4   ;...
            0                       0                           0                           (7*PHI^2+14*PHI+8)*L^2/4    ];

    % construct full matrix
    M_basic_trans_x = triu(M_basic_trans_x) + triu(M_basic_trans_x,1)';

    % definition of the upper triangle element mass matrix for rotational
    % inertia
    M_basic_rot_x   = rho*Ir / (30*(1+PHI)^2*L)* ...
        [   36      -(15*PHI-3)*L           -36             -(15*PHI-3)*L           ;...
            0       (10*PHI^2+5*PHI+4)*L^2  (15*PHI-3)*L    (5*PHI^2-5*PHI-1)*L^2   ;...
            0       0                       36              (15*PHI-3)*L            ;...
            0       0                       0               (10*PHI^2+5*PHI+4)*L^2  ];
    % construct full matrix
    M_basic_rot_x = triu(M_basic_rot_x) + triu(M_basic_rot_x,1)';
    M_basic_x   = M_basic_rot_x+M_basic_trans_x;
else
    msg = 'Wrong type defined, choose between TIM or EB as "type"';
    error(msg)
end

% create element matrices in XZ PLANE
if strcmp(typ,'EB')
    % definition of the euler bernoulli beam element stiffness matrix nodal
    % movement in in xz plane
    K_basic_x = E*Ix/(L^3) * ...
        [   12      6*L     -12     6*L     ; ...
            6*L     4*L^2   -6*L    2*L^2   ; ...
            -12     -6*L    12      -6*L    ; ...
            6*L     2*L^2   -6*L    4*L^2   ];
    % mass matrix
    M_basic_x = rho*A*L/420* ...
        [   156     22*L    54      -13*L   ;...
            22*L    4*L^2   13*L    -3*L^2  ;...
            54      13*L    156     22*L    ;...
            -13*L   -3*L^2  -22*L   4*L^2   ];
elseif strcmp(typ,'TIM')
    % definition of the upper triangle element stiffness matrix from Ref 1
    % in xz plane
    K_basic_x = E*Ix/((1+PHI)*(L^3)) * ...
        [   12      6*L             -12     6*L             ; ...
            0       (4+PHI)*L^2     -6*L    (2-PHI)*L^2     ; ...
            0       0               12      -6*L            ; ...
            0       0               0       (4+PHI)*L^2     ];
    % construct full matrix
    K_basic_x = triu(K_basic_x) + triu(K_basic_x,1)';
    % definition of the upper triangle element mass matrix for translational
    % inertia
    M_basic_trans_x = rho*A*L/(210*(1+PHI)^2)* ...
        [   (70*PHI^2+147*PHI+78)   (35*PHI^2+77*PHI+44)*L/4    (35*PHI^2+63*PHI+27)        -(35*PHI^2+63*PHI+26)*L/4   ;...
            0                       (7*PHI^2+14*PHI+8)*L^2/4    (35*PHI^2+63*PHI+26)*L/4    -(7*PHI^2+14*PHI+6)*L^2/4   ;...
            0                       0                           (70*PHI^2+147*PHI+78)       -(35*PHI^2+77*PHI+44)*L/4   ;...
            0                       0                           0                           (7*PHI^2+14*PHI+8)*L^2/4    ];

    % construct full matrix
    M_basic_trans_x = triu(M_basic_trans_x) + triu(M_basic_trans_x,1)';

    % definition of the upper triangle element mass matrix for rotational
    % inertia
    M_basic_rot_x   = rho*Ir / (30*(1+PHI)^2*L)* ...
        [   36      -(15*PHI-3)*L           -36             -(15*PHI-3)*L           ;...
            0       (10*PHI^2+5*PHI+4)*L^2  (15*PHI-3)*L    (5*PHI^2-5*PHI-1)*L^2   ;...
            0       0                       36              (15*PHI-3)*L            ;...
            0       0                       0               (10*PHI^2+5*PHI+4)*L^2  ];
    % construct full matrix
    M_basic_rot_x = triu(M_basic_rot_x) + triu(M_basic_rot_x,1)';
    M_basic_x   = M_basic_rot_x+M_basic_trans_x;
else
    msg = 'Wrong type defined, choose between TIM or EB as "type"';
    error(msg)
end

% create element matrices in XZ PLANE
if strcmp(typ,'EB')
    % definition of the euler bernoulli beam element stiffness matrix nodal
    % movement in in yx plane
    K_basic_y = E*Iy/(L^3) * ...
        [   12      6*L     -12     6*L     ; ...
            6*L     4*L^2   -6*L    2*L^2   ; ...
            -12     -6*L    12      -6*L    ; ...
            6*L     2*L^2   -6*L    4*L^2   ];
    % mass matrix
    M_basic_y = rho*A*L/420* ...
        [   156     22*L    54      -13*L   ;...
            22*L    4*L^2   13*L    -3*L^2  ;...
            54      13*L    156     22*L    ;...
            -13*L   -3*L^2  -22*L   4*L^2   ];
elseif strcmp(typ,'TIM')
    % definition of the upper triangle element stiffness matrix from Ref 1
    % in yx plane
    K_basic_y = E*Iy/((1+PHI)*(L^3)) * ...
        [   12      6*L             -12     6*L             ; ...
            0       (4+PHI)*L^2     -6*L    (2-PHI)*L^2     ; ...
            0       0               12      -6*L            ; ...
            0       0               0       (4+PHI)*L^2     ];
    % construct full matrix
    K_basic_y = triu(K_basic_y) + triu(K_basic_y,1)';
    % definition of the upper triangle element mass matrix for translational
    % inertia
    M_basic_trans_y = rho*A*L/(210*(1+PHI)^2)* ...
        [   (70*PHI^2+147*PHI+78)   (35*PHI^2+77*PHI+44)*L/4    (35*PHI^2+63*PHI+27)        -(35*PHI^2+63*PHI+26)*L/4   ;...
            0                       (7*PHI^2+14*PHI+8)*L^2/4    (35*PHI^2+63*PHI+26)*L/4    -(7*PHI^2+14*PHI+6)*L^2/4   ;...
            0                       0                           (70*PHI^2+147*PHI+78)       -(35*PHI^2+77*PHI+44)*L/4   ;...
            0                       0                           0                           (7*PHI^2+14*PHI+8)*L^2/4    ];

    % construct full matrix
    M_basic_trans_y = triu(M_basic_trans_y) + triu(M_basic_trans_y,1)';

    % definition of the upper triangle element mass matrix for rotational
    % inertia
    M_basic_rot_y   = rho*Ir / (30*(1+PHI)^2*L)* ...
        [   36      -(15*PHI-3)*L           -36             -(15*PHI-3)*L           ;...
            0       (10*PHI^2+5*PHI+4)*L^2  (15*PHI-3)*L    (5*PHI^2-5*PHI-1)*L^2   ;...
            0       0                       36              (15*PHI-3)*L            ;...
            0       0                       0               (10*PHI^2+5*PHI+4)*L^2  ];
    % construct full matrix
    M_basic_rot_y = triu(M_basic_rot_y) + triu(M_basic_rot_y,1)';
    M_basic_y   = M_basic_rot_y+M_basic_trans_y;
end

%%% transformation of the element to the global coordinates
% definition of the basic transformation matrix (NO TRANSFORMATION FOR
% ANGLES)
R = [ cos(dir) sin(dir) ;...
     -sin(dir) cos(dir) ];
% assembly of the element transformation matrix
G0 = zeros(2,2);
G1 = ones(2,2);
T = [ 1 sin(dir) 0 0;...
      0 cos(dir) 0 0;...
      0 0       sin(dir) 0;...
      0 0      -cos(dir) 1];

% transform the partial element matrices
K_el_x  = T'*K_basic_x*T;
M_el_x  = T'*M_basic_x*T;
K_el_y  = T'*K_basic_y*T;
M_el_y  = T'*M_basic_y*T;

%%% create element matrix for a 2D element
%%% displacement vector: first 4 entries are the displacements at node 1
%%%                       last 4 entries are the displacements at node 2

% define K and M matrix (empty for now)
K_el = zeros(8,8);
M_el = K_el;

% get the partitioned matrices to fill up the element matrix correctly
K_part_x = mat2cell(K_el_x,[2 2],[2 2]);
K_part_y = mat2cell(K_el_y,[2 2],[2 2]);
M_part_x = mat2cell(M_el_x,[2 2],[2 2]);
M_part_y = mat2cell(M_el_y,[2 2],[2 2]);

% fill up the element matrix K
K_el(1:2,1:2) = K_part_x{1,1};
K_el(3:4,3:4) = K_part_y{1,1};
K_el(5:6,5:6) = K_part_x{2,2};
K_el(7:8,7:8) = K_part_y{2,2};
K_el(5:6,1:2) = K_part_x{2,1};
K_el(1:2,5:6) = K_part_x{1,2};
K_el(3:4,7:8) = K_part_y{1,2};
K_el(7:8,3:4) = K_part_y{2,1};
% fill up element matrix M
M_el(1:2,1:2) = M_part_x{1,1};
M_el(3:4,3:4) = M_part_y{1,1};
M_el(5:6,5:6) = M_part_x{2,2};
M_el(7:8,7:8) = M_part_y{2,2};
M_el(5:6,1:2) = M_part_x{2,1};
M_el(1:2,5:6) = M_part_x{1,2};
M_el(3:4,7:8) = M_part_y{1,2};
M_el(7:8,3:4) = M_part_y{2,1};



%%% some checks on the matrices
%%% check if element matrices are symmetric and print warning if not UNNNECCESARY DUE TO ALTERNATIVE WAY OF CONSTRUCTING ELEMENT MATRICES
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



    
    
    
    