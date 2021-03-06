function [K_el , M_el] = element_matrix_3D(M_prop,Sec_prop,Length,Angle,N_el)
% Generation of a element matrix in the global coordinate system for 
% a 3 dimensional Timoshenko beam.
% 
% The matrices are from:
% Generating a Parametric Finite Element Model of a 3D Cantilever Timoshenko Beam Using Matlab
%
% Arthur Schout
% 07/12/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input = material & section properties, element length & orientation
% Output = Stiffness and mass matrix in the global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:    M_prop  = [ E , nu , m ];
%           Sec_prop= [ b , h   ];
%           length  = L
%           angle   =  angle in radians between the element and the x-axis
%           for now we remain in plane so only one angle is given
%           type  = 'string' 
% output:   K_el = element stiffness matrix in global coordinates
%           M_el = element mass matrix in global coordinates
% 
% constants:    k = [-] shear factor of the area, dependent on area
%                       geometry. Determines how much of the area is
%                       effective in shear.
%
% The crossection is assumed to be a square with height h and width b.
% This assumption affects the definition of:
%           1. Moment of inertia
%           2. Polar moment of inertia
%           3. Shear factor
%           4. Torsion constant
%           
% TO BE IMPLEMENTED, PROPORTIONAL DAMPING

%========= constants ============
k = 5/6;                        % [-] shear factor for square cross section

%========= parameters ===========
E = M_prop(1);                  % [N/m^2] youngs modulus
nu = M_prop(2);                 % [-] Poissons ratio
rho = 7800;                     % [kg/m^3] mass of the element
b = Sec_prop(1);                % [m] width of the beam
h = Sec_prop(2);                % [m] height of the beam
G = E/(2*(1+nu));               % [N/m^2] shear modulus
L = Length;                     % [m] length of the element
A = b*h;                        % [m^2] Area of the crossection
As_y = k*A;                     % [m^2] Area effective in shear y 
As_z = k*A;                     % [m^2] Area effective in shear z
Iy = b^3*h/12;                  % [m^4] Moment of inertia in bending in y
Iz = b*h^3/12;                  % [m^4] Moment of inertia in bending in z
Ip = 1/12*b*h*(h^2+b^2);        % [m^4] polar moment of inertia
It = min([h b])^3*max([h b])/7; % [m^4] torsion constant
Phi_y = 12*E*Iz/(G*As_y*L^2);   % [-] ratio between shear and bending in y
Phi_z = 12*E*Iy/(G*As_z*L^2);   % [-] rario between shear and bending in z
% disp(L)
% m = A*7800*L;                   

%========= Element matrix definition ========

% stiffness matrix 
% partition 1,1
K11 = zeros(6,6);
K11(1,1) = E*A/L;
K11(2,2) = 12*E*Iz/(L^3*(1+Phi_y));
K11(3,3) = 12*E*Iy/(L^3*(1+Phi_z));
K11(4,4) = G*It/L;
K11(5,5) = (4+Phi_z)*E*Iy/(L*(1+Phi_z));
K11(6,6) = (4+Phi_y)*E*Iz/(L*(1+Phi_y));
K11(3,5) = -6*E*Iy/(L^2*(1+Phi_z));
K11(5,3) = K11(3,5);
K11(2,6) = 6*E*Iz/(L^2*(1+Phi_y));
K11(6,2) = K11(2,6);
% partition 2,2
K22 = -K11 + 2*diag(diag(K11));
% partition 1,2
K21 = K11 - 2*diag(diag(K11));
K21(5,5) = (2-Phi_z)*E*Iy/(L*(1+Phi_z));
K21(6,6) = (2-Phi_y)*E*Iz/(L*(1+Phi_y));
K21(3,5) = -K21(5,3); 
K21(2,6) = -K21(6,2); 
K12 = K21';
% assembling the partitions
K_el = [K11, K12; K21 K22]; 

% mass matrix 
% partition 1,1
M11 = zeros(6,6);
M11(1,1) = 1/3;
M11(2,2) = 13/35+6*Iz/(5*A*L^2);
M11(3,3) = 13/35+6*Iy/(5*A*L^2);
M11(4,4) = Ip/3*A;
M11(5,5) = L^2/105+2*Iy/(15*A);
M11(6,6) = L^2/105+2*Iz/(15*A);
M11(5,3) = -11*L/210-Iy/(10*A*L);
M11(3,5) = M11(5,3);
M11(6,2) = 11*L/210+Iz/(10*A*L);
M11(2,6) = M11(6,2);
% partition 2,2
M22 = -M11 + 2*diag(diag(M11));
% partition 2,1
M21 = zeros(6,6);
M21(1,1) = 1/6;
M21(2,2) = 9/70-6*Iz/(5*A*L^2);
M21(3,3) = 9/70-6*Iy/(5*A*L^2);
M21(4,4) = Ip/(6*A);
M21(5,5) = -L^2/140-Iy/(30*A);
M21(6,6) = -L^2/140-Iz/(30*A);
M21(6,2) = -13*L/420-Iz/(10*A*L);
M21(2,6) = -M21(6,2);
M21(5,3) = 13*L/420-Iy/(10*A*L);
M21(3,5) = -M21(5,3);
% partition 1,2
M12 = M21';
% assembling the full matrix
M_el = rho*A*L*[M11, M12; M21, M22];

%======== Transformation of the matrices =========
% for now the only transformation that is needed is a rotation of pi/2
% about the z-axis, so it is hard coded into the transformation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the transformation matrix is given for a general rotation 
% however the input is just one angle, to solve the problem at hand
% a rotation about the z-axis suffices

% see if a rotation is needed
if Angle ~=0
    % find the transformation matrix
    T = zeros(size(M_el));
    T_11 = T_matrix(-pi); % <<<< general transformation set to pi/2
    for i=1:4
        a = 1+3*(i-1):3*i;
        T(a,a) = T_11;
    end
    % transform the element matrices into the global coordinates
    K_el = T'*K_el*T;
    M_el = T'*M_el*T;
end

%========= remove small elements from the matrices
remove_small = 1;
small_threshold = 0.1; %threshodl for which small elements should be removed

if remove_small ==1
    K_el = K_el.*(abs(K_el)>small_threshold); % remove small elements from K
    M_el = M_el.*(abs(M_el)>small_threshold); % remove small elements from M
end

%========= Checks of the matrices ==================
%%% some checks on the matrices
%%% check if element matrices are symmetric and print warning if not 
check =0;
if check ==1
    sym_K_el = isequal(K_el,K_el');
    sym_M_el = isequal(M_el,M_el');
    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
    fprintf('!!! ELEMENT: %G \n',N_el)
    
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
    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
end    