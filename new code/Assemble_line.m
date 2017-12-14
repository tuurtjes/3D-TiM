%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the Global stiffness and mass matrix from the 
% Arthur Schout
% m files that need to be executed before running this file
% - mesh_v02.m
% - 
% Functions that are needed in this file
% - element_matrix_3D.m -> generates a element matrix in global
% coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants should be implemented inside the loop and to get element
% dependent properties
E = 200e9;
nu = 0.3;
%m = 10000*3;
b = 1;
h = 1;
M_prop  = [ E , nu ];% , m ];
Sec_prop= [ b , h ];
%Angle = 0;

% create empty global matrices depending on the number of nodes and NDOF
K_global = zeros(N_nodes*6);
M_global = K_global;

for i = 1:N_elements
    % input for the function 
    Angle = Element.angle(i);   % Angle of the element between the x-axis in the x-y plane
    L = Element.length(i);      % Length of an element
    %%%%%% different element dependent properties added here %%%%%%%
    % get the element matrix for an element, different for each element
    [ K_el, M_el ] = element_matrix_3D(M_prop,Sec_prop,L,Angle,i);
    % partition the matrix for assembly in the global matrix
    a=1+6*(i-1):6*(i+1);
    K_global(a,a) = K_global(a,a) + K_el;
    M_global(a,a) = M_global(a,a) + M_el;
end

