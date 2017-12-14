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
    K_part = mat2cell(K_el,[6 6],[6 6]);
    M_part = mat2cell(M_el,[6 6],[6 6]);
    % get the nodes that connect the element
    node_1 = Element.connect(i,1);
    node_2 = Element.connect(i,2);
    % find where the element matrices should be added to the global matrix
    R1 = ((node_1-1)*6+1):node_1*6;     % find the row numbers for node 1
    R2 = ((node_2-1)*6+1):node_2*6;     % find the row numbers for node 2
    C1 = R1;                            % column numbers for node 1
    C2 = R2;                            % column numbers for node 2
    % assign partitioned element matrices to the global matrix for the
    % stiffness matrix
    K_global(R1,C1) = K_global(R1,C1) + K_part{1,1};
    K_global(R2,C1) = K_global(R2,C1) + K_part{2,1};
    K_global(R1,C2) = K_global(R1,C2) + K_part{1,2};
    K_global(R2,C2) = K_global(R2,C2) + K_part{2,2};
    % mass matrix
    M_global(R1,C1) = M_global(R1,C1) + M_part{1,1};
    M_global(R2,C1) = M_global(R2,C1) + M_part{2,1};
    M_global(R1,C2) = M_global(R1,C2) + M_part{1,2};
    M_global(R2,C2) = M_global(R2,C2) + M_part{2,2};
end

