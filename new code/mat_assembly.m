%%% Assembly of the global matrices

% Material properties
E = 210e9;           % [Pa] Youngs modulus
% Beam dimensions
h = 0.1;            % height
w = 0.1;            % width

A = h*w;            % Area of the beam
m = 1;
G = 70e9;
rho = 8e3;
Ix = 1/12*h^3*w;     % moment of inertia
Iy = 1/12*w^3*h;  
nu = 0.33;           % poissons ratio

M_P = [ m , E , G , rho , nu];
Sec_prop = [ A , Ix , Iy  ];

% define the global mass and stiffness matrix
[N_nodes , ~ ] = size(Node.x);
M_g = zeros(N_nodes*4);
K_g = M_g;

for i = 1:N_elements
    %%% location dependent material properties some if loops to be
    %%% integrated later
    % get element matrices
    [K_el , M_el] = element_matrix_2D(M_P,Sec_prop,Element.length(i),Element.angle(i),'TIM');
    %%% fill up the global matrix
    % find nodes that are connected to the element
    n1 = Element.connect(i,1);
    n2 = Element.connect(i,2);
    % split the element matrix up into four parts
    K_el_part = mat2cell(K_el,[4 4],[4 4]);
    M_el_part = mat2cell(M_el,[4 4],[4 4]);
    % fill global matrices with the parts of the element matrices
    K_g(n1*4-3:n1*4,n1*4-3:n1*4) = K_g(n1*4-3:n1*4,n1*4-3:n1*4) + K_el_part{1,1};
    K_g(n2*4-3:n2*4,n1*4-3:n1*4) = K_g(n2*4-3:n2*4,n1*4-3:n1*4) + K_el_part{2,1};
    K_g(n1*4-3:n1*4,n2*4-3:n2*4) = K_g(n1*4-3:n1*4,n2*4-3:n2*4) + K_el_part{1,2};
    K_g(n2*4-3:n2*4,n2*4-3:n2*4) = K_g(n2*4-3:n2*4,n2*4-3:n2*4) + K_el_part{2,2};
    % mass matrix
    M_g(n1*4-3:n1*4,n1*4-3:n1*4) = M_g(n1*4-3:n1*4,n1*4-3:n1*4) + M_el_part{1,1};
    M_g(n2*4-3:n2*4,n1*4-3:n1*4) = M_g(n2*4-3:n2*4,n1*4-3:n1*4) + M_el_part{2,1};
    M_g(n1*4-3:n1*4,n2*4-3:n2*4) = M_g(n1*4-3:n1*4,n2*4-3:n2*4) + M_el_part{1,2};
    M_g(n2*4-3:n2*4,n2*4-3:n2*4) = M_g(n2*4-3:n2*4,n2*4-3:n2*4) + M_el_part{2,2};
    
end

%% eigen value solver

% n = 10
% 
% [U , D] = eig(K_g(1:n*4,1:n*4),M_g(1:n*4,1:n*4));
% d = diag(D)



[U , D] = eig(K_g,M_g);
%[u,d] = eigs(K_g,M_g,N_nodes*4,'LM');
DD = diag(D)
%dd = diag(d)

%% static deflection

% constrained at node 1

K_g_static = K_g;
K_g_static(1:4,:) = 0;
K_g_static(:,1:4) = 0;
K_g_static(1:4,1:4) = eye(4);

f = zeros(N_nodes*4,1);
f(end-3,1) = 1000;

Displ_nodes = K_g_static\f;


