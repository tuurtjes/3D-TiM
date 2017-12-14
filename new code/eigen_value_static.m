
%find the eigen vectors
[U,D] = eig(K_global,M_global);

%%

K_g_static = K_global;
% constrained at first node
K_g_static(1:6,:) = 0;
K_g_static(:,1:6) = 0;
K_g_static(1:6,1:6) = eye(6);

f = zeros(N_nodes*6,1);
f(end-4,1) = -1e6;

Displ_nodes = K_g_static\f;

%% relative nodal displacements

nodal_displacements  = [ Displ_nodes(7:end)-Displ_nodes(1:end-6)] ;

for i =1:N_nodes-1
    a = [(i-1)*6+1:i*6];
    f_el{i} = K_el*[ zeros(6,1) ; nodal_displacements(a) ];
end


plot( 1:N_nodes, [0 ; node_rel_y]);

%% find normal stress in the elements

Rot_z_el = Displ_nodes(6+6:6:end) - Displ_nodes(6:6:end-6);
S_2 = E*(Rot_z_el)/2*b/4;

