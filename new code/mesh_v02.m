%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation script
% Arthur Schout
% 27/11/2017
% V0.2 06/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this script will generate two structures: element & nodes
% Element: CoG, element length, nodal connectivity, element orientation
% Nodes : position of the nodes
% two scalars : N_elements, N_nodes

%% Definition of nodal coordinates

% specify the y-coordinates of the nodes in 
B(1,1) = -61.9;
B(2,1) = -47.7;
% B(3,1) = -29.9;
% B(4,1) = -9.9;
% B(5,1) = 9.9;
% B(6,1) = 29.9;
% B(7,1) = 47.4;
% B(8,1) = 61.9;

% specify the x-location of the nodes
L(1,1) = -5.7;
L(2,1) = 15.5;
L(3,1) = 38;
L(4,1) = 65.5;
L(5,1) = 93;
L(6,1) = 120.5;
L(7,1) = 148;
L(8,1) = 175.5;
L(9,1) = 203;
L(10,1) = 228;
L(11,1) = 247.8;
L(12,1) = 283;
L(13,1) = 310.5;
L(14,1) = 335.5;
L(15,1) = 376.2;

% test mesh locations
% n_1=10;
% L = [0:1:n_1]';
% % 
% B(1,1) = 0;
% B(2,1) = 1;


% define the bounds of the bow gap
X_bound = 247.8;    % x location of the bow hole
Y_bound = 29.9;     % +/- y location of the bow hole 

%% Generate the nodes on the coordinates

% the matrix that is generated containting the nodal coordinates: Nodal is
% nx3, with the row number indicating the nodal number and the columns
% containing the coordinates: x y z

% get the lengths of the vectors
[Length, ~ ] = size(L);
[Width , ~ ] = size(B);

% creation of the matrix with the nodal coordinates 
Nodes = zeros((Length)*(Width),3);
% loop to create the nodal coordinates information for nodes x-dir, y=0
for i = 1:Length
    Nodes(i:Length:end,1) = L(i,1); % write x coordinates 
end

% loop to create the y-coordinates of the nodes
for i = 1:Width
    Nodes((i-1)*Length+1:i*Length,2) = B(i,1) * ones(Length,1); %write y-coordinates
end

%% Remove the nodes in between the bows

% remove the nodes in the bow gap
[ nodes_rem.y , ~ ] = find(abs(Nodes(:,2)) < Y_bound ); % find the row numbers x coordinates of the nodes
[ nodes_rem.x , ~ ] = find(Nodes(:,1) > X_bound); % find the row numbers y coordinate of the nodes
nodes_rem.xy = intersect(nodes_rem.x,nodes_rem.y); % find which nodes are in the bow gap

Nodes(nodes_rem.xy,:) = [ ]; 

%% calculate the element connectivity matrix
%{
% start in the corner then look for which nodes are closest in x and y
% direction. Note down these nodes: each node that is closest needs to be
% assigned an element number together with the nodal connectivity
% later the element can also be assigned: a coordinate & property that can
% be used for location dependent element properties and for location dependent updating

% the elements will be made in x and y direction consecutively
% to be implemented a length check on the elements that will be created
% an estimate of the elements can be used to create a matrix of a certain
% size to prevent the constant changing of the size of a matrix
%}

% the number of elements is not known yet
N_elements = 0;

% create elements in x-direction
for y = 1:length(B)
    % find the nodes on a y-coordinate
    [ node_nr , ~ ] = find(Nodes(:,2) == B(y,1));
    % generate element connectivity matrix
    for i = 1:length(node_nr)-1 % were looking in x-direction first, influences the definition of the location of the element
        loc_node_1 = Nodes(node_nr(i),:); % get the coordinates of node 1
        loc_node_2 = Nodes(node_nr(i+1),:);% get the coordinates of node 2
        Element_node_1(N_elements+1,:) = loc_node_1; % write the location of node 1
        Element_node_2(N_elements+1,:) = loc_node_2; % write the location of node 2
        Element_cog = (loc_node_1 + loc_node_2)/2; % find the "COG" of the element
        Element_loc(N_elements+1,:) = Element_cog; % write the COG location to the right matrix
        Element_len(N_elements+1,1) = norm(loc_node_1 - loc_node_2); % find the length of the element
        Element_conectivity(N_elements+1,:) = [ node_nr(i) node_nr(i+1) ]; % assign element connectivity
        N_elements = N_elements + 1; % number of element created
    end
    clear node_nr % removing node nr from the workspace because the number of nodes on a diffente x coordinate may change
end

% create elements in y-direction

for x = 1:length(L)
    % find the nodes on a y-coordinate
    [ node_nr , ~ ] = find(Nodes(:,1) == L(x,1));
    % generate element connectivity matrix
    for i = 1:length(node_nr)-1 % were looking in y-direction, influences the definition of the location of the element
        loc_node_1 = Nodes(node_nr(i),:); % get the coordinates of node 1
        loc_node_2 = Nodes(node_nr(i+1),:);% get the coordinates of node 2
        Element_node_1(N_elements+1,:) = loc_node_1; % write the location of node 1
        Element_node_2(N_elements+1,:) = loc_node_2; % write the location of node 2
        Element_cog = (loc_node_1 + loc_node_2)/2; % find the "COG" of the element
        Element_loc(N_elements+1,:) = Element_cog; % write the COG location to the right matrix
        Element_len(N_elements+1,1) = norm(loc_node_1 - loc_node_2); % find the length of the element
        Element_conectivity(N_elements+1,:) = [ node_nr(i) node_nr(i+1) ]; % assign element connectivity
        N_elements = N_elements + 1; % number of element created
        
    end
    clear node_nr % removin node nr from the workspace because the number of nodes on a diffente x coordinate may change
end

%% Remove elements in between the bows

% nodal connectivity matrix containing nodal connectivity,element cog location, and
% element length
Element_LOC = [ Element_conectivity Element_loc Element_len Element_node_1 Element_node_2 ];
% find the element in between the bows that don't belong
[ Element_rem.x , ~ ] = find(Element_LOC(:,3) > X_bound); % filter the x-coordinate and save row numbers
[ Element_rem.y , ~ ] = find(abs(Element_LOC(:,4)) < Y_bound ); % filter the y_coordinate and save row numbers
Element_rem.xy = intersect(Element_rem.x,Element_rem.y); % find which elements are in the bow gap

Element_LOC(Element_rem.xy,:) = []; % remove elements in between the bows

%% BETTER VARIABLE DEFINITION FOR FURTHER ANALYSIS

% element information
Element.cog.x = Element_LOC(:,3);
Element.cog.y = Element_LOC(:,4);
Element.cog.z = Element_LOC(:,5);
Element.connect = Element_LOC(:,1:2);
Element.length = Element_LOC(:,6);
Element.node1.x = Element_LOC(:,7);
Element.node1.y = Element_LOC(:,8);
Element.node1.z = Element_LOC(:,9);
Element.node2.x = Element_LOC(:,10);
Element.node2.y = Element_LOC(:,11);
Element.node2.z = Element_LOC(:,12);

% number of elements
N_elements = length(Element.length);

% location of the nodes
Node.x = Nodes(:,1);
Node.y = Nodes(:,2);
Node.z = Nodes(:,3);

% number of nodes
N_nodes = length(Node.x);

%% Find the orientation of each element

% determine the angle of a element to the reference, in this case the
% x-axis is taken as the reference vector
% the angle is given in raidans

ref_vec = [1 0 0]; % reference vector normalised

% vector of element 1 
Element.dir = [Element.node2.x - Element.node1.x ...
            Element.node2.y - Element.node1.y ...
            Element.node2.z - Element.node1.z ];

for  i = 1:N_elements
    Element.angle(i,1) = acos(dot(ref_vec,Element.dir(i,:))/(norm(Element.dir(i,:)))); % find the angel in radians of element i
    %Element.angle(i,2) = atan2(norm(cross(Element.(i,:),ref_vec),dot(Element.dir(i,:),ref_vec));
end

%% clearing up the workspace variables

% removing variables not used in further analysis
rem_var_mesh = 1;

var_rem = {'A','B','Length','Width','x','y','X_bound','Y_bound','i','L','nodes_rem'};

if  rem_var_mesh == 1
    clear(var_rem{:})
    clear -regexp ^Element_ ^loc_
    clear rem_var_mesh
    clear var_rem
end

%% plotting the nodes
% 
% figure(1)
% %plot the nodes
% scatter3(Node.x,Node.y,Node.z)
% 
% hold on
% % plot elements as points
% scatter3(Element.cog.x,Element.cog.y,Element.cog.z)
% title('nodes and element centres')
% legend('nodes','element centre')

