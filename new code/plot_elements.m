%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the elements
%%% Arthur Schout
%%% 03/12/2017
%%% v_02 - 07/12/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aspect_ratio = (max(Node.x)-min(Node.x))/(max(Node.y)-min(Node.y));

figure(1)
clf
hold on
title('Mesh of the vessel')
xlabel('[m]')
ylabel('[m]')
ax = gca;
ax.Color ='white';
ax.PlotBoxAspectRatioMode = 'manual';
%ax.PlotBoxAspectRatio = [aspect_ratio 1 1];
grid on

% plot the elements
for i = 1:N_elements
    X_el = [ Element.node1.x(i) ; Element.node2.x(i) ];
    Y_el = [ Element.node1.y(i) ; Element.node2.y(i) ];
    Z_el = [ Element.node1.z(i) ; Element.node2.z(i) ];
    plot3(X_el,Y_el,Z_el,'-k')
end

% plot node numbers
for i = 1:N_elements
    X_el = Element.cog.x(i);
    Y_el = Element.cog.y(i);
    Z_el = Element.cog.z(i);
    text(X_el,Y_el,Z_el,num2str(i),'Color','black');
end

% plot nodes
scatter3(Node.x,Node.y,Node.z,'r')

% plot node numbers
for i = 1:N_nodes
    text(Node.x(i)+1,Node.y(i)+1,Node.z(i)+0.1,num2str(i),'Color','red');
end
