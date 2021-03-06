%%% plot the mode shapes

% input: the eigen vector, the frequency 
% output a plot of the mode shape

n = 0;
% nodal displacements
for j = 1:6
    j = j + n;
    
    Eig = U(:,j);
    fac = 1;
    
    %Eig = Eig.*(abs(Eig)>1e-6);
    
    Node.defx = Node.x + fac*Eig(1:6:end,1);
    Node.defy = Node.y + fac*Eig(2:6:end,1);
    Node.defz = Node.z + fac*Eig(3:6:end,1);
         
    figure(n/6+1)
    % aspect_ratio = (max(Node.x)-min(Node.x))/(max(Node.y)-min(Node.y));
    subplot(2,3,j-n)

    %clf
    hold on
    title(['Mode shape ',num2str(j),', \omega ^2 = ',num2str(D(j,j))])
    xlabel('[m]')
    ylabel('[m]')
    ax = gca;
    ax.Color ='white';
    ax.PlotBoxAspectRatioMode = 'manual';
    %ax.PlotBoxAspectRatio = [aspect_ratio 1 1];
    grid on
    view(45,35)
    for i = 1:N_elements
        % find the nodes which are connected to the elements
        node_1 = Element.connect(i,1);
        node_2 = Element.connect(i,2);
        % get the coordinates of the nodes in the deformed state
        X_el_def = [ Node.defx(node_1) ; Node.defx(node_2) ];
        Y_el_def = [ Node.defy(node_1) ; Node.defy(node_2) ];
        Z_el_def = [ Node.defz(node_1) ; Node.defz(node_2) ];
        % get the coordinates of the nodes in the undeformed state
        X_el = [ Element.node1.x(i) ; Element.node2.x(i) ];
        Y_el = [ Element.node1.y(i) ; Element.node2.y(i) ];
        Z_el = [ Element.node1.z(i) ; Element.node2.z(i) ];
        % plot the deformed shape
        plot3(X_el_def,Y_el_def,Z_el_def,'-k')
        hold on
        grid on
        % plot the undeformed shape
        plot3(X_el,Y_el,Z_el,'--b');
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
    %legend('deformed','undeformed')
    
end
