function Solution8_44
%
%          Example FEM code using elastic Reissner-Mindlin plate elements
%
%
%       This code solves problem 8.44 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
    close all
    rad = 1.;    % Plate radius
    nr = 15;   nq = 30; % No. nodes in radial and circumferential directions
    pmag = 1;    % Pressure
    Tmag = 20;    % In-plane tension magnitude (applies biaxial membrane stress)
    h = 0.01;    % Plate thickness
    E = 1/h^3;   % Plate modulus
    nu = 0.3;    % Plate poisson ratio
    mu = E/(2*(1+nu));   % Shear modulus
    beta = 5/3;  % Shear coefficient
    modifyshearstiffness = false;  % Set to true to add the MIN3 shear correction

    
  % Generate the mesh
    coord = zeros(nr*nq,2);

       fixnodes = zeros(nq,3);

    for j = 1:nq
        q = (j-1)*2*pi/nq;
        for i = 1:nr
            rr = (rad/nr)*i/nr + rad*(1-1/nr)*(i/nr);
            coord((j-1)*nr+i+1,1) = rr*cos(q);
            coord((j-1)*nr+i+1,2) = rr*sin(q);
        end
        fixnodes(j,1) = j*nr+1;
        fixnodes(j,2) = 1;
    end



    [nfix,~] = size(fixnodes);

    connect = delaunay(coord);

    [nelem,~] = size(connect);
    [nnode,~] = size(coord);
    
    T = zeros(nelem,3);
    p3 = ones(nelem,1)*pmag;
    T(:,1) = Tmag;
    T(:,2) = Tmag;


    Stif = zeros(3*nnode,3*nnode);
    resid = zeros(3*nnode,1);


    for lmn=1:nelem
    %
    %   Set up the stiffness for the current element
    %
        a = connect(lmn,1);
        b = connect(lmn,2);
        c = connect(lmn,3);
        [k,r]=elstif(coord(a,:),coord(b,:),coord(c,:),E,nu,mu,beta,h,p3(lmn),T(lmn,:),modifyshearstiffness);
    %
    %   Add the current element stiffness and force vector to the global matrices
    %
        for i = 1 : 3
            for ii = 1 : 3
                rw =3*(connect(lmn,i)-1)+ii;
                resid(rw) = resid(rw) + r(3*(i-1)+ii);
                for j = 1 : 3
                    for jj = 1 : 3
                        cl = 3*(connect(lmn,j)-1)+jj;
                        Stif(rw,cl) = Stif(rw,cl) + k(3*(i-1)+ii,3*(j-1)+jj);
                    end
                end
            end
        end
    end
%   Impose boundary constraints
    for i=1:nfix
        rw=3*(fixnodes(i,1)-1)+fixnodes(i,2);
        Stif(rw,:)=0;
        Stif(rw,rw)=1.0;
        resid(rw)=fixnodes(i,3);
    end
    %
    % ================== Solve the FEM equations ===================
    %
    u=Stif\resid;

    close all
%   Plot the initial mesh
    figure0 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);     
    trisurf(connect,coord(:,1),coord(:,2),zeros(nnode,1));
    hold on
    axis square
    axis equal
    set(gca,'XColor', 'none','YColor','none')
%   Plot the deformed mesh
    trisurf(connect,coord(:,1),coord(:,2),u(1:3:3*nnode-2));
    
 %  Plot a graph comparing the FEA and analytical solutions   
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]); 
    axes2 = axes('Parent',figure1);
    hold(axes2,'on');
    set(axes2,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On');

   uplot = u(1:3:3*nr+1)*Tmag/(pmag*rad^2);
    plot(coord(1:nr+1,1),uplot,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0]);
    hold on 

    u_exact = pmag*(rad^2-coord(1:nr+1,1).^2)/(4*Tmag);
    u_exact = u_exact*Tmag/(pmag*rad^2);
    q_exact = pmag*coord(1:nr+1,1)/(2*Tmag);
    q_exact = q_exact*Tmag/(pmag*rad);

    plot(coord(1:nr+1,1),u_exact,'LineWidth',2,'LineStyle','-','Color',[0 0 0]);

    xlabel('$r/R$','FontName','Times','Interpreter','latex');
    ylabel('Deflection and rotation','FontName','Times','Interpreter','latex');
    title('Pressurized circular membrane');
    
    qplot = u(3:3:3*nr+3)*Tmag/(pmag*rad);
    plot(coord(1:nr+1,1),qplot,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
    plot(coord(1:nr+1,1),q_exact,'LineWidth',2,'LineStyle','--','Color',[0 0 0])
        
 
    if (modifyshearstiffness)
        linex = [0.1,0.25];
        liney = [0.72,0.72];
        plot(linex,liney,'LineWidth',2,'LineStyle','-','Color',[0 0 0])
        linex = [0.1,0.25];
        liney = [0.68,0.68];
        plot(linex,liney,'LineWidth',2,'LineStyle','--','Color',[0 0 0])
        linex = [0.1,0.25];
        liney = [0.6,0.6];
        plot(linex,liney,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
        % Create textbox
        annotation(figure1,'textbox',...
            [0.344642857142853 0.778571428571434 0.267857142857143 0.0870952380952407],...
            'String','Exact Solution',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % Create textbox
        annotation(figure1,'textbox',...
            [0.355357142857141 0.678571428571438 0.110714285714287 0.0870952380952407],...
            'String','FEA',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');
        % Create textarrow
        annotation(figure1,'textarrow',[0.583928571428567 0.676785714285709],...
            [0.204761904761908 0.26666666666667],'String','$T_0 u/(p R^2)$',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman');

        % Create textarrow
        annotation(figure1,'textarrow',[0.712499999999993 0.596428571428571],...
            [0.430952380952388 0.595238095238095],'String','$T_0 \theta/(p R)$',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman');
    else
        linex = [0.1,0.25];
        liney = [0.45,0.45];
        plot(linex,liney,'LineWidth',2,'LineStyle','-','Color',[0 0 0])
        linex = [0.1,0.25];
        liney = [0.4,0.4];
        plot(linex,liney,'LineWidth',2,'LineStyle','--','Color',[0 0 0])
        linex = [0.1,0.25];
        liney = [0.35,0.35];
        plot(linex,liney,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
        % Create textbox
        annotation(figure1,'textbox',...
            [0.344642857142853 0.778571428571434 0.267857142857143 0.0870952380952407],...
            'String','Exact Solution',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % Create textbox
        annotation(figure1,'textbox',...
            [0.355357142857141 0.678571428571438 0.110714285714287 0.0870952380952407],...
            'String','FEA',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman',...
            'FitBoxToText','off',...
            'EdgeColor','none');

        % Create textarrow
        annotation(figure1,'textarrow',[0.694642857142853 0.637499999999994],...
            [0.659523809523813 0.695238095238102],'String','$T_0 \theta/(p R)$',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman');

        % Create textarrow
        annotation(figure1,'textarrow',[0.646428571428568 0.73928571428571],...
            [0.335714285714288 0.39761904761905],'String','$T_0 u/(p R^2)$',...
            'Interpreter','latex',...
            'FontSize',15,...
            'FontName','Times New Roman');

    end

end

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function [kel,rel] = elstif(coorda,coordb,coordc,E,nu,mu,beta,h,p3,T,modifyshearstiffness)

kk = zeros(9,9);
kg = zeros(9,9);
kt = zeros(9,9);
rel = zeros(9,1);
%  Gaussian integration points and weights
n_int_pts = 3;

w(1:3) =  1/6;
xi = [0.5,0.0;...
    0.0,0.5;...
    0.5,0.5;];

 %   Define the element stiffness


Dg = 2*beta*h*mu*eye(2);
Dk = E*h^3/(12*(1-nu^2))*...
    [1,nu,0;...
     nu,1,0;...
     0,0,(1-nu)/2];
Tmx = [T(1),T(3);...
       T(3),T(2)];
for i = 1:n_int_pts

   [N,dNdx,L,dLdx,Area] = shape(xi(i,:),coorda,coordb,coordc);
 
    Bg = [dLdx(1,1),dNdx(1,1),dNdx(2,1)+L(1),dLdx(2,1),dNdx(3,1),dNdx(4,1)+L(2),dLdx(3,1),dNdx(5,1),dNdx(6,1)+L(3);...
          dLdx(1,2),dNdx(1,2)-L(1),dNdx(2,2),dLdx(2,2),dNdx(3,2)-L(2),dNdx(4,2),dLdx(3,2),dNdx(5,2)-L(3),dNdx(6,2)];

    Bk = zeros(3,9);
    Bk(1,3) = Bk(1,3) + dLdx(1,1);
    Bk(1,6) = Bk(1,6) + dLdx(2,1);
    Bk(1,9) = Bk(1,9) + dLdx(3,1);
    Bk(2,2) = Bk(2,2) - dLdx(1,2);
    Bk(2,5) = Bk(2,5) - dLdx(2,2);
    Bk(2,8) = Bk(2,8) - dLdx(3,2);
    Bk(3,:) = [0,-dLdx(1,1),dLdx(1,2),0,-dLdx(2,1),dLdx(2,2),0,-dLdx(3,1),dLdx(3,2)];

    Bq = [dLdx(1,1),dNdx(1,1),dNdx(2,1),dLdx(2,1),dNdx(3,1),dNdx(4,1),dLdx(3,1),dNdx(5,1),dNdx(6,1);...
          dLdx(1,2),dNdx(1,2),dNdx(2,2),dLdx(2,2),dNdx(3,2),dNdx(4,2),dLdx(3,2),dNdx(5,2),dNdx(6,2)];   
%     
      
    kg = kg + transpose(Bg)*Dg*Bg*2*w(i)*Area;
    kk = kk + transpose(Bk)*Dk*Bk*2*w(i)*Area;
    kt = kt + transpose(Bq)*Tmx*(Bq)*w(i)*2*Area;    
    rel = rel + N(1,:)*p3*2*w(i)*Area;
    
end
  kgsum = kg(2,2)+kg(3,3)+kg(5,5)+kg(6,6)+kg(8,8)+kg(9,9);
  kksum = kk(2,2)+kk(3,3)+kk(5,5)+kk(6,6)+kk(8,8)+kk(9,9); 
  alpha = kgsum/kksum;
  if (modifyshearstiffness)
      betastar = 1/(1+alpha/2);
  else
      betastar = 1;
  end
  kel = kk + betastar*kg + kt;

 end
function [Mout,Vout,xout] = internalmoments(coorda,coordb,coordc,dofa,dofb,dofc,E,nu,mu,beta,h,T,smoothdata)

%  Gaussian integration points and weights
n_int_pts = 3;

w(1:3) =  1/6;
xi = [0.5,0.0;...
    0.0,0.5;...
    0.5,0.5;];

%     Define the element stiffness
%

Dg = 2*beta*h*mu*eye(2);
Dk = E*h^3/(12*(1-nu^2))*...
    [1,nu,0;...
     nu,1,0;...
     0,0,(1-nu)/2];
Tmx = [T(1),T(3);...
       T(3),T(2)];
xM = zeros(2,n_int_pts);
M = zeros(3,n_int_pts);
V = zeros(2,n_int_pts);
for i = 1:n_int_pts
    xM(:,i) = xi(i,1)*coorda + xi(i,2)*coordb + (1-xi(i,1)-xi(i,2))*coordc; 
   [~,dNdx,L,dLdx,~] = shape(xi(i,:),coorda,coordb,coordc);
    Bg = [dLdx(1,1),dNdx(1,1),dNdx(2,1)+L(1),dLdx(2,1),dNdx(3,1),dNdx(4,1)+L(2),dLdx(3,1),dNdx(5,1),dNdx(6,1)+L(3);...
          dLdx(1,2),dNdx(1,2)-L(1),dNdx(2,2),dLdx(2,2),dNdx(3,2)-L(2),dNdx(4,2),dLdx(3,2),dNdx(5,2)-L(3),dNdx(6,2)];

    Bk = zeros(3,9);
    Bk(1,3) = Bk(1,3) + dLdx(1,1);
    Bk(1,6) = Bk(1,6) + dLdx(2,1);
    Bk(1,9) = Bk(1,9) + dLdx(3,1);
    Bk(2,2) = Bk(2,2) - dLdx(1,2);
    Bk(2,5) = Bk(2,5) - dLdx(2,2);
    Bk(2,8) = Bk(2,8) - dLdx(3,2);
    Bk(3,:) = [0,-dLdx(1,1),dLdx(1,2),0,-dLdx(2,1),dLdx(2,2),0,-dLdx(3,1),dLdx(3,2)]; 

    Bq = [dLdx(1,1),dNdx(1,1),dNdx(2,1),dLdx(2,1),dNdx(3,1),dNdx(4,1),dLdx(3,1),dNdx(5,1),dNdx(6,1);...
          dLdx(1,2),dNdx(1,2),dNdx(2,2),dLdx(2,2),dNdx(3,2),dNdx(4,2),dLdx(3,2),dNdx(5,2),dNdx(6,2)];       
    
    %     
    V(:,i) = Dg*Bg*[dofa;dofb;dofc] + Tmx*Bq*[dofa;dofb;dofc];
    M(:,i) = Dk*Bk*[dofa;dofb;dofc];
    
end
    if (smoothdata)
        Mout = sum(M,2)/n_int_pts;
        Vout = sum(V,2)/n_int_pts;
        xout = sum(xM,2)/n_int_pts;
    else
        Mout = M;
        Vout = V;
        xout = xM;
    end

 end
function [N,dNdx,L,dLdx,Area] = shape(xi,coorda,coordb,coordc)

x1 = coorda(1); y1 = coorda(2);
x2 = coordb(1); y2 = coordb(2);
x3 = coordc(1); y3 = coordc(2);

x = xi(1)*x1 + xi(2)*x2 + (1-xi(1)-xi(2))*x3;
y = xi(1)*y1 + xi(2)*y2 + (1-xi(1)-xi(2))*y3;

b = [y2-y3,y3-y1,y1-y2];
c = [x3-x2,x1-x3,x2-x1];

Delta = (b(1)*c(2)-b(2)*c(1));
L1 = ((x2-x)*(y3-y)-(x3-x)*(y2-y));
L1 = L1/Delta;
L2 = ((x3-x)*(y1-y)-(x1-x)*(y3-y));
L2 = L2/Delta;
L3 = ((x1-x)*(y2-y)-(x2-x)*(y1-y));
L3 = L3/Delta;

L = [L1,L2,L3];

dLdx = [b(1),c(1);...
    b(2),c(2);...
    b(3),c(3)]/Delta;

N = [L1,(b(2)*L3*L1-b(3)*L1*L2)/2,(c(2)*L3*L1-c(3)*L1*L2)/2,...
    L2,(b(3)*L1*L2-b(1)*L2*L3)/2,(c(3)*L1*L2-c(1)*L2*L3)/2,...
    L3,(b(1)*L2*L3-b(2)*L3*L1)/2,(c(1)*L2*L3-c(2)*L3*L1)/2;...
    0,L1,0,0,L2,0,0,L3,0;...
    0,0,L1,0,0,L2,0,0,L3];

dNdL = 0.5*[b(2)*L3-b(3)*L2, -b(3)*L1,  b(2)*L1;...
            c(2)*L3-c(3)*L2, -c(3)*L1,  c(2)*L1;...
            b(3)*L2, b(3)*L1-b(1)*L3,   -b(1)*L2;...
            c(3)*L2, c(3)*L1-c(1)*L3,   -c(1)*L2;...
           -b(2)*L3, b(1)*L3, b(1)*L2-b(2)*L1;...
           -c(2)*L3, c(1)*L3, c(1)*L2-c(2)*L1];
   
         
dNdx = dNdL*dLdx;

 Area = Delta/2;

end
%
%  =================== Function to extract variables from input file ======
%
function [nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,E,mu,e1dir,inertias,areas,...
                        nfix,fixnodes,nload,loads] = read_file(cellarray) 
%
%  Extract no. dof, no. nodes and nodal coordinates
%
n_total_nodes = 0;
n_total_elements = 0;
nbeams=str2num(cellarray{1}{2});
dum = 3;
for beam=1:nbeams
    dum = dum + 2;
    nnode=str2num(cellarray{1}{dum})

    dum=dum + 2;
    for i=1:nnode
      for k = 1:3
         coord(i,n_total_nodes + k) = str2num(cellarray{1}{dum});
         dum=dum+1;
      end
    end
    first_node(beam) = n_total_nodes + 1;
    n_total_nodes = n_total_nodes + nnode;
    %
%   Extract no. elements and connectivity
%
    dum=dum + 1;
    nelem=str2num(cellarray{1}{dum})

    dum = dum + 2;
    for i = 1 : nelem
      for j = 1 : 2
        connect(n_total_elements+i,j) = str2num(cellarray{1}{dum});
        dum=dum+1;
      end
    end
    first_element(beam) = n_total_elements + 1;
    n_total_elements = n_total_elements + nelem;    
    
% Extract no. element props and property values
   dum = dum + 1;
   E(beam) = str2num(cellarray{1}{dum});
   mu(beam) = str2num(cellarray{1}{dum+1});
   dum = dum + 3;
   for i=1:3
     e1dir(beam,i) = str2num(cellarray{1}{dum})
     dum = dum + 1;
   end
   dum = dum + 1;
   for i = 1:4
       inertias(beam,i) = str2num(cellarray{1}{dum});
       dum = dum + 1;
   end
   dum = dum + 1;
   areas(beam)  = str2num(cellarray{1}{dum});
   dum = dum + 1;
end
%
%   Extract no. nodes with prescribed displacements and the prescribed displacements
%
dum = dum + 1;
nfix=str2num(cellarray{1}{dum})
dum = dum + 4;
fixnodes = zeros(nfix,3);
for i = 1 : nfix
    for j = 1 : 3
        fixnodes(i,j) = str2num(cellarray{1}{dum});
        dum=dum+1;
    end
end
%
%   Extract no. loaded element faces, with the loads
%
dum = dum + 1;
nload=str2num(cellarray{1}{dum});
dum=dum + 3;
loads = zeros(nload,3);
for i = 1 : nload
    for j=1:7
        loads(i,j)=str2num(cellarray{1}{dum});
        dum=dum+1;
    end
end

end
function plot_beam(nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,u)
hold on

    for i = 1 : nelem
        xvals = [coord(connect(i,1),1),coord(connect(i,2),1)];
        yvals = [coord(connect(i,1),2),coord(connect(i,2),2)];
        zvals = [coord(connect(i,1),3),coord(connect(i,2),3)];
        plot3(xvals,yvals,zvals,'LineWidth',2,'Color',[0,0,0]);
    end

    scatter3(coord(:,1),coord(:,2),coord(:,3),'MarkerFaceColor',[1,0,0]);
        
    




end
