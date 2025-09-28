function Solution8_43
%
%          Example FEM code using elastic Kirchoff plate elements
%
%
%       This code solves problem 8.43 from the text
%       A.F. Bower 'Solved Problems in Mechanics of Solids'
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
%
    rad = 1.;    % Plate radius
    nr = 30;   nq = 60; % No. nodes in radial and circumferential directions
    pmag = 1;    % Pressure
    Tmag = 20;    % In-plane tension magnitude (applies biaxial membrane stress)
    h = 0.01;    % Plate thickness
    E = 1/h^3;   % Plate modulus
    nu = 0.3;    % Plate poisson ratio
    
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
    Stif = zeros(3*nnode,3*nnode);
    resid = zeros(3*nnode,1);
    T(:,1) = Tmag;
    T(:,2) = Tmag;


    for lmn=1:nelem
    %
    %   Set up the stiffness for the current element
    %
        a = connect(lmn,1);
        b = connect(lmn,2);
        c = connect(lmn,3);
        [k,r]=elstif(coord(a,:),coord(b,:),coord(c,:),E,nu,h,p3(lmn),T(lmn,:));
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
    title('Pressurized circular membrane')

    qplot = u(3:3:3*nr+3)*Tmag/(pmag*rad);

    plot(coord(1:nr+1,1),qplot,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
    

    plot(coord(1:nr+1,1),q_exact,'LineWidth',2,'LineStyle','--','Color',[0 0 0])
    
    

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

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function [kel,rel] = elstif(coorda,coordb,coordc,E,nu,h,p3,T)

    kel = zeros(9,9);
    rel = zeros(9,1);
    %  Three point Gaussian integration points and weights
    n_int_pts = 3;
    w(1:3) =  1/6;
    xi = [0.5,0.0;...
        0.0,0.5;...
        0.5,0.5;];

    D = [1,nu,0;...
        nu,1,0;...
        0,0,(1-nu)/2];
    D = D*E*h^3/(12*(1-nu^2));

    Tmx = [T(1),T(3);...
        T(3),T(2)];


    for i = 1:n_int_pts
        [N,dNdx,d2Ndx2,Area] = shape(xi(i,:),coorda,coordb,coordc);
        Bk = -transpose(d2Ndx2);
        Bk(3,:) = 2*Bk(3,:);
        Bq = transpose(dNdx);

        kel = kel + (transpose(Bk)*D*Bk + transpose(Bq)*Tmx*Bq)*w(i)*2*Area;
        rel = rel + transpose(N)*p3*w(i)*2*Area;

    end

end
function [Mout,xMout] = internalmoments(coorda,coordb,coordc,dofa,dofb,dofc,E,nu,h,smoothdata)

    
    %  Three point Gaussian integration points and weights
    n_int_pts = 3;
    w(1:3) =  1/6;
    xi = [0.5,0.0;...
        0.0,0.5;...
        0.5,0.5;];

    D = [1,nu,0;...
        nu,1,0;...
        0,0,(1-nu)/2];
    D = D*E*h^3/(12*(1-nu^2));

    xM = zeros(2,n_int_pts);
    M = zeros(3,n_int_pts);
    for i = 1:n_int_pts
        xM(:,i) = xi(i,1)*coorda + xi(i,2)*coordb + (1-xi(i,1)-xi(i,2))*coordc;
        [~,~,d2Ndx2,~] = shape(xi(i,:),coorda,coordb,coordc);
        Bk = -transpose(d2Ndx2);
        Bk(3,:) = 2*Bk(3,:);
        M(:,i) = D*Bk*[dofa;dofb;dofc];
    end
    if (smoothdata)
        Mout = sum(M,2)/3;
        xMout = sum(xM,2)/3;
    else
        Mout = M;
        xMout = xM;
    end
end
function [N,dNdx,d2Ndx2,Area] = shape(xi,coorda,coordb,coordc)

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

    d1 = sqrt((x3-x2)^2 + (y3-y2)^2);
    d2 = sqrt((x3-x1)^2 + (y3-y1)^2);
    d3 = sqrt((x2-x1)^2 + (y2-y1)^2);

    mu1 = (d3^2-d2^2)/d1^2;
    mu2 = (d1^2-d3^2)/d2^2;
    mu3 = (d2^2-d1^2)/d3^2;


    g = [L1,L2,L3,...
        L1*L2,L2*L3,L3*L1,...
        L1^2*L2 + L1*L2*L3*(3*(1-mu3)*L1 - (1+3*mu3)*L2 + (1+3*mu3)*L3)/2,...
        L2^2*L3 + L1*L2*L3*(3*(1-mu1)*L2 - (1+3*mu1)*L3 + (1+3*mu1)*L1)/2,...
        L3^2*L1 + L1*L2*L3*(3*(1-mu2)*L3 - (1+3*mu2)*L1 + (1+3*mu2)*L2)/2];

    dgdL = [eye(3);...
        L2,L1,0;...
        0,L3,L2;...
        L3,0,L1;...
        2*L1*L2,L1^2,0;...
        0,2*L2*L3,L2^2;...
        L3^2,0,2*L3*L1];

    dgdL = dgdL + ...
        1.5*[zeros(6,3);...
        2*(1-mu3)*L1*L2*L3, (1-mu3)*L1*L1*L3, (1-mu3)*L1*L1*L2;...
        (1-mu1)*L2*L2*L3, 2*(1-mu1)*L1*L2*L3, (1-mu1)*L1*L2*L2;...
        (1-mu2)*L2*L3*L3, (1-mu2)*L1*L3*L3, 2*(1-mu2)*L1*L2*L3];

    dgdL = dgdL + ...
        0.5*[zeros(6,3);...;
        (1+3*mu3)*(L3-L2)*L2*L3,  (1+3*mu3)*(L3-2*L2)*L1*L3, (1+3*mu3)*(2*L3-L2)*L1*L2;...
        (1+3*mu1)*(2*L1-L3)*L2*L3,(1+3*mu1)*(L1-L3)*L1*L3,(1+3*mu1)*(L1-2*L3)*L1*L2;...
        (1+3*mu2)*(L2-2*L1)*L2*L3,(1+3*mu2)*(2*L2-L1)*L1*L3,(1+3*mu2)*(L2-L1)*L1*L2];

    dLdx = [b(1),c(1);...
        b(2),c(2);...
        b(3),c(3)]/Delta;

    dgdx = dgdL*dLdx;

    d2gdL2 = zeros(6,6);
    d2gdL2(1,4) = 1; d2gdL2(2,6) = 1; d2gdL2(3,5) = 1;
    d2gdL2(4,1) = 2*L2; d2gdL2(4,4) = 2*L1;
    d2gdL2(5,2) = 2*L3; d2gdL2(5,6) = 2*L2;
    d2gdL2(6,3) = 2*L1; d2gdL2(6,5) = 2*L3;

    d2gdL2 = d2gdL2 + ...
        1.5*[zeros(3,6);...
        2*(1-mu3)*L2*L3,0,0,2*(1-mu3)*L1*L3,2*(1-mu3)*L1*L2,(1-mu3)*L1^2;...
        0,2*(1-mu1)*L1*L3,0,2*(1-mu1)*L2*L3,(1-mu1)*L2*L2,2*(1-mu1)*L1*L2;...
        0,0,2*(1-mu2)*L1*L2,(1-mu2)*L3*L3,2*(1-mu2)*L2*L3,2*(1-mu2)*L1*L3];

    d2gdL2 = d2gdL2 + ...
        1.5*[zeros(3,6);...
        2*(1-mu3)*L2*L3,0,0,2*(1-mu3)*L1*L3,2*(1-mu3)*L1*L2,(1-mu3)*L1^2;...
        0,2*(1-mu1)*L1*L3,0,2*(1-mu1)*L2*L3,(1-mu1)*L2*L2,2*(1-mu1)*L1*L2;...
        0,0,2*(1-mu2)*L1*L2,(1-mu2)*L3*L3,2*(1-mu2)*L2*L3,2*(1-mu2)*L1*L3];

    d2gdL2 = d2gdL2 + ...
        0.5*[zeros(3,6);...
        0,-2*(1+3*mu3)*L1*L3,2*(1+3*mu3)*L1*L2,(1+3*mu3)*(L3-2*L2)*L3,(1+3*mu3)*(2*L3-L2)*L2,2*(1+3*mu3)*(L3-L2)*L1;...
        2*(1+3*mu1)*L2*L3,0,-2*(1+3*mu1)*L1*L2,(1+3*mu1)*(2*L1-L3)*L3,2*(1+3*mu1)*(L1-L3)*L2,(1+3*mu1)*(L1-2*L3)*L1;...
        -2*(1+3*mu2)*L2*L3,2*(1+3*mu2)*L1*L3,0,2*(1+3*mu2)*(L2-L1)*L3,(1+3*mu3)*(L2-2*L1)*L2,(1+3*mu3)*(2*L2-L1)*L1];

    d2Ldx2 = [b(1)*b(1), c(1)*c(1), b(1)*c(1);...
        b(2)*b(2), c(2)*c(2), b(2)*c(2);...
        b(3)*b(3), c(3)*c(3), b(3)*c(3);...
        2*b(1)*b(2), 2*c(1)*c(2), b(1)*c(2) + b(2)*c(1);...
        2*b(1)*b(3), 2*c(1)*c(3), b(1)*c(3) + b(3)*c(1);...
        2*b(2)*b(3), 2*c(2)*c(3), b(2)*c(3) + b(3)*c(2)]/Delta^2;

    d2gdx2 = d2gdL2*d2Ldx2;

    N = zeros(9,1);
    d2Ndx2 = zeros(9,3);

    for i = 1:3
        j = i+1; if (j>3) j = j-3; end
        k = j+1; if (k>3) k = k-3; end

        nn = [g(i)-g(i+3)+g(k+3)+2*(g(i+6)-g(k+6));...
            -b(j)*(g(k+6)-g(k+3)) - b(k)*g(i+6);...
            -c(j)*(g(k+6)-g(k+3)) - c(k)*g(i+6)];
        dnn = [dgdx(i,:)-dgdx(i+3,:)+dgdx(k+3,:)+2*(dgdx(i+6,:)-dgdx(k+6,:));...
            -b(j)*(dgdx(k+6,:)-dgdx(k+3,:)) - b(k)*dgdx(i+6,:);...
            -c(j)*(dgdx(k+6,:)-dgdx(k+3,:)) - c(k)*dgdx(i+6,:)];
        d2nn = [-d2gdx2(i,:)+d2gdx2(k,:)+2*(d2gdx2(i+3,:)-d2gdx2(k+3,:));...
            -b(j)*(d2gdx2(k+3,:)-d2gdx2(k,:)) - b(k)*d2gdx2(i+3,:);...
            -c(j)*(d2gdx2(k+3,:)-d2gdx2(k,:)) - c(k)*d2gdx2(i+3,:)];

        N(3*(i-1)+1:3*i) = nn;
        dNdx(3*(i-1)+1:3*i,:) = dnn;
        d2Ndx2(3*(i-1)+1:3*i,:) = d2nn;
    end

    Area = Delta/2;

end
