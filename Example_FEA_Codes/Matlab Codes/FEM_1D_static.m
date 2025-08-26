function FEM_1D_static
%
%              EXAMPLE STATIC FEM CODE FOR A SIMPLE 1D PROBLEM
%    Calculates static displacement for a 1D bar subject to uniform body
%    force, constrained so u=0 at x=0 and subjected to prescribed traction
%    at x=L
%
%       This code is an example from the text
%       A.F. Bower 'Applied Mechanics of Solids' (2nd ed.) 
%       CRC press, Baton Rouge, 2026
%
%   It was downloaded from
%   https://github.com/albower/Applied_Mechanics_of_Solids
%
%%  Define values for parameters below
%
%   Length and x-sect area of bar
    Length = 5.;
    A = 1.;
%       Material props
    mu = 50.; %Shear modulus
    nu = 0.3; % Poissons ratio
    %
    const = 2*mu*A*(1-nu)/(1-2*nu);
    %       Loading
    bodyforce = 10.;
    traction = 2.;
%    total no. elements, no. nodes on each element (2 for linear, 3 for quadratic elemnts), total no. nodes;
    L = 10;
    Ne = 3;
    nnodes = (Ne-1)*L+1;
%
%%     Set up some data structures storing the mesh
%
    coords = zeros(1,nnodes);
    for i= 1 : nnodes
        coords(i) = Length*(i-1)/(nnodes-1);
    end
%
%   Element connectivity (speficies node numbers on each element)
%
    connect = zeros(Ne,L);
    for lmn=1:L
        if (Ne==3)
            connect(1,lmn) = 2*lmn-1;
            connect(3,lmn) = 2*lmn;
            connect(2,lmn) = 2*lmn+1;
        elseif (Ne == 2)
            connect(1,lmn) = lmn;
            connect(2,lmn) = lmn+1;
        end 
    end
%
%    Integration points and weights for 2 point integration
%
    npoints = Ne-1;
    if (npoints == 2)
        w = [1,1];
        xi = [-0.5773502692,0.5773502692];
    elseif (npoints == 1)
        w = [2.,0.];
        xi = [0.,0.];
    end
%
%%     Assemble the global stiffness and force vector
%
    K = zeros(nnodes,nnodes);
    F = zeros(nnodes,1);
    %
    for lmn = 1 : L
    %
    %       Extract the coords of each node on the current element
    %
        lmncoords = zeros(Ne,1);
        for a = 1 : Ne
            lmncoords(a) = coords(connect(a,lmn));
        end
    %
    %      For the current element, loop over integration points and assemble element stiffness
    %
        kel = zeros(Ne,Ne);
        fel = zeros(Ne,1);
        %
        for II = 1 : npoints
        %
        %        Compute N and dN/dxi at the current integration point
        %
            N = zeros(1,Ne);
            dNdxi = zeros(1,Ne);
            if (Ne == 3)
                N(1) = -0.5*xi(II)*(1.-xi(II));
                N(2) = 0.5*xi(II)*(1.+xi(II));
                N(3) = (1.-xi(II)^2);
                dNdxi(1) = -0.5+xi(II);
                dNdxi(2) = 0.5+xi(II);
                dNdxi(3) = -2.*xi(II);
            elseif (Ne == 2)
                N(1) = 0.5*(1.-xi(II));
                N(2) = 0.5*(1.+xi(II));
                dNdxi(1) = -0.5;
                dNdxi(2) = 0.5;
            end
        %
        %        Compute dx/dxi, J and dN/dx
        %
            dxdxi = dot(dNdxi,lmncoords);
            J = abs(dxdxi);
            dNdx = dNdxi/dxdxi;
        %
        %         Add contribution to element stiffness and force vector from current integration pt
        %
            for a = 1 : Ne
                fel(a) = fel(a) + w(II)*bodyforce*J*N(a);
                for b = 1 : Ne
                    kel(a,b) = kel(a,b) + const*w(II)*J*dNdx(a)*dNdx(b);
                end
            end
            %
        end
    %
    %       Add the stiffness and residual from the current element into global matrices
    %
        for a = 1 : Ne
            rw = connect(a,lmn);
            F(rw) = F(rw) + fel(a);
            for b = 1 : Ne
                cl = connect(b,lmn);
                K(rw,cl) = K(rw,cl) + kel(a,b);
            end
        end
    end
%     Add the extra forcing term from the traction at x=L
%
    F(nnodes) = F(nnodes) + traction;
%      Modify FEM equations to enforce displacement boundary condition
%      To do this we simply replace the equation for the first node with u=0
%
    for a = 1 : nnodes
        K(1,a) = 0.;
    end
    K(1,1) = 1.;
    F(1) = 0.;
%
%%    Solve the equations
%
    u = K\F;
%
%%    Plot the displacement field
%
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes1 = axes('Parent',figure1,...
        'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes1,'on');
    set(axes1,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On');
    hold on
    plot(coords,u,'r - s','LineWidth',2)
    xlabel('Distance $x$','Interpreter','latex','FontSize',16)
    ylabel('Displacement $u$','Interpreter','latex','FontSize',16)
    title('Displacement of 1D bar','FontSize',16)
end
  