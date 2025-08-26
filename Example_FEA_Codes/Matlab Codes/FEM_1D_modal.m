function FEM_1D_modal
%                  EXAMPLE MODAL DYNAMIC FEM CODE FOR 1D PROBLEM
%    Calculates displacement as function of time for a 1D bar subject to
%    traction on its right end and fixed on its left end.   Time
%    integration is accomplished using modal dynamics.
%
%       This code is an example from the text
%       A.F. Bower 'Applied Mechanics of Solids' (2nd ed.) 
%       CRC press, Baton Rouge, 2009 
%
%   It was downloaded from
%   https://github.com/albower/Applied_Mechanics_of_Solids
%
%
        close all
%%  Define values for parameters below
%
%    Length and x-sect area of bar
        Length = 5.;
        A = 1.; 
%       Material props - rho is mass density, mu is shear modulus and nu is
%       Poisson's ratio
        rho = 100.;
        mu = 50.;
        nu = 0.;
%
        const = 2*mu*A*(1-nu)/(1-2*nu);
%       Loading
        bodyforce = 0.;
        traction = 10.;
%    total no. elements, no. nodes on each element (2 for linear, 3 for quadratic elemnts), total no. nodes; 
        L = 31;
        Ne = 2;
        nnodes = (Ne-1)*L+1;
%
%    Initial conditions
%
        un = zeros(nnodes,1);
        vn = zeros(nnodes,1);
%
%   Time stepping
         dt = 0.1;
         nsteps = 800;
         lumpedmass = false;
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
       npointsK = Ne-1;
       npointsM = Ne;
       if (npointsK == 2)
         wK = [1,1];
         xiK = [-0.5773502692,0.5773502692];
       elseif (npointsK == 1)
         wK = [2.,0.];
         xiK = [0.,0.];
       end

     if (lumpedmass)
         if (npointsM == 3)
           wM = [2/3.,2/3.,2/3.];
           xiM = [-1.,0.,1.];
         elseif (npointsM == 2)
           wM = [1,1];
           xiM = [-1.,1.];
         end
      else
         if (npointsM == 3)
           wM = [0.55555555,0.8888888888,0.5555555555];
           xiM = [-0.7745967,0.,0.7745967];
         elseif (npointsM == 2)
           wM = [1,1];
           xiM = [-0.5773502692,0.5773502692];
         end
     end

%          
%
%%     Assemble the global stiffness, mass and force vector
%
       M = zeros(nnodes,nnodes);
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
         mel = zeros(Ne,Ne);
         kel = zeros(Ne,Ne);
         fel = zeros(Ne,1);
%
         for II = 1 : npointsK
%
%        Compute N and dN/dxi at the current integration point
%
           N = zeros(1,Ne);
           dNdxi = zeros(1,Ne);
           if (Ne == 3) 
             N(1) = -0.5*xiK(II)*(1.-xiK(II));
             N(2) = 0.5*xiK(II)*(1.+xiK(II));
             N(3) = (1.-xiK(II)^2);
             dNdxi(1) = -0.5+xiK(II);
             dNdxi(2) = 0.5+xiK(II);
             dNdxi(3) = -2.*xiK(II);
           elseif (Ne == 2)
             N(1) = 0.5*(1.-xiK(II));
             N(2) = 0.5*(1.+xiK(II));
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
             fel(a) = fel(a) + wK(II)*bodyforce*J*N(a);
             for b = 1 : Ne 
               kel(a,b) = kel(a,b) + const*wK(II)*J*dNdx(a)*dNdx(b);
             end
          end
% 
         end
%
%        Assemble the mass matrix
%
         for II = 1 : npointsM
%
%        Compute N and dN/dxi at the current integration point
%
           if (Ne == 3) 
             N(1) = -0.5*xiM(II)*(1.-xiM(II));
             N(2) = 0.5*xiM(II)*(1.+xiM(II));
             N(3) = (1.-xiM(II)^2);
             dNdxi(1) = -0.5+xiM(II);
             dNdxi(2) = 0.5+xiM(II);
             dNdxi(3) = -2.*xiM(II);
           elseif (Ne == 2)
             N(1) = 0.5*(1.-xiM(II));
             N(2) = 0.5*(1.+xiM(II));
             dNdxi(1) = -0.5;
             dNdxi(2) = 0.5;
           end
% 
%        Compute dx/dxi, J and dN/dx
%
            dxdxi = dot(dNdxi,lmncoords);
            J = abs(dxdxi);
%  
%         Add contribution to element stiffness and force vector from current integration pt
%
          for a = 1 : Ne
             for b = 1 : Ne 
               mel(a,b) = mel(a,b) + rho*wM(II)*J*N(a)*N(b);
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
             M(rw,cl) = M(rw,cl) + mel(a,b);
             K(rw,cl) = K(rw,cl) + kel(a,b);
           end
         end
       end

%     Add the extra forcing term from the traction at x=L
%
       F(nnodes) = F(nnodes) + traction;
%       
%      Modify FEM equations to enforce displacement boundary condition
%      To do this we simply replace the equation for the first node with u=0
%      and set the first column of M and K to zero
%
       for a = 1 : nnodes
         M(1,a) = 0.;
         M(a,1) = 0;
         K(1,a) = 0.;
         K(a,1) = 0;
       end
       M(1,1) = 1;
       F(1) = 0.;
%
%%     Perform the modal decomposition
%
%      Matlab has lots of matrix functions we can use
%
       if (lumpedmass)
         rootM = sqrt(M);
         inverserootM = zeros(nnodes);
         for i = 1 : nnodes
           inverserootM(i,i) = 1/rootM(i,i);
         end
       else
         rootM = sqrtm(M);
         inverserootM = inv(rootM);
       end
       
       H = inverserootM*(K*inverserootM);
%      svd is the singular value decomposition we need
       [Q,Lambda,~] = svd(H);
%
%      Lambda has eigenvalues in decreasing order. The last one is a rigid
%      body mode
       
       u = inverserootM*Q(:,nnodes-1);
       plot(coords,u)
%
       Fmodal = transpose(Q)*(inverserootM*F);
       
       nprint = 4;
       u = zeros(nnodes,1);
       figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
       axes1 = axes('Parent',figure1,...
           'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
       hold(axes1,'off');
       set(axes1,'FontName','Times','FontSize',16,...
           'LineWidth',2,'XGrid','on','YGrid','on',...
           'Box','On');
       plot(coords,u);
       axis([0,Length,0,1]);
       xlabel('Position $x$','Interpreter','latex','FontSize',16);
       ylabel('Displacement $u$','Interpreter','latex','FontSize',16);
       count = 1;
       vxy = zeros(2,nsteps);
       for i = 1 : nsteps
         t = (i-1)*dt;
         for n = 1 : nnodes
            if (Lambda(n,n) >0) 
            u(n) = Fmodal(n)*(1.-cos(sqrt(Lambda(n,n))*t))/Lambda(n,n);
            else
              u(n) = 0.;
            end
         end
         u = inverserootM*(Q*u);

%      Save the end displacement for an x-y plot
         vxy(1,i) = t;
         vxy(2,i) = u(nnodes);
%
%      Save plots of displacements as a function of x for an animation
%
         if (count==nprint)
%            Replace the next line with plot(coords,vn1) to plot the
%            velocity
             plot(coords,u,'Color','r','LineWidth',2);
             axis([0,Length,0,1]);
             xlabel('Position $x$','Interpreter','latex','FontSize',16)
             ylabel('Displacement $u$,','Interpreter','latex','FontSize',16')
             pause(0.1)
             count=0;
         end
         count = count + 1;
       end
       figure2 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
       axes2 = axes('Parent',figure2,...
           'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
       hold(axes2,'on');
       set(axes2,'FontName','Times','FontSize',16,...
           'LineWidth',2,'XGrid','on','YGrid','on',...
           'Box','On');
       plot(vxy(1,:),vxy(2,:),'r - ','LineWidth',2)
       xlabel('Time $t$','Interpreter','latex','FontSize',16)
       ylabel('Displacement $u$','Interpreter','latex','FontSize',16)
       title('Displacement of 1D bar','FontSize',16)
end