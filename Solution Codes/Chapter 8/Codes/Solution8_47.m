function Solution8_47

%
%
%          Example 2D and 3D linear elastic FEM code, extended to
%          include cohesive zone elements.
%
%       This code solves problem 8.47 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
%
%%
%        Variables read from input file;
%        nprops              No. material parameters
%        materialprops(i)    List of material parameters
%        ncoord              No. spatial coords (2 for 2D, 3 for 3D)
%        ndof                No. degrees of freedom per node (2 for 2D, 3 for 3D)
%                            (here ndof=ncoord, but the program allows them to be different
%                            to allow extension to plate & beam elements with C^1 continuity)
%        nelnodes               No. nodes
%        coords(i,j)         ith coord of jth node, for i=1..ncoord; j=1..nelnodes
%        nelem               No. elements
%        n_solid_elem        No. standard continuum elements (the remainder
%                            are cohesives zones
%        maxnodes            Max no. nodes on any one element (used for array dimensioning)
%        nelnodes(i)         No. nodes on the ith element
%        elident(i)          An integer identifier for the ith element.
%                            For 2D elements elident=0 means plane stress
%                                            elident=1 means plane strain
%        connect(i,j)        List of nodes on the jth element
%        nfix                Total no. prescribed displacements
%        fixnodes(i,j)       List of prescribed displacements at nodes
%                            fixnodes(1,j) Node number
%                            fixnodes(2,j) Displacement component number (1, 2 or 3)
%                            fixnodes(3,j) Value of the displacement
%        ndload              Total no. element faces subjected to tractions
%        dloads(i,j)         List of element tractions
%                            dloads(1,j) Element number
%                            dloads(2,j) face number
%                            dloads(3,j), dloads(4,j), dloads(5,j) Components of traction
%                            (assumed uniform)
% 
%
%% ==================== Read data from the input file ===========================
%
   filename = '../Input Files/Solution8_47.txt';
   while (~isfile(filename))
        f = msgbox({'The input file was not found';...
                    'Please use the browser to select the file' });
        uiwait(f);
        [filename,location] = uigetfile('*.txt');
        filename = strcat(location,filename);
    end
    
    infile=fopen(filename,'r');
%   Element data will be written to the file specified below
    outfile=fopen('FEM_results.txt','w');  

    [nprops,materialprops,ncoord,ndof,nnode,coords,n_solid_elem,nelem,maxnodes,connect,nelnodes,elident,nfix,fixnodes,ndload,dloads] = read_input_file(infile);

    if (nelem>n_solid_elem)
       npoints = numberofintegrationpoints(ncoord-1,nelnodes(n_solid_elem+1)/2);
       Damage0 = zeros(npoints,nelem-n_solid_elem);  %Cohesive zone damage at start of step
       Damage1 = Damage0;  %CZ damage at end of step
    end
    fclose(infile);

    close all

%
%============================ MAIN FEM ANALYSIS PROCEDURE ========================
%
%   w           Nodal displacements.  Let w_i^a be ith displacement component
%               at jth node.  Then dofs contain (w_1^1, w_2^1, w_1^2, w_2^2....) for 2D
%               and (w_1^1, w_2^1, w_3^1, w_1^2, w_2^2, w_3^2....) for 3D
%   dw          Correction to nodal displacements.  Let w_i^a be ith displacement component
%               at jth node.  Then dofs contain (w_1^1, w_2^1, w_1^2, w_2^2....) for 2D
%               and (w_1^1, w_2^1, w_3^1, w_1^2, w_2^2, w_3^2....) for 3D
%   K           Global stiffness matrix.  Stored as [K_1111 K_1112  K_1121  K_1122...
%                                                    K_1211 K_1212  K_1221  K_1222...
%                                                    K_2111 K_2112  K_2121  K_2122...]
%               for 2D problem and similarly for 3D problem
%   F           Force vector.  Currently only includes contribution from tractions
%               acting on element faces (i.e. body forces are neglected)
%   R           Volume contribution to residual
%   b           RHS of equation system
%
%
%  Here we specify how the Newton Raphson iteration should run
%  Load is applied in nsteps increments;
%  tol is the tolerance used in checking Newton-Raphson convergence
%  maxit is the max no. Newton-Raphson iterations
%  relax is the relaxation factor (Set to 1 unless big convergence problems)
%
  nsteps = 100;
  dt = 1/nsteps;
  tol = 0.0001;
  maxit = 30;
  relax = 1.;

  forcevdisp = zeros(2,2*nsteps+2,2);


    w = zeros(nnode*ndof,1);
    u = w;
    
    for cases = 1:2
        if (cases==2)
            materialprops(3) = 6;
            Damage0 = zeros(npoints,nelem-n_solid_elem);
        end
    for step = 1 : 2*nsteps+2

        loadfactor = step/nsteps;
        if (step>nsteps)   % Unloading
            loadfactor = 2-step/nsteps;
        end

        err1 = 1.;
        nit = 0;

        fprintf(1,'\n Step %f Load %f\n',step,loadfactor);

        while ((err1>tol) && (nit<maxit))          % Newton Raphson loop
            nit = nit + 1;

            [R,K,Damage1] = globalmatrices(dt,ncoord,ndof,nnode,coords,n_solid_elem, ...
                nelem,maxnodes,elident,nelnodes,connect,materialprops,u,w,Damage0);
            F = globaltraction(ncoord,ndof,nnode,ndload,coords, ...
                nelnodes,elident,connect,dloads,w);

            b = loadfactor*F - R;

            %          Fix constrained nodes.
            for n = 1 : nfix
                rw = ndof*(fixnodes(1,n)-1) + fixnodes(2,n);
                for cl = 1 : ndof*nnode
                    K(rw,cl) = 0;
                end
                K(rw,rw) = 1.;
                b(rw) = loadfactor*fixnodes(3,n)-w(rw)-u(rw);
            end
            %
            %        Solve for the correction
            %
            dw = K\b;
            %
            %        Check convergence
            %

            w = w + relax*dw;
            wnorm = dot(w,w);
            err1 = dot(dw,dw);
            err2 = dot(b,b);
            err1 = sqrt(err1/wnorm);
            err2 = sqrt(err2)/(ndof*nnode);
            fprintf(1,'Iteration number %d Correction %f Residual %f tolerance %f\n',nit,err1,err2,tol);
        end
%
%      Store traction and displacement for plotting later
    %
        u = u + w;
        Damage0 = Damage1;
        s33 = (u(end)-u(27))*materialprops(1);
        forcevdisp(2,step+1,cases) = s33/materialprops(5) ; 
        forcevdisp(1,step+1,cases) = u(end)*materialprops(3)/materialprops(5);  

    end
    end
    
    close all
    
    %
    % Create a plot of the deformed mesh
    defcoords = zeros(ndof,nnode);
    scalefactor = 1.0;
    for i = 1:nnode
        for j = 1:ndof
            defcoords(j,i) = coords(j,i) + scalefactor*u(ndof*(i-1)+j);
        end
        defcoords(1,i) = defcoords(1,i) + 1.5;
    end
    figure
    % plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g');
    plotmesh(defcoords,ncoord,nnode,connect,nelem,elident,nelnodes,'k');

    view(gca,[-25.3109079280651 7.11656316575341]);
    set(gca,'XColor', 'none','YColor','none','ZColor','none')
    
    
    
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes2 = axes('Parent',figure1);
    hold(axes2,'on');
    set(axes2,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On');
    plot(forcevdisp(1,:,1),forcevdisp(2,:,1),'LineWidth',2,'Color',[0,0,0],'DisplayName','$u_{max} k_0/\sigma_{max}=2.5$');
    plot(forcevdisp(1,:,2),forcevdisp(2,:,2),'LineWidth',2,'LineStyle','--','Color',[0,0,0],'DisplayName','$u_{max} k_0/\sigma_{max}=1.5$');
    ylabel('Stress $\sigma/\sigma_{max}$','FontName','Times','Interpreter','latex');
    xlabel('Displacement $u k_0/\sigma_{max}$','FontName','Times','Interpreter','latex');
    hold on
           legend1 = legend(axes2,'show');
       set(legend1,...
           'Position',[0.482578202828691 0.23571428571429 0.359236709232172 0.206190474146893],...
           'Interpreter','latex',...
           'FontSize',14,...
           'FontName','Times New Roman');

  
    
end


function [rel,kel] = elforcnstif(ncoord,ndof,nelnodes,elident,coord,materialprops,displacement)
%
%  Assemble the element stiffness matrix and residual force vector
%
%    Arguments;
%
%      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
%      ndof               No. degrees of freedom per node (often ndof = ncoord)
%      nelnodes           No. nodes on the element
%      elident            Element identifier (not used here - for future enhancements!)
%      coords(i,a)        ith coord of ath node
%      materialprops      Material properties passed on:constitutive procedures
%      displacement(i,a)  ith displacement component at ath node
%
%   Local variables
%      npoints            No. integration points
%      xi(i,inpt)         ith local coord of integration point no. intpt
%      w(intpt)           weight for integration point no. intpt
%      N(a)               Shape function associated with ath node on element
%      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
%      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
%      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
%      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
%      det                Determinant of jacobian
%      strain(i,j)        strain_ij components
%      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
%      kel(row,col)       Rows && cols of element stiffness
%
%
   E = materialprops(1);
   nu = materialprops(2);

   npoints = numberofintegrationpoints(ncoord,nelnodes,elident);
   rel = zeros(ndof*nelnodes,1);
   kel = zeros(ndof*nelnodes,ndof*nelnodes);

   xilist = integrationpoints(ncoord,nelnodes,npoints,elident);
   w = integrationweights(ncoord,nelnodes,npoints,elident);
   %
   uvec = reshape(displacement,1,ndof*nelnodes);
 
      
   for intpt = 1 : npoints

       xi(1:ncoord) = xilist(1:ncoord,intpt);
 
      dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi);

      dxdxi = coord*dNdxi;
      dxidx = inv(dxdxi);
      dt = det(dxdxi);
      
      dNdx = dNdxi*dxidx;


      if (ndof==2)
        B = zeros(3,2*nelnodes);
        B(1,1:2:2*nelnodes-1) = dNdx(1:nelnodes,1);
        B(2,2:2:2*nelnodes) = dNdx(1:nelnodes,2);
        B(3,1:2:2*nelnodes-1) = dNdx(1:nelnodes,2);
        B(3,2:2:2*nelnodes) = dNdx(1:nelnodes,1);
        D = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];  % Default is plane stress
        if (elident==1)
           D = E/(1+nu)/(1-2*nu)*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
        end       
      else
        B = zeros(6,3*nelnodes);
        B(1,1:3:3*nelnodes-2) = dNdx(1:nelnodes,1);
        B(2,2:3:3*nelnodes-1) = dNdx(1:nelnodes,2);
        B(3,3:3:3*nelnodes)   = dNdx(1:nelnodes,3);
        B(4,1:3:3*nelnodes-2) = dNdx(1:nelnodes,2);
        B(4,2:3:3*nelnodes-1) = dNdx(1:nelnodes,1);
        B(5,1:3:3*nelnodes-2) = dNdx(1:nelnodes,3);
        B(5,3:3:3*nelnodes)   = dNdx(1:nelnodes,1);
        B(6,2:3:3*nelnodes-1) = dNdx(1:nelnodes,3);
        B(6,3:3:3*nelnodes)   = dNdx(1:nelnodes,2);       
        D = eye(6)*(1-2*nu)/2;
        D(1:3,1:3) = nu;
        D([1,8,15]) = (1-nu);
        D = E/(1+nu)/(1-2*nu)*D;
      end           
 
      strainvec = B*uvec';
      stressvec = D*strainvec;   

      rel = rel + transpose(B)*stressvec*w(intpt)*dt;
      kel = kel + (transpose(B)*D*B)*w(intpt)*dt;
   end

end

function [r,kel,Damage1] = elforcnstif_cz(dt,ncoord,ndof,nnode,elident,coord,materialprops,displacement,displacement_increment,Damage0)

  npoints = numberofintegrationpoints(ncoord-1,nnode/2);
  xi = zeros(ncoord-1,1);
  r = zeros(ndof*nnode,1);
  kel = zeros(ndof*nnode,ndof*nnode);
  B = zeros(ndof,ndof*nnode);
  Damage1 = zeros(npoints,1);
   
  xilist = integrationpoints(ncoord-1,nnode/2,npoints);
  w = integrationweights(ncoord-1,nnode/2,npoints);

  for intpt = 1:npoints

    for i = 1:ncoord-1
      xi(i) = xilist(i,intpt);
    end

    N = shapefunctions(nnode/2,ncoord-1,elident,xi);
    dNdxi = shapefunctionderivs(nnode/2,ncoord-1,elident,xi);
%
%     Compute the jacobian matrix && its determinant
%
    dxdxi = coord(:,1:nnode/2)*dNdxi;
    
    if (ncoord == 2)
      nn = [dxdxi(2,1),-dxdxi(1,1)];
      e1 = dxdxi(:,1)/norm(dxdxi(:,1));
      J = norm(nn);
      nn = nn/J;
      for a = 1:nnode/2
         B(1,2*a-1:2*a) = -N(a)*nn(:);
         B(2,2*a-1:2*a) = -N(a)*e1(:);
         B(1,nnode + 2*a-1:nnode+2*a) = N(a)*nn(:);
         B(2,nnode + 2*a-1:nnode+2*a) = N(a)*e1(:);
      end
    elseif (ncoord == 3) 
      nn = [dxdxi(2,1)*dxdxi(3,2)-dxdxi(2,2)*dxdxi(3,1), ...
            dxdxi(1,1)*dxdxi(3,2)-dxdxi(1,2)*dxdxi(3,1), ...
          + dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)];
      J = norm(nn);
      nn = nn/J;
      e1 = dxdxi(:,1)/norm(dxdxi(:,1));
      e2 = cross(nn,e1);
      for a = 1:nnode/2
         B(1,3*a-2:3*a) = -N(a)*nn(:);
         B(2,3*a-2:3*a) = -N(a)*e1(:);
         B(3,3*a-2:3*a)= -N(a)*e2(:);
         B(1,3*nnode/2 + 3*a-2:3*nnode/2+3*a) = N(a)*nn(:);
         B(2,3*nnode/2 + 3*a-2:3*nnode/2+3*a) = N(a)*e1(:);
         B(3,3*nnode/2 + 3*a-2:3*nnode/2+3*a) = N(a)*e2(:);
     end
    end

    Delta0 = B*displacement(:);
    dDelta = B*displacement_increment(:);
    Delta1 = Delta0+dDelta;
    [Tract,D,Damage1(intpt)] = interface_law(dt,materialprops,Delta1,Delta0,Damage0(intpt));
    
    r = r + transpose(B)*Tract*w(intpt)*J;
    kel = kel + transpose(B)*D*B*w(intpt)*J;
  end
end

function [T,Dmat,D1] = interface_law(dt,materialprops,Delta1,Delta0,D0)

    % dt - time increment
    % Delta1 - interface displacements at end of step
    % Delta0 - interface displacements at start of step
    % D0 damage at start of step
    % T - interface traction
    % Dmat - tangent stiffness

    k0 = materialprops(3);
    k1 = materialprops(4);
    smax = materialprops(5);  
    eta = materialprops(6);    
    
    lam0 = norm(Delta0);
    lam1 = norm(Delta1);
    ndim = length(Delta0);  % No. dimensions (2 or 3)
    
    if (lam1==0 && D0==0)
        T = zeros(ndim,1);
        Dmat = (k0 + eta/dt)*eye(ndim);
        D1 = D0;
        return
    end
   
    
    if (Delta1(1)>0)
        if (D0>1) D0=1; end

        T = k0*(1-D0)*Delta1 + eta*(Delta1-Delta0)/dt; % Elastic predictor
        tau = norm(T);
        lam1 = norm(Delta1);
        if (tau<2*(1-D0)*smax/(2-D0) || lam1<lam0 || D0==1)
            Dmat = (k0*(1-D0) + eta/dt)*eye(ndim);
            D1 = D0;
        else
            D1 = (lam1*D0 + 2*(lam1-lam0))/(2*lam1-lam0);
            dDdlam = lam0*(2-D0)/(2*lam1-lam0)^2;
            if (D1>1)
                D1 = 1;
                dDdlam = 0;
            end
            T = k0*(1-D1)*Delta1 + eta*(Delta1-Delta0)/dt;
            Dmat = (k0*(1-D1) + eta/dt)*eye(3) - k0*dDdlam*tensorprod(Delta1,Delta1)/lam1;
        end
        return
    end

    
    if (D0>1) D0=1; end
    T = k0*(1-D0)*Delta1 + eta*(Delta1-Delta0)/dt; % Elastic predictor
    T(1) = k1*Delta1(1);
    tau = norm(T(2:ndim));
    lam1 = norm(Delta1(2:ndim));
    if (tau<2*(1-D0)*smax/(2-D0) || lam1<lam0 || D0==1)
        Dmat = (k0*(1-D0) + eta/dt)*eye(ndim);
        Dmat(1,1) = k1;
        D1 = D0;
    else
        D1 = (lam1*D0 + 2*(lam1-lam0))/(2*lam1-lam0);
        dDdlam = (D0*(lam1-lam0) + 2*lam1)/(2*lam1-lam0)^2;
        if (D1>1)
            D1 = 1;
            dDdlam = 0;
        end
        T = k0*(1-D1)*Delta1 + eta*(Delta1-Delta0)/dt;
        T(1) = k1*Delta1;
        Dmat = (k0*(1-D1) + eta/dt)*eye(ndim);
        Dmat(1,1) = k1;
        if (ndim>2)
            Dmat(2,2) = Dmat(2,2) - k0*dDdlam*Delta1(2)^2/lam1;
            Dmat(3,3) = Dmat(3,3) - k0*dDdlam*Delta1(3)^2/lam1;
            Dmat(2,3) = -k0*dDdlam*Delta1(2)*Delta1(3)/lam1;
            Dmat(3,2) = Dmat(2,3);
        end
    end



end


function M = tensorprod(a,b)
  M = zeros(length(a),length(b));
  for i = 1:length(a)
      M(i,:) = a(i)*b(:);
  end

end

function r = eldload(ncoord,ndof,nfacenodes,elident,coords,traction)

  npoints = numberofintegrationpoints(ncoord-1,nfacenodes);
  xi = zeros(ncoord-1,1);
  dxdxi = zeros(ncoord,ncoord-1);
  r = zeros(ndof*nfacenodes,1);
   
  xilist = integrationpoints(ncoord-1,nfacenodes,npoints);
  w = integrationweights(ncoord-1,nfacenodes,npoints);

  for intpt = 1:npoints

    for i = 1:ncoord-1
      xi(i) = xilist(i,intpt);
    end

    N = shapefunctions(nfacenodes,ncoord-1,elident,xi);
    dNdxi = shapefunctionderivs(nfacenodes,ncoord-1,elident,xi);
%
%     Compute the jacobian matrix && its determinant
%
    for i = 1:ncoord
      for j = 1:ncoord-1
        dxdxi(i,j) = 0.;
        for a = 1:nfacenodes
          dxdxi(i,j) = dxdxi(i,j) + coords(i,a)*dNdxi(a,j);
        end
      end
    end
    if (ncoord == 2) 
      dt = sqrt(dxdxi(1,1)^2+dxdxi(2,1)^2);
    elseif (ncoord == 3) 
      dt = sqrt( ((dxdxi(2,1)*dxdxi(3,2))-(dxdxi(2,2)*dxdxi(3,1)))^2 ...
          + ((dxdxi(1,1)*dxdxi(3,2))-(dxdxi(1,2)*dxdxi(3,1)))^2 ...
          + ((dxdxi(1,1)*dxdxi(2,2))-(dxdxi(1,2)*dxdxi(2,1)))^2 );
    end
   for a = 1:nfacenodes
      for i = 1:ndof
        row = ndof*(a-1)+i;
        r(row) = r(row) + N(a)*traction(i)*w(intpt)*dt;
      end
    end
  end
end


function [resid,Stif,Damage1] = globalmatrices(dt,ncoord,ndof,nnode,coords,n_solid_elem,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs,ddofs,Damage0)
%
%   Assemble the global stiffness matrix
%
   resid = zeros(ndof*nnode,1);
   Stif = zeros(ndof*nnode,ndof*nnode);
   lmncoord = zeros(ncoord,maxnodes);
   lmndof = zeros(ndof,maxnodes);
   dlmndof = lmndof;
   Damage1 = Damage0;
%
%   Loop over all the elements
%
   for lmn = 1:nelem
%
%   Extract coords of nodes, DOF for the current element
%
      for a = 1:nelnodes(lmn)
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
          lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
          dlmndof(i,a) = ddofs(ndof*(connect(a,lmn)-1)+i);
        end
      end
    n = nelnodes(lmn);
    ident = elident(lmn);
    if (lmn<n_solid_elem+1)
       [rel,kel] = elforcnstif(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof+dlmndof);
    else
       [rel,kel,Damage1(:,lmn-n_solid_elem)] = elforcnstif_cz(dt,ncoord,ndof,n,ident,lmncoord,materialprops,lmndof,dlmndof,Damage0(:,lmn-n_solid_elem));
    end
            %
%   Add the current element stiffness:the global stiffness
%
    for a = 1:nelnodes(lmn)
      for i = 1:ndof
        rw = ndof*(connect(a,lmn)-1)+i;
        resid(rw) = resid(rw) + rel(ndof*(a-1)+i);  
        for b = 1:nelnodes(lmn)
          for k = 1:ndof
            cl = ndof*(connect(b,lmn)-1)+k;
            Stif(rw,cl) = Stif(rw,cl) + kel(ndof*(a-1)+i,ndof*(b-1)+k);
          end
        end
      end
    end
   end
end

function r = globaltraction(ncoord,ndof,nnodes,ndload,coords,nelnodes,elident,connect,dloads,dofs)

   r = zeros(ndof*nnodes,1);
   traction = zeros(ndof,1);

   for load = 1:ndload
%
%     Extract the coords of the nodes on the appropriate element face
%
      lmn = dloads(1,load);
      face = dloads(2,load);
      n = nelnodes(lmn);
      ident = elident(lmn);
      nfnodes = nfacenodes(ncoord,n,ident,face); 
      nodelist = facenodes(ncoord,n,ident,face);     
      lmncoord = zeros(ncoord,nfnodes);
      for a = 1:nfnodes
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(nodelist(a),dloads(1,load)));
        end
        for i = 1:ndof
          lmndof(i,a) = dofs(ndof*(connect(nodelist(a),dloads(1,load))-1)+i);
        end
      end
%
%    Compute the element load vector
%
     for i = 1:ndof
       traction(i) = dloads(i+2,load);
     end

     rel = eldload(ncoord,ndof,nfnodes,ident,lmncoord,traction);
%
%    Assemble the element load vector into global vector
%
     for a = 1:nfnodes
       for i = 1:ndof
         rw = (connect(nodelist(a),dloads(1,load))-1)*ndof+i;
         r(rw) = r(rw) + rel((a-1)*ndof+i);
       end
     end

   end
end    


function n = numberofintegrationpoints(ncoord,nelnodes,elident)
  
   if (ncoord == 1) 
     n = nelnodes;   
   elseif (ncoord == 2) 
     if (nelnodes == 3)
         n = 1;
     end
     if (nelnodes == 6)
         n = 3;
     end
     if (nelnodes == 4)
         n = 4;
     end
     if (nelnodes == 8)
         n = 9;
     end
   elseif (ncoord == 3) 
     if (nelnodes == 4)
         n = 1 ;
     end
     if (nelnodes == 10)
         n = 4;
     end
     if (nelnodes == 8)
         n = 8;
     end
     if (nelnodes == 20)
         n = 27;
     end
   end
end   

function xi = integrationpoints(ncoord,nelnodes,npoints,elident)

   xi = zeros(ncoord,npoints);
%
%  1D elements
%
   if (ncoord == 1) 
     if (npoints==1) 
       xi(1,1) = 0.;
     elseif (npoints == 2) 
       xi(1,1) = -0.5773502692;
       xi(1,2) = -xi(1,1);
     elseif (npoints == 3) 
       xi(1,1) = -0.7745966692;
       xi(1,2) = 0.0;
       xi(1,3) = -xi(1,1);
     end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 || nelnodes == 6 ) 
       if (npoints == 1) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
       elseif (npoints == 3) 
         xi(1,1) = 0.6;
         xi(2,1) = 0.2;
         xi(1,2) = 0.2;
         xi(2,2) = 0.6;
         xi(1,3) = 0.2;
         xi(2,3) = 0.2;
       elseif (npoints == 4) 
         xi(1,1) = 1./3.;
         xi(2,1) = 1./3.;
         xi(1,2) = 0.6;
         xi(2,2) = 0.2;
         xi(1,3) = 0.2;
         xi(2,3) = 0.6;
         xi(1,4) = 0.2;
         xi(2,4) = 0.2;
       end
%
%    Rectangular element
%                  
     elseif ( nelnodes==4 || nelnodes==8 ) 

       if (npoints == 1) 
         xi(1,1) = 0.;
         xi(2,1) = 0.;
       elseif (npoints == 4) 
         xi(1,1) = -0.5773502692;
         xi(2,1) = xi(1,1);
         xi(1,2) = -xi(1,1);
         xi(2,2) = xi(1,1);
         xi(1,3) = xi(1,1);
         xi(2,3) = -xi(1,1);
         xi(1,4) = -xi(1,1);
         xi(2,4) = -xi(1,1);
       elseif (npoints == 9) 
         xi(1,1) = -0.7745966692;
         xi(2,1) = xi(1,1);
         xi(1,2) = 0.0;
         xi(2,2) = xi(1,1);
         xi(1,3) = -xi(1,1);
         xi(2,3) = xi(1,1);
         xi(1,4) = xi(1,1);
         xi(2,4) = 0.0;
         xi(1,5) = 0.0;
         xi(2,5) = 0.0;
         xi(1,6) = -xi(1,1);
         xi(2,6) = 0.0;
         xi(1,7) = xi(1,1);
         xi(2,7) = -xi(1,1);
         xi(1,8) = 0.;
         xi(2,8) = -xi(1,1);
         xi(1,9) = -xi(1,1);
         xi(2,9) = -xi(1,1);
       end
     end
%
%   3D elements
%
   elseif (ncoord == 3) 
%
%  3D elements
%
     if (nelnodes == 4 || nelnodes==10 ) 
       if (npoints == 1) 
         xi(1,1) = 0.25;
         xi(2,1) = 0.25;
         xi(3,1) = 0.25;
       elseif (npoints == 4) 
         xi(1,1) = 0.58541020;
         xi(2,1) = 0.13819660;
         xi(3,1) = xi(2,1);
         xi(1,2) = xi(2,1);
         xi(2,2) = xi(1,1);
         xi(3,2) = xi(2,1);
         xi(1,3) = xi(2,1);
         xi(2,3) = xi(2,1);
         xi(3,3) = xi(1,1);
         xi(1,4) = xi(2,1);
         xi(2,4) = xi(2,1);
         xi(3,4) = xi(2,1);
       end
     elseif ( nelnodes==8 || nelnodes==20 ) 
       if (npoints == 1) 
         xi(1,1) = 0.;
         xi(2,1) = 0.;
         xi(3,1) = 0.;
       elseif (npoints == 8) 
         x1D = [-0.5773502692,0.5773502692];
         for k = 1:2
           for j = 1:2 
             for i = 1:2
               n = 4*(k-1) + 2*(j-1) + i;
               xi(1,n) = x1D(i);
               xi(2,n) = x1D(j);
               xi(3,n) = x1D(k);
             end
           end
         end
       elseif (npoints == 27) 
         x1D = [-0.7745966692,0.,0.7745966692];
         for k = 1:3
           for j = 1:3
             for i = 1:3
               n = 9*(k-1) + 3*(j-1) + i;
               xi(1,n) = x1D(i);
               xi(2,n) = x1D(j);
               xi(3,n) = x1D(k);
             end
           end
         end
       end
     end
   end
end

function w = integrationweights(ncoord,nelnodes,npoints,elident)

   w = zeros(npoints,1);

%
%  1D elements
%
   if (ncoord == 1) 
     if (npoints == 1)
       w(1) = 2.;
     elseif (npoints == 2) 
       w = [1.,1.];
     elseif (npoints == 3) 
       w = [0.555555555,0.888888888,0.555555555];
     end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 || nelnodes == 6 ) 
       if (npoints == 1) 
         w(1) = 0.5;
       elseif (npoints == 3) 
         w(1) = 1./6.;
         w(2) = 1./6.;
         w(3) = 1./6.;
       elseif (npoints == 4) 
         w = [-27./96.,25./96.,25/96.,25/96.];
       end
%
%    Rectangular element
%                  
     elseif ( nelnodes==4 || nelnodes==8 ) 

       if (npoints == 1) 
         w(1) = 4.;
       elseif (npoints == 4) 
         w = [1.,1.,1.,1.];
       elseif (npoints == 9 ) 
         w1D = [0.555555555,0.888888888,0.55555555555];
         for j = 1:3
           for i = 1:3
             n = 3*(j-1)+i;
             w(n) = w1D(i)*w1D(j);
           end
         end    
       end
     end 

   elseif (ncoord == 3) 
%
%  3D elements
%
     if (nelnodes == 4 || nelnodes==10 ) 
       if (npoints == 1) 
         w(1) = 1./6.;
       elseif (npoints == 4) 
         w = [1./24.,1./24.,1./24.,1./24.];
       end
     elseif ( nelnodes==8 || nelnodes==20 ) 
       if (npoints == 1) 
         w(1) = 8.;
       elseif (npoints == 8) 
         w = [1.,1.,1.,1.,1.,1.,1.,1.];
       elseif (npoints == 27) 
         w1D = [0.555555555,0.888888888,0.55555555555];
         for k = 1:3
           for j = 1:3
             for i = 1:3
               n = 9*(k-1)+3*(j-1)+i;
               w(n) = w1D(i)*w1D(j)*w1D(k);
             end
           end    
         end
       end
     end
   end
end

function N = shapefunctions(nelnodes,ncoord,elident,xi)
 

   N = zeros(nelnodes,1);
%
%  1D elements
%
  if (ncoord == 1) 
    if (nelnodes==2) 
      N(1) = 0.5*(1.+xi(1));
      N(2) = 0.5*(1.-xi(1));
    elseif (nelnodes == 3) 
      N(1) = -0.5*xi(1)*(1.-xi(1));
      N(2) =  0.5*xi(1)*(1.+xi(1));
      N(3) = (1.-xi(1))*(1.+xi(1));
    end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 ) 
       N(1) = xi(1);
       N(2) = xi(2);
       N(3) = 1.-xi(1)-xi(2);               
     elseif ( nelnodes == 6 ) 
       xi3 = 1.-xi(1)-xi(2);
       N(1) = (2.*xi(1)-1.)*xi(1);
       N(2) = (2.*xi(2)-1.)*xi(2);
       N(3) = (2.*xi3-1.)*xi3;
       N(4) = 4.*xi(1)*xi(2);
       N(5) = 4.*xi(2)*xi3;
       N(6) = 4.*xi3*xi(1);
%
%    Rectangular element
%                  
     elseif ( nelnodes == 4 ) 
       N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
       N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
       N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
       N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
     elseif (nelnodes == 8) 
       N(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
       N(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
       N(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
       N(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
       N(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
       N(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
       N(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
       N(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
     end
%
   elseif (ncoord==3) 

     if (nelnodes == 4) 
       N(1) = xi(1);
       N(2) = xi(2);
       N(3) = xi(3);
       N(4) = 1.-xi(1)-xi(2)-xi(3);
     elseif (nelnodes == 10) 
       xi4 = 1.-xi(1)-xi(2)-xi(3);
       N(1) = (2.*xi(1)-1.)*xi(1);
       N(2) = (2.*xi(2)-1.)*xi(2);
       N(3) = (2.*xi(3)-1.)*xi(3);
       N(4) = (2.*xi4-1.)*xi4;
       N(5) = 4.*xi(1)*xi(2);
       N(6) = 4.*xi(2)*xi(3);
       N(7) = 4.*xi(3)*xi(1);
       N(8) = 4.*xi(1)*xi4;
       N(9) = 4.*xi(2)*xi4;
       N(10) = 4.*xi(3)*xi4;
     elseif (nelnodes == 8) 
       N(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.;
       N(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.;
       N(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.;
       N(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.;
       N(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.;
       N(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.;
       N(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.;
       N(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.;
     elseif (nelnodes == 20) 
       N(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.;
       N(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.;
       N(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.;
       N(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.;
       N(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.;
       N(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.;
       N(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.;
       N(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.;
       N(9)  = (1.-xi(1)^2)*(1.-xi(2))*(1.-xi(3))/4.;
       N(10) = (1.+xi(1))*(1.-xi(2)^2)*(1.-xi(3))/4.;
       N(11) = (1.-xi(1)^2)*(1.+xi(2))*(1.-xi(3))/4.;
       N(12) = (1.-xi(1))*(1.-xi(2)^2)*(1.-xi(3))/4.;
       N(13) = (1.-xi(1)^2)*(1.-xi(2))*(1.+xi(3))/4.;
       N(14) = (1.+xi(1))*(1.-xi(2)^2)*(1.+xi(3))/4.;
       N(15) = (1.-xi(1)^2)*(1.+xi(2))*(1.+xi(3))/4.;
       N(16) = (1.-xi(1))*(1.-xi(2)^2)*(1.+xi(3))/4.;
       N(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)^2)/4.;
       N(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)^2)/4.;
       N(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)^2)/4.;
       N(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)^2)/4.;
     end
   end

end

function dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi)

  dNdxi = zeros(nelnodes,ncoord);
%
% 1D elements
%
  if (ncoord == 1) 
    if (nelnodes==2) 
      dNdxi(1,1) = 0.5;
      dNdxi(2,1) = -0.5;
    elseif (nelnodes == 3) 
      dNdxi(1,1) = -0.5+xi(1);
      dNdxi(2,1) =  0.5+xi(1);
      dNdxi(3,1) = -2.*xi(1);
    end
%
%  2D elements
%
   elseif (ncoord == 2) 
%
%    Triangular element
%
     if ( nelnodes == 3 ) 
       dNdxi(1,1) = 1.;
       dNdxi(2,2) = 1.;
       dNdxi(3,1) = -1.;
       dNdxi(3,2) = -1.;               
     elseif ( nelnodes == 6 ) 
       z = 1.D0 - xi(1) - xi(2);
       dzdp = -1.D0;
       dzdq = -1.D0;
       dNdxi(1, 1) = 4.D0*xi(1) - 1.D0;
       dNdxi(3, 1) = 4.D0*z*dzdp - dzdp;
       dNdxi(4, 1) = 4.D0*xi(2);
       dNdxi(5, 1) = 4.D0*xi(2)*dzdp;
       dNdxi(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp;
       dNdxi(1, 2) = 0.D0;
       dNdxi(2, 2) = 4.D0*xi(2) - 1.D0;
       dNdxi(3, 2) = 4.D0*z*dzdq - dzdq;
       dNdxi(4, 2) = 4.D0*xi(1);
       dNdxi(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq;
       dNdxi(6, 2) = 4.D0*xi(1)*dzdq;
%
%    Rectangular element
%                  
     elseif ( nelnodes == 4 ) 
       dNdxi(1,1) = -0.25*(1.-xi(2));
       dNdxi(1,2) = -0.25*(1.-xi(1));
       dNdxi(2,1) = 0.25*(1.-xi(2));
       dNdxi(2,2) = -0.25*(1.+xi(1));
       dNdxi(3,1) = 0.25*(1.+xi(2));
       dNdxi(3,2) = 0.25*(1.+xi(1));
       dNdxi(4,1) = -0.25*(1.+xi(2));
       dNdxi(4,2) = 0.25*(1.-xi(1));
     elseif (nelnodes == 8) 
       dNdxi(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
       dNdxi(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
       dNdxi(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
       dNdxi(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
       dNdxi(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
       dNdxi(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
       dNdxi(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
       dNdxi(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
       dNdxi(5,1) = -xi(1)*(1.-xi(2));
       dNdxi(5,2) = -0.5*(1.-xi(1)*xi(1));
       dNdxi(6,1) = 0.5*(1.-xi(2)*xi(2));
       dNdxi(6,2) = -(1.+xi(1))*xi(2);
       dNdxi(7,1) = -xi(1)*(1.+xi(2));
       dNdxi(7,2) = 0.5*(1.-xi(1)*xi(1));
       dNdxi(8,1) = -0.5*(1.-xi(2)*xi(2));
       dNdxi(8,2) = -(1.-xi(1))*xi(2);
      end
%
%    3D elements
%
   elseif (ncoord==3) 

     if (nelnodes == 4) 
       dNdxi(1,1) = 1.;
       dNdxi(2,2) = 1.;
       dNdxi(3,3) = 1.;
       dNdxi(4,1) = -1.;
       dNdxi(4,2) = -1.;
       dNdxi(4,3) = -1.;
     elseif (nelnodes == 10) 
       xi4 = 1.-xi(1)-xi(2)-xi(3);
       dNdxi(1,1) = (4.*xi(1)-1.);
       dNdxi(2,2) = (4.*xi(2)-1.);
       dNdxi(3,3) = (4.*xi(3)-1.);
       dNdxi(4,1) = -(4.*xi4-1.);
       dNdxi(4,2) = -(4.*xi4-1.);
       dNdxi(4,3) = -(4.*xi4-1.);
       dNdxi(5,1) = 4.*xi(2);
       dNdxi(5,2) = 4.*xi(1);
       dNdxi(6,2) = 4.*xi(3);
       dNdxi(6,3) = 4.*xi(2);
       dNdxi(7,1) = 4.*xi(3);
       dNdxi(7,3) = 4.*xi(1); 
       dNdxi(8,1) = 4.*(xi4-xi(1));
       dNdxi(8,2) = -4.*xi(1);
       dNdxi(8,3) = -4.*xi(1);
       dNdxi(9,1) = -4.*xi(2);
       dNdxi(9,2) = 4.*(xi4-xi(2));
       dNdxi(9,3) = -4.*xi(2);
       dNdxi(10,1) = -4.*xi(3)*xi4;
       dNdxi(10,2) = -4.*xi(3);
       dNdxi(10,3) = 4.*(xi4-xi(3));
     elseif (nelnodes == 8) 
       dNdxi(1,1) = -(1.-xi(2))*(1.-xi(3))/8.;
       dNdxi(1,2) = -(1.-xi(1))*(1.-xi(3))/8.;
       dNdxi(1,3) = -(1.-xi(1))*(1.-xi(2))/8.;
       dNdxi(2,1) = (1.-xi(2))*(1.-xi(3))/8.;
       dNdxi(2,2) = -(1.+xi(1))*(1.-xi(3))/8.;
       dNdxi(2,3) = -(1.+xi(1))*(1.-xi(2))/8.;
       dNdxi(3,1) = (1.+xi(2))*(1.-xi(3))/8.;
       dNdxi(3,2) = (1.+xi(1))*(1.-xi(3))/8.;
       dNdxi(3,3) = -(1.+xi(1))*(1.+xi(2))/8.;
       dNdxi(4,1) = -(1.+xi(2))*(1.-xi(3))/8.;
       dNdxi(4,2) = (1.-xi(1))*(1.-xi(3))/8.;
       dNdxi(4,3) = -(1.-xi(1))*(1.+xi(2))/8.;
       dNdxi(5,1) = -(1.-xi(2))*(1.+xi(3))/8.;
       dNdxi(5,2) = -(1.-xi(1))*(1.+xi(3))/8.;
       dNdxi(5,3) = (1.-xi(1))*(1.-xi(2))/8.;
       dNdxi(6,1) = (1.-xi(2))*(1.+xi(3))/8.;
       dNdxi(6,2) = -(1.+xi(1))*(1.+xi(3))/8.;
       dNdxi(6,3) = (1.+xi(1))*(1.-xi(2))/8.;
       dNdxi(7,1) = (1.+xi(2))*(1.+xi(3))/8.;
       dNdxi(7,2) = (1.+xi(1))*(1.+xi(3))/8.;
       dNdxi(7,3) = (1.+xi(1))*(1.+xi(2))/8.;
       dNdxi(8,1) = -(1.+xi(2))*(1.+xi(3))/8.;
       dNdxi(8,2) = (1.-xi(1))*(1.+xi(3))/8.;
       dNdxi(8,3) = (1.-xi(1))*(1.+xi(2))/8.;
     elseif (nelnodes == 20) 
       dNdxi(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;

       dNdxi(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;
       dNdxi(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.;

       dNdxi(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)-(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;

       dNdxi(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.;
       dNdxi(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)-(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)+(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)-(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)+(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.;
       dNdxi(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)+(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)-(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)+(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.;
       dNdxi(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.;
       dNdxi(9,2)  = -(1.-xi(1)^2)*(1.-xi(3))/4.;
       dNdxi(9,3)  = -(1.-xi(1)^2)*(1.-xi(2))/4.;
       dNdxi(10,1)  = (1.-xi(2)^2)*(1.-xi(3))/4.;
       dNdxi(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.;
       dNdxi(10,3)  = -(1.-xi(2)^2)*(1.+xi(1))/4.;
       dNdxi(11,1)  = -2.*xi(1)*(1.+xi(2))*(1.-xi(3))/4.;
       dNdxi(11,2)  = (1.-xi(1)^2)*(1.-xi(3))/4.;
       dNdxi(11,3)  = -(1.-xi(1)^2)*(1.+xi(2))/4.;
       dNdxi(12,1)  = -(1.-xi(2)^2)*(1.-xi(3))/4.;
       dNdxi(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.;
       dNdxi(12,3)  = -(1.-xi(2)^2)*(1.-xi(1))/4.;
       dNdxi(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.;
       dNdxi(13,2)  = -(1.-xi(1)^2)*(1.+xi(3))/4.;
       dNdxi(13,3)  = (1.-xi(1)^2)*(1.-xi(2))/4.;
       dNdxi(14,1)  = (1.-xi(2)^2)*(1.+xi(3))/4.;
       dNdxi(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.;
       dNdxi(14,3)  = (1.-xi(2)^2)*(1.+xi(1))/4.;
       dNdxi(15,1)  = -2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.;
       dNdxi(15,2)  =  (1.-xi(1)^2)*(1.+xi(3))/4.;
       dNdxi(15,3)  = (1.-xi(1)^2)*(1.+xi(2))/4.;
       dNdxi(16,1)  = -(1.-xi(2)^2)*(1.+xi(3))/4.;
       dNdxi(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.;
       dNdxi(16,3)  = (1.-xi(2)^2)*(1.-xi(1))/4.;
       dNdxi(17,1) = -(1.-xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(17,2) = -(1.-xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.;
       dNdxi(18,1) = (1.-xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(18,2) = -(1.+xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.;
       dNdxi(19,1) = (1.+xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(19,2) = (1.+xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.;
       dNdxi(20,1) = -(1.+xi(2))*(1.-xi(3)^2)/4.;
       dNdxi(20,2) = (1.-xi(1))*(1.-xi(3)^2)/4.;
       dNdxi(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.;
     end
   end

end

function n = nfacenodes(ncoord,nelnodes,elident,face)
   if (ncoord == 2) 
     if (nelnodes == 3 || nelnodes == 4)
         n = 2;
     elseif (nelnodes == 6 || nelnodes == 8)
         n=3;
     end
   elseif (ncoord == 3) 
     if (nelnodes == 4)
         n = 3;
     elseif (nelnodes == 10)
         n = 6;
     elseif (nelnodes == 8)
         n = 4;
     elseif (nelnodes == 20)
         n = 8;
     end
   end
end

function list = facenodes(ncoord,nelnodes,elident,face)

   i3 = [2,3,1];
   i4 = [2,3,4,1]; 

   list = zeros(nfacenodes(ncoord,nelnodes,elident,face),1);

   if (ncoord == 2) 
     if (nelnodes == 3) 
       list(1) = face;
       list(2) = i3(face);
     elseif (nelnodes == 6) 
       list(1) = face;
       list(2) = i3(face);
       list(3) = face+3;
     elseif (nelnodes==4) 
       list(1) = face;
       list(2) = i4(face);
     elseif (nelnodes==8) 
       list(1) = face;
       list(2) = i4(face);
       list(3) = face+4;
     end
   elseif (ncoord == 3) 
     if (nelnodes==4) 
       if   (face == 1)
           list = [1,2,3];
       elseif (face == 2)
           list = [1,4,2];
       elseif (face == 3)
           list = [2,4,3];
       elseif (face == 4)
           list = [3,4,1];
       end
     elseif (nelnodes == 10) 
       if   (face == 1)
           list = [1,2,3,5,6,7];
       elseif (face == 2)
           list = [1,4,2,8,9,5];
       elseif (face == 3)
           list = [2,4,3,9,10,6];
       elseif (face == 4)
           list = [3,4,1,10,8,7];
       end
     elseif (nelnodes == 8) 
       if   (face == 1)
           list = [1,2,3,4];
       elseif (face == 2)
           list = [5,8,7,6];
       elseif (face == 3)
           list = [1,5,6,2];
       elseif (face == 4)
           list = [2,3,7,6];
       elseif (face == 5)
           list = [3,7,8,4];
       elseif (face == 6)
           list = [4,8,5,1];
       end
     elseif (nelnodes == 20)  
       if   (face == 1)
           list = [1,2,3,4,9,10,11,12];
       elseif (face == 2)
           list = [5,8,7,6,16,15,14,13];
       elseif (face == 3)
           list = [1,5,6,2,17,13,18,9];
       elseif (face == 4)
           list = [2,6,7,3,18,14,19,10];
       elseif (face == 5)
           list = [3,7,8,4,19,15,20,11];
       elseif (face == 6)
           list = [4,8,5,1,20,16,17,12];
       end
     end
   end
end


function print_results(outfile, ...
    nprops,materialprops,ncoord,ndof,nnode,coords, ...
    nelem,maxnodes,connect,nelnodes,elident, ...
    nfix,fixnodes,ndload,dloads,dofs)
%      Print nodal displacements, element strains && stresses:a file




    fprintf(outfile,'Nodal Displacements: \n');
    if (ndof == 2)
        fprintf(outfile,' Node      Coords         u1       u2 \n');
        for i = 1:nnode
            fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f\n', ...
                i,coords(1,i),coords(2,i),dofs(2*i-1),dofs(2*i));
        end
    elseif (ndof == 3)
        fprintf(outfile,' Node            Coords            u1       u2       u3 \n');
        for i = 1:nnode
            fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', ...
                i,coords(1,i),coords(2,i),coords(3,i),dofs(3*i-2),dofs(3*i-1),dofs(3*i));
        end
    end

    fprintf(outfile,'\n\n Strains and Stresses \n');


    lmncoord = zeros(ncoord,maxnodes);
    displacement = zeros(ndof,maxnodes);

    E = materialprops(1);
    nu = materialprops(2);
    
    %
    %   Loop over all the elements
    %
    for lmn = 1:nelem

        fprintf(outfile,' \n Element; %d ',lmn);
        if (ncoord == 2)
            fprintf(outfile,'  \n int pt    Coords          e_11      e_22     e_12      s_11       s_22      s_12 \n');

        elseif (ncoord == 3)
            fprintf(outfile,'\n int pt         Coords            e_11      e_22     e_33      e_12       e_13      e_23      s_11      s_22      s_33      s_12      s_13      s_23 \n');
        end
        %
        %   Extract coords of nodes, DOF for the current element
        %
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                lmncoord(i,a) = coords(i,connect(a,lmn));
            end
            for i = 1:ndof
                displacement(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
            end
        end
        n = nelnodes(lmn);

        npoints = numberofintegrationpoints(ncoord,n);

        xi = zeros(ncoord,1);
        xilist = integrationpoints(ncoord,n,npoints);

        uvec = reshape(displacement,1,ndof*nelnodes(lmn));

        for intpt = 1 : npoints

            xi(1:ncoord) = xilist(1:ncoord,intpt);

            N = shapefunctions(n,ncoord,elident(lmn),xi); 
            dNdxi = shapefunctionderivs(nelnodes(lmn),ncoord,elident(lmn),xi);
            x = lmncoord*N;
            dxdxi = lmncoord*dNdxi;
            dxidx = inv(dxdxi);

            dNdx = dNdxi*dxidx;


           if (ndof==2)
                B = zeros(3,2*nelnodes(lmn));
                B(1,1:2:2*nelnodes(lmn)-1) = dNdx(1:nelnodes(lmn),1);
                B(2,2:2:2*nelnodes(lmn)) = dNdx(1:nelnodes(lmn),2);
                B(3,1:2:2*nelnodes(lmn)-1) = dNdx(1:nelnodes(lmn),2);
                B(3,2:2:2*nelnodes(lmn)) = dNdx(1:nelnodes(lmn),1);
                D = E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];  % Default is plane stress
                if (elident(lmn)==1)
                   D = E/(1+nu)/(1-2*nu)*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
                end
            else
                B = zeros(6,3*nelnodes(lmn));
                B(1,1:3:3*nelnodes(lmn)-2) = dNdx(1:nelnodes(lmn),1);
                B(2,2:3:3*nelnodes(lmn)-1) = dNdx(1:nelnodes(lmn),2);
                B(3,3:3:3*nelnodes(lmn))   = dNdx(1:nelnodes(lmn),3);
                B(4,1:3:3*nelnodes(lmn)-2) = dNdx(1:nelnodes(lmn),2);
                B(4,2:3:3*nelnodes(lmn)-1) = dNdx(1:nelnodes(lmn),1);
                B(5,1:3:3*nelnodes(lmn)-2) = dNdx(1:nelnodes(lmn),3);
                B(5,3:3:3*nelnodes(lmn))   = dNdx(1:nelnodes(lmn),1);
                B(6,2:3:3*nelnodes(lmn)-1) = dNdx(1:nelnodes(lmn),3);
                B(6,3:3:3*nelnodes(lmn))   = dNdx(1:nelnodes(lmn),2);
                D = eye(6)*(1-2*nu)/2;
                D(1:3,1:3) = nu;
                D([1,8,15]) = (1-nu);
                D = E/(1+nu)/(1-2*nu)*D;
            end
            e = B*lmndof;
            s = D*e;
            if (ncoord == 2)
                fprintf(outfile,'%5d %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n', ...
                    intpt,x(1),x(2),e(1),e(2),e(3)/2,s(1),s(2),s(3));
            elseif (ncoord == 3)
                fprintf(outfile,'%5d %7.4f %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n',...
                    intpt,x(1),x(2),x(3), ...
                    e(1:3),e(4:6)/2,s(1:6));
            end

        end
    end
end



%==============================Function read_input_file =================
 function [nprops,materialprops,ncoord,ndof,nnode,coords,n_solid_elem,nelem,maxnodes, ...
                connect,nelnodes,elident,nfix,fixnodes,ndload,dloads] = read_input_file(infile) 

 cellarray=textscan(infile,'%s');

%  Total no. material parameters, && list of parameters
%
   nprops = str2num(cellarray{1}{2});
   materialprops = zeros(nprops,1);
   cellno = 2;
   for i = 1:nprops
     cellno = cellno + 2;  
     materialprops(i) = str2num(cellarray{1}{cellno});
   end
%
%    no. coords (1:3), no. DOF, no. nodes && nodal coordinates
%
   cellno = cellno + 2;
   ncoord = str2num(cellarray{1}{cellno});
   cellno = cellno + 2;
   ndof = str2num(cellarray{1}{cellno});
   cellno = cellno + 2;
   nnode = str2num(cellarray{1}{cellno});
 
   coords = zeros(ncoord,nnode);
   cellno = cellno + 1;
   for i = 1 : nnode
     for j = 1 : ncoord
       cellno = cellno + 1;
       coords(j,i) = str2num(cellarray{1}{cellno});
     end
   end
%
%    No. solid elements && connectivity
%
   cellno = cellno + 2;
   nelem = str2num(cellarray{1}{cellno});
   cellno = cellno + 2;
   maxnodes = str2num(cellarray{1}{cellno});
   connect = zeros(maxnodes,nelem);
   nelnodes = zeros(nelem,1);
   elident = zeros(nelem,1);
   cellno = cellno + 3;
   for i = 1 : nelem
     cellno = cellno + 1;
     elident(i) = str2num(cellarray{1}{cellno});
     cellno = cellno + 1;
     nelnodes(i) = str2num(cellarray{1}{cellno});
     for j = 1 : nelnodes(i)
       cellno = cellno + 1;
       connect(j,i) = str2num(cellarray{1}{cellno});
     end
   end
%
%    No. CZ elements && connectivity
%
   cellno = cellno + 2;
   n_solid_elem = nelem;
   nelem = nelem + str2num(cellarray{1}{cellno});
   cellno = cellno + 3;
   for i = n_solid_elem + 1 : nelem
     cellno = cellno + 1;
     elident(i) = str2num(cellarray{1}{cellno});
     cellno = cellno + 1;
     nelnodes(i) = str2num(cellarray{1}{cellno});
     for j = 1 : nelnodes(i)
       cellno = cellno + 1;
       connect(j,i) = str2num(cellarray{1}{cellno});
     end
   end   
   
%
%    No. nodes with prescribed displacements, with the prescribed displacements
% 
   cellno = cellno + 2;
   nfix = str2num(cellarray{1}{cellno});
   cellno = cellno + 3;
   fixnodes = zeros(3,nfix);
   for i = 1 : nfix
      cellno = cellno + 1;
      fixnodes(1,i) = str2num(cellarray{1}{cellno});
      cellno = cellno + 1;
      fixnodes(2,i) = str2num(cellarray{1}{cellno});
      cellno = cellno + 1;
      fixnodes(3,i) = str2num(cellarray{1}{cellno});
   end
%
%    No. loaded element faces, with the loads
%
    cellno = cellno + 2;
    ndload = str2num(cellarray{1}{cellno});
    cellno = cellno + 3;
    dloads = zeros(2+ndof,ndload);
    for i = 1 : ndload
       cellno = cellno + 1;
       dloads(1,i) = str2num(cellarray{1}{cellno});
       cellno = cellno + 1;
       dloads(2,i) = str2num(cellarray{1}{cellno});
       for j = 1 : ndof
         cellno = cellno + 1;
         dloads(j+2,i) = str2num(cellarray{1}{cellno});
       end
    end
    
 end
 
 function plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,color)
% Function to plot a mesh.  
    f2D_3 = [1,2,3];
    f2D_4 = [1,2,3,4];
    f2D_6 = [1,4,2,5,3,6];
    f2D_8 = [1,5,2,6,3,7,4,8];
    f3D_4 = [[1,2,3];[1,4,2];[2,4,3];[3,4,1]];
    f3D_10 = [[1,5,2,6,3,7];[1,8,4,9,2,5];[2,9,4,10,3,6];[3,10,4,8,1,7]];
    f3D_8 = [[1,2,3,4];[5,8,7,6];[1,5,6,2];[2,3,7,6];[3,7,8,4];[4,8,5,1]];
    f3D_20 = [[1,9,2,10,3,11,4,12];[5,16,8,15,7,14,6,13];
              [1,17,5,13,6,18,2,9];[2,18,6,14,7,19,3,10];
              [3,19,7,15,8,20,4,11];[4,20,8,16,5,17,1,12]];

   hold on
   if (ncoord==2)  % Plot a 2D mesh
       for lmn = 1:nelem
           if (elident(lmn)==1)
           for i = 1:nelnodes(lmn)
               x(i,1:2) = coords(1:2,connect(i,lmn));
           end
 %          scatter(x(:,1),x(:,2),'MarkerFaceColor','r');
           if (nelnodes(lmn)==3) 
               patch('Vertices',x,'Faces',f2D_3,'FaceColor','none','EdgeColor',color);
           elseif (nelnodes(lmn)==4)
               patch('Vertices',x,'Faces',f2D_4,'FaceColor','none','EdgeColor',color,'LineWidth',2,'Marker','o','MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor',[0,0,0]);
           elseif (nelnodes(lmn)==6) 
               patch('Vertices',x,'Faces',f2D_6,'FaceColor','none','EdgeColor',color);
           elseif (nelnodes(lmn)==8 || nelnodes(lmn)==9)
               patch('Vertices',x,'Faces',f2D_8,'FaceColor','none','EdgeColor',color);
           end
           end
       end
   elseif (ncoord==3) % Plot a 3D mesh
       for lmn = 1:nelem
           if (elident(lmn)==1)
           for i = 1:nelnodes(lmn)
               x(i,1:3) = coords(1:3,connect(i,lmn));
           end
  %         scatter3(x(:,1),x(:,2),x(:,3),'MarkerFaceColor','r');
           if (nelnodes(lmn)==4) 
               patch('Vertices',x,'Faces',f3D_4,'FaceColor','none','EdgeColor',color);
           elseif (nelnodes(lmn)==10)
               patch('Vertices',x,'Faces',f3D_10,'FaceColor','none','EdgeColor',color);
           elseif (nelnodes(lmn)==8) 
               patch('Vertices',x,'Faces',f3D_8,'FaceColor','none','EdgeColor',color,'LineWidth',2,'Marker','o','MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor',[0,0,0]);
           elseif (nelnodes(lmn)==20)
               patch('Vertices',x,'Faces',f3D_20,'FaceColor','none','EdgeColor',color);
           end
           end
       end    
   end
   axis equal
   hold off
  end
function check_stiffness(outfile,lmn,ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs)
%
%   Assemble the global stiffness matrix
%

   lmncoord = zeros(ncoord,maxnodes);
   lmndof = zeros(ndof,maxnodes);
   kel_num = zeros(ndof*maxnodes,ndof*maxnodes);
   %
%
%   Extract coords of nodes, DOF for the current element
%
      for a = 1:nelnodes(lmn)
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
          lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
      end
      lmndof1 = lmndof;
    n = nelnodes(lmn);
    ident = elident(lmn);
    for a = 1:nelnodes(lmn)
        for i = 1:ndof
          lmndof(i,a) = lmndof(i,a) - 1.e-08;  
          [rel0,kel] = elforcnstif(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof);
          lmndof(i,a) = lmndof(i,a) + 2.e-08;
          [rel1,kel1] = elforcnstif(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof);
          kel_num(ndof*(a-1)+i,:) = (rel1(:)-rel0(:))/2.e-08;
        end
    end

    fprintf(outfile,'Coded stiffness \n ');
    for a = 1:nelnodes(lmn)
      for i = 1:ndof
%         fprintf(outfile,'\n\n Node %d DOF %d \n',a,i); 
        rw = ndof*(a-1)+i;
        for b = 1:nelnodes(lmn)
          for k = 1:ndof
            cl = ndof*(b-1)+k;
              fprintf(outfile,'%d ',kel(rw,cl));
          end
        end
        fprintf(outfile,'\n');
      end
    end
    
%
    fprintf(outfile,'\n \n Numerical stiffness \n');
    for a = 1:nelnodes(lmn)
      for i = 1:ndof
%         fprintf(outfile,'\n\n Node %d DOF %d \n',a,i); 
        rw = ndof*(a-1)+i;
        for b = 1:nelnodes(lmn)
          for k = 1:ndof
            cl = ndof*(b-1)+k;
              fprintf(outfile,'%d ',kel_num(rw,cl));
          end
        end
        fprintf(outfile,'\n');
      end
    end    
end

