function fem_conststrain_triangles
%
%
%       Simple FEA program using axisymmetric linear triangle elements
%
%       The code is ste up to calculate the displacements and stresses
%       in an internally pressurized thick-walled circular cylinder
%       constrained to deform in plane strain
%
%       This code solves problem 7.31 from the text
%       A.F. Bower, 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%   It was downloaded from
%   https://github.com/albower/Applied_Mechanics_of_Solids
%
%% =================== Generate mesh and BCs ===================

aa = 1;   % Inner radius
bb = 3;   % Outer radius
p = 1.;   % Pressure

E = 100;  % Young's modulus
nu = 0.3; % Poisson's ratio
nelem = 30;  % No-elements in radial direction
nnode = nelem + 2; % No. nodes
% Generate the nodal coordinates
coord(1:nnode/2,1) = [aa:2*(bb-aa)/nelem:bb];
coord(1:nnode/2,2) = 0.;
coord(nnode/2+1:nnode,1) = [bb:-2*(bb-aa)/nelem:aa];
coord(nnode/2+1:nnode,2) = 1;
% Generate the element connectivity
connect(1:2:nelem-1,1) = 1:nnode/2-1;
connect(1:2:nelem-1,2) = 2:nnode/2;
connect(1:2:nelem-1,3) = nnode:-1:nnode/2+2;
connect(2:2:nelem,1) = 2:nnode/2;
connect(2:2:nelem,2) = nnode-1:-1:nnode/2+1;
connect(2:2:nelem,3) = nnode:-1:nnode/2+2;
% Constrain the axial displacements of all the nodes
nfix = nnode;
fixnodes(1:nnode,1) = 1:nnode;
fixnodes(1:nnode,2) = 2;
fixnodes(1:nnode,3) = 0.;

% Set up the distributed loads acting on the interior wall
ndload = 1;
dloads(1,1) = 1;
dloads(1,2) = 3;
dloads(1,3) = p;
dloads(1,4) = 0.0; 

% Calculate the exact displacements at each node 
ur_exact(1:nnode/2) = (1+nu)*aa^2*bb^2*p/(bb^2-aa^2)/E.*(1./coord(1:nnode/2,1) + (1-2*nu).*coord(1:nnode/2,1)./bb^2);

%
% Plot the undeformed mesh as a check
close all
figure;
triplot(connect,coord(:,1),coord(:,2),'k');
axis equal
xlim([0 3])
title('Undeformed Mesh')
%
% =============  Define the D matrix ===========================
%
Dmat = [[1-nu,nu,nu,0];[nu,1-nu,nu,0];[nu,nu,1-nu,0];[0,0,0,(1-2*nu)/2]]*E/((1+nu)*(1-2*nu));
%
%% ===========  Assemble the global stiffness matrix ====================
%
Stif=zeros(2*nnode,2*nnode);
for lmn=1:nelem    % Loop over all the elements
    %
    %   Set up the stiffness for the current element
    %
    a = connect(lmn,1);
    b = connect(lmn,2);
    c = connect(lmn,3);
    k=elstif(coord(a,1),coord(a,2),coord(b,1),coord(b,2),coord(c,1),coord(c,2),Dmat);
    %
    %   Add the current element stiffness to the global stiffness
    %
    for i = 1 : 3
        for ii = 1 : 2
            for j = 1 : 3
                for jj = 1 : 2
                    rw = 2*(connect(lmn,i)-1)+ii;
                    cl = 2*(connect(lmn,j)-1)+jj;
                    Stif(rw,cl) = Stif(rw,cl) + k(2*(i-1)+ii,2*(j-1)+jj);
                end
            end
        end
    end
end
%
% ==================== Assemble global force vector ============
%
%   Define the force
%
resid=zeros(2*nnode);
pointer=[2,3,1];
for i=1:ndload   % Loop over elements with loaded faces
    lmn=dloads(i,1);
    face=dloads(i,2);
    a=connect(lmn,face);
    b=connect(lmn,pointer(face));
    r=elresid(coord(a,1),coord(a,2),coord(b,1),coord(b,2),dloads(i,3),dloads(i,4));
    resid(2*a-1)=resid(2*a-1)+r(1);
    resid(2*a)=resid(2*a)+r(2);
    resid(2*b-1)=resid(2*b-1)+r(3);
    resid(2*b)=resid(2*b)+r(4);
end
%
%   Modify the global stiffness and residual to include constraints
%
for i=1:nfix
    rw=2*(fixnodes(i,1)-1)+fixnodes(i,2);
    for j=1:2*nnode
        Stif(rw,j)=0;
    end
    Stif(rw,rw)=1.0;
    resid(rw)=fixnodes(i,3);
end
%
%% ================== Solve the FEM equations ===================
%
u=Stif\resid;
%
%% ================== Post Processing ============================
%  Compare solution to analytical solution
%
figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]); 
axes3 = axes('Parent',figure1,...
    'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')

plot(coord(1:nnode/2,1),u(1:2:nnode-1),'DisplayName','FEA Solution','MarkerFaceColor',[1 1 1],...
    'MarkerSize',5,'Color','k',...
    'Marker','o',...
    'LineStyle','none')
plot(coord(1:nnode/2,1),ur_exact,'DisplayName','Exact solution','LineWidth',2,'Color','k')
 ylabel({'Displacement '},'Interpreter','latex','FontSize',16);
 xlabel({'Position $r$'},'Interpreter','latex','FontSize',16);
   legend1 = legend(axes3,'show');
 set(legend1,...
    'Position',[0.596130964913888 0.649206355805441 0.278869035086111 0.260714278334663],...
    'Interpreter','tex',...
    'FontSize',14,...
    'FontName','Times New Roman');
%
% Compute stresses for plotting
%
for lmn = 1 : nelem   % Loop over all the elements
    a = connect(lmn,1);
    b = connect(lmn,2);
    c = connect(lmn,3);
    xa = coord(a,1);
    ya = coord(a,2);
    xb = coord(b,1);
    yb = coord(b,2);
    xc = coord(c,1);
    yc = coord(c,2);
    uxa = u(2*a-1);
    uya = u(2*a);
    uxb = u(2*b-1);
    uyb = u(2*b);
    uxc = u(2*c-1);
    uyc = u(2*c);
    strain = elstrain(xa,ya,xb,yb,xc,yc,uxa,uya,uxb,uyb,uxc,uyc);
    stress = Dmat*strain;
    srr_fea(lmn) = stress(1);
    szz_fea(lmn) = stress(2);
    sqq_fea(lmn) = stress(3);
    rplot(lmn) = (xa+xb+xc)/3.;
    [srr_exact(lmn),sqq_exact(lmn),szz_exact(lmn)] = exact_stresses(rplot(lmn),p,aa,bb,nu);
end

figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]); 
axes3 = axes('Parent',figure1,...
    'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')


plot(rplot,srr_fea,'DisplayName','\sigma_r_r (FEA)','MarkerFaceColor','k',...
    'MarkerEdgeColor','k',...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none');
hold on
plot(rplot,srr_exact,'DisplayName','\sigma_r_r (Exact)','LineWidth',2,...
    'Color','k');
plot(rplot,sqq_fea,'DisplayName','\sigma_\theta_\theta (FEA)',...
    'MarkerFaceColor',[1 1 1],...
    'MarkerSize',4,...
    'Marker','o',...
    'LineStyle','none',...
    'Color','k');
plot(rplot,sqq_exact,'DisplayName','\sigma_\theta_\theta (Exact)','LineWidth',2,...
    'Color','k','LineStyle','--');

 ylabel({'Stress '},'Interpreter','latex','FontSize',16);
 xlabel({'Position $r$'},'Interpreter','latex','FontSize',16);

   legend1 = legend(axes3,'show');
 set(legend1,...
    'Position',[0.612202393485316 0.613492070658819 0.278869035086111 0.298809515294575],...
    'Interpreter','tex',...
    'FontSize',14,...
    'FontName','Times New Roman');



end

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function kel = elstif(xa,ya,xb,yb,xc,yc,Dmat)
        %
        %     Define the B matrix
        %
        x = (xa+xb+xc)/3;
        y = (ya+yb+yc)/3;
        Na = ( (y-yb)*(xc-xb) - (x-xb)*(yc-yb) )/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        Nb = ( (y-yc)*(xa-xc) - (x-xc)*(ya-yc) )/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        Nc = ( (y-ya)*(xb-xa) - (x-xa)*(yb-ya) )/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        nax = -(yc-yb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nay =  (xc-xb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nbx = -(ya-yc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        nby =  (xa-xc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        ncx = -(yb-ya)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        ncy =  (xb-xa)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        area = (1/2)*abs( (xb-xa)*(yc-ya) - (xc-xa)*(yb-ya) );
        Bmat = [[nax,  0,nbx,  0,ncx,  0];...
            [0,nay,  0,nby,  0,ncy];...
            [Na/x,  0, Nb/x, 0,Nc/x, 0];...
            [nay,nax,nby,nbx,ncy,ncx]];
        %
        %     Define the element stiffness
        %
        kel = 2*pi*x*area*transpose(Bmat)*Dmat*Bmat;
end


%====================== ELEMENT FORCE VECTOR ==============
%
function rel = elresid (xa,ya,xb,yb,tx,ty)
        length = 2*pi*sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb));
        A = length*(xa+2*xb)/6;
        B = length*(2*xa+xb)/6;
        rel = [A*tx,A*ty,B*tx,B*ty];
end
%
%
%    Function to calculate the element strains
%
    function strain = elstrain (xa,ya,xb,yb,xc,yc,uax,uay,ubx,uby,ucx,ucy)
        %
        %   B matrix
        %
        x = (xa+xb+xc)/3;
        y = (ya+yb+yc)/3;
        Na = ( (y-yb)*(xc-xb) - (x-xb)*(yc-yb) )/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        Nb = ( (y-yc)*(xa-xc) - (x-xc)*(ya-yc) )/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        Nc = ( (y-ya)*(xb-xa) - (x-xa)*(yb-ya) )/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        nax = -(yc-yb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nay =  (xc-xb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nbx = -(ya-yc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        nby =  (xa-xc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        ncx = -(yb-ya)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        ncy =  (xb-xa)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        Bmat = [[nax,  0,nbx,  0,ncx,  0];...
            [0,nay,  0,nby,  0,ncy];...
            [Na/x,  0, Nb/x, 0,Nc/x, 0];...
            [nay,nax,nby,nbx,ncy,ncx]];
        %
        %     Element displacement vector
        %
        uel = [uax;uay;ubx;uby;ucx;ucy];
        %
        %     Element strains
        %
        strain = Bmat*uel;

    end
    function [srr,sqq,szz] = exact_stresses(r,pa,a,b,nu)
    
        srr = pa*a^2/(b^2-a^2)*(1-b^2/r^2);
        sqq = pa*a^2/(b^2-a^2)*(1 + b^2/r^2);
        szz = 2*nu*pa*a^2/(b^2-a^2);
    

    end
