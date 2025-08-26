function FEM_truss_plastic
%
%          Example FEM code using elastic-plastic truss elements
%          (illustrates Newton-Raphson and relaxation iteration)
%
%       The code reads an input file that defines member properties,
%       geometry and boundary conditions.  It is set up to run with
%       the example file 'FEM_truss_plastic_input' and will plot
%       the deformed shape of a simple 2-member truss for progressively
%       increasing force applied to the unconstrained joint, as well as
%       a graph showing the joint displacement as a function of load
%
%       This code is an example from the text
%       A.F. Bower 'Applied Mechanics of Solids' (2nd ed.)
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%%
% ================= Read data from the input file ==================
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
filename = '../Input Files/FEM_truss_plastic_input.txt';
while (~isfile(filename))
    f = msgbox({'The input file was not found';...
        'Please use the browser to select the file' });
    uiwait(f);
    [filename,location] = uigetfile('*.txt');
    filename = strcat(location,filename);
end

infile=fopen(filename,'r');
cellarray=textscan(infile,'%s');
[ndof,nnode,coord,nelem,connect,nprops,elprops,nfix,fixnodes,nload,loads,nsteps,tol] = read_file(cellarray);

fclose(infile);
%%
%  Specify solver parameters
%
% Use Relaxation or Newton (Newton will fail at the point of plastic collapse)
% load step
solvertype = 'Relaxation';
relax_factor = 0.001;
maxit_newton = 20;
maxit_relaxation = 10000;
%
% Plot the undeformed mesh
%
close all
figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

set(axes1,'DataAspectRatio',[1 1 1],'FontName','Times','FontSize',16,...
    'LineWidth',1,'PlotBoxAspectRatio',[1 1 1.66666666666667],'XGrid','on',...
    'XLimitMethod','tight','YGrid','on','YLimitMethod','tight','ZLimitMethod',...
    'tight','Box','On');
u1 = zeros(nsteps+1,1);
u2 = zeros(nsteps+1,1);
Pmag = zeros(nsteps+1,1);
loadfactor = 0;
loadscalefactor = 0.075;
plot_truss(ndof,nnode,coord,nelem,connect,nload,loads,loadscalefactor,loadfactor)
xlim([0 3.0])
ylim([0 3.0])
axis square
%
%   Unstretched member lengths
%
L0 = zeros(nelem,1);
for lmn = 1 : nelem
    coorda = coord(connect(lmn,1),:);
    coordb = coord(connect(lmn,2),:);
    L0(lmn) = sqrt(dot(coorda-coordb,coorda-coordb));
end
%
%   Initial plastic stretch and accumulated plastic strain
%
lamp(1:nelem) = 1.;
epsp(1:nelem) = 0.;
lampnew = lamp;
epspnew = epsp;

undef_coord = coord;   % Store original joint positions

for step = 1:nsteps
    
    loadfactor = step/nsteps;
    err = 1;
    w = coord;
    
    disp( ' ')
    disp(['Starting step number ',num2str(step),' of ',num2str(nsteps),' total steps'])
    disp(['Load Factor ',num2str(loadfactor)])
    
    
    if (strcmp(solvertype,'Newton'))
        %
        %================ Newton - Raphson iterations
        %
        nit = 0;
        
        disp('Starting Newton-Raphson solver')
        
        while(err>tol)
            %
            %===========  Assemble the global stiffness matrix and residual force vector ====================
            %
            Stif=zeros(ndof*nnode,ndof*nnode);
            resid = zeros(ndof*nnode,1);
            for lmn=1:nelem    % Loop over all the members
                %
                %   Set up the stiffness for the current element
                %
                a = connect(lmn,1);
                b = connect(lmn,2);
                [k,r,lampnew(lmn),epspnew(lmn)]=elstif(w(a,:),w(b,:),elprops(lmn,:),L0(lmn),lamp(lmn),epsp(lmn));
                
                %
                %   Add the current element stiffness to the global stiffness
                %
                for i = 1 : 2
                    for ii = 1 : ndof
                        rw = ndof*(connect(lmn,i)-1)+ii;
                        resid(rw) = resid(rw) + r(ndof*(i-1)+ii);
                        for j = 1 : 2
                            for jj = 1 : ndof
                                cl = ndof*(connect(lmn,j)-1)+jj;
                                Stif(rw,cl) = Stif(rw,cl) + k(ndof*(i-1)+ii,ndof*(j-1)+jj);
                            end
                        end
                    end
                end
            end
            
            %
            % ==================== Add external forces ============
            %
            for i=1:nload   % Loop over joints with external forces
                node=loads(i,1);
                resid(ndof*(node-1)+1)=resid(ndof*(node-1)+1) + loadfactor*loads(i,2);
                resid(ndof*(node-1)+2)=resid(ndof*(node-1)+2)+ loadfactor*loads(i,3);
                if (ndof>2)
                    resid(ndof*(node-1)+3)=resid(ndof*(node-1)+3)+ loadfactor*loads(i,4);
                end
            end
            %
            %   Modify the global stiffness and residual to include constraints
            %
            for i=1:nfix
                rw=ndof*(fixnodes(i,1)-1)+fixnodes(i,2);
                for j=1:2*nnode
                    Stif(rw,j)=0;
                end
                Stif(rw,rw)=1.0;
                resid(rw)=fixnodes(i,3);
            end
            %
            % ================== Solve the FEM equations ===================
            %
            dw=Stif\resid;
            for i = 1 : nnode
                for j = 1 : ndof
                    w(i,j) = w(i,j) + dw(ndof*(i-1)+j);
                end
            end
            err = sqrt(dot(resid,resid));
            
            nit = nit + 1;
            if (nit>maxit_newton)
                disp([' Newton solver failed to converge after ',num2str(nit),' iterations'])
                disp(' The analysis has terminated');
                return
            end
            
        end
        
        disp([' Newton solver converged after ',num2str(nit),' iterations'])
        
    elseif (strcmp(solvertype,'Relaxation'))
        
        
        nit = 0;
        disp('Starting Relaxation solver')
        
        while(err>tol)
            %
            %   Apply constraints
            %
            for i=1:nfix
                w(fixnodes(i,1),fixnodes(i,2)) = coord(fixnodes(i,1),fixnodes(i,2))+fixnodes(i,3);
            end
            
            %
            %===========  Assemble the global stiffness matrix and residual force vector ====================
            %
            resid = zeros(ndof*nnode,1);
            for lmn=1:nelem    % Loop over all the members
                %
                %   Set up the stiffness for the current element
                %
                a = connect(lmn,1);
                b = connect(lmn,2);
                [~,r]=elstif(w(a,:),w(b,:),elprops(lmn,:),L0(lmn),lamp(lmn),epsp(lmn));
                %
                %   Add the current element stiffness to the global stiffness
                %
                for i = 1 : 2
                    for ii = 1 : ndof
                        rw = ndof*(connect(lmn,i)-1)+ii;
                        resid(rw) = resid(rw) + r(ndof*(i-1)+ii);
                    end
                end
            end
            
            %
            % ==================== Add external forces ============
            %
            for i=1:nload   % Loop over joints with external forces
                node=loads(i,1);
                resid(ndof*(node-1)+1)=resid(ndof*(node-1)+1) + loadfactor*loads(i,2);
                resid(ndof*(node-1)+2)=resid(ndof*(node-1)+2)+ loadfactor*loads(i,3);
                if (ndof>2)
                    resid(ndof*(node-1)+3)=resid(ndof*(node-1)+3)+ loadfactor*loads(i,4);
                end
            end
            %
            %   Correct the residual for constrained nodes
            %
            for i=1:nfix
                rw=ndof*(fixnodes(i,1)-1)+fixnodes(i,2);
                resid(rw)=0;
            end
            dw = relax_factor*resid;
            for i = 1 : nnode
                for j = 1 : ndof
                    w(i,j) = w(i,j) + dw(ndof*(i-1)+j);
                end
            end
            err = sqrt(dot(resid,resid));
            
            nit = nit + 1;
            if (nit>maxit_relaxation)
                disp([' Relaxation failed to converge after ',num2str(nit),' iterations'])
                disp(' The analysis has terminated');
                return
            end
            
        end
        
        disp([' Relaxation solver converged after ',num2str(nit),' iterations'])
        
    end
    
    %  Update coords, plastic stretch, and accumulated plastic strain
    coord = w;
    lamp = lampnew;
    epsp = epspnew;
    plot_truss(ndof,nnode,coord,nelem,connect,nload,loads,loadscalefactor,loadfactor)
    
    u = coord - undef_coord;
    u1(step+1) = u(2,1);
    u2(step+1) = u(2,2);
    Pmag(step+1) = loadfactor*sqrt(0.5^2 + 5^2)/5;
    
    
end


x1 = [1.75,2.250];
x2 = x1 + loadscalefactor*[5.,0.];
h = annotation('arrow');
set(h,'parent', gca, ...
    'position', [x1(1),x1(2),x2(1)-x1(1),x2(2)-x1(2)], ...
    'HeadLength', 5, 'HeadWidth', 5, 'HeadStyle', 'cback3', ...
    'Color',[0,0,0],'LineWidth',1 );


annotation(figure1,'textbox',...
    [0.510714285714285 0.604761904761907 0.246428571428571 0.0942380952380986],...
    'String','$P/(AY_0) = 1$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor','none');

title('Elastic-Plastic Truss Demo','FontSize',15)
figure2 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');

set(axes2,'FontName','Times','FontSize',16,...
    'LineWidth',2,'XGrid','on',...
    'XLimitMethod','tight','YGrid','on','YLimitMethod','tight','ZLimitMethod',...
    'tight','Box','On');
plot(Pmag,u1,'Linewidth',2,'Color',[0,0,0])
plot(Pmag,u2,'Linewidth',2,'Color',[0,0,0])
ylim([-1 2.75])

xlabel('$\vert P \vert/(A Y_0)$','FontName','Times','Interpreter','latex');
ylabel('$u/H$','FontName','Times','Interpreter','latex');

% Create textarrow
annotation(figure2,'textarrow',[0.760714285714284 0.83392857142857],...
    [0.247619047619049 0.295238095238097],'String','$u_1/H$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FontName','Times New Roman');

% Create textarrow
annotation(figure2,'textarrow',[0.7625 0.832142857142857],...
    [0.633333333333334 0.583333333333334],'String','$u_2/H$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FontName','Times New Roman');

title('Elastic-Plastic Truss Demo','FontSize',15)
end

function [stress,tangent,lampnew,epspnew] = stressvstretch(lam,lamp,epsp,E,Y0,m)

tol = 1.e-04;

spred = E*(log(lam/lamp));
if (abs(spred)<Y0*(1+epsp)^m)
    stress = spred;
    tangent  = E/lam;
    lampnew = lamp;
    epspnew = epsp;
else
    lampnew = lamp;
    epspnew = epsp;
    err = 1.0;
    while (err>tol)
        if (spred>0)
            f = Y0*(1+epspnew)^m - E*log(lam/lampnew);
            depsdlamp = 1/lampnew;
            dfdlamp = m*Y0*(1+epspnew)^(m-1)*depsdlamp + E/lampnew;
        else
            f = -Y0*(1+epspnew)^m - E*log(lam/lampnew);
            depsdlamp = -1/lampnew;
            dfdlamp = -m*Y0*(1+epspnew)^(m-1)*depsdlamp + E/lampnew;
        end
        lampnew = lampnew - f/dfdlamp;
        epspnew = epsp + abs(log(lampnew/lamp));
        err = abs(f);
    end
    stress = E*log(lam/lampnew);
    tangent = (E/lam)*(m*Y0*(1+epspnew)^(m-1)/(E+m*Y0*(1+epspnew)^(m-1)));
end


end

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function [kel,rel,lampnew,epspnew] = elstif(coorda,coordb,props,L0,lamp,epsp)
%
len = sqrt(dot(coorda-coordb,coorda-coordb));
lam = len/L0;
%
%     Define the element stiffness
%
[stress,tangent,lampnew,epspnew] = stressvstretch(lam,lamp,epsp,props(2),props(3),props(4));
if (length(coorda)==2)
    gstif = [1,0,-1,0;...
        0,1,0,-1;...
        -1,0,1,0;...
        0,-1,0,1];
else
    gstif = [1,0,0,-1,0,0;...
        0,1,0,0,-1,0;...
        0,0,1,0,0,-1;...
        -1,0,0,1,0,0;...
        0,-1,0,0,1,0;
        0,0,-1,0,0,1];
end

vec = [coorda-coordb,coordb-coorda];
kel = zeros(length(vec),length(vec));
for i = 1 : length(vec)
    for j = 1 : length(vec)
        kel(i,j) = (tangent - 2*stress/lam)*vec(i)*vec(j)/(L0^2*lam^3)...
            + (stress/lam^2)*gstif(i,j);
    end
end
rel = -(props(1)/L0)*(stress/lam^2)*vec';
kel = (props(1)/L0)*kel;
end
%
%  =================== Function to extract variables from input file ======
%
function [ndof,nnode,coord,nelem,connect,nprops,elprops,nfix,fixnodes,nload,loads,nsteps,tol] = read_file(cellarray)
%
%  Extract no. dof, no. nodes and nodal coordinates
%
ndof=str2num(cellarray{1}{2});
nnode=str2num(cellarray{1}{4});
dum=6;
coord=zeros(nnode,ndof);
for i=1:nnode
    coord(i,1) = str2num(cellarray{1}{dum});
    dum=dum+1;
    coord(i,2) = str2num(cellarray{1}{dum});
    dum=dum+1;
    if (ndof>2) coord(i,3) = str2num(cellarray{1}{dum}); dum = dum+1; end
end
%
%   Extract no. elements and connectivity
%
dum=dum + 1;
nelem=str2num(cellarray{1}{dum});
connect = zeros(nelem,3);
dum = dum + 2;
for i = 1 : nelem
    for j = 1 : 2
        connect(i,j) = str2num(cellarray{1}{dum});
        dum=dum+1;
    end
end
% Extract no. element props and property values
dum = dum + 1;
nprops=str2num(cellarray{1}{dum})
dum = dum + 1;
elprops = zeros(nelem,nprops);
for i = 1 : nelem
    for j = 1 : nprops
        elprops(i,j) = str2num(cellarray{1}{dum});
        dum = dum + 1;
    end
end
%  Extract no. load steps
dum = dum + 1;
nsteps = str2num(cellarray{1}{dum});
dum = dum + 2;
tol = str2num(cellarray{1}{dum});
%
%   Extract no. nodes with prescribed displacements and the prescribed displacements
%
dum = dum + 2;
nfix=str2num(cellarray{1}{dum})
dum = dum + 4;
fixnodes = zeros(nfix,3);
for i = 1 : nfix
    for j = 1 : 3
        fixnodes(i,j) = str2num(cellarray{1}{dum});
        dum=dum+1;
    end
end
fixnodes
%
%   Extract no. loaded element faces, with the loads
%
dum = dum + 1;
nload=str2num(cellarray{1}{dum})
dum=dum + 3;
loads = zeros(nload,1+ndof);
for i = 1 : nload
    for j=1:1+ndof
        loads(i,j)=str2num(cellarray{1}{dum});
        dum=dum+1;
    end
end

end
function plot_truss(ndof,nnode,coord,nelem,connect,nload,loads,loadscalefactor,loadfactor)
hold on
if (ndof==3)
    
    for i = 1 : nelem
        xvals = [coord(connect(i,1),1),coord(connect(i,2),1)];
        yvals = [coord(connect(i,1),2),coord(connect(i,2),2)];
        zvals = [coord(connect(i,1),3),coord(connect(i,2),3)];
        plot3(xvals,yvals,zvals,'LineWidth',2,'Color',[0,0,0]);
    end
    
    scatter3(coord(:,1),coord(:,2),coord(:,3),'MarkerFaceColor',[1,0,0]);
    
    
else
    for i = 1 : nelem
        xvals = [coord(connect(i,1),1),coord(connect(i,2),1)];
        yvals = [coord(connect(i,1),2),coord(connect(i,2),2)];
        plot(xvals,yvals,'LineWidth',2,'Color',[0,0,0]);
    end
    
    scatter(coord(:,1),coord(:,2),'MarkerFaceColor',[1,0,0]);
    
    if (loadfactor>0)
        for i = 1:nload
            node = loads(i,1);
            x1 = coord(node,:);
            x2 = x1 + loadscalefactor*loadfactor*loads(i,2:3);
            h = annotation('arrow');
            set(h,'parent', gca, ...
                'position', [x1(1),x1(2),x2(1)-x1(1),x2(2)-x1(2)], ...
                'HeadLength', 5, 'HeadWidth', 5, 'HeadStyle', 'cback3', ...
                'Color',[0,0,0],'LineWidth',1 );
        end
    end
    
end
hold off


end