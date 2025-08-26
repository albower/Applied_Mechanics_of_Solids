function fem_beam_Euler
%
%          Example FEM code using elastic Euler beam elements
%
%       The code reads an input file that defines the geometry and
%       cross-section properties of one or more beams. The code will
%       plot the shape of the deformed beams, and print the nodal
%       displacements and rotations to a file.
%
%       If the code is run with the demonstration input file
%       'FEM_beam_euler_input.txt' it will plot a graph comparing the
%       predicted deflection of a cantilever beam with the exact
%       (Euler-Bernoulli) solution.
%
%       This code is an example from the text
%       A.F. Bower 'Applied Mechanics of Solids' (2nd ed.)
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids


% ================= Read data from the input file ==================
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
filename = '../Input Files/FEM_beam_euler_input.txt';
while (~isfile(filename))
    f = msgbox({'The input file was not found';...
        'Please use the browser to select the file' });
    uiwait(f);
    [filename,location] = uigetfile('*.txt');
    filename = strcat(location,filename);
end

infile=fopen(filename,'r');
outfile=fopen('FEM_results.txt','w');


% Read the input file into a cell array
cellarray=textscan(infile,'%s');
[nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,E,mu,e1dir,inertias,areas,...
    nfix,fixnodes,nload,loads] = read_file(cellarray) ;
fclose(infile);

plotexactsol = false;
if (strcmp(filename,'../Input Files/FEM_beam_euler_input.txt')) plotexactsol = true; end

% Plot the initial mesh
close all
figure
plot_beam(nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,[0,1,0]);
%
%  Assemble the stiffness matrix
%
Stif=zeros(6*n_total_nodes,6*n_total_nodes);
resid = zeros(6*n_total_nodes,1);
for beam = 1:nbeams    % Loop over all the beams
    if (beam<nbeams)
        end_element = first_element(beam+1);
    else
        end_element = n_total_elements;
    end
    for lmn=first_element:end_element  % Loop over elements on current beam
        %
        %   Set up the stiffness for the current element
        %
        a = connect(lmn,1);
        b = connect(lmn,2);
        [k,r]=elstif(coord(a,:),coord(b,:),E(beam),mu(beam),e1dir(beam,:),inertias(beam,:),areas(beam));
        
        %
        %   Add the current element stiffness to the global stiffness
        %
        for i = 1 : 2
            for ii = 1 : 6
                rw =6*(connect(lmn,i)-1)+ii;
                resid(rw) = resid(rw) + r(6*(i-1)+ii);
                for j = 1 : 2
                    for jj = 1 : 6
                        cl = 6*(connect(lmn,j)-1)+jj;
                        Stif(rw,cl) = Stif(rw,cl) + k(6*(i-1)+ii,6*(j-1)+jj);
                    end
                end
            end
        end
    end
end
%
%  Assemble the force vector
for i=1:nload   % Loop over nodes with external forces and moments
    node=loads(i,1);
    for k = 1:6
        resid(6*(node-1)+k)=resid(6*(node-1)+k) + loads(i,1+k);
    end
end
%
%   Modify the global stiffness and residual to include constraints
%
for i=1:nfix
    rw=6*(fixnodes(i,1)-1)+fixnodes(i,2);
    for j=1:6*n_total_nodes
        Stif(rw,j)=0;
    end
    Stif(rw,rw)=1.0;
    resid(rw)=fixnodes(i,3);
end
%
% Solve the equations
%
u=Stif\resid;

% Plot the displaced mesh
displaced_coords(:,1) = coord(:,1) + u(1:6:end-5);
displaced_coords(:,2) = coord(:,2) + u(2:6:end-4);
displaced_coords(:,3) = coord(:,3) + u(3:6:end-3);
plot_beam(nbeams,n_total_nodes,first_node,displaced_coords,n_total_elements,first_element,connect,[0,0,0]);
axis equal
title('Rotate view to see deflection')

%
% Print displacements and rotations to an output file
fprintf(outfile,'%s\n','Nodal Displacements (u) and rotations (q):');
fprintf(outfile,'%s\n',' Node    u1       u2        u3        q1        q2        q3');
for i = 1 : n_total_nodes
    fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
        i,u(6*i-5),u(6*i-4),u(6*i-3),u(6*i-2),u(6*i-1),u(6*i));
end
fclose(outfile);


% Plot the exact and FEA solution on the same graph if requested
if (plotexactsol)
    x = coord(:,3);
    disp = u(1:6:end-5)*E(1)*inertias(1,3)/(loads(2)*x(end)^3);
    rot = u(5:6:end-1)*E(1)*inertias(1,3)/(loads(2)*x(end)^2);
    x = x/x(end);
    
    uexact = x.^2.*(3-x)/6;
    rexact = x.*(2-x)/2;
    
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes2 = axes('Parent',figure1);
    hold(axes2,'on');
    set(axes2,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On');
    plot(x,disp,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
    plot(x,rot,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
    plot(x,uexact,'LineWidth',2,'LineStyle','-','Color',[0 0 0])
    plot(x,rexact,'LineWidth',2,'LineStyle','--','Color',[0 0 0])
    xlabel('$x_3/L$','FontName','Times','Interpreter','latex');
    ylabel('Deflection and rotation','FontName','Times','Interpreter','latex');
    
    linex = [0.05,0.25];
    liney = [0.55,0.55];
    plot(linex,liney,'LineWidth',2,'LineStyle','-','Color',[0 0 0])
    linex = [0.05,0.25];
    liney = [0.5,0.5];
    plot(linex,liney,'LineWidth',2,'LineStyle','--','Color',[0 0 0])
    linex = [0.05,0.25];
    liney = [0.45,0.45];
    plot(linex,liney,'Marker','o','LineStyle','none','MarkerEdgeColor',[0 0 0],'Color',[0 0 0])
    
    annotation(figure1,'textbox',...
        [0.351785714285714 0.769047619047619 0.267857142857143 0.0870952380952403],...
        'String','Exact Solution',...
        'Interpreter','latex',...
        'FontSize',15,...
        'FontName','Times New Roman',...
        'FitBoxToText','off',...
        'EdgeColor','none');
    
    annotation(figure1,'textbox',...
        [0.355357142857142 0.67857142857143 0.1125 0.0870952380952403],...
        'String','FEA',...
        'Interpreter','latex',...
        'FontSize',15,...
        'FontName','Times New Roman',...
        'FitBoxToText','off',...
        'EdgeColor','none');
    
    annotation(figure1,'textarrow',[0.549999999999997 0.64285714285714],...
        [0.434714285714291 0.380952380952386],'String','$E I u/P L^3$',...
        'Interpreter','latex',...
        'FontSize',15,...
        'FontName','Times New Roman');
    
    annotation(figure1,'textarrow',[0.632142857142853 0.564285714285714],...
        [0.618047619047632 0.657142857142859],'String','$E I \theta/P L^2$',...
        'Interpreter','latex',...
        'FontSize',15,...
        'FontName','Times New Roman');
    
    title('Euler Beam Demo','FontSize',15)
end


end

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function [kel,rel] = elstif(coorda,coordb,E,mu,e1dir,inertias,area)

        kel = zeros(12,12);
        rel = zeros(12,1);
%  Two point Gaussian integration points and weights        
        n_int_pts = 2;
        w(1:2) = 1;
        xi(1) = -1./sqrt(3);
        xi(2) = 1./sqrt(3);
        len = sqrt(dot(coorda-coordb,coorda-coordb));
        e3dir = (coordb-coorda)/len;
        e1dir = e1dir - dot(e1dir,e3dir)*e3dir;        
        e1dir = e1dir/norm(e1dir);
        e2dir = cross(e3dir,e1dir);
        R = [e1dir,0,0,0,0,0,0,0,0,0;
             e2dir,0,0,0,0,0,0,0,0,0;
             e3dir,0,0,0,0,0,0,0,0,0;
            0,0,0,e1dir,0,0,0,0,0,0;
            0,0,0,e2dir,0,0,0,0,0,0;
            0,0,0,e3dir,0,0,0,0,0,0;
            0,0,0,0,0,0,e1dir,0,0,0;
            0,0,0,0,0,0,e2dir,0,0,0;
            0,0,0,0,0,0,e3dir,0,0,0;
            0,0,0,0,0,0,0,0,0,e1dir;
            0,0,0,0,0,0,0,0,0,e2dir;
            0,0,0,0,0,0,0,0,0,e3dir];
        %
        %     Define the element stiffness
        %
        for i = 1:n_int_pts
             dxidx = 2./len;
             dN1 = -0.5*dxidx; dN2 = 0.5*dxidx;
             ddMu1 = 3*xi(i)*dxidx^2/2;  ddMu2 = -3*xi(i)*dxidx^2/2;
             ddMq1 = (3*xi(i)-1)*dxidx^2/4;  ddMq2 = (3*xi(i)+1)*dxidx^2/4;
            
            B = [0, -ddMu1, 0, ddMq1, 0,0,0,-ddMu2,0,ddMq2,0,0;...
                 ddMu1,0,0,0,ddMq1,0,ddMu2,0,0,0,ddMq2,0;...
                 0,0,0,0,0,dN1,0,0,0,0,0,dN2;...
                 0,0,dN1,0,0,0,0,0,dN2,0,0,0];
             
             D = [E*inertias(1),-E*inertias(2), 0, 0;...
                 -E*inertias(2),E*inertias(3), 0,0;...
                 0,0,mu*inertias(4),0;...
                 0,0,0,E*area];
             kel = kel + transpose(R)*transpose(B)*D*B*R*w(i)*len/2;
        end

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
    nbeams=str2double(cellarray{1}{2});
    dum = 3;
    for beam=1:nbeams
        dum = dum + 2;
        nnode=str2double(cellarray{1}{dum});

        dum=dum + 2;
        for i=1:nnode
          for k = 1:3
             coord(i,n_total_nodes + k) = str2double(cellarray{1}{dum});
             dum=dum+1;
          end
        end
        first_node(beam) = n_total_nodes + 1;
        n_total_nodes = n_total_nodes + nnode;
        %
    %   Extract no. elements and connectivity
    %
        dum=dum + 1;
        nelem=str2double(cellarray{1}{dum});

        dum = dum + 2;
        for i = 1 : nelem
          for j = 1 : 2
            connect(n_total_elements+i,j) = str2double(cellarray{1}{dum});
            dum=dum+1;
          end
        end
        first_element(beam) = n_total_elements + 1;
        n_total_elements = n_total_elements + nelem;    

    % Extract no. element props and property values
       dum = dum + 1;
       E(beam) = str2double(cellarray{1}{dum});
       mu(beam) = str2double(cellarray{1}{dum+1});
       dum = dum + 3;
       for i=1:3
         e1dir(beam,i) = str2double(cellarray{1}{dum});
         dum = dum + 1;
       end
       dum = dum + 1;
       for i = 1:4
           inertias(beam,i) = str2double(cellarray{1}{dum});
           dum = dum + 1;
       end
       dum = dum + 1;
       areas(beam)  = str2double(cellarray{1}{dum});
       dum = dum + 1;
    end
    %
    %   Extract no. nodes with prescribed displacements and the prescribed displacements
    %
    dum = dum + 1;
    nfix=str2double(cellarray{1}{dum});
    dum = dum + 4;
    fixnodes = zeros(nfix,3);
    for i = 1 : nfix
        for j = 1 : 3
            fixnodes(i,j) = str2double(cellarray{1}{dum});
            dum=dum+1;
        end
    end
    %
    %   Extract no. loaded element faces, with the loads
    %
    dum = dum + 1;
    nload=str2double(cellarray{1}{dum});
    dum=dum + 3;
    loads = zeros(nload,3);
    for i = 1 : nload
        for j=1:7
            loads(i,j)=str2double(cellarray{1}{dum});
            dum=dum+1;
        end
    end

end
function plot_beam(nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,color)
hold on

    for i = 1 : n_total_elements
        xvals = [coord(connect(i,1),1),coord(connect(i,2),1)];
        yvals = [coord(connect(i,1),2),coord(connect(i,2),2)];
        zvals = [coord(connect(i,1),3),coord(connect(i,2),3)];
        plot3(xvals,yvals,zvals,'LineWidth',2,'Color',color);
    end

    scatter3(coord(:,1),coord(:,2),coord(:,3),'MarkerFaceColor',[1,0,0]);

end