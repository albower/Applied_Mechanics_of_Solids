function Solution8_37
%
%          Example FEM code using elastic Euler beam elements
%
%       This code solves problem 8.36(a) and 8.36(c) from the text
%       A.F. Bower 'Solved Problems in Mechanics of Solids'
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids


% ================= Read data from the input file ==================
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
    filename = '../Input Files/Solution8_37.txt';
    while (~isfile(filename))
        f = msgbox({'The input file was not found';...
            'Please use the browser to select the file' });
        uiwait(f);
        [filename,location] = uigetfile('*.txt');
        filename = strcat(location,filename);
    end

    infile=fopen(filename,'r');

% The moments and internal forces are returned
% as components in the e1,e2,e3 basis aligned with the beam, with 
% e3 parallel to beam axis, and e1 transverse to the
% beam with orientation specified by the user.    
%   Plot M_2 and T_1 for beam with length 10 and end force 1
    Mcomptoplot = 2;
    Tcomptoplot = 1;
    L = 10;
    P = 1;
    
    % Read the input file into a cell array
    cellarray=textscan(infile,'%s');
    [nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,E,mu,rho,e1dir,inertias,areas,...
                            nfix,fixnodes,nload,loads] = read_file(cellarray) ;
    fclose(infile);

    close all

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
    scale = 1;
    displaced_coords(:,1) = coord(:,1) + scale*u(1:6:end-5);
    displaced_coords(:,2) = coord(:,2) + scale*u(2:6:end-4);
    displaced_coords(:,3) = coord(:,3) + scale*u(3:6:end-3);
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes3 = axes('Parent',figure1,...
        'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')
    plot_beam(nbeams,n_total_nodes,first_node,displaced_coords,n_total_elements,first_element,connect,[0,0,0]);
    axis equal
    title('Deflection of cantilever beam')
    annotation(figure1,'textbox',...
        [0.345642857142856 0.266666666666669 0.382928571428572 0.0690476190476213],...
        'String',{['Displacement Scale Factor: ',num2str(scale)]},...
        'FitBoxToText','on',...
        'BackgroundColor',[1 1 1]);
    
% Print displacements and rotations to an output file
    
    outfile=fopen('FEM_results.txt','w');
    
    fprintf(outfile,'%s\n','Nodal Displacements (u) and rotations (q):');
    fprintf(outfile,'%s\n',' Node    u1       u2        u3        q1        q2        q3');
    for i = 1 : n_total_nodes
        fprintf(outfile,'%3d %8.4g %8.4g %8.4g %8.4g %8.4g %8.4g\n',...
            i,u(6*i-5),u(6*i-4),u(6*i-3),u(6*i-2),u(6*i-1),u(6*i));
    end

% Print internal forces to the output file
    
    for beam = 1:nbeams
        fprintf(outfile,'\n\n%s%d\n','Beam number',beam);
        if (beam<nbeams)
            end_element = first_element(beam+1);
        else
            end_element = n_total_elements;
        end
        
        
        [~,~,x]=internalforces(coord(connect(1,1),:),coord(connect(1,2),:),u(1:6),u(7:12),E(beam),mu(beam),e1dir(beam,:),inertias(beam,:),areas(beam));
        [~,nintpts] = size(x);
        nintpts = nintpts-1;
        xprevi = coord(connect(1,1),:)';
        xprevc = coord(connect(1,1),:)';
        arcleni = zeros(n_total_elements*nintpts,1); % Arc length at int. pts.
        arclenc = zeros(n_total_elements,1); % Arc length at element center
        Tplot = arclenc;
        Mplot = arcleni;
        
        for lmn=first_element:end_element
            
            a = connect(lmn,1);
            b = connect(lmn,2);
            dofa = u(6*a-5:6*a);
            dofb = u(6*b-5:6*b);
            
            [T,M,x]=internalforces(coord(a,:),coord(b,:),dofa,dofb,E(beam),mu(beam),e1dir(beam,:),inertias(beam,:),areas(beam));
            
            fprintf(outfile,'\n%s\n','Moments (M):');
            fprintf(outfile,'%s\n',' Element Int Pt     x1       x2       x3        M1        M2        M3');
            
            for i = 1:nintpts
                fprintf(outfile,'%8d %6d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
                    lmn,i,x(1:3,i),M(1:3,i));
                if ((lmn-1)*nintpts+i)>1
                    arcleni((lmn-1)*nintpts+i) = arcleni((lmn-1)*nintpts+i-1) + norm(x(:,i) - xprevi);
                    xprevi = x(:,i);
                else
                    arcleni((lmn-1)*nintpts+i) =  norm(x(:,i) - xprevi);
                    xprevi = x(:,i);
                end
                Mplot((lmn-1)*nintpts+i) = M(Mcomptoplot,i);
            end
            fprintf(outfile,'\n%s\n','Forces (M):');
            fprintf(outfile,'%s\n',' Element    x1       x2       x3        T1        T2        T3');
            fprintf(outfile,'%8d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
                lmn,x(1:3,nintpts+1),T(1:3));
            if (lmn)>1
                arclenc(lmn) = arclenc(lmn-1) + norm(x(:,nintpts+1) - xprevc);
                xprevc = x(:,nintpts+1);
            else
                arclenc(lmn) =  norm(x(:,nintpts+1) - xprevc);
                xprevc = x(:,nintpts+1);
            end
            Tplot(lmn) = T(Tcomptoplot);
        end
    end
    fclose(outfile);
    
    
    
    figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
    axes3 = axes('Parent',figure1,...
        'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
    hold(axes3,'on');
    set(axes3,'FontName','Times','FontSize',16,...
        'LineWidth',2,'XGrid','on','YGrid','on',...
        'Box','On')
    plot(arcleni/L,Mplot/(P*L),'LineStyle','none','Color','k','Marker','o','Displayname','Moment (FEA)')
    hold on
    plot(arclenc/L,Tplot/P,'LineStyle','none','Color','k','MarkerFaceColor','k','Marker','o','Displayname','Shear force (FEA)')
    plot([0,1],[1,1],'LineStyle','--','Color','k','LineWidth',2,'Displayname','Shear force (exact)')
    plot([0,1],[1,0],'LineStyle','-','Color','k','LineWidth',2,'Displayname','Moment (exact)')
    xlabel({'Position $x_1/L$'},'Interpreter','latex','FontSize',16);
    ylabel({'$M/(PL)$ or $T/P$'},'Interpreter','latex','FontSize',16);
    title('Forces in Euler beam')
    legend1 = legend(axes3,'show');
    set(legend1,...
        'Position',[0.196067879785483 0.222301592935656 0.387923776771482 0.226428565979004],...
        'Interpreter','latex',...
        'FontSize',14,...
        'FontName','Times New Roman');
    
    
    
    



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
        e1dir = e1dir/norm(e1dir);
        e2dir = -cross(e1dir,e3dir);
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
function mel = elmass(coorda,coordb,rho,e1dir,inertias,area)

        mel = zeros(12,12);
%  Two point Gaussian integration points and weights        
        n_int_pts = 2;
        w(1:2) = 1;
        xi(1) = -1./sqrt(3);
        xi(2) = 1./sqrt(3);
        len = sqrt(dot(coorda-coordb,coorda-coordb));
        e3dir = (coordb-coorda)/len;
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
            N1 = 0.5*(1-xi(i)); N2 = 0.5*(1+xi(i));
            Mu1 = (xi(i)-1)^2*(xi(i)+2)/4; Mu2 = -(xi(i)+1)^2*(xi(i)-2)/4;
            Mq1 = (xi(i)-1)^2*(xi(i)+1)/8; Mq2 = (xi(i)-1)^2*(xi(i)-1)/8;
            dMu1 = 3*(xi(i)-1)*(xi(i)+1)/(2*len);  dMu2 = -3*(xi(i)+1)*(xi(i)-1)/(2*len);
            dMq1 = (3*xi(i)^2 - 2*xi(i) - 1)/(4*len);   dMq2 = (3*xi(i)^2 + 2*xi(i) - 1)/(4*len);
           
            N = [Mu1, 0, 0, 0, Mq1,0,Mu2,0,0,0,Mq2,0;...
                 0,Mu1,0,-Mq1,0,0,0,Mu2,0,-Mq2,0,0;...
                 0,0,N1,0,0,0,0,0,N2,0,0,0;...
                 0,-dMu1,0,dMq1,0,0,0,-dMu2,0,dMq2,0,0;...
                 dMu1,0,0,0,dMq1,0,dMu2,0,0,0,dMq2,0;...
                 0,0,0,0,0,N1,0,0,0,0,0,N2];

             mu = rho*area*eye(6);
             mu(4,4) = rho*inertias(1); 
             mu(5,5) = rho*inertias(3);
             mu(6,6) = rho*inertias(4);
             mu(4,5) = -rho*inertias(2); mu(5,4) = mu(4,5);
             mel = mel + transpose(R)*transpose(N)*mu*N*R*w(i)*len/2;
        end

end
function [T,M,x] = internalforces(coorda,coordb,dofa,dofb,E,mu,e1dir,inertias,area)
        n_int_pts = 2;
        xi(1) = -1./sqrt(3);
        xi(2) = 1./sqrt(3);

        len = sqrt(dot(coorda-coordb,coorda-coordb));
        e3dir = (coordb-coorda)/len;
        e1dir = e1dir - dot(e1dir,e3dir)*e3dir;
        e1dir = e1dir/norm(e1dir);
        e2dir = -cross(e1dir,e3dir);
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
        T = zeros(3,1);
        M = T;
        x = zeros(3,n_int_pts+1);
        avT3 = 0;
        for i = 1:n_int_pts

            N1 = (1-xi(i))/2; N2 = (1+xi(i))/2;
            x(:,i) = coorda*N1 + coordb*N2;
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
              intFs = D*B*R*[dofa;dofb];
              M(1:3,i) = intFs(1:3);
              avT3 = avT3 + intFs(4);
        end
        avT3 = avT3/n_int_pts;
        dist = norm(x(:,end-1)-x(:,1));
        T(2) = (M(1,end)-M(1,1))/dist;
        T(1) = -(M(2,end)-M(2,1))/dist;
        T(3) = avT3;
        x(:,n_int_pts+1) = 0.5*(x(:,1)+x(:,end-1));
  
end
%
%  =================== Function to extract variables from input file ======
%
function [nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,E,mu,rho,e1dir,inertias,areas,...
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
       rho(beam) = str2double(cellarray{1}{dum+2});
       dum = dum + 4;
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