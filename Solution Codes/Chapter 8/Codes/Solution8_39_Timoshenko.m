function Solution8_39_Timoshenko
%
%          Example FEM code using elastic Timoshenko beam elements
%
%       This code solves problem 8.39 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids


% ================= Read data from the input file ==================
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
    filename = '../Input Files/Solution8_39_Timoshenko.txt';
    while (~isfile(filename))
        f = msgbox({'The input file was not found';...
            'Please use the browser to select the file' });
        uiwait(f);
        [filename,location] = uigetfile('*.txt');
        filename = strcat(location,filename);
    end

    infile=fopen(filename,'r');
    
    % lumpedmass = 0 : full mass matrix
    %              1 : row sum method
    %              2 : diagonal scaling (not recommended!)
    lumpedmass = 0;    
    n_modes=5;
% Read the input file into a cell array
    cellarray=textscan(infile,'%s');
    [nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,E,mu,rho,e1dir,inertias,areas,betas,...
                            nfix,fixnodes,nload,loads] = read_file(cellarray) ;
    fclose(infile);

    sfactor = E(1)*inertias(1,1)/(betas(1)*mu(1)*areas(1)*coord(end,1));
    
    disp(['Shear factor for Timoshenko beam: ',num2str(sfactor)]);
    
    close all
%
%  Assemble the stiffness matrix
%
    Stif=zeros(6*n_total_nodes,6*n_total_nodes);
    resid = zeros(6*n_total_nodes,1);
    Mass = Stif;
    for beam = 1:nbeams
        if (beam<nbeams)
            end_element = first_element(beam+1);
        else
            end_element = n_total_elements;
        end
        for lmn=first_element:end_element    % Loop over all the members
            %
            %   Set up the stiffness for the current element
            %
            a = connect(lmn,1);
            b = connect(lmn,2);
            [k,r]=elstif(coord(a,:),coord(b,:),E(beam),mu(beam),e1dir(beam,:),inertias(beam,:),areas(beam),betas(beam));

            mel = elmass(coord(a,:),coord(b,:),rho(beam),e1dir(beam,:),inertias(beam,:),areas(beam),lumpedmass);

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

                           Mass(rw,cl) = Mass(rw,cl) + mel(6*(i-1)+ii,6*(j-1)+jj);

                        end
                    end
                end
            end
        end
    end

   %   Modify the global stiffness and mass to include constraints
    %
    for i=1:nfix
        rw=6*(fixnodes(i,1)-1)+fixnodes(i,2);
        Stif(rw,:)=0;
        Stif(:,rw)=0;
        Mass(:,rw) = 0;
        Mass(rw,:) = 0;
        Mass(rw,rw)=1.0;
    end

    if (lumpedmass>0)
     rootM = sqrt(Mass);
     inverserootM = zeros(6*n_total_nodes);
     for i = 1 : 6*n_total_nodes
       inverserootM(i,i) = 1/rootM(i,i);
     end
    else
     rootM = sqrtm(Mass);
     inverserootM = inv(rootM);
    end

    H = inverserootM*(Stif*inverserootM);
    %      svd is the singular value decomposition we need
    [Q,Lambda,~] = svd(H);
    
    %
    %      Lambda has eigenvalues in decreasing order. The last nfix
    %      frequencies will be zero and should be ignored.
    disp('Normalized Natural Frequencies');
    freqs = sqrt(diag(Lambda));
    % No. zero frequency modes is 6(3D), or # constrained DOFs
    nzfmodes = 6;    
    if nfix>nzfmodes
        nzfmodes = nfix;
    end
    L = coord(end,1);    
    for m = 1:length(freqs)-nzfmodes
         u = inverserootM*Q(:,6*n_total_nodes+1-nzfmodes-m);
         u1 = u(1:6:end-5);
         u2 = u(2:6:end-4);
         u3 = u(3:6:end-3);
         q1 = u(4:6:end-2);
         if (norm(u2)+norm(u3)>0.0001)
             wnnorm = freqs(6*n_total_nodes+1-nzfmodes-m)*sqrt(rho(beam)*L^4*areas(beam)/(E(beam)*inertias(beam,1)));
             disp([num2str(wnnorm),'(transverse)'])
         elseif (norm(u1)>0.0001)
             wnnorm = freqs(6*n_total_nodes+1-nzfmodes-m)*sqrt(rho(beam)*L^2/E(beam));
             disp([num2str(wnnorm),'rad/s (axial)'])   
         elseif (norm(q1)>0.0001)
             wnnorm = freqs(6*n_total_nodes+1-nzfmodes-m)*sqrt(rho(beam)*L^2/mu(beam));                 
             disp([num2str(wnnorm),'rad/s (torsional)'])
         end
    end
 %
 % Plot modes with distinct frequencies
    plottedfreq = 0;
    m_plotted = 0;
    m = 0;
    while m_plotted < n_modes
        m = m + 1;
        if (freqs(6*n_total_nodes+1-nzfmodes-m)-plottedfreq<1.e-05) continue; end
        m_plotted = m_plotted + 1;
        u = inverserootM*Q(:,6*n_total_nodes+1-nzfmodes-m);
        scale = 15;
        displaced_coords(:,1) = coord(:,1) + scale*u(1:6:end-5);
        displaced_coords(:,2) = coord(:,2) + scale*u(2:6:end-4);
        displaced_coords(:,3) = coord(:,3) + scale*u(3:6:end-3);
        u1 = u(1:6:end-5);
        u2 = u(2:6:end-4);
        u3 = u(3:6:end-3);
        q1 = u(4:6:end-2);
        if (norm(u2)+norm(u3)>0.0001)
            wnnorm = freqs(6*n_total_nodes+1-nzfmodes-m)*sqrt(rho(beam)*L^4*areas(beam)/(E(beam)*inertias(beam,1)));
            modestr = ' (Transverse) ';
        elseif (norm(u1)>0.0001)
            wnnorm = freqs(6*n_total_nodes+1-nzfmodes-m)*sqrt(rho(beam)*L^2/E(beam));
            modestr = ' (Axial) ';
        elseif (norm(q1)>0.0001)
            wnnorm = freqs(6*n_total_nodes+1-nzfmodes-m)*sqrt(rho(beam)*L^2/mu(beam));
            modestr = ' (Torsional) ';
        end
        
        figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
        axes3 = axes('Parent',figure1,...
            'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
        hold(axes3,'on');
        set(axes3,'FontName','Times','FontSize',16,...
            'LineWidth',2,'XGrid','on','YGrid','on',...
            'Box','On')
        plot_beam(nbeams,n_total_nodes,first_node,displaced_coords,n_total_elements,first_element,connect,[0,0,0]);
        axis square
        axis equal
        axis off
        title([' Normalized Frequency: ',num2str(wnnorm),' rad/s',modestr])
        plottedfreq = freqs(6*n_total_nodes+1-nzfmodes-m);
    end
    
end

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function [kel,rel] = elstif(coorda,coordb,E,mu,e1dir,inertias,area,beta)

        kel = zeros(12,12);
        rel = zeros(12,1);
%  One point Gaussian integration points and weights        
        n_int_pts = 1;
        w(1) = 2.;
        xi(1) = 0.;
%         w(1:2) = 1;
%         xi(1) = -1./sqrt(3);
%         xi(2) = 1./sqrt(3);
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
            
            N1 = (1-xi(i))/2; N2 = (1+xi(i))/2;
           
            B = [-eye(6)/len,eye(6)/len];
            B(1,5) = -N1;   B(1,11) = -N2;
            B(2,4) = N1; B(2,10) = N2;
            
            D = eye(6);
            D(1,1) = beta*area*mu;
            D(2,2) = beta*area*mu;
            D(3,3) = E*area;
            D(4,4) = E*inertias(1);
            D(4,5) = -E*inertias(2);
            D(5,4) = -E*inertias(2);
            D(5,5) = E*inertias(3);
            D(6,6) = mu*inertias(4);
            
             kel = kel + transpose(R)*transpose(B)*D*B*R*w(i)*len/2;
        end

end
function mel = elmass(coorda,coordb,rho,e1dir,inertias,area,lumpedmass)

        mel = zeros(12,12);
%  One point Gaussian integration points and weights        
        n_int_pts = 1;
        w(1) = 2.;
        xi(1) = 0.;
%         w(1:2) = 1;
%         xi(1) = -1./sqrt(3);
%         xi(2) = 1./sqrt(3);
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
            
            N1 = (1-xi(i))/2; N2 = (1+xi(i))/2;
           
            N = [N1*eye(6),N2*eye(6)];
            
             mu = rho*area*eye(6);
             mu(4,4) = rho*inertias(1); 
             mu(5,5) = rho*inertias(3);
             mu(6,6) = rho*inertias(4);
             mu(4,5) = -rho*inertias(2); mu(5,4) = mu(4,5);
             mel = mel + transpose(R)*transpose(N)*mu*N*R*w(i)*len/2;
        end
        
        
        if (lumpedmass==1)  % Row sum mass lumping
            mel_new = eye(12).*sum(mel,2);
            mel = mel_new;
        elseif (lumpedmass==2)  % Scaled diag mass lumping
            totalmass = rho*area*len;
            totaldiag = sum(diag(mel))/3;  % Mass is assigned to each DOF
            factor = totalmass/totaldiag;
            mel_new = zeros(12);
            for i = 1:12
                mel_new(i,i) = mel(i,i)*factor;
            end
            mel = mel_new;
        end

end
function [T,M,x] = internalforces(coorda,coordb,dofa,dofb,E,mu,e1dir,inertias,area,beta)
% Calculate internal foreces and coordinates of integration points
%  One point Gaussian integration points and weights        
        n_int_pts = 1;
        xi(1) = 0.;

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
        T = zeros(3,n_int_pts);
        M = T;
        x = T;
        for i = 1:n_int_pts
            
            N1 = (1-xi(i))/2; N2 = (1+xi(i))/2;
            x(:,i) = coorda*N1 + coordb*N2;
            B = [-eye(6)/len,eye(6)/len];
            B(1,5) = -N1;   B(1,11) = -N2;
            B(2,4) = N1; B(2,10) = N2;
            
            D = eye(6);
            D(1,1) = beta*area*mu;
            D(2,2) = beta*area*mu;
            D(3,3) = E*area;
            D(4,4) = E*inertias(1);
            D(4,5) = -E*inertias(2);
            D(5,4) = -E*inertias(2);
            D(5,5) = E*inertias(3);
            D(6,6) = mu*inertias(4);
            intFs = D*B*R*[dofa;dofb];
            T(1:3,i) = intFs(1:3);
            M(1:3,i) = intFs(4:6);
        end


end
%
%  =================== Function to extract variables from input file ======
%
function [nbeams,n_total_nodes,first_node,coord,n_total_elements,first_element,connect,E,mu,rho,e1dir,inertias,areas,betas,...
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
        dum = dum + 2;
        betas(beam) = str2double(cellarray{1}{dum});
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
