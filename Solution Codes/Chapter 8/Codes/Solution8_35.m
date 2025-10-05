function Solution8_35
%
%          Example FEM code using linear elastic truss elements
%
%       This code solves problem 8.34 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
% ================= Read data from the input file ==================
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
    filename = '../Input Files/Solution8_35.txt';
        while (~isfile(filename))
            f = msgbox({'The input file was not found';...
                        'Please use the browser to select the file' });
            uiwait(f);
            [filename,location] = uigetfile('*.txt');
            filename = strcat(location,filename);
        end

        infile=fopen(filename,'r');


    % Read the input file into a cell array
    cellarray=textscan(infile,'%s');
    [ndof,nnode,coord,nelem,connect,nprops,elprops,nfix,fixnodes,nload,loads] = read_file(cellarray);

    fclose(infile);

    close all

    lumpedmass = false;
    n_modes = 3;

%
%
%===========  Assemble the global stiffness and mass matrix ====================
%
    Stif=zeros(ndof*nnode,ndof*nnode);
    M = Stif;
    for lmn=1:nelem    % Loop over all the members
%
%   Set up the stiffness for the current element
%
        a = connect(lmn,1);
        b = connect(lmn,2);
        [kel,mel]=elstifandmass(coord(a,:),coord(b,:),elprops(lmn,:),lumpedmass);
%
%   Add the current element stiffness to the global stiffness
%
        for i = 1 : 2
            for ii = 1 : ndof
                rw = ndof*(connect(lmn,i)-1)+ii;
                for j = 1 : 2
                    for jj = 1 : ndof
                        cl = ndof*(connect(lmn,j)-1)+jj;
                        Stif(rw,cl) = Stif(rw,cl) + kel(ndof*(i-1)+ii,ndof*(j-1)+jj);
                        M(rw,cl) = M(rw,cl) + mel(ndof*(i-1)+ii,ndof*(j-1)+jj);                    
                    end
                end
            end
        end
    end

%
%  Modify the global stiffness and residual to include constraints
%
    for i=1:nfix
        rw=ndof*(fixnodes(i,1)-1)+fixnodes(i,2);
        for j=1:ndof*nnode
            Stif(rw,j)=0;
            M(rw,j) = 0;
            M(j,rw) = 0;
            Stif(j,rw) = 0;
        end
        M(rw,rw)=1.0;
    end

%
    if (lumpedmass)
     rootM = sqrt(M);
     inverserootM = zeros(nnode);
     for i = 1 : ndof*nnode
       inverserootM(i,i) = 1/rootM(i,i);
     end
    else
     rootM = sqrtm(M);
     inverserootM = inv(rootM);
    end

    H = inverserootM*(Stif*inverserootM);
%      svd is the singular value decomposition we need
    [Q,Lambda,~] = svd(H);
%
%      Lambda has eigenvalues in decreasing order. The last nfix
%      frequencies will be zero and should be ignored.
    disp('Natural Frequencies');
    freqs = sqrt(diag(Lambda));
    disp(freqs(1:6))
    disp(freqs(7:end))

    scale = 5;
% No. zero frequency modes is 3 (2D), 6(3D), or # constrained DOFs
    if (ndof==2)
        nzfmodes = 3;
    else
        nzfmodes = 6;
    end
    if nfix>nzfmodes
        nzfmodes = nfix;
    end
    for m = 1:n_modes
        u = inverserootM*Q(:,ndof*nnode-nzfmodes-m+1);
        for i = 1 : nnode
            x = coord(i,1);
            y = coord(i,2);
            coord(i,1) = x + scale*u(ndof*(i-1)+1);
            coord(i,2) = y + scale*u(ndof*(i-1)+2);
            if (ndof>2)
                z = coord(i,3);
                coord(i,3) = z + scale*u(ndof*(i-1)+3);
            end
        end
        pause(1);
        figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
        axes3 = axes('Parent',figure1,...
            'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
        hold(axes3,'on');
        set(axes3,'FontName','Times','FontSize',16,...
            'LineWidth',2,'XGrid','on','YGrid','on',...
            'Box','On')
        plot_truss(ndof,nnode,coord,nelem,connect)
        titstr = ['Mode: ',num2str(m),' Frequency: ',num2str(freqs(ndof*nnode-nzfmodes-m+1),'%.2f'),' rad/s'];
        title(titstr);
        axis square
        axis equal
        axis off
    end






 end

%
%================= ELEMENT STIFFNESS AND MASS MATRIX  ===================
%
function [kel,mel] = elstifandmass(coorda,coordb,props,lumpedmass)

        len = sqrt(dot(coorda-coordb,coorda-coordb));
        %
        %     Define the element stiffness
        %
        vec = [coorda-coordb,coordb-coorda];
        for i = 1 : length(vec)
            for j = 1 : length(vec)
                kel(i,j) = vec(i)*vec(j);
            end
        end
                
        kel = (props(1)*props(2)/len^3)*kel;
        if (lumpedmass)
            mel = 0.5*eye(4)*props(2)*props(3)*len;
        else
           mel = [2,0,1,0;0,2,0,1;1,0,2,0;0,1,0,2]*props(2)*props(3)*len/(6);
        end
end

%
%  =================== Function to extract variables from input file ======
%
function [ndof,nnode,coord,nelem,connect,nprops,elprops,nfix,fixnodes,nload,loads] = read_file(cellarray) 
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
%
%   Extract no. nodes with prescribed displacements and the prescribed displacements
%
dum = dum + 1;
nfix=str2num(cellarray{1}{dum});
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
loads = zeros(nload,1+ndof);
for i = 1 : nload
    for j=1:1+ndof
        loads(i,j)=str2num(cellarray{1}{dum});
        dum=dum+1;
    end
end

end
function plot_truss(ndof,nnode,coord,nelem,connect)
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
    
    
end
hold off



end
