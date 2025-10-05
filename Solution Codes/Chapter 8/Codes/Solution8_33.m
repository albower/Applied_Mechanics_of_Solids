function Solution8_33
%
%          Example FEM code using linear elastic truss elements
%
%       This code solves problem 8.33 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
%%
% ================= Read data from the input file ==================
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
filename = '../Input Files/Solution8_33.txt';
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
[ndof,nnode,coord,nelem,connect,nprops,elprops,nfix,fixnodes,nload,loads] = read_file(cellarray);

fclose(infile);

close all

%
%%
%===========  Assemble the global stiffness matrix ====================
%
Stif=zeros(ndof*nnode,ndof*nnode);
for lmn=1:nelem    % Loop over all the members
    %
    %   Set up the stiffness for the current element
    %
    a = connect(lmn,1);
    b = connect(lmn,2);
    k=elstif(coord(a,:),coord(b,:),elprops(lmn,:));
    %
    %   Add the current element stiffness to the global stiffness
    %
    for i = 1 : 2
        for ii = 1 : ndof
            for j = 1 : 2
                for jj = 1 : ndof
                    rw = ndof*(connect(lmn,i)-1)+ii;
                    cl = ndof*(connect(lmn,j)-1)+jj;
                    Stif(rw,cl) = Stif(rw,cl) + k(ndof*(i-1)+ii,ndof*(j-1)+jj);
                end
            end
        end
    end
end
%
% ==================== Assemble global force vector ============
%
resid=zeros(ndof*nnode);
for i=1:nload   % Loop over joints with external forces
    node=loads(i,1);
    resid(ndof*(node-1)+1)=resid(ndof*(node-1)+1) + loads(i,2);
    resid(ndof*(node-1)+2)=resid(ndof*(node-1)+2)+ loads(i,3);
    if (ndof>2)
      resid(ndof*(node-1)+3)=resid(ndof*(node-1)+3)+ loads(i,4);
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
%%
% ================== Solve the FEM equations ===================
%
u=Stif\resid;

%%
% =================== Print the results to a file ===============
%
fprintf(outfile,'%s\n','Nodal Displacements:');
if (ndof==2)
    fprintf(outfile,'%s\n',' Node    u1       u2');
    for i = 1 : nnode
        fprintf(outfile,'%3d %8.4g %8.4g\n',i,u(2*i-1),u(2*i));
    end
    fprintf(outfile,'\n\n%s\n','Member forces:');
    fprintf(outfile,'%s\n',' Element    T ');
    for lmn = 1 : nelem
        a = connect(lmn,1);
        b = connect(lmn,2);
        nvec = coord(b,:)-coord(a,:);
        len = norm(nvec);
        e = dot(nvec,u(2*b-1:2*b)-u(2*a-1:2*a))/len^2;
        T = elprops(lmn,1)*elprops(lmn,2)*e;
        fprintf(outfile,'%3d %8.4g \n',lmn,T);
    end
elseif (ndof==3)
    fprintf(outfile,'%s\n',' Node    u1       u2       u3');
    for i = 1 : nnode
        fprintf(outfile,'%3d %8.4g %8.4g %8.4g\n',i,u(3*i-2),u(3*i-1),u(3*i));
    end
    fprintf(outfile,'\n\n%s\n','Member forces:');
    fprintf(outfile,'%s\n',' Element    T ');
    for lmn = 1 : nelem
        a = connect(lmn,1);
        b = connect(lmn,2);
        nvec = coord(b,:)-coord(a,:);
        len = norm(nvec);
        e = dot(nvec,u(3*b-2:3*b)-u(3*a-2:3*a))/len^2;
        T = elprops(lmn,1)*elprops(lmn,2)*e;
        fprintf(outfile,'%3d %8.4g \n',lmn,T);
    end
    
    
end
fclose(outfile);
figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
axes2 = gca;
set(axes2,'FontName','Times','FontSize',16,...
    'LineWidth',2,'XGrid','off','YGrid','off',...
    'Box','On');
hold on
plot_truss(ndof,nnode,coord,nelem,connect,'--')
title('Static deflection of truss')
%
% Plot the displaced mesh
%
scale = 100;
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

plot_truss(ndof,nnode,coord,nelem,connect,'-')
annotation(figure1,'textbox',...
    [0.167857142857143 0.133333333333333 0.535714285714286 0.0751904761904761],...
    'String',{'Deformation Scale Factor: 100'},'FontName','Times','FontSize',16);
axis square
axis equal
 end

%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function kel = elstif(coorda,coordb,props)

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
function plot_truss(ndof,nnode,coord,nelem,connect,ls)
hold on
if (ndof==3)

    for i = 1 : nelem
        xvals = [coord(connect(i,1),1),coord(connect(i,2),1)];
        yvals = [coord(connect(i,1),2),coord(connect(i,2),2)];
        zvals = [coord(connect(i,1),3),coord(connect(i,2),3)];
        plot3(xvals,yvals,zvals,'LineWidth',2,'Color',[0,0,0],'LineStyle',ls);
    end

    scatter3(coord(:,1),coord(:,2),coord(:,3),'MarkerFaceColor',[1,0,0]);
        
    
else
    for i = 1 : nelem
        xvals = [coord(connect(i,1),1),coord(connect(i,2),1)];
        yvals = [coord(connect(i,1),2),coord(connect(i,2),2)];
        plot(xvals,yvals,'LineWidth',2,'Color',[0,0,0],'LineStyle',ls);
    end
    
    scatter(coord(:,1),coord(:,2),'MarkerFaceColor',[1,0,0]);
    
    
end
hold off



end
