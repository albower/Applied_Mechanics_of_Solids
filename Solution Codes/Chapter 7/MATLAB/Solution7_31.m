function Solution7_31
%
%       Simple FEA program using constant strain triangle elements
%
%       The code reads an input file that defines material properties,
%       mesh geometry and boundary conditions.  It produces a text
%       file that lists the strains and stresses at the centroid
%       of each element; plots the displaced mesh; and produces a
%       mesh plot with the elements colored by the magnitude of the
%       von Mises stress
%
%       This code solves problem 7.31 from the text
%       A.F. Bower, 'Solved Problems in Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%   It was downloaded from
%   https://github.com/albower/Applied_Mechanics_of_Solids
%
%
%% ====== Read data from the input file and open output file===============
%
% Change the name of the file below to point to your input file
% or leave blank to have user select the file
%
%filename = '';
filename = 'Solution7_31.txt';
    while (~isfile(filename))
        f = msgbox({'The input file was not found';...
                    'Please use the browser to select the file' });
        uiwait(f);
        [filename,location] = uigetfile('*.txt');
        filename = strcat(location,filename);
    end
    
infile=fopen(filename,'r');
% Element data will be written to the file specified below
outfile=fopen('Solution7_31_results.txt','w');
%
%% Extract data from the file
[E,nu,pstrain,nnode,coord,nelem,connect,nfix,fixnodes,ndload,dloads] = read_file(infile);
fclose(infile);

close all
%%
% =============  Define the D matrix  =======================
%
if (pstrain) %Plane strain
    Dmat = [[1-nu,nu,0];[nu,1-nu,0];[0,0,(1-2*nu)/2]]*E/((1+nu)*(1-2*nu));
else  % Plane stress
    Dmat = [[1,nu,0];[nu,1,0];[0,0,(1-nu)/2]]*E/((1-nu^2));
end
%%
%===========  Assemble the global stiffness matrix ====================
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
%%
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
%%
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
%%
% ================== Solve the FEM equations ===================
%

u=Stif\resid;
%%
% =================== Print the results to a file ===============
%
fprintf(outfile,'%s\n','Nodal Displacements:');
fprintf(outfile,'%s\n',' Node    u1       u2');
for i = 1 : nnode
    fprintf(outfile,'%3d %8.4f %8.4f\n',i,u(2*i-1),u(2*i));
end
fprintf(outfile,'\n %s\n','Strains and Stresses');
fprintf(outfile,'%s\n',' Element   e_11      e_22     e_12      s_11       s_22      s_12');
smax = 0;
splot = zeros(nelem,1);
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
    fprintf(outfile,'%5d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n',lmn,strain(1),strain(2),strain(3)/2,stress(1),stress(2),stress(3));
    splot(lmn) = sqrt((stress(1)-stress(2))^2 + 2*stress(3)^2);
    if (splot(lmn)>smax); smax = splot(lmn); end
end
fclose(outfile);
%%
% Plot the displaced mesh
%
for i = 1 : nnode
    x = coord(i,1);
    y = coord(i,2);
    coord(i,1) = x + u(2*(i-1)+1);
    coord(i,2) = y + u(2*(i-1)+2);
end

figure
triplot(connect,coord(:,1),coord(:,2),'g');
hold on;
triplot(connect,coord(:,1),coord(:,2),'r');
axis equal
%%
% Plot mesh with triangles colored according to Von Mises stress
 f2D_3 = [1,2,3];
 figure
 hold on
 for lmn=1:nelem
    for i = 1:3
      x(i,1:2) = coord(connect(lmn,i),1:2);
    end
    ss = splot(lmn)/smax;
    rval = 2*ss-1;
    if (rval<0); rval = 0;end
    if (rval>1); rval = 1; end
    gval = 4*(1-ss)*ss;
    if (gval<0); gval = 0;end
    if (gval>1); gval = 1; end
    bval = 1-2*ss;
    if (bval<0); bval = 0;end
    if (bval>1); bval = 1; end
    color = [rval,gval,bval];
    patch('Vertices',x,'Faces',f2D_3,'FaceColor',color,'EdgeColor','r');
 end
axis equal

end

%%
%================= ELEMENT STIFFNESS MATRIX ===================
%
function kel = elstif(xa,ya,xb,yb,xc,yc,Dmat)
        %
        %     Define the B matrix
        %
        nax = -(yc-yb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nay =  (xc-xb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nbx = -(ya-yc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        nby =  (xa-xc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        ncx = -(yb-ya)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        ncy =  (xb-xa)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        area = (1/2)*abs( (xb-xa)*(yc-ya) - (xc-xa)*(yb-ya) );
        Bmat = [[nax,  0,nbx,  0,ncx,  0]; [0,nay,  0,nby,  0,ncy];[nay,nax,nby,nbx,ncy,ncx]];
        %
        %     Define the element stiffness
        %
        kel = area*transpose(Bmat)*Dmat*Bmat;
end

%%
%====================== ELEMENT FORCE VECTOR ==============
%
function rel = elresid (xa,ya,xb,yb,tx,ty)
        length = sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb));
        rel = [tx,ty,tx,ty]*length/2;
end
%%
%
%    Function to calculate the element strains
%
    function strain = elstrain (xa,ya,xb,yb,xc,yc,uax,uay,ubx,uby,ucx,ucy)
        %
        %   B matrix
        %
        nax = -(yc-yb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nay =  (xc-xb)/( (ya-yb)*(xc-xb) - (xa-xb)*(yc-yb) );
        nbx = -(ya-yc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        nby =  (xa-xc)/( (yb-yc)*(xa-xc) - (xb-xc)*(ya-yc) );
        ncx = -(yb-ya)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        ncy =  (xb-xa)/( (yc-ya)*(xb-xa) - (xc-xa)*(yb-ya) );
        Bmat = [[nax,  0,nbx,  0,ncx,  0]; [ 0,nay,  0,nby,  0,ncy];[nay,nax,nby,nbx,ncy,ncx]];
        %
        %     Element displacement vector
        %
        uel = [uax;uay;ubx;uby;ucx;ucy];
        %
        %     Element strains
        %
        strain = Bmat*uel;

    end
%%
%  =================== Function to extract variables from input file ======
%
function [E,nu,pstrain,nnode,coord,nelem,connect,nfix,fixnodes,ndload,dloads] = read_file(infile) 
cellarray=textscan(infile,'%s');
%
%  Extract the material properties
%
E=str2double(cellarray{1}{3});
nu=str2double(cellarray{1}{5});
pstrain = false;   % Default element is plane stress
if(strcmp(cellarray{1}{7},'Plane_Strain'))
    pstrain = true;
end
%
%  Extract no. nodes and nodal coordinates
%
nnode=str2double(cellarray{1}{9});
dum=11;
coord=zeros(nnode,2);
for i=1:nnode
    coord(i,1) = str2double(cellarray{1}{dum});
    dum=dum+1;
    coord(i,2) = str2double(cellarray{1}{dum});
    dum=dum+1;
end
%
%   Extract no. elements and connectivity
%
dum=dum + 1;
nelem=str2double(cellarray{1}{dum});
connect = zeros(nelem,3);
dum = dum + 2;
for i = 1 : nelem
    for j = 1 : 3
        connect(i,j) = str2double(cellarray{1}{dum});
        dum=dum+1;
    end
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
ndload=str2double(cellarray{1}{dum});
dum=dum + 4;
dloads = zeros(ndload,4);
for i = 1 : ndload
    for j=1:4
        dloads(i,j)=str2double(cellarray{1}{dum});
        dum=dum+1;
    end
end

end