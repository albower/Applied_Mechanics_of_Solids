function Solution8_17
%
%       This code solves problem 8.17 the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
       n_modes = 3; % No. of non-rigid body modes to be plotted
%       
%  Define values for parameters below
%
       Length = 10;  % String length
       T = 1;  % Tension
       m = 1;  % Mass per unit length
       L = 119; %    total no. elements,
       Ne = 2; % no. nodes on each element (2 for linear, 3 for quadratic elemnts)
       nnodes = (Ne-1)*L+1;
       
       lumpedmass = false;
       
%
%     Set up some data structures storing the mesh
%
       coords = zeros(1,nnodes);
       for i= 1 : nnodes
           coords(i) = Length*(i-1)/(nnodes-1);
       end
       %
       connect = zeros(Ne,L);
       for lmn=1:L
           connect(1,lmn) = lmn;
           connect(2,lmn) = lmn+1;
       end
%
%
%     Assemble the global stiffness and force vector
%
       M = zeros(nnodes,nnodes);
       K = zeros(nnodes,nnodes);
       %
       for lmn = 1 : L
%
%      Element stiffness and mass
%
           lel = abs(coords(connect(2,lmn))-coords(connect(1,lmn)));
           if (lumpedmass)
               mel = m*lel/2*[1,0;0,1];
           else
               mel = m*lel/6*[2,1;1,2];
           end
           kel = T/lel*[1,-1;-1,1];
%
%       Add the stiffness and residual from the current element into global matrices
%
           for a = 1 : Ne
               rw = connect(a,lmn);
               for b = 1 : Ne
                   cl = connect(b,lmn);
                   M(rw,cl) = M(rw,cl) + mel(a,b);
                   K(rw,cl) = K(rw,cl) + kel(a,b);
               end
           end
       end

%      Modify FEM equations to enforce displacement boundary condition at
%      the ends
%
       M(1,1:nnodes) = 0.;
       M(1:nnodes,1) = 0;
       K(1,1:nnodes) = 0.;
       K(1:nnodes,1) = 0;
       M(nnodes,1:nnodes) = 0.;
       K(nnodes,1:nnodes) = 0.;
       M(1:nnodes,nnodes) = 0;
       K(1:nnodes,nnodes) = 0;
       M(1,1) = 1.;
       M(nnodes,nnodes) = 1.0;
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
%       
%      svd is the singular value decomposition we need
       [Q,Lambda,~] = svd(H);
%
%      Lambda has eigenvalues in decreasing order. The last two are rigid
%      body mode
%
       disp('Lowest 10 Natural Frequencies:');
       freqs = sqrt(diag(Lambda));
       disp(freqs(end-2:-1:end-11))
       
       close all
       
       figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
       axes3 = axes('Parent',figure1,...
           'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
       hold(axes3,'on');
       set(axes3,'FontName','Times','FontSize',16,...
           'LineWidth',2,'XGrid','on','YGrid','on',...
           'Box','On')
       
       styles = ["-","--","-."];
       
       for i = 1:n_modes
           u = inverserootM*Q(:,nnodes-1-i);
           plot(coords,u,'LineStyle',styles(i),'Color','k','LineWidth',2,'Displayname',['$\omega= $',num2str(freqs(nnodes-1-i))])
       end
       xlabel({'Time $t$'},'Interpreter','latex','FontSize',16);
       ylabel({'Displacement $u$'},'Interpreter','latex','FontSize',16);
       title('Mode shapes for a vibrating string')
       legend1 = legend(axes3,'show');
       set(legend1,...
           'Position',[0.174393114889557 0.207142857142859 0.264892599396158 0.206190474146893],...
           'Interpreter','latex',...
           'FontSize',14,...
           'FontName','Times New Roman');

  

end
