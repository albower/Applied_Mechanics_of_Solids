function Solution8_16
%
%       This code solves problem 8.16 the text
%       A.F. Bower 'Solved Problems in Mechanics of Solids'  
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%       
%  Define values for parameters below
%
        Length = 10;  % String length
        T = 1;  % Tension
        m = 1;  % Mass per unit length
        L = 119; %    total no. elements,
        Ne = 2; % no. nodes on each element (2 for linear, 3 for quadratic elemnts)
        nnodes = (Ne-1)*L+1;
   
        cL = sqrt(T/m); % Wave speed
        beta1 = 0.5;  %Newmark parameter
        beta2 = 0.;   %Newmark parameter
        dt = 0.5*(Length/L)/cL; % Timestep size (1/2 critical )
        lumpedmass = false;
        nsteps = floor(20/dt);  % Period of fundamental mode is 20
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
       
%      Modify FEM equations to enforce displacement boundary condition at the ends
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
       
       
%      Initial displacement and velocity (fundamental mode)
       u0 = sin(pi*coords/Length)';
       v0 = zeros(nnodes,1);
       
% plotvars stores the time and central displacement of the string
       plotvars = integrate_vibrating_string(nsteps,dt,beta1,beta2,lumpedmass,u0,v0,M,K);
%
% Plot the center displacement as a function of time
%
       figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
       axes3 = axes('Parent',figure1,...
           'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
       hold(axes3,'on');
       set(axes3,'FontName','Times','FontSize',16,...
           'LineWidth',2,'XGrid','on','YGrid','on',...
           'Box','On')
       plot(plotvars(1,:),plotvars(2,:),'LineStyle','-','Color','k','LineWidth',2,'Displayname','n=1')
       
%      Initial displacement (3rd mode)
       u0 = sin(3*pi*coords/Length)';
       plotvars = integrate_vibrating_string(nsteps,dt,beta1,beta2,lumpedmass,u0,v0,M,K);
       plot(plotvars(1,:),plotvars(2,:),'LineStyle','--','Color','k','LineWidth',2,'Displayname','n=3')
       
       xlabel({'Time $t$'},'Interpreter','latex','FontSize',16);
       ylabel({'Displacement $u$'},'Interpreter','latex','FontSize',16);
       title('Central displacement of vibrating string')
       legend1 = legend(axes3,'show');
       set(legend1,...
           'Position',[0.215773822056742 0.22341270423126 0.162797606514686 0.140873010054454],...
           'Interpreter','latex',...
           'FontSize',14,...
           'FontName','Times New Roman');
     
end

    function plotvars = integrate_vibrating_string(nsteps,dt,beta1,beta2,lumpedmass,u0,v0,M,K)
%
%    Newmark time stepping scheme (problem is linear so all matrices are constant)
%
%    Initial accelerations are -M^-1(Ku + F)
%
       nnodes = length(u0);

       un = u0;
       vn = v0;
       an1 = zeros(nnodes,1);

     MK = M + 0.5*beta2*dt*dt*K;
     an = M\(-K*un);

     plotvars = zeros(2,nsteps);    % This stores time and end velocity for plotting

     for n = 1:nsteps
       if (lumpedmass>0 && beta2<1.e-08)  % Explicit dynamics with lumped mass
         rn =  - K*(un +dt*vn + 0.5*(1.-beta2)*dt*dt*an );
         for i = 1 : nnodes
           an1(i) = rn(i)/MK(i,i);
         end
       else         
         an1 = MK\(-K*(un +dt*vn + 0.5*(1.-beta2)*dt*dt*an ) );
       end

       vn1 = (vn + dt*(1.-beta1)*an + dt*beta1*an1);
       un1 = (un + dt*vn + (1.-beta2)*0.5*dt*dt*an + 0.5*beta2*dt*dt*an1 );  
       un = un1;
       vn = vn1;
       an = an1;
       plotvars(1,n) = n*dt;
       plotvars(2,n) = un(nnodes/2);
% 
     end
    end