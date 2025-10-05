function Solution8_20
%
%       This code solves problem 8.19 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids' 
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
   Nel = 10;  % No. elements
   L = 10;    % Bar length
   coords = 0:L/Nel:L;
   tstar = 0; % Traction on left end (use 10 for case (i))
   b = 1;     % Body force (use 0 for case(i))
   
%   materialprops = [5,4,0.01];  %[s0, n, e0]   (case(i))
   materialprops = [5,5,0.01];  %[s0, n, e0]    (case(ii))
   
   F = b*L/Nel*ones(Nel+1,1);   % The external force vector can be calculated explicitly
   F(1) = F(1) + tstar;
   
   nsteps = 15;
   tol = 0.0001;
   maxit = 30;
   relax = 1.;
   
   w = zeros(Nel+1,1);
   for step = 1 : 15
       
       loadfactor = step/nsteps;
       
       err1 = 1.;
       nit = 0;
       
       fprintf(1,'\n Step %f Load %f\n',step,loadfactor);
       
       while ((err1>tol) && (nit<maxit))          % Newton Raphson loop
           nit = nit + 1;
           
           [R,K] = globalmatrices(coords,materialprops,w);
           
           rhs = loadfactor*F - R;
           
           %          Fix constrained node.
           K(Nel+1,:) = 0.0;
           K(Nel+1,Nel+1) = 1.0;
           rhs(Nel+1) = 0;
           %
           %        Solve for the correction
           %
           dw = K\rhs;
           %
           %        Check convergence
           %
           w = w + relax*dw;
           wnorm = dot(w,w);
           err1 = dot(dw,dw);
           err2 = dot(b,b);
           err1 = sqrt(err1/wnorm);
           err2 = sqrt(err2)/(Nel+1);
           fprintf(1,'Iteration number %d Correction %f Residual %f tolerance %f\n',nit,err1,err2,tol);
           
       end
       
   end
   
   if (b>0)
       close all
       
       figure1 = figure('NumberTitle','off','Name','Figure','Color',[1 1 1]);
       axes3 = axes('Parent',figure1,...
           'Position',[0.175 0.207142857142857 0.73 0.717857142857145]);
       hold(axes3,'on');
       set(axes3,'FontName','Times','FontSize',16,...
           'LineWidth',2,'XGrid','on','YGrid','on',...
           'Box','On')
       
       plot(coords,w,'LineStyle','-','Color','k','LineWidth',2)
       xlabel({'Position $x_1$'},'Interpreter','latex','FontSize',16);
       ylabel({'Displacement $u_1$'},'Interpreter','latex','FontSize',16);
       title('Displacement of 1D bar')
   else
       disp([' Average strain: ',num2str(-w(1)/L)])
       
   end
   
end
function [R,K] = globalmatrices(coords,materialprops,w)
    s0 = materialprops(1); n = materialprops(2); e0 = materialprops(3);
    Nel = length(coords)-1;
    R = 0*coords';
    K = zeros(Nel+1);
    for el = 1:Nel
        Lel = coords(el+1)-coords(el); % Element length
        e = (w(el+1)-w(el))/Lel;        % Strain
        signe = sign(e);
        if (e==0)  signe = 1; end
        e = abs(e);
        if (e>e0)
            s = signe*s0*(e/e0)^(1/n)
            dsde = signe*s/(n*e);
        else
            f = sqrt( (1+n^2)/(n-1)^2 - (n/(n-1) - e/e0)^2 );
            s = signe*s0*(f-1/(n-1));
            dsde = s0*(n/(n-1) - e/e0)/(f*e0);
        end
        R(el) = R(el) - s;
        R(el+1) = R(el+1) + s;
        K(el,el) = K(el,el) + dsde/Lel;
        K(el,el+1) = K(el,el+1) - dsde/Lel;
        K(el+1,el) = K(el+1,el) - dsde/Lel;
        K(el+1,el+1) = K(el+1,el+1) + dsde/Lel;
    end

end
