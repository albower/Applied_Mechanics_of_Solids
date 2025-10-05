function Solution8_19
%
%       This code solves problem 8.19 from the text
%       A.F. Bower 'Solved Problems in Applied Mechanics of Solids' 
%       CRC press, Baton Rouge, 2026
%
%       It was downloaded from
%       https://github.com/albower/Applied_Mechanics_of_Solids
%
  tol = 1.e-05;   % Tolerance for Newton solver
  err = 1;        % Current error
  nit = 0;        % Iteration number
  
  x = 0; y = 0;   % Initial guess
  c = 1;          % x^2+y^2+1
  w = [x;y];
  r = [-5;3]; 
  while (err>tol)
      K = [3*c+2*x^2, 2*x*y; 2*x*y, 3*c + 2*y^2]/(3*c^(2/3));
      w = w - K\r;
      x = w(1); y = w(2);
      c = (x^2+y^2+1);
      r = [x*c^(1/3)-5;y*c^(1/3)+3];
      err = norm(r);
      nit = nit + 1;
      disp(' ');
      disp(['Iteration: ',num2str(nit)]);
      disp(['   Current solution ',num2str(w')]);
      disp(['   Current residual ',num2str(r')]);
      disp(['   Error ',num2str(err)]);
  end
  disp(' ');
  disp(['Solution: x = ',num2str(x),' y = ',num2str(y)]) 




end
