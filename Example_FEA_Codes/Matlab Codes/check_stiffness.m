function check_stiffness(outfile,dt,lmn,ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,stress0,eplas0,dofs)
%
%   This is an auxiliary file that can be used to check that the
%   stiffness matrix and residual force vector are consistent.
%   The code compares the coded stiffness matrix with an approximate
%   numerical derivative of the residual force vector.
%   The two should agree to 3-4 decimal places if the stiffness and
%   residual force are coded correctly
%
%   Assemble the global stiffness matrix
%

   lmncoord = zeros(ncoord,maxnodes);
   lmndof = zeros(ndof,maxnodes);
   kel_num = zeros(ndof*maxnodes,ndof*maxnodes);
   %
%   Loop over all the elements
%
%
%   Extract coords of nodes, DOF for the current element
%
      for a = 1:nelnodes(lmn)
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
          lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
      end

    n = nelnodes(lmn);
    ident = elident(lmn);
    [rel0,kel,~,~] = elforcnstif(dt,ncoord,ndof,n,ident,lmncoord,materialprops,stress0(:,:,:,lmn),eplas0(:,lmn),lmndof);
    for a = 1:nelnodes(lmn)
        for i = 1:ndof
            lmndof1 = lmndof;
            lmndof1(i,a) = lmndof(i,a) + 1.e-07;
          [rel1,~,~,~] = elforcnstif(dt,ncoord,ndof,n,ident,lmncoord,materialprops,stress0(:,:,:,lmn),eplas0(:,lmn),lmndof1);
          kel_num(:,ndof*(a-1)+i) = (rel1(:)-rel0(:))/1.e-07;
        end
    end

    fprintf(outfile,'Coded stiffness \n ');
    for a = 1:nelnodes(lmn)
      for i = 1:ndof
%         fprintf(outfile,'\n\n Node %d DOF %d \n',a,i); 
        rw = ndof*(a-1)+i;
        for b = 1:nelnodes(lmn)
          for k = 1:ndof
            cl = ndof*(b-1)+k;
              fprintf(outfile,'%d ',kel(rw,cl));
          end
        end
        fprintf(outfile,'\n');
      end
    end
    
%
    fprintf(outfile,'\n \n Numerical stiffness \n');
    for a = 1:nelnodes(lmn)
      for i = 1:ndof
%         fprintf(outfile,'\n\n Node %d DOF %d \n',a,i); 
        rw = ndof*(a-1)+i;
        for b = 1:nelnodes(lmn)
          for k = 1:ndof
            cl = ndof*(b-1)+k;
              fprintf(outfile,'%d ',kel_num(rw,cl));
          end
        end
        fprintf(outfile,'\n');
      end
    end    
    

end