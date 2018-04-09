%-----------------------------------------------------------------------
function [gs] = egauss(xc,conn,gauss,numcell)
% routine to set up gauss points, jacobian, and weights
index=0;
one = ones(1,4);
psiJ = [-1,+1,+1,-1]; etaJ = [-1,-1,+1,+1];
l = size(gauss);
l = l(2);
for e=1:numcell
% DETERMINE NODES IN EACH CELL
   for j = 1:4
     je=conn(j,e);xe(j)=xc(1,je);ye(j)=xc(2,je);
   end
   for i=1:l
     for j=1:l
       	index = index+1;
        eta=gauss(1,i);psi=gauss(1,j);
        N = .25*(one+psi*psiJ).*(one+eta*etaJ);
        NJpsi=.25*psiJ.*(one+eta*etaJ);
        NJeta=.25*etaJ.*(one+psi*psiJ);
        xpsi=NJpsi*xe';ypsi=NJpsi*ye';xeta=NJeta*xe';yeta=NJeta*ye';
        jcob=xpsi*yeta-xeta*ypsi;
        xq = N*xe';yq = N*ye';
        gs(1,index) = xq;
	gs(2,index) = yq;
	gs(3,index) = gauss(2,i)*gauss(2,j);
	gs(4,index) = jcob;
     end
   end
end
