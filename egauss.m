function [gs] = egauss(xc,conn,gauss,numcell)

% °/0 routine to set up gauss points, jacobian, and weights

index=0;

one = ones(1,4);

psiJ = [-1, + 1, + 1,-1] ; etaJ = [-1 ,-1 , + 1 ,+1] ;
1 = size(gauss);
1 = 1(2);
for e=l:numcell
% °/0 DETERMINE NODES IN EACH CELL
for j = 1:4

j e=conn(j,e);xe(j)=xc(l,je);ye(j)=xc(2,je);
end
for i=l:1

for j = 1:1

index = index+1;

eta=gauss(1,i);psi=gauss(l,j);

N = .25*(one+psi*psiJ).*(one+eta*etaJ);

38

NJpsi=.25*psiJ.*(one+eta*etaJ);
NJeta=.25*etaJ.*(one+psi*psiJ);

xpsi=NJpsi*xe;;ypsi=NJpsi*ye;;xeta=NJeta*xe;;yeta=NJeta*ye;;
j cob=xpsi*yeta-xeta*ypsi;
xq = N*xe;;yq = N*ye;;
gs(l,index) = xq;
gs(2,index) = yq;

gs(3,index) = gauss(2,i)*gauss(2,j);
gs(4,index) = jcob;
end
end

end
