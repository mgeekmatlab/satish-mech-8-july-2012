% Appendix B: Two-Dimensional EFG Program

% °/o 2D MATLAB EFG CODE - SQUARE DOMAIN

% °/0 John Dolbow 2/17/97

clear;
clc
l=1
% °/0 DEFINE BOUNDARIES/PARAMETERS
Lb = 48;
D=12;

young = 30e6;nu=0.3;
P=1000;

% °/0 PLANE STRESS DMATRIX

Dmat = (young/(1-nu^2))* [1 nu 0;nu 1 0;0 0 (1-nu)/2];

% °/0 SET UP NODAL COORDINATES
ndivl = 10;ndivw = 4;

[x,connl,conn2,numnod] = mesh2(Lb,D,ndivl,ndivw);

% °/o DETERMINE DOMAINS OF INFLUENCE - UNIFORM NODAL SPACING

dmax=3.5;

xspac = Lb/ndivl;

yspac = D/ndivw;

dm(1,1:numnod)=dmax*xspac*ones(1,numnod);
dm(2,1:numnod)=dmax*yspac*ones(1,numnod);

% °/0 SET UP QUADRATURE CELLS
ndivlq = 10;ndivwq = 4;

[xc,conn,numcell,numq] = mesh2(Lb,D,ndivlq,ndivwq);

% °/0 SET UP GAUSS POINTS, WEIGHTS, AND JACOBIAN FOR EACH CELL
quado = 4;

[gauss] = gauss2(quado);
numq2 = numcell*quado^2;
gs = zeros(4,numq2);

[gs] = egauss(xc,conn,gauss,numcell);

% °/0 LOOP OVER GAUSS POINTS TO ASSEMBLE DISCRETE EQUATIONS
k = sparse(numnod*2,numnod*2);
for gg=gs

gpos = gg(l:2);
weight = gg(3);
jac = gg(4);

33

% °/o DETERMINE NODES IN NEIGHBORHOOD OF GAUSS POINT
v = domain(gpos,x,dm,numnod);
L = length(v);
en = zeros(l,2*L);

[phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
Bmat=zeros(3,2*L);
for j=l:L

Bmat(1:3,(2*j-1):2*j) = [dphix(j) 0;0 dphiy(j);dphiy(j) dphix(j)];
end

for i=l:L

en(2*i-l) = 2*v(i)-l;
en(2*i) = 2*v(i);
end

k(en,en) = k(en,en)+sparse((weight*jac)*Bmat*Dmat*Bmat);

end

% °/0 DETERMINE NODES ON BOUNDARY, SET UP BC;S

indl = 0;ind2 = 0;

for j=l:numnod

if(x(l,j)==0.0)

indl=indl+l;

nnu(l,indl) = x(l,j);

nnu(2,indl) = x(2,j);

end

if(x(l,j)==Lb)

ind2=ind2+l;

nt(1,ind2) = x(l,j);

nt(2,ind2) = x(2,j);

end

end

lthu = length(nnu);
ltht = length(nt);
ubar = zeros(lthu*2,1);
f = zeros(numnod*2,1);

% °/0SET UP GAUSS POINTS ALONG TRACTION BOUNDARY
ind=0;

gauss=gauss2(quado) ;

for i=l:(ltht-1)

ycen=(nt(2,i+l)+nt(2,i))/2;

jcob = abs((nt(2,i+1)-nt(2,i))/2);

for j=l:quado

mark(j) = ycen-gauss(l,j)*jcob;

% 34

ind = ind+1;

gst(1,ind)=nt(1,i);

gst(2,ind)=mark(j);

gst(3,ind)=gauss(2,j);

gst(4,ind)=jcob;

end

end

% °/0SET UP GAUSS POINTS ALONG DISPLACEMENT BOUNDARY
gsu=gst;

gsu(l,1:ind)=zeros(l,ind);
qk = zeros(1,2*lthu);

% 0/oINTEGRATE FORCES ALONG BOUNDARY
Imo = (1/12)*D^3;
for gt=gst

gpos = gt (1:2) ;
weight=gt(3);
jac = gt(4);

v = domain(gpos,x5dm5numnod);
L = length(v);
en = zeros(1,2*L);
force=zeros(1,2*L);

[phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
tx=0;

ty= -(P/(2*Imo))*((D^2)/4-gpos(2,ir2) );
for i=l:L

en(2*i-l) = 2*v(i)-l;
en(2*i) = 2*v(i);
force(2*i-l)=tx*phi(i);
force(2*i) = ty*phi(i);
end

f(en) = f(en) + jac*weight*force5;
end

% °/0 INTEGRATE G MATRIX AND Q VECTOR ALONG DISPLACEMENT BOUNDARY
GG = zeros(numnod*2,lthu*2);
indl=0;ind2=0;
for i=l:(lthu-1)
indl=indl+l;
ml = indl; m2 = ml+1;
yl = nnu(2,ml); y2 = nnu(2,m2);
len = yl-y2;

% 35

for j=l:quado

ind2=ind2+l;

gpos = gsu(l:2,ind2);

weight = gsu(3,j);

jac = gsu(4,j);

facll = (-P*nnu(2,ml))/(6*young*Imo);
fac2 = P/(6*young*Imo);
xpl = gpos(l,l);
ypl = gpos(2,1);

uxexl = (6*Lb-3*xpl)*xpl + (2+nu)*(ypl^2 - (D/2)^2);
uxexl = uxexl*facll;

uyexl = 3*nu*ypl^2*(Lb-xpl)+0.25*(4+5*nu)*xpl*D^2+(3*Lb-xpl)*xpl^2;
uyexl = uyexl*fac2;
N1 = (gpos(2,l)-y2)/len; N2 = 1-N1;
qk(2*ml-l) = qk(2*ml-l)-weight*jac*Nl*uxexl;
qk(2*ml) = qk(2*ml) - weight*jac*Nl*uyexl;
qk(2*m2-l) = qk(2*m2-l) -weight*jac*N2*uxexl;
qk(2*m2) = qk(2*m2) - weight*jac*N2*uyexl;
v = domain(gpos,x,dm,numnod);
[phi5dphix5dphiy] = shape(gpos,dmax,x,v,dm);
L = length(v);
for n=l:L

G1 = -weight*jac*phi(n)* [N1 0;0 Nl];
G2 = -weight*jac*phi(n)* [N2 0;0 N2];
cl=2*v(n)-l;c2=2*v(n);c3=2*ml-l;c4=2*ml;
c5=2*m2-l;c6=2*m2;

GG(cl:c2,c3:c4)=GG(cl:c2,c3:c4)+ G1;
GG(cl:c2,c5:c6)=GG(cl:c2,c5:c6)+G2;
end

end

end

% °/0 ENFORCE BC; S USING LAGRANGE MULTIPLIERS
f = [f;zeros(lthu*2,l)];
f(numnod*2+l:numnod*2+2*lthu,1) = -qk;;
m = sparse([k GG;GG; zeros(lthu*2)]);
d=m\f;

u=d(1:2*numnod);
for i=l:numnod
u2(l,i) = u(2*i-l);
u2(2,i) = u(2*i);
end

% 36

% °/o SOLVE FOR OUTPUT VARIABLES - DISPLACEMENTS
displ=zeros(1,2*numnod);
ind=0;
for gg=x

ind = ind+1;

gpos = gg(1:2);

v = domain(gpos,x,dm,numnod);

[phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);

displ(2*ind-l) = phi*u2(l,v);;

displ(2*ind) = phi*u2(2,v);;

end

% °/0 SOLVE FOR STRESSES AT GAUSS POINTS
ind = 0;
enorm=0;
for gg=gs

ind = ind+1;
gpos = gg(l:2);
weight = gg(3);
jac = gg(4);

v = domain(gpos,x,dm,numnod);
L = length(v);
en = zeros(1,2*L);

[phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
Bmat=zeros(3,2*L);
for j=l:L

Bmat(1:3,(2*j-1):2*j) = [dphix(j) 0;0 dphiy(j);dphiy(j) dphix(j)];
end

for i=l:L

en(2*i-l) = 2*v(i)-l;
en(2*i) = 2*v(i);
end

stress(l:3,ind) = Dmat*Bmat*u(en);

stressex(l,ind) = (1/Imo)*P*(Lb-gpos(1,1))*gpos(2,1);
stressex(2,ind) = 0;

stressex(3,ind) = -0.5*(P/Imo)*(0.25*D^2 - gpos(2,1)^2);
err = stress(1:3,ind)-stressex(l:3,ind);
err2 = weight*jac*(0.5*(inv(Dmat)*err)*(err));
enorm = enorm + err2;

end

enorm=sqrt(enorm);
% °/0 MESH GENERATION PROGRAM FOR 2D BEAM IN BENDING


