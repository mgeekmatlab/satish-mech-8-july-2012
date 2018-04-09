% 2D MATLAB EFG CODE - SQUARE DOMAIN
% John Dolbow 12/17/96
clear;
clc
% DEFINE BOUNDARIES/PARAMETERS
Lb = 48;
D = 12;
young = 30e6;
nu=0.3;
P=1000;

% PLANE STRESS DMATRIX
Dmat = (young/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2];

% SET UP NODAL COORDINATES
ndivl = 12;
ndivw = 6;
[x,conn1,conn2,numnod] = mesh2(Lb,D,ndivl,ndivw);

% DETERMINE DOMAINS OF INFLUENCE - UNIFORM NODAL SPACING
dmax=3.5;
xspac = Lb/ndivl;
yspac = D/ndivw;
dm = zeros(2,numnod);
dm(1,1:numnod)=dmax*xspac*ones(1,numnod);
dm(2,1:numnod)=dmax*yspac*ones(1,numnod);

% SET UP QUADRATURE CELLS
ndivlq = 10;
ndivwq = 4;
[xc,conn,numcell,numq] = mesh2(Lb,D,ndivlq,ndivwq);

% SET UP GAUSS POINTS, WEIGHTS, AND JACOBIAN FOR EACH CELL
quado = 4;
[gauss] = pgauss(quado);
numq2 = numcell*quado^2;
gs = zeros(4,numq2);
[gs] = egauss(xc,conn,gauss,numcell);

% LOOP OVER GAUSS POINTS TO ASSEMBLE DISCRETE EQUATIONS
k = sparse(numnod*2,numnod*2);
for gg=gs
   gpos = gg(1:2);
   weight = gg(3);
   jac = gg(4);
   
% DETERMINE NODES IN NEIGHBORHOOD OF GAUSS POINT
   v = domain(gpos,x,dm,numnod);
   L = length(v);
   en = zeros(1,2*L);
   [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
   Bmat=zeros(3,2*L);
   for j=1:L
   Bmat(1:3,(2*j-1):2*j) = [dphix(j) 0;0 dphiy(j);dphiy(j) dphix(j)];
   end
   for i=1:L
   en(2*i-1) = 2*v(i)-1;
   en(2*i) = 2*v(i);
   end
   k(en,en) = k(en,en)+sparse((weight*jac)*Bmat'*Dmat*Bmat);
end

% DETERMINE NODES ON BOUNDARY, SET UP BC'S
ind1 = 0;ind2 = 0;
for j=1:numnod
   if(x(1,j)==0.0)
	ind1=ind1+1;
	nnu(1,ind1) = x(1,j);
	nnu(2,ind1) = x(2,j);
   end
   if(x(1,j)==Lb)
	ind2=ind2+1;
	nt(1,ind2) = x(1,j);
	nt(2,ind2) = x(2,j);
   end
end
lthu = length(nnu);
ltht = length(nt);
ubar = zeros(lthu*2,1); 
f = zeros(numnod*2,1);

%SET UP GAUSS POINTS ALONG TRACTION BOUNDARY
ind=0;
gauss=pgauss(quado);
for i=1:(ltht-1)
   ycen=(nt(2,i+1)+nt(2,i))/2;
   jcob = abs((nt(2,i+1)-nt(2,i))/2);
   for j=1:quado
	mark(j) = ycen-gauss(1,j)*jcob;
	ind = ind+1;
	gst(1,ind)=nt(1,i);
	gst(2,ind)=mark(j);
	gst(3,ind)=gauss(2,j);
	gst(4,ind)=jcob;
   end
end

%SET UP GAUSS POINTS ALONG DISPLACEMENT BOUNDARY
gsu=gst;
gsu(1,1:ind)=zeros(1,ind);
qk = zeros(1,2*lthu);

%INTEGRATE FORCES ALONG BOUNDARY
Imo = (1/12)*D^3;
for gt=gst
   gpos = gt(1:2);
   weight=gt(3);
   jac = gt(4);
   v = domain(gpos,x,dm,numnod);
   L = length(v);
   en = zeros(1,2*L);
   force=zeros(1,2*L);
   [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
   tx=0;
   ty= -(P/(2*Imo))*((D^2)/4-gpos(2,1)^2);
   for i=1:L
      en(2*i-1) = 2*v(i)-1;
      en(2*i) = 2*v(i);
      force(2*i-1)=tx*phi(i);
      force(2*i) = ty*phi(i);
   end
   f(en) = f(en) + jac*weight*force';
end

% INTEGRATE G MATRIX AND Q VECTOR ALONG DISPLACEMENT BOUNDARY
GG = zeros(numnod*2,lthu*2);
ind1=0; ind2=0;
for i=1:(lthu-1)
   ind1=ind1+1;
   m1 = ind1; m2 = m1+1;
   y1 = nnu(2,m1);  y2 = nnu(2,m2);
   len = y1-y2;
   for j=1:quado
   ind2=ind2+1;
   gpos = gsu(1:2,ind2);
   weight = gsu(3,j);
   jac = gsu(4,j);
   fac11 = (-P*nnu(2,m1))/(6*young*Imo);
   fac2 = P/(6*young*Imo);
   xp1 = gpos(1,1);
   yp1 = gpos(2,1);
   uxex1 = (6*Lb-3*xp1)*xp1 + (2+nu)*(yp1^2 - (D/2)^2);
   uxex1 = uxex1*fac11;
   uyex1 = 3*nu*yp1^2*(Lb-xp1)+0.25*(4+5*nu)*xp1*D^2+(3*Lb-xp1)*xp1^2;
   uyex1 = uyex1*fac2;   
   N1 = (gpos(2,1)-y2)/len; N2 = 1-N1;
   qk(2*m1-1) = qk(2*m1-1)-weight*jac*N1*uxex1;
   qk(2*m1) = qk(2*m1) - weight*jac*N1*uyex1;
   qk(2*m2-1) = qk(2*m2-1) -weight*jac*N2*uxex1;
   qk(2*m2) = qk(2*m2) - weight*jac*N2*uyex1;
   v = domain(gpos,x,dm,numnod);
   [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
   L = length(v);
      for n=1:L
      G1 = -weight*jac*phi(n)*[N1 0;0 N1];
      G2 = -weight*jac*phi(n)*[N2 0;0 N2];
      c1=2*v(n)-1;c2=2*v(n);c3=2*m1-1;c4=2*m1;
      c5=2*m2-1;c6=2*m2;
      GG(c1:c2,c3:c4)=GG(c1:c2,c3:c4)+ G1;
      GG(c1:c2,c5:c6)=GG(c1:c2,c5:c6)+G2;
      end
   end
end

% ENFORCE BC'S USING LAGRANGE MULTIPLIERS 
f = [f;zeros(lthu*2,1)];
f(numnod*2+1:numnod*2+2*lthu,1) = -qk';
m = sparse([k GG;GG' zeros(lthu*2)]);
d=m\f;
u=d(1:2*numnod);
for i=1:numnod
   u2(1,i) = u(2*i-1);
   u2(2,i) = u(2*i);
end

% SOLVE FOR OUTPUT VARIABLES - DISPLACEMENTS
displ=zeros(1,2*numnod);
ind=0;
for gg=x
   ind = ind+1;
   gpos = gg(1:2);
   v = domain(gpos,x,dm,numnod);
   [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
   displ(2*ind-1) = phi*u2(1,v)';
   displ(2*ind) = phi*u2(2,v)';
end

%SOLVE FOR STRESSES AT GAUSS POINTS
ind = 0;
enorm = 0;
for gg=gs
   ind = ind+1;
   gpos = gg(1:2);
   weight = gg(3);
   jac = gg(4);   
   v = domain(gpos,x,dm,numnod);
   L = length(v);
   en = zeros(1,2*L);
   [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm);
   Bmat=zeros(3,2*L);
   for j=1:L
   Bmat(1:3,(2*j-1):2*j) = [dphix(j) 0;0 dphiy(j);dphiy(j) dphix(j)];
   end
   for i=1:L
   en(2*i-1) = 2*v(i)-1;
   en(2*i) = 2*v(i);
   end
   stress(1:3,ind) = Dmat*Bmat*u(en);
   stressex(1,ind) = (1/Imo)*P*(Lb-gpos(1,1))*gpos(2,1);
   stressex(2,ind) = 0;
   stressex(3,ind) = -0.5*(P/Imo)*(0.25*D^2 - gpos(2,1)^2);
   err = stress(1:3,ind)-stressex(1:3,ind);
   err2 = weight*jac*(0.5*(inv(Dmat)*err)'*(err));
   enorm = enorm + err2;
end
enorm = sqrt(enorm)