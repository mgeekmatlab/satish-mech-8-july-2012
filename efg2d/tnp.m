% Appendix A: One-Dimensional EFG Program

% ONE-DIMENSIONAL EFG PROGRAM

% SET UP NODAL COORDINATES ALONG BAR, DETERMINE NUMBER OF CELLS
x = [0.0: 0.1:1.0] ;
nnodes = length(x);
ncells = nnodes-1;

% °/„ SET PARAMETERS FOR WEIGHT FUNCTION, MATERIAL PROPERITES
dmax = 2.0; %RATIO OF DMI TO CI
E=1.0;

% °/„ DETERMINE DMI FOR EACH NODE

dm = dmax*(x(2)-x(l))*ones(l,nnodes);

% °/,SET UP GAUSS POINTS, WEIGHTS, AND JACOBIAN FOR EACH CELL
gg = zeros(1,ncells);
jac = (x(2)-x(l))/2;
weight = 2;

gg = -.05:.1:0.95; gg(l) = 0.0;

% INITIALIZE MATRICES
k = zeros(nnodes);
f = zeros(nnodes,1);
GG = zeros(nnodes,1);

% LOOP OVER GAUSS POINTS
for j = 1:length(gg)

xg = gg(j);

% °/„ DETERMINE DISTANCE BETWEEN NODES AND GAUSS POINT
dif = xg*ones(l,nnodes)-x;

% °/„ SET UP WEIGHTS W AND DW FOR EACH NODE
clear w dw
for i=l:nnodes
drdx = sign(dif(i))/dm(i);
r = abs(dif(i))/dm(i);
if r<=0.5

w(i) = (2/3) - 4*r*r + 4*r-3;
dw(i) = (-8*r + 12*r-2)*drdx;
elseif r<=1.0

w(i) = (4/3)-4*r+4*r*r -(4/3)*r-3;

% 31

dw(i) = (-4 + 8*r-4*r^2)*drdx;
elseif r>l.0
w(i) = 0.0;
dw(i) = 0.0;

end
end

% °/0SET UP SHAPE FUNCTIONS AND DERIVATIVES
won = ones(1,nnodes);
p = [won;x];
B = p .* [w; w] ;
pp = zeros(2);
A = zeros(2);
dA = zeros(2);
for i=l:nnodes
pp = p(1:2,i)*p(l:2,i);;
A = A+w(l,i)*pp;
dA = dA+dw(l,i)*pp;
end

Ainv = inv(A);
pg = [1 xg] ;
phi = pg*Ainv*B;
db = p.* [dw;dw];
da = -Ainv*(dA*Ainv);

dphi = [0 1]*Ainv*B+pg*(da*B+Ainv*db);

% °/0ASSEMBLE DISCRETE EQUATIONS
if j == 1
GG(1:3,1) = -phi(1:3);

else if j>l
k = k+(weight*E*jac)*(dphi.*dphi);
fbody = xg;

f = f+(weight*fbody*jac)*phi;;
end

end
end

% °/0 ENFORCE BOUNDARY CONDITION USING LAGRANGE MULTIPLIERS
q=[0];

m = [k GG;GG; zeros(l)];
% °/0 SOLVE FOR NODAL PARAMETERS
d = m\[f > q] > ;
u = d(l:nnodes)

% 32


