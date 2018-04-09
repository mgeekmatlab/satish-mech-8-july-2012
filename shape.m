function [phi,dphix,dphiy] = shape(gpo s,dmax}x,v,dm)

% °/0 EFG shape function and it's derivative with linear base

L = length(v);

won = ones(1,L);

nv = x(l:2,v);

p = [won;nv];

dif = gpos*won-nv;

t = dm(l:2,v)/dmax;

% °/0 WEIGHTS—W and dw are vectors

[w,dwdx,dwdy] = cubwgt(dif,t,v,dmax,dm);

% 39

B = p . * [w; w; w] ;
pp = zeros(3);
aa = zeros(3);
daax = zeros(3);
daay = zeros(3);
for i=l:L

pp = p(1:3,i)*p(l:3,i);;
aa = aa+w(l,i)*pp;
daax = daax+dwdx(l,i)*pp;
daay = daay+dwdy(l,i)*pp;

end

Pg = [1 gpos;];
[L,U,PERM] = lu(aa);
for i=l:3
if i==l

C = PERM*pg;;
elseif i==2

C = PERM*( [0 1 0]' - daax*gam(1:3,1));
elseif i==3

C = PERM*( [0 0 1]; - daay*gam(1:3,1));

end

D1 = C(l) ;

D2 = (C(2) - L(2,1)*D1);

D3 = (C(3) - L(3,1)*D1-L(3,2)*D2);

gam(3,i) = D3/U(3,3);

gam(2,i) = (D2 - U(2,3)*gam(3,i))/(U(2,2));

gam(l,i) = (D1 - U(1,2)*gam(2,i)-U(l,3)*gam(3,i))/U(l,1);

end

phi = gam(l:3,1);*B;

dbx = p.*[dwdx;dwdx;dwdx];

dby = p.* [dwdy;dwdy;dwdy];

dphix = gam(l:3,2);*B + gam(l:3,1);*dbx;

dphiy = gam(l:3,3);*B + gam(l:3,1);*dby;

