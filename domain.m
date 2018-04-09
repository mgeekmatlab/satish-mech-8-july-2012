function v = domain(gpos,x,dm,nnodes)

°/0DETERMINES NODES IN DOMAIN OF INFLUENCE OF POINT GPOS

dif = gpos*ones(l5nnodes)-x;

a = dm-abs(dif);

b = [a;zeros(15nnodes)] ;

c = all(b>=-100*eps);

v = find(c);