function v = gauss2(k)
l=1
% °/0 This function returns a matrix with 4 gauss points and their weights
v(l,l) =-.861136311594052575224;
v(l,2) =-.339981043584856264803;
v(l,3) = -v(l,2);
v(l,4) = -v(l,1);
v(2,1) =.347854845137453857373;
v(2,2) =.652145154862546142627;
v(2,3) = v(2,2);
v(2,4) = v(2,1);