function cc = sy(il,iu,bb,dd,aa,cc)
% SY - Solves a tridiagonal system of equations using Thomas' algorithm
% cc = sy(bb,dd,aa,cc)
%
%  il = subscript of first equation
%  iu = subscript of last equation
%  bb = coefficient behind diagonal
%  dd = coefficient on diagonal
%  aa = coefficient ahead of diagonal
%  cc = element of constant vector on input, solution vector on return
%
% Anderson, Tannehill, and Pletcher (1984) pp. 549-550

% csherwood@usgs.gov
% last revised December 21, 2004

%...establish upper triangular matrix
lp = il+1;
for i = lp:iu
    r = bb(i)/dd(i-1);
    dd(i) = dd(i) -r*aa(i-1);
    cc(i) = cc(i) -r*cc(i-1);
end
%...back substitution
cc(iu) = cc(iu) ./dd(iu) ;
for i = lp:iu
    j = iu-i+il;
    cc(j) = (cc(j)-aa(j).*cc(j+1)) ./dd(j);
end
%...solution stored in cc
return