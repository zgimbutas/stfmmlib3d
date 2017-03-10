function y0=st3dmultfmmflat_dir_s0(A,x0);
%
%     The exterior Stokes solver with Dirichlet boundary conditions,
%     on a flat triangulated surface. 
%
%     Integral representation is 
%
%     u = S_0 sigma
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%     S_0 sol = f
%
%
'st3dmultfmmflat_dir_s0'

x=x0.';

x = reshape(x,3,A.ntri);

ifcharge = 1;
ifdipole = 0;
ifpot = 1;
ifgrad = 0;
[U]=stfmm3dtria(A.iprec,...
     A.ntri,A.triangles,A.trianorm,A.source,...
     ifcharge,x,ifdipole,x,ifpot,ifgrad);

y = U.pot;

y = reshape(y,1,3*A.ntri);

y0=y.';

toc
