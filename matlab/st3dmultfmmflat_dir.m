function y0=st3dmultfmmflat_dir(A,x0)
%
%     The exterior Stokes solver with Dirichlet boundary conditions,
%     on a flat triangulated surface. 
%
%     Integral representation is 
%
%     u = (beta S_0 + D_0) sigma
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%     (4*pi/2) sol + (beta S_0 + D_0) sol = f
%
%
'st3dmultfmmflat_dir'

x=x0.';

x = reshape(x,3,A.ntri);

charge = A.beta * x;
dipstr = x;

ifcharge = 1;
ifdipole = 1;
ifpot = 1;
ifgrad = 0;
[U]=stfmm3dtria(A.iprec,...
     A.ntri,A.triangles,A.trianorm,A.source,...
     ifcharge,charge,ifdipole,dipstr,ifpot,ifgrad);

y = U.pot + ((4*pi)/2) * x;

y = reshape(y,1,3*A.ntri);

y0=y.';

toc
