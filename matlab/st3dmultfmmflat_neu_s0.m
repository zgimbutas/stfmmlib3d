function y0=st3dmultfmmflat_neu_s0(A,x0)
%
%     The exterior Stokes solver with Neumann boundary conditions,
%     on a flat triangulated surface. 
%
%     Integral representation is 
%
%     u = S_0 sigma
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%     -(4*pi/2) sol + S_0' sol = f
%
%
'sl3dmultfmmflat_neu_s0'

x=x0.';

x = reshape(x,3,A.ntri);
%
%     Compute S_0' x 
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
ifgrad = 1;
[S]=stfmm3dtria(A.iprec,...
     A.ntri,A.triangles,A.trianorm,A.source,...
     ifcharge,x,ifdipole,x,ifpot,ifgrad);

S.potn = st3dtraction(S.pre,S.grad,A.trianorm);

%
%     Initialize y = -(4*pi/2)x  and increment by S_0'x 
%
y = S.potn - ((4*pi)/2)*x;


y = reshape(y,1,3*A.ntri);

y0=y.';

toc
