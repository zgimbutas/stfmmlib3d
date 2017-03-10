function y0=st3dmultfmmflat_neu(A,x0)
%
%     The exterior Stokes solver with Neumann boundary conditions,
%     on a flat triangulated surface. 
%
%     Integral representation is 
%
%     u = S_0 sigma + alpha * D_0 S_0 sigma
%
%     We do not normalize the Green's function by 4 pi 
%     Integral equation is 
%
%     -(4*pi/2) sol + (S_0' + alpha D_0' S_0) sol = f
%
%
'sl3dmultfmmflat_neu'

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

if( A.alpha ~= 0 )
%
%     Now compute S_0' S_0' x (needed in Calderon identity)
%
ifcharge = 1;
ifdipole = 0;
ifpot = 1;
ifgrad = 1;
[C]=stfmm3dtria(A.iprec,...
       A.ntri,A.triangles,A.trianorm,A.source,...
       ifcharge,S.potn,ifdipole,x,ifpot,ifgrad);
C.potn = st3dtraction(C.pre,C.grad,A.trianorm);

%
%     Add Calderon contributions to y:
%     y = y+alpha*(D'_0 S_0)x = y+alpha*(-(4*pi)**2/4 + S_0' S_0')x
%
y = y + A.alpha * (-((4*pi)/2)^2*x+C.potn);
end


y = reshape(y,1,3*A.ntri);

y0=y.';

toc
