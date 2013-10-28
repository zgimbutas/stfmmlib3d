function [U]=stfmm3dpart_matlab2(iprec,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg)
%STFMM3DPART Stokes particle target FMM in R^3.
%
% Stokes FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interaction) and interactions with targets.
%
% No slip (zero-velocity) boundary condition at z=0
%
% [U]=STFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,SIGMA_DV,IFPOT,IFGRAD);
%
% [U]=STFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,SIGMA_DV,IFPOT,IFGRAD,...
%         NTARGET,TARGET,IFPOTTARG,IFGRADTARG);
%
%
% This subroutine evaluates the Stokes potential and gradient due
% to a collection of Stokes single and double forces. We use
%
%
%       ifsingle=1, stokeslet, f = sigma_sl
%       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
%       p = [r_j / r^3] f_j
%
%       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
%       u_i = [3 r_i r_j r_k / r^5] n_k g_j
%       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
%
%       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
%       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
%       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
%
%       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
%       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
%       p = 0
%
%       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
%                   g = sigma_dl, n = sigma_dv
%       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
%             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
%       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
%
% for the free-space Green's function, without the (1/4 pi) scaling.  
%
%
% Input parameters:
% 
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
%
% nsource - number of sources
% source - double (3,nsource): source locations
% ifsingle - single force computation flag
%
%         0 => do not compute
%         1 => include Stokes single force contribution
% 
% sigma_sl - double (3,nsource): single force strengths
% ifdouble - double force computation flag
%
%         0 => do not compute
%         1 => include Stokes double force contribution
%         2 => include Stokes stresslet contribution
%         3 => include Stokes rotlet contribution
%         4 => include Stokes doublet contribution
% 
% sigma_dl - double (3,nsource): double force strengths
% sigma_dv - double (3,nsource): double force orientation vectors 
%
% ifpot - velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgrad - gradient computation flag, 
%         1 => compute the gradient, otherwise no
%
% ntarget - number of targets
% target - double (3,ntarget): target locations
%
% ifpottarg - target velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgradtarg - target gradient computation flag, 
%         1 => compute the gradient, otherwise no
%
%
% Output parameters: 
%
% U.pot - double (3,nsource) - velocity field at source locations
% U.pre - double (nsource) - pressure at source locations
% U.grad - double (3,3,nsource) - gradient at source locations
% U.pottarg - double (3,ntarget) - velocity field at targets
% U.pretarg - double (ntarget) - pressure at targets
% U.gradtarg - double (3,3,ntarget) - gradient at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%


if( nargin == 8 ) 
  ifpot = 1;
  ifgrad = 1;
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 10 ) 
  ntarget = 0;
  target = zeros(3,ntarget);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 12 ) 
  ifpottarg = 1;
  ifgradtarg = 1;
end

ifsingle = double(ifsingle); ifdouble = double(ifdouble);
ifpot = double(ifpot); ifgrad = double(ifgrad); 
ifpottarg = double(ifpottarg); ifgradtarg = double(ifgradtarg); 

pot=zeros(3,nsource);
pre=zeros(1,nsource);
grad=zeros(3,3,nsource);
if( ntarget > 0 ),
pottarg=zeros(3,ntarget);
pretarg=zeros(1,ntarget);
gradtarg=zeros(3,3,ntarget);
else
pottarg=zeros(3,1);
pretarg=zeros(1,1);
gradtarg=zeros(3,3,1);
end

ier=0;

if_use_fmm = 1;



ifpot0 = 0; iffld0 = 0; ifhess0 = 0;
ifpottarg0 = 0; iffldtarg0 = 0; ifhesstarg0 = 0;

if( ifpot == 1 ), ifpot0 = 1; iffld0 = 1; end;
if( ifgrad == 1 ), iffld0 = 1; ifhess0 = 1; end;

if( ifpottarg == 1 ), ifpottarg0 = 1; iffldtarg0 = 1; end;
if( ifgradtarg == 1 ), iffldtarg0 = 1; ifhesstarg0 = 1; end;



for j=1:3

ifcharge = 0;
ifdipole = 0;

ifcharge = 0;
charge = zeros(1,nsource);

ifdipole = 0;
dipstr = zeros(1,nsource);
dipvec = zeros(3,nsource);

ifquad  = 0;
quadstr = zeros(1,nsource);
quadvec = zeros(6,nsource);

if( ifsingle == 1 ),

  ifcharge = 1;
  charge = sigma_sl(j,:)/2;

end

if( ifdouble == 1 || ifdouble == 2 || ifdouble == 4 ),

  ifdipole = 1;
  dipstr = ones(1,nsource);
  dipvec(1,:) = sigma_dl(j,:).*sigma_dv(1,:) + sigma_dl(1,:).*sigma_dv(j,:);
  dipvec(2,:) = sigma_dl(j,:).*sigma_dv(2,:) + sigma_dl(2,:).*sigma_dv(j,:);
  dipvec(3,:) = sigma_dl(j,:).*sigma_dv(3,:) + sigma_dl(3,:).*sigma_dv(j,:);
  dipvec = dipvec/2;

end

if( if_use_fmm == 1 ),
H=lfmm3dpart(iprec,nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot0,iffld0,...
        ntarget,target,ifpottarg0,iffldtarg0);
else
H=l3dpartdirect(nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot0,iffld0,...
        ntarget,target,ifpottarg0,iffldtarg0);
end
if( ifpot0 == 1 ), H.pot = real(H.pot); end;
if( iffld0 == 1 ), H.fld = real(H.fld); end;
if( ifhess0 == 1 ), H.hess = real(H.hess); end;
if( ifpottarg0 == 1 ), H.pottarg = real(H.pottarg); end;
if( iffldtarg0 == 1 ), H.fldtarg = real(H.fldtarg); end;
if( ifhesstarg0 == 1 ), H.hesstarg = real(H.hesstarg); end;

if( ifpot == 1 ),
   pot(j,:) = pot(j,:) + H.pot;
   pot(1,:) = pot(1,:) + source(j,:).*H.fld(1,:);
   pot(2,:) = pot(2,:) + source(j,:).*H.fld(2,:);
   pot(3,:) = pot(3,:) + source(j,:).*H.fld(3,:);
   pre = pre + 2*H.fld(j,:);
end

if( ifgrad == 1 ),

F.grad = zeros(3,3,nsource);
F.grad(1,1,:) = H.hess(1,:).*source(j,:);
F.grad(2,2,:) = H.hess(2,:).*source(j,:);
F.grad(3,3,:) = H.hess(3,:).*source(j,:);
F.grad(1,2,:) = H.hess(4,:).*source(j,:);
F.grad(1,3,:) = H.hess(5,:).*source(j,:);
F.grad(2,3,:) = H.hess(6,:).*source(j,:);
F.grad(2,1,:) = H.hess(4,:).*source(j,:);
F.grad(3,1,:) = H.hess(5,:).*source(j,:);
F.grad(3,2,:) = H.hess(6,:).*source(j,:);

F.grad(1:3,j,:) = reshape(F.grad(1:3,j,:),3,nsource) - H.fld;
F.grad(j,1:3,:) = reshape(F.grad(j,1:3,:),3,nsource) + H.fld;

grad = grad - F.grad;

end

if( ifpottarg == 1 ),
   pottarg(j,:) = pottarg(j,:) + H.pottarg;
   pottarg(1,:) = pottarg(1,:) + target(j,:).*H.fldtarg(1,:);
   pottarg(2,:) = pottarg(2,:) + target(j,:).*H.fldtarg(2,:);
   pottarg(3,:) = pottarg(3,:) + target(j,:).*H.fldtarg(3,:);
   pretarg = pretarg + 2*H.fldtarg(j,:);
end

if( ifgradtarg == 1 ),

F.gradtarg = zeros(3,3,ntarget);
F.gradtarg(1,1,:) = H.hesstarg(1,:).*target(j,:);
F.gradtarg(2,2,:) = H.hesstarg(2,:).*target(j,:);
F.gradtarg(3,3,:) = H.hesstarg(3,:).*target(j,:);
F.gradtarg(1,2,:) = H.hesstarg(4,:).*target(j,:);
F.gradtarg(1,3,:) = H.hesstarg(5,:).*target(j,:);
F.gradtarg(2,3,:) = H.hesstarg(6,:).*target(j,:);
F.gradtarg(2,1,:) = H.hesstarg(4,:).*target(j,:);
F.gradtarg(3,1,:) = H.hesstarg(5,:).*target(j,:);
F.gradtarg(3,2,:) = H.hesstarg(6,:).*target(j,:);

F.gradtarg(1:3,j,:) = reshape(F.gradtarg(1:3,j,:),3,ntarget) - H.fldtarg;
F.gradtarg(j,1:3,:) = reshape(F.gradtarg(j,1:3,:),3,ntarget) + H.fldtarg;

gradtarg = gradtarg - F.gradtarg;

end

end


ifcharge = 0;
ifdipole = 0;

ifcharge = 0;
charge = zeros(1,nsource);

ifdipole = 0;
dipstr = zeros(1,nsource);
dipvec = zeros(3,nsource);

ifquad  = 0;
quadstr = zeros(1,nsource);
quadvec = zeros(6,nsource);

if( ifsingle == 1 ),

  ifcharge = 1;
  charge = (sigma_sl(1,:).*source(1,:) + ...
            sigma_sl(2,:).*source(2,:) + ...
            sigma_sl(3,:).*source(3,:))/2;

end

if( ifdouble == 2 || ifdouble == 4 ),

  ifcharge = 1;
  charge = charge + ...
           (sigma_dl(1,:).*sigma_dv(1,:) + ...
            sigma_dl(2,:).*sigma_dv(2,:) + ...
            sigma_dl(3,:).*sigma_dv(3,:));

end

if( ifdouble == 1 || ifdouble == 2 || ifdouble == 4 ),

  ifdipole = 1;
  dipstr = ones(1,nsource);
  ul = (sigma_dl(1,:).*source(1,:) + ...
        sigma_dl(2,:).*source(2,:) + ...
        sigma_dl(3,:).*source(3,:));
  uv = (sigma_dv(1,:).*source(1,:) + ...
        sigma_dv(2,:).*source(2,:) + ...
        sigma_dv(3,:).*source(3,:));
  dipvec(1,:) = sigma_dv(1,:).*ul + sigma_dl(1,:).*uv;
  dipvec(2,:) = sigma_dv(2,:).*ul + sigma_dl(2,:).*uv;
  dipvec(3,:) = sigma_dv(3,:).*ul + sigma_dl(3,:).*uv;
  dipvec = dipvec/2;

end

if( if_use_fmm == 1 ),
H=lfmm3dpart(iprec,nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot0,iffld0,...
        ntarget,target,ifpottarg0,iffldtarg0);
else
H=l3dpartdirect(nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot0,iffld0,...
        ntarget,target,ifpottarg0,iffldtarg0);
end
if( ifpot0 == 1 ), H.pot = real(H.pot); end;
if( iffld0 == 1 ), H.fld = real(H.fld); end;
if( ifhess0 == 1 ), H.hess = real(H.hess); end;
if( ifpottarg0 == 1 ), H.pottarg = real(H.pottarg); end;
if( iffldtarg0 == 1 ), H.fldtarg = real(H.fldtarg); end;
if( ifhesstarg0 == 1 ), H.hesstarg = real(H.hesstarg); end;


if( ifpot == 1 ),
   pot = pot - H.fld;
end

if( ifgrad == 1 ),

F.grad = zeros(3,3,nsource);
F.grad(1,1,:) = H.hess(1,:);
F.grad(2,2,:) = H.hess(2,:);
F.grad(3,3,:) = H.hess(3,:);
F.grad(1,2,:) = H.hess(4,:);
F.grad(1,3,:) = H.hess(5,:);
F.grad(2,3,:) = H.hess(6,:);
F.grad(2,1,:) = H.hess(4,:);
F.grad(3,1,:) = H.hess(5,:);
F.grad(3,2,:) = H.hess(6,:);

grad = grad + F.grad;

end

if( ifpottarg == 1 ),
   pottarg = pottarg - H.fldtarg;
end

if( ifgradtarg == 1 ),

F.gradtarg = zeros(3,3,ntarget);
F.gradtarg(1,1,:) = H.hesstarg(1,:);
F.gradtarg(2,2,:) = H.hesstarg(2,:);
F.gradtarg(3,3,:) = H.hesstarg(3,:);
F.gradtarg(1,2,:) = H.hesstarg(4,:);
F.gradtarg(1,3,:) = H.hesstarg(5,:);
F.gradtarg(2,3,:) = H.hesstarg(6,:);
F.gradtarg(2,1,:) = H.hesstarg(4,:);
F.gradtarg(3,1,:) = H.hesstarg(5,:);
F.gradtarg(3,2,:) = H.hesstarg(6,:);

gradtarg = gradtarg + F.gradtarg;

end



if( ifdouble == 3 || ifdouble == 4 ),

for j = 1:3,

ifcharge = 0;
ifdipole = 0;

ifcharge = 0;
charge = zeros(1,nsource);

ifdipole = 0;
dipstr = zeros(1,nsource);
dipvec = zeros(3,nsource);

ifquad  = 0;
quadstr = zeros(1,nsource);
quadvec = zeros(6,nsource);

ifdipole = 1;
dipstr = ones(1,nsource);
dipvec(1,:) = sigma_dv(1,:).*sigma_dl(j,:) - sigma_dl(1,:).*sigma_dv(j,:);
dipvec(2,:) = sigma_dv(2,:).*sigma_dl(j,:) - sigma_dl(2,:).*sigma_dv(j,:);
dipvec(3,:) = sigma_dv(3,:).*sigma_dl(j,:) - sigma_dl(3,:).*sigma_dv(j,:);

if( if_use_fmm == 1 ),
H=lfmm3dpartquad(iprec,nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifpot0,iffld0,ifhess0,...
        ntarget,target,ifpottarg0,iffldtarg0,ifhesstarg0);
else
H=l3dpartquaddirect(nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifpot0,iffld0,ifhess0,...
        ntarget,target,ifpottarg0,iffldtarg0,ifhesstarg0);
end
if( ifpot0 == 1 ), H.pot = real(H.pot); end;
if( iffld0 == 1 ), H.fld = real(H.fld); end;
if( ifhess0 == 1 ), H.hess = real(H.hess); end;
if( ifpottarg0 == 1 ), H.pottarg = real(H.pottarg); end;
if( iffldtarg0 == 1 ), H.fldtarg = real(H.fldtarg); end;
if( ifhesstarg0 == 1 ), H.hesstarg = real(H.hesstarg); end;


if( ifpot == 1 ),
   pot(j,:) = pot(j,:) + H.pot;
end

if( ifgrad == 1 ),
   grad(j,1:3,:) = reshape(grad(j,1:3,:),3,nsource) - H.fld;
end

if( ifpottarg == 1 ),
   pottarg(j,:) = pottarg(j,:) + H.pottarg;
end

if( ifgradtarg == 1 ),
   gradtarg(j,1:3,:) = reshape(gradtarg(j,1:3,:),3,ntarget) - H.fldtarg;
end


end

end


if( ifpot == 1 ), U.pot = pot; end
if( ifpot == 1 ), U.pre = pre; end
if( ifgrad == 1 ), U.grad = grad; end
if( ifpottarg == 1 ), U.pottarg = pottarg; end
if( ifpottarg == 1 ), U.pretarg = pretarg; end
if( ifgradtarg == 1 ), U.gradtarg = gradtarg; end


