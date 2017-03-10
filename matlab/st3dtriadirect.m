function [U]=st3dtriadirect(nsource,triaflat,trianorm,source,ifsingle,sigma_sl,ifdouble,sigma_dl,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg)
%STFMM3DTRIADIRECT Stokes triangle interactions in R^3, direct algorithm.
%
%
% [U]=ST3DTRIADIRECT(NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,IFPOT,IFGRAD);
%
% [U]=ST3DTRIADIRECT(NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFSINGLE,SIGMA_SL,IFDOUBLE,SIGMA_DL,IFPOT,IFGRAD,...
%         NTARGET,TARGET,IFPOTTARG,IFGRADTARG);
%
%
% This subroutine evaluates the Stokes potential (velocity/pressure)
% and velocity gradient due
% to a collection of flat triangles with constant single and/or
% double layer densities. We use 
%
%       \delta u = \grad p, div u = 0, mu = 1.
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
% for the Green's function, without the (1/4 pi) scaling.  
% Self-interactions are included.
%
% It is capable of evaluating the layer potentials either on or 
% off the surface (or both).            
%
%
% Input parameters:
% 
% nsource - number of triangles
% triaflat - double (3,3,ntriangles): array of triangle vertex coordinates
% trianorm - double (3,ntriangles): triangle normals
% source - double (3,ntriangles): triangle centroids
% ifsingle - single layer computation flag
%
%         0 => do not compute
%         1 => include Stokes SLP contribution
% 
% sigma_sl - double (3,ntriangles): piecewise constant SLP 
%                                   (single force) strength 
% ifdouble - double layer computation flag
%
%         0 => do not compute
%         1 => include Stokes DLP contribution
%         2 => include Stokes stresslet contribution
%         3 => include Stokes rotlet contribution
%         4 => include Stokes doublet contribution
% 
% sigma_dl - double (3,ntriangles): piecewise constant DLP 
%                                   (double force) strength 
%
%     In the present version, double force orientation vector is assumed to
%     BE SET EQUAL to the triangle normal. 
%
% ifpot - velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgrad - velocity gradient computation flag, 
%         1 => compute the velocity gradient, otherwise no
%
% ntarget - number of targets
% target - double (3,ntarget): target locations
%
% ifpottarg - target velocity field/pressure computation flag, 
%         1 => compute the velocity field/pressure, otherwise no
% ifgradtarg - target velocity gradient computation flag, 
%         1 => compute the velocity gradient, otherwise no
%
%
% Output parameters: 
%
% U.pot - double (3,nsource) - velocity field at triangle centroids
% U.pre - double (nsource) - pressure at triangle centroids
% U.grad - double (3,3,nsource) - velocity gradient at triangle centroids
% U.pottarg - double (3,ntarget) - velocity field at targets
% U.pretarg - double (ntarget) - pressure at targets
% U.gradtarg - double (3,3,ntarget) - velocity gradient at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%

if( nargin == 8 ) 
  ifpot = 1;
  ifgrad = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  ifgradtarg = 0;
end

if( nargin == 10 ) 
  ntarget = 0;
  target = zeros(3,1);
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


mex_id_ = 'st3dtriadirect(i double[], i double[], i int[x], i double[xx], i int[x], i double[], i int[x], i double[], i int[x], io double[], io double[], i int[x], io double[], i int[x], i double[], i int[x], io double[], io double[], i int[x], io double[])';
[pot, pre, grad, pottarg, pretarg, gradtarg] = stfmm3d_r2012b(mex_id_, triaflat, trianorm, nsource, source, ifsingle, sigma_sl, ifdouble, sigma_dl, ifpot, pot, pre, ifgrad, grad, ntarget, target, ifpottarg, pottarg, pretarg, ifgradtarg, gradtarg, 1, 3, nsource, 1, 1, 1, 1, 1, 1, 1);


if( ifpot == 1 ), U.pot=pot; end
if( ifpot == 1 ), U.pre=pre; end
if( ifgrad == 1 ), U.grad=reshape(grad,3,3,nsource); end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( ifpottarg == 1 ), U.pretarg=pretarg; end
if( ifgradtarg == 1 ), U.gradtarg=reshape(gradtarg,3,3,ntarget); end
U.ier=ier;


