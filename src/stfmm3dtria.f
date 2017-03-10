cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains the FMM routines for Stokes layer
c        potentials in free space in R^3. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      stfmm3dtriaself - evaluates the Stokes potential ON SURFACE due
c         to a collection of flat triangles with piecewise constant
c         single and/or double layer densities using the Fast Multipole
c         Method.
c
c      stfmm3dtriatarg - evaluates the Stokes potential ON OR OFF
c         SURFACE due to a collection of flat triangles with constant
c         single and/or double layer densities using the Fast Multipole
c         Method.
c
c      st3dtriadirect - evaluates the Stokes potential ON OR
c         OFF SURFACE due to a collection of flat triangles with
c         constant single and/or double layer densities using the direct
c         O(N^2) algorithm.
c
c
c      supst3triadirectself - evaluates the Stokes potential at a
c         single surface point (triangle centroid) due to a collection
c         of flat triangles with piecewise constant single layer density
c         BY DIRECT CALCULATION.
c
c      supst3triadirecttarg - evaluates the Stokes potential at an
c         (OFF SURFACE) target due to a collection of flat triangles
c         with piecewise constant single layer density BY DIRECT
c         CALCULATION.
c
c      stpst3triadirectself - evaluates the Stokes potential at a
c         single surface point (triangle centroid) due to a collection
c         of flat triangles with piecewise constant double layer density
c         BY DIRECT CALCULATION.
c     
c      stpst3triadirecttarg - evaluates the Stokes potential at an
c         (OFF SURFACE) target due to a collection of flat triangles
c         with piecewise constant double layer density BY DIRECT
c         CALCULATION.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
      subroutine stfmm3dtriaself
     $     (ier,iprec,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifpot,pot,pre,ifgrad,grad)
c
c     This subroutine evaluates the Stokes potential 
c     due to a collection of flat triangles with constant
c     single and/or double layer densities using the Fast Multipole Method.
c
c     Unfortunately, in different communities the source strengths
c     are referred to by different names. Here we will use the 
c     mathematical conventions of single and double layer potentials.
c     
c   
c     The single layer Green's function maps 
c     traction (surface force density) to velocity/gradient.
c     (In some classical literature, this is called the <<single force>>.)
c     We will use SIGMA_SL to describe the source density.
c
c     (The resulting velocity is continuous across the surface.)
c
c     The double layer Green's function maps 
c     generalized slip/jump in velocity to velocity/gradient
c     (In some classical literature, this is called the <<double force>>.)
c     We will use SIGMA_DL to describe the source density.
c
c     (The resulting velocity is NOT continuous across the surface.)
c
c     \delta u = \grad p, div u = 0, mu = 1.
c
c     Free space Stokes Green's functions:
c
c       ifsingle=1, stokeslet, f = sigma_sl
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
c       u_i = [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
c       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 0
c
c       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
c                   g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
c             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c     without the (1/4 pi) scaling.
c
c
c     INPUT:
c
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of single layer source 
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of double layer source
c
c     iprec:  FMM precision flag
c
c             -2 => tolerance =.5d0
c             -1 => tolerance =.5d-1
c             0 => tolerance =.5d-2
c             1 => tolerance =.5d-3
c             2 => tolerance =.5d-6
c             3 => tolerance =.5d-9
c             4 => tolerance =.5d-12
c             5 => tolerance =.5d-15
c
c     These errors were established for the electrostatic case and are
c     approximately valid for velocity. For grad, there is
c     approximately one digit of loss in the above table.
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     pre(3,nparts) = pressure at source locations
c     grad(3,3,nparts) = gradient at source locations
c
c     The main FMM routine permits both evaluation on surface
c     and at a collection of off-surface targets. 
c     This subroutine is used to simplify the user interface 
c     (by setting the number of targets to zero) and calling the more 
c     general FMM.
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)
        real *8 pot(3,nparts),pre(nparts),grad(3,3,nparts)
        integer nparts,ntargs

        ntargs=0
        ifpottarg=0
        ifgradtarg=0
        call stfmm3dtriatarg
     $     (ier,iprec,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,POTtarg,PREtarg,
     $     ifgradtarg,GRADtarg)

        return
        end
c
c
c
c
c*********************************
      subroutine stfmm3dtriatarg
     $     (ier,iprec,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
c       This is the principal subroutine for evaluating Stokes
c       layer potentials on (flat) triangulated surfaces.  It permits
c       the evaluation of a single layer potential with piecewise
c       constant density defined by the (force) vector sigma_sl and a
c       dipole layer potential with piecewise constant density vector
c       sigma_dl and dipole orientation defined by trinorm.
c
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutines 
c       stfmm3dtriatargmain_fast  for the far field interactions and 
c       and stfmm3dtriatarg0  for the near field interactions.
c 
c     INPUT:
c
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     pre(3,nparts) = pressure at source locations
c     grad(3,3,nparts) = gradient at source locations
c     pottarg(3,ntargs) = velocity at target locations
c     pretarg(3,ntargs) = pressure at target locations
c     gradtarg(3,3,ntargs) = gradient at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts)
        real *8 pot(3,nparts),pre(nparts),grad(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 pottarg(3,ntargs),pretarg(ntargs),
     $     gradtarg(3,3,ntargs)
        integer nparts,ntargs
        real *8, allocatable :: w(:)
c       
        ier=0
        lused=0
c
c       ... allocate work arrays
c
        icharge=lused+1
        lcharge=2*nparts *3
        lused=lused+lcharge

        idipstr=lused+1
        ldipstr=2*nparts *3
        lused=lused+ldipstr

        idipvec=lused+1
        ldipvec=3*nparts *3
        lused=lused+ldipvec

        icpot=lused+1
        lcpot=2*nparts
        lused=lused+lcpot

        icfld=lused+1
        lcfld=2*3*nparts
        lused=lused+lcfld

        ichess=lused+1
        lchess=2*6*nparts
        lused=lused+lchess

        ihessmatr=lused+1
        lhessmatr=3*3*nparts
        lused=lused+lhessmatr

        icpottarg=lused+1
        lcpottarg=2*ntargs
        lused=lused+lcpottarg

        icfldtarg=lused+1
        lcfldtarg=2*3*ntargs
        lused=lused+lcfldtarg

        ichesstarg=lused+1
        lchesstarg=2*6*ntargs
        lused=lused+lchesstarg

        ihessmatrtarg=lused+1
        lhessmatrtarg=3*3*ntargs
        lused=lused+lhessmatrtarg
c
        allocate( w(lused), stat=ier)
        if( ier .ne. 0 ) return
c
c     call FMM to account for all far field interactions.
c     
c     The subroutine stfmm3dtriatargmain_fast makes 4 calls
c     to a scalar (Laplace FMM) according to the canonical
c     decomposition of the Stokes Green fucntions.
c
c
        call stfmm3dtriatargmain_fast
     $     (ier,iprec,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,trinorm,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg,
     $     w(icharge),w(idipstr),w(idipvec),
     $     w(icpot),w(icfld),w(ichess),w(ihessmatr),
     $     w(icpottarg),w(icfldtarg),w(ichesstarg),w(ihessmatrtarg))
c
c     reconstruct FMM data structure and account for  all local 
c     interactions using quadrature routines for piecewise
c     constant densities - no interactions are saved in the 
c     present version.
c     
        call stfmm3dtriatarg0
     $     (ier,iprec,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,trinorm,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
        return
        end
c
c
c
c
c
c*********************************
      subroutine stfmm3dtriatargmain_fast
     $     (ier,iprec,triangle,trinorm,nparts,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg,
     $     charge,dipstr,dipvec,cpot,cfld,chess,hessmatr,
     $     cpottarg,cfldtarg,chesstarg,hessmatrtarg)
c
c     FMM calculation subroutine for Stokes N-body problem
c
c     4 Laplace FMM calls, uses linear charge and dipole densities
c
c     3 Laplace FMM calls for rotlet
c     7 Laplace FMM calls for doublet
c
c
c     INPUT:
c
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nparts = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c             -2 => tolerance =.5d0
c             -1 => tolerance =.5d-1
c             0 => tolerance =.5d-2
c             1 => tolerance =.5d-3
c             2 => tolerance =.5d-6
c             3 => tolerance =.5d-9
c             4 => tolerance =.5d-12
c             5 => tolerance =.5d-15
c
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     pre(3,nparts) = pressure at source locations
c     grad(3,3,nparts) = gradient at source locations
c     pottarg(3,ntargs) = velocity at target locations
c     pretarg(3,ntargs) = pressure at target locations
c     gradtarg(3,3,ntargs) = gradient at target locations
c
c
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,nparts),trinorm(3,nparts)
        real *8 source(3,nparts)
        real *8 sigma_sl(3,nparts)
        real *8 sigma_dl(3,nparts),sigma_dv(3,nparts)
        real *8 pot(3,nparts),pre(nparts),grad(3,3,nparts)
        real *8 target(3,ntargs)
        real *8 pottarg(3,ntargs),pretarg(ntargs),
     $     gradtarg(3,3,ntargs)
        integer nparts,ntargs
c       
        complex *16 charge(3,1)
        complex *16 dipstr(3,1)
        real *8 dipvec(3,3,1)
c
        complex *16 cpot(1)
        complex *16 cfld(3,1)
        complex *16 chess(6,1) 
        real *8 hessmatr(3,3,1)

        complex *16 cpottarg(1)
        complex *16 cfldtarg(3,1)
        complex *16 chesstarg(6,1) 
        real *8 hessmatrtarg(3,3,1)
c
c
        do k=1,nparts
c
        if( ifpot .eq. 1 ) then
           pot(1,k) = 0.0d0
           pot(2,k) = 0.0d0
           pot(3,k) = 0.0d0
           pre(k) = 0.0d0
        endif
C
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
           grad(i,j,k) = 0.0d0
           hessmatr(i,j,k) = 0.0d0
        enddo
        enddo
        endif        
c
        enddo
c
c
        do k=1,ntargs
c
        if( ifpottarg .eq. 1 ) then
           pottarg(1,k) = 0.0d0
           pottarg(2,k) = 0.0d0
           pottarg(3,k) = 0.0d0
           pretarg(k) = 0.0d0
        endif
C
        if( ifgradtarg .eq. 1 ) then
        do i=1,3
        do j=1,3
           gradtarg(i,j,k) = 0.0d0
           hessmatrtarg(i,j,k) = 0.0d0
        enddo
        enddo
        endif
c
        enddo
c
c
c
        ifpot0=0
        iffld0=0
        ifhess0=0
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
        ifpot0=1
        iffld0=1
        endif
        if( ifgrad .eq. 1 ) ifhess0=1

        ifpottarg0=0
        iffldtarg0=0
        ifhesstarg0=0
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
        ifpottarg0=1
        iffldtarg0=1
        endif
        if( ifgradtarg .eq. 1 ) ifhesstarg0=1
c
c
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
c       Note, that the constructed charge and dipole vectors are constant
c       functions on the triangle,  we call lfmm3dtrilhesstftarg routine 
c       to simplify the things (we can use lfmm3dtrilhesstftarg here)
c       
        do j = 1,3

           ifcharge=0
           ifdipole=0

           do m=1,3
           do k = 1,nparts
              charge(m,k) = 0
              dipstr(m,k) = 0
              dipvec(1,m,k) = 0
              dipvec(2,m,k) = 0
              dipvec(3,m,k) = 0
              if( ifsingle .eq. 1 ) then
                charge(m,k) = sigma_sl(j,k)/2
                ifcharge=1
              endif
              if( ifdouble .eq. 1 .or. ifdouble .eq. 2 
     $           .or. ifdouble .eq. 4 ) then
                DIPSTR(m,k) = 1
                dipvec(1,m,k)=sigma_dv(1,k)*sigma_dl(j,k)
                dipvec(2,m,k)=sigma_dv(2,k)*sigma_dl(j,k)
                dipvec(3,m,k)=sigma_dv(3,k)*sigma_dl(j,k)
                dipvec(1,m,k)=dipvec(1,m,k)+sigma_dl(1,k)*sigma_dv(j,k)
                dipvec(2,m,k)=dipvec(2,m,k)+sigma_dl(2,k)*sigma_dv(j,k)
                dipvec(3,m,k)=dipvec(3,m,k)+sigma_dl(3,k)*sigma_dv(j,k)
                dipvec(1,m,k)=dipvec(1,m,k)/2
                dipvec(2,m,k)=dipvec(2,m,k)/2
                dipvec(3,m,k)=dipvec(3,m,k)/2
                ifdipole=1
              endif
           enddo
           enddo
c
           if( ifcharge .ne. 0 .or. ifdipole .ne. 0 ) then

           call lfmm3dtrilhesstftarg(ier,iprec,
     $        nparts,triangle,trinorm,source,
     $        ifcharge,charge,ifdipole,dipstr,dipvec,
     $        ifpot0,cpot,iffld0,cfld,ifhess0,chess,
     $        ntargs,target,ifpottarg0,cpottarg,iffldtarg0,cfldtarg,
     $        ifhesstarg0,chesstarg)
c
           call stfmm3dlap1(nparts,j,cpot,cfld,chess,
     $        source,ifpot,pot,pre,ifgrad,grad)
           call stfmm3dlap1(ntargs,j,cpottarg,cfldtarg,chesstarg,
     $        target,ifpottarg,pottarg,pretarg,
     $        ifgradtarg,gradtarg)

           endif

        enddo
c
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c
c       Note, that the constructed charge and dipole vectors are linear
c       functions on the triangle, we must call lfmm3dtrilhesstftarg routine
c       here (a significant loss of accuracy in near field will occur if
c       we use lfmm3dtriahesstarg instead)
c       
        ifcharge=0
        ifdipole=0
c        
        do m=1,3
        do k = 1,nparts
          charge(m,k) = 0
          dipstr(m,k) = 0
          dipvec(1,m,k) = 0
          dipvec(2,m,k) = 0
          dipvec(3,m,k) = 0
          if( ifsingle .eq. 1 ) then
          charge(m,k) = 
     $      (sigma_sl(1,k)*triangle(1,m,k)+
     $       sigma_sl(2,k)*triangle(2,m,k)+
     $       sigma_sl(3,k)*triangle(3,m,k))/2
          ifcharge = 1
          endif
          if( ifdouble .eq. 2 .or. ifdouble .eq. 4 ) then
          charge(m,k) = charge(m,k) + 
     $        (sigma_dl(1,k)*sigma_dv(1,k) + 
     1         sigma_dl(2,k)*sigma_dv(2,k) + 
     2         sigma_dl(3,k)*sigma_dv(3,k))
          ifcharge = 1
          endif
          if( ifdouble .eq. 1 .or. ifdouble. eq. 2 
     $       .or. ifdouble .eq. 4 ) then
          dipstr(m,k) = 1 
          dipvec(1,m,k) = sigma_dv(1,k)*
     $        (sigma_dl(1,k)*triangle(1,m,k) + 
     1         sigma_dl(2,k)*triangle(2,m,k) + 
     2         sigma_dl(3,k)*triangle(3,m,k) )
          dipvec(2,m,k) = sigma_dv(2,k)*
     $        (sigma_dl(1,k)*triangle(1,m,k) + 
     1         sigma_dl(2,k)*triangle(2,m,k) + 
     2         sigma_dl(3,k)*triangle(3,m,k) )
          dipvec(3,m,k) = sigma_dv(3,k)*
     $        (sigma_dl(1,k)*triangle(1,m,k) + 
     1         sigma_dl(2,k)*triangle(2,m,k) + 
     2         sigma_dl(3,k)*triangle(3,m,k) )
          dipvec(1,m,k) = dipvec(1,m,k) + sigma_dl(1,k)*
     $        (sigma_dv(1,k)*triangle(1,m,k) + 
     1         sigma_dv(2,k)*triangle(2,m,k) + 
     2         sigma_dv(3,k)*triangle(3,m,k))
          dipvec(2,m,k) = dipvec(2,m,k) + sigma_dl(2,k)*
     $        (sigma_dv(1,k)*triangle(1,m,k) + 
     1         sigma_dv(2,k)*triangle(2,m,k) + 
     2         sigma_dv(3,k)*triangle(3,m,k))
          dipvec(3,m,k) = dipvec(3,m,k) + sigma_dl(3,k)*
     $        (sigma_dv(1,k)*triangle(1,m,k) + 
     1         sigma_dv(2,k)*triangle(2,m,k) + 
     2         sigma_dv(3,k)*triangle(3,m,k))
          dipvec(1,m,k) = dipvec(1,m,k)/2
          dipvec(2,m,k) = dipvec(2,m,k)/2
          dipvec(3,m,k) = dipvec(3,m,k)/2
          ifdipole = 1          
          endif
        enddo
        enddo
c
        if( ifcharge .ne. 0 .or. ifdipole .ne. 0 ) then

        call lfmm3dtrilhesstftarg(ier,iprec,
     $     nparts,triangle,trinorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot0,cpot,iffld0,cfld,ifhess0,chess,
     $     ntargs,target,ifpottarg0,cpottarg,iffldtarg0,cfldtarg,
     $     ifhesstarg0,chesstarg)
c
        call stfmm3dlap2(nparts,cpot,cfld,chess,
     $     ifpot,pot,ifgrad,grad)
        call stfmm3dlap2(ntargs,cpottarg,cfldtarg,chesstarg,
     $     ifpottarg,pottarg,ifgradtarg,gradtarg)
c
        endif
c
c
        if( ifdouble .eq. 3 .or. ifdouble .eq. 4 ) then
c
c       ... rotlet part
c
c       Combine dipoles linearly. It is possible to do so, since both
c       dipstr and dipvec are real numbers in this calculation (in
c       general case, one would have to introduce complex dipvec
c       vectors, and rewrite the underlying FMM). 
c        
        do j = 1,3
c
        ifcharge=0
        ifdipole=0
c        
        do m=1,3
        do k = 1,nparts
          charge(m,k) = 0
          dipstr(m,k) = 0
          dipvec(1,m,k) = 0
          dipvec(2,m,k) = 0
          dipvec(3,m,k) = 0
          if( ifdouble .eq. 3 .or. ifdouble .eq. 4 ) then
          dipstr(m,k) = 1
          dipvec(1,m,k) = sigma_dv(1,k)*sigma_dl(j,k)
          dipvec(2,m,k) = sigma_dv(2,k)*sigma_dl(j,k)
          dipvec(3,m,k) = sigma_dv(3,k)*sigma_dl(j,k)
          dipvec(1,m,k) = dipvec(1,m,k)-sigma_dl(1,k)*sigma_dv(j,k)
          dipvec(2,m,k) = dipvec(2,m,k)-sigma_dl(2,k)*sigma_dv(j,k)
          dipvec(3,m,k) = dipvec(3,m,k)-sigma_dl(3,k)*sigma_dv(j,k)
          ifdipole = 1
          endif
        enddo
        enddo

        if( ifcharge .ne. 0 .or. ifdipole .ne. 0 ) then

        ifhess0=0
        ifhesstarg0=0
        call lfmm3dtrilhesstftarg(ier,iprec,
     $     nparts,triangle,trinorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot0,cpot,iffld0,cfld,ifhess0,chess,
     $     ntargs,target,ifpottarg0,cpottarg,iffldtarg0,cfldtarg,
     $     ifhesstarg0,chesstarg)

        call stfmm3dlap3(nparts,j,cpot,cfld,
     $     ifpot,pot,ifgrad,grad)
        call stfmm3dlap3(ntargs,j,cpottarg,cfldtarg,
     $     ifpottarg,pottarg,ifgradtarg,gradtarg)

        endif

        enddo

        endif

        return
        end
C
c
c
c
c
c
c
c
        subroutine stfmm3dtriatarg0(ier,iprec,
     $     triangle,trinorm,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntarget,target,
     $     ifpottarg,pottarg,pretarg,ifgradtarg,gradtarg)
c
c     FMM calculation subroutine for Stokes N-body problem
c
c     Direct evaluation routine, for local interactions only
c
c     Constant densities on flat triangles.
c     Note that currently, SIGMA_DV must be the same as TRINORM.
c
c     INPUT:
c
c     triangle(3,3,nparts) = array of triangles in standard format
c     trinorm(3,nparts)    = array of triangle normals
c     nsource = number of sources
c     source(3,nparts) = source locations
c     ifsingle = single layer computation flag  
c     sigma_sl(3,nparts) = vector strength of nth charge (single layer)
c     ifdouble = double layer computation flag  
c     sigma_dl(3,nparts) = vector strength of nth dipole (double layer)
c     sigma_dv(3,nparts) = dipole orientation vectors (double layer)
c     target(3,ntargs) = evaluation target points
c
c     iprec:  FMM precision flag
c
c     OUTPUT:
c
c     pot(3,nparts) = velocity at source locations
c     pre(3,nparts) = pressure at source locations
c     grad(3,3,nparts) = gradient at source locations
c     pottarg(3,ntargs) = velocity at target locations
c     pretarg(3,ntargs) = pressure at target locations
c     gradtarg(3,3,ntargs) = gradient at target locations
c
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,1),trinorm(3,1)
        real *8 source(3,1)
        real *8 sigma_sl(3,1)
        real *8 sigma_dl(3,1)
        real *8 sigma_dv(3,1)
        real *8 pot(3,1)
        real *8 pre(1)
        real *8 grad(3,3,1)
        real *8 target(3,1)
        real *8 pottarg(3,1)
        real *8 pretarg(1)
        real *8 gradtarg(3,3,1)
c
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
c
        real *8 timeinfo(10)
c       
        real *8 center(3)
c       
        integer laddr(2,200)
        real *8 scale(0:200)
        real *8 bsize(0:200)
        integer nterms(0:200)
c       
        integer box(20)
        real *8 center0(3),corners0(3,8)
c       
        integer box1(20)
        real *8 center1(3),corners1(3,8)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
c
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c       ... build the oct-tree
c       
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c       
        call prin2('epsfmm=*',epsfmm,1)
c
        if( iprec .eq. -2 ) nbox=8/3
        if( iprec .eq. -1 ) nbox=15/3
        if( iprec .eq. 0 ) nbox=30/3
        if( iprec .eq. 1 ) nbox=60/3
        if( iprec .eq. 2 ) nbox=120/3
        if( iprec .eq. 3 ) nbox=240/3
        if( iprec .eq. 4 ) nbox=480/3
        if( iprec .eq. 5 ) nbox=700/3
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        call prinf('nbox=*',nbox,1)
c
c
c     create oct-tree data structure
c
        ntot = 100*(nsource+ntarget)+10000
        if( iprec .eq. -2 ) ntot = ntot * 1.5*1.5*1.5
        if( iprec .eq. -1 ) ntot = ntot * 1.5*1.5
        do ii = 1,10
           allocate (wlists(ntot))
           call lfmm3dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
              goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return          
        endif
c
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
c       
        call prin2('scale=*',scale,nlev+1)
c       
c       
c       carve up workspace further
c
c     itriaflatsort is pointer for sorted triangle coordinates
c     itrianormsort is pointer for sorted triangle normals
c     isourcesort is pointer for sorted triangle centroids
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
c
        itrianglesort = lused7 
        ltrianglesort = 3*3*nsource
        itrinormsort = itrianglesort + ltrianglesort
        ltrinormsort = 3*nsource

        isourcesort = itrinormsort + ltrinormsort 
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget

        isigma_slsort = itargetsort+ltargetsort
        if (ifsingle.eq.1) then
          lsigma_slsort = 3*nsource
        else
          lsigma_slsort = 3
        endif
        isigma_dlsort = isigma_slsort+lsigma_slsort
        if (ifdouble.ge.1) then
          lsigma_dlsort = 3*nsource
        else
          lsigma_dlsort = 3
        endif
        isigma_dvsort = isigma_dlsort+lsigma_dlsort
        if (ifdouble.ge.1) then
          lsigma_dvsort = 3*nsource
        else
          lsigma_dvsort = 3
        endif

        lused7 = isigma_dvsort+lsigma_dvsort
c
c
c       ... allocate the potential and field arrays
c
c
        ipot = lused7 
        if( ifpot .eq. 1) then
        lpot = 2*(3*nsource)
        else
        lpot=6
        endif
        lused7=lused7+lpot
c      
        ipre = lused7 
        if( ifpot .eq. 1) then
        lpre = 2*(nsource)
        else
        lpre=2
        endif
        lused7=lused7+lpre
c      
        igrad = lused7
        if( ifgrad .eq. 1) then
        lgrad = 2*(3*3*nsource)
        else
        lgrad= 2*3*3
        endif
        lused7=lused7+lgrad
c      
        ipottarg = lused7
        if( ifpottarg .eq. 1) then
        lpottarg = 2*(3*ntarget)
        else
        lpottarg=6
        endif
        lused7=lused7+lpottarg
c      
        ipretarg = lused7
        if( ifpottarg .eq. 1) then
        lpretarg = 2*(ntarget)
        else
        lpretarg=2
        endif
        lused7=lused7+lpretarg
c      
        igradtarg = lused7
        if( ifgradtarg .eq. 1) then
        lgradtarg = 2*(3*3*ntarget)
        else
        lgradtarg= 2*3*3
        endif
        lused7=lused7+lgradtarg
c      
c
        call prinf(' lused7 is *',lused7,1)
c
c   
c       ... allocate temporary arrays
c
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                   lused7 is *',lused7,1)
           ier = 8
           return          
        endif
c

        call l3dreordertarg
     $     (nsource,source,wlists(iisource),w(isourcesort))
        if( ifsingle .eq. 1 ) then
        call l3dreordertarg
     $     (nsource,sigma_sl,wlists(iisource),w(isigma_slsort))
        endif
        if( ifdouble .ge. 1 ) then
        call l3dreordertarg
     $     (nsource,sigma_dl,wlists(iisource),w(isigma_dlsort))
        call l3dreordertarg
     $     (nsource,sigma_dv,wlists(iisource),w(isigma_dvsort))
        endif
        
        call l3dreordertria(nsource,wlists(iisource),
     $     triangle,w(itrianglesort),trinorm,w(itrinormsort))
c
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
c
c
c
        call stfmm3dtriatarg0_evalloc(ier,iprec,
     $     w(itrianglesort),w(itrinormsort),
     $     nsource,w(isourcesort),
     $     ifsingle,w(isigma_slsort),
     $     ifdouble,w(isigma_dlsort),w(isigma_dvsort),
     $     ifpot,w(ipot),w(ipre),ifgrad,w(igrad),
     $     ntarget,w(itargetsort),
     $     ifpottarg,w(ipottarg),w(ipretarg),
     $     ifgradtarg,w(igradtarg),
     $     nboxes,laddr,nlev,wlists(iwlists),lwlists)
c
        call prinf('lwlists=*',lwlists,1)
        call prinf('lused total =*',lused7,1)
c       
        call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l3dptsort(nsource,wlists(iisource),w(ipot),pot)
        if(ifpot .eq. 1) 
     $     call l3dftsort(nsource,wlists(iisource),w(ipre),pre)
        if(ifgrad .eq. 1) 
     $     call l3dstsort(nsource,wlists(iisource),w(igrad),grad)
c
        if(ifpottarg .eq. 1 )
     $     call l3dptsort(ntarget,wlists(iitarget),
     $     w(ipottarg),pottarg)
        if(ifpottarg .eq. 1 )
     $     call l3dftsort(ntarget,wlists(iitarget),
     $     w(ipretarg),pretarg)
        if(ifgradtarg .eq. 1) 
     $     call l3dstsort(ntarget,wlists(iitarget),
     $     w(igradtarg),gradtarg)
c       
        return
        end
c
c
c
c
c
        subroutine l3dptsort(n,isource,potsort,pot)
        implicit real *8 (a-h,o-z)
        integer isource(1)
        real *8 pot(3,1),potsort(3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        do m=1,3
ccc        pot(m,isource(i))=potsort(m,i)
        pot(m,isource(i))=pot(m,isource(i))+potsort(m,i)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l3dftsort(n,isource,presort,pre)
        implicit real *8 (a-h,o-z)
        integer isource(1)
        real *8 pre(1),presort(1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
ccc        pre(isource(i))=presort(i)
        pre(isource(i))=pre(isource(i))+presort(i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine l3dstsort(n,isource,gradsort,grad)
        implicit real *8 (a-h,o-z)
        integer isource(1)
        real *8 grad(3,3,1),gradsort(3,3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        do j=1,3
        do m=1,3
ccc        grad(m,j,isource(i))=gradsort(m,j,i)
        grad(m,j,isource(i))=
     $     grad(m,j,isource(i))+gradsort(m,j,i)
        enddo
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine stfmm3dtriatarg0_evalloc(ier,iprec,
     $     triangle,trinorm,
     $     nsource,source,
     $     ifsingle,sigma_sl,
     $     ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntarget,target,
     $     ifpottarg,pottarg,pretarg,ifgradtarg,gradtarg,
     $     nboxes,laddr,nlev,wlists,lwlists)
        implicit real *8 (a-h,o-z)
        real *8 triangle(3,3,1),trinorm(3,1)
        real *8 source(3,1)
        real *8 sigma_sl(3,1)
        real *8 sigma_dl(3,1)
        real *8 sigma_dv(3,1)
        real *8 pot(3,1)
        real *8 pre(1)
        real *8 grad(3,3,1)
        real *8 target(3,1)
        real *8 pottarg(3,1)
        real *8 pretarg(1)
        real *8 gradtarg(3,3,1)
c
        integer laddr(2,200)
        integer list(10 000)
c
        real *8 timeinfo(10)
c
        real *8 wlists(1)
c
        integer box(20)
        real *8 center0(3),corners0(3,8)
        integer box1(20)
        real *8 center1(3),corners1(3,8)
c
ccc        save
c
c     
c       ... set the velocity/pressure and gradient to zero
c
        do i=1,nsource
        if (ifpot .eq. 1) then
        pot(1,i)=0
        pot(2,i)=0
        pot(3,i)=0
        pre(i)=0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,i)=0
        grad(2,1,i)=0
        grad(3,1,i)=0
        grad(1,2,i)=0
        grad(2,2,i)=0
        grad(3,2,i)=0
        grad(1,3,i)=0
        grad(2,3,i)=0
        grad(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if (ifpottarg .eq. 1) then
        pottarg(1,i)=0
        pottarg(2,i)=0
        pottarg(3,i)=0
        pretarg(i)=0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,i)=0
        gradtarg(2,1,i)=0
        gradtarg(3,1,i)=0
        gradtarg(1,2,i)=0
        gradtarg(2,2,i)=0
        gradtarg(3,2,i)=0
        gradtarg(1,3,i)=0
        gradtarg(2,3,i)=0
        gradtarg(3,3,i)=0
        endif
        enddo
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
        call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        ifprint=0
        if (ifprint .eq. 1) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .eq. 1) then
               call prinf('npts=*',npts,1)
            endif
        endif
c
c
        if (nkids .eq. 0 ) then
c
c       ... evaluate self interactions
c
        call stfmm3dtria_direct_self(box,
     $     triangle,trinorm,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .eq. 1) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c    
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c
               call stfmm3dtria_direct(box1,box,
     $            triangle,trinorm,source,
     $            ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $            ifpot,pot,pre,ifgrad,grad,
     $            target,ifpottarg,pottarg,pretarg,
     $            ifgradtarg,gradtarg)
c
 6203           continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
c
ccc        call prin2('inside fmm, pot=*',pot,3*nsource)
ccc        call prin2('inside fmm, pottarg=*',pottarg,3*ntarget)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        call prin2('timeinfo=*',timeinfo,8)
c       
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
c       
        return
        end
c
c
c
c
c
        subroutine stfmm3dtria_direct(box,box1,
     $     triangle,trinorm,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        real *8 triangle(3,3,1),trinorm(3,1),source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8 pot0(3),grad0(3,3)
c
c       ... sources
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(14),box1(14)+box1(15)-1
        if (ifsingle .eq. 1 ) then
        call supst3triadirecttarg
     $     (box(15),triangle(1,1,box(14)),
     $     sigma_sl(1,box(14)),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call stpst3triadirecttarg
     $     (ifdouble,box(15),triangle(1,1,box(14)),
     $     sigma_dl(1,box(14)),sigma_dv(1,box(14)),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(16),box1(16)+box1(17)-1
        if (ifsingle .eq. 1 ) then
        call supst3triadirecttarg
     $     (box(15),triangle(1,1,box(14)),
     $     sigma_sl(1,box(14)),
     1     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call stpst3triadirecttarg
     $     (ifdouble,box(15),triangle(1,1,box(14)),
     $     sigma_dl(1,box(14)),sigma_dv(1,box(14)),
     1     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine stfmm3dtria_direct_self(box,
     $     triangle,trinorm,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,
     $     ifpot,pot,pre,ifgrad,grad,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        real *8 triangle(3,3,1),trinorm(3,1),source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8 pot0(3),grad0(3,3)
c
c       ... sources
c
        ione=1
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(14),box(14)+box(15)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifsingle .eq. 1 ) then
        if( i .eq. j ) then
        call supst3triadirectself
     $     (ione,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        else
        call supst3triadirecttarg
     $     (ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        endif
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        if( i .eq. j ) then 
        call stpst3triadirectself
     $     (ifdouble,ione,ione,triangle(1,1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        else
        call stpst3triadirecttarg
     $     (ifdouble,ione,triangle(1,1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        endif
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c       
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,pot0,pre0,grad0)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(16),box(16)+box(17)-1
        if (ifsingle .eq. 1 ) then
        call supst3triadirecttarg
     $     (box(15),triangle(1,1,box(14)),
     $     sigma_sl(1,box(14)),
     1     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call stpst3triadirecttarg
     $     (ifdouble,box(15),triangle(1,1,box(14)),
     $     sigma_dl(1,box(14)),sigma_dv(1,box(14)),
     1     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine st3dtriadirect(
     $     triangle,trinorm,nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifpot,pot,pre,ifgrad,grad,ntarget,
     $     target,ifpottarg,pottarg,pretarg,
     $     ifgradtarg,gradtarg)
        implicit real *8 (a-h,o-z)
c
c
c     Stokes interactions in R^3: evaluate all pairwise triangle
c     interactions and interactions with targets using the direct
c     O(N^2) algorithm.
c
c     \delta u = \grad p, div u = 0, mu = 1.
c
c     Free space Stokes Green's functions:
c
c       ifsingle=1, stokeslet, f = sigma_sl
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       ifdouble=1, double layer stresslet (type 1), g = sigma_dl, n = sigma_dv
c       u_i = [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=2, symmetric stresslet (type 2), g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c       ifdouble=3, rotlet, g = sigma_dl, n = sigma_dv
c       u_i = [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 0
c
c       ifdouble=4, doublet = symmetric stresslet (type 2) + rotlet, 
c                   g = sigma_dl, n = sigma_dv
c       u_i = [-r_i /r^3] n_j g_j + [3 r_i r_j r_k / r^5] n_k g_j 
c             + [r_j n_j /r^3] g_i - [r_j g_j/ r^3] n_i
c       p = 2 [-n_j g_j / r^3 + 3 r_k n_k r_j g_j / r^5 ]
c
c     without the (1/4 pi) scaling.
c
c
c       INPUT:
c
c       triangle(3,3,nsource) - array of triangles in standard format
c       trinorm(3,nsource) - array of triangle normals
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifpot - velocity computation flag
c       ifgrad - gradient computation flag
c       ifpottarg - target velocity computation flag
c       ifgradtarg - target gradient computation flag
c
c       OUTPUT:
c
c       pot(3,nsource) - velocity at source locations
c       pre(3,nsource) - pressure at source locations
c       grad(3,3,nsource) - gradient at source locations
c       pottarg(3,ntarget) - velocity at target locations
c       pretarg(3,ntarget) - pressure at target locations
c       gradtarg(3,3,ntarget) - gradient at target locations
c
c
        real *8 triangle(3,3,1),trinorm(3,1),source(3,1)
        real *8 sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        real *8 target(3,1)
c
        real *8 pot(3,1),pre(1),grad(3,3,1)
        real *8 pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        real *8 pot0(3),grad0(3,3)
c
c
c     NOTE: In the present version, dipole vectors SIGMA_DV must be SET EQUAL
c     to the triangle normal for stp3d direct routines
c
        do i=1,nsource
        if( ifpot .eq. 1) then
           pot(1,i)=0
           pot(2,i)=0
           pot(3,i)=0
           pre(i)=0
        endif
        if( ifgrad .eq. 1) then
           grad(1,1,i)=0
           grad(2,1,i)=0
           grad(3,1,i)=0
           grad(1,2,i)=0
           grad(2,2,i)=0
           grad(3,2,i)=0
           grad(1,3,i)=0
           grad(2,3,i)=0
           grad(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) then
           pottarg(1,i)=0
           pottarg(2,i)=0
           pottarg(3,i)=0
           pretarg(i)=0
        endif
        if( ifgradtarg .eq. 1) then
           gradtarg(1,1,i)=0
           gradtarg(2,1,i)=0
           gradtarg(3,1,i)=0
           gradtarg(1,2,i)=0
           gradtarg(2,2,i)=0
           gradtarg(3,2,i)=0
           gradtarg(1,3,i)=0
           gradtarg(2,3,i)=0
           gradtarg(3,3,i)=0
        endif
        enddo
c
c       ... sources
c
        ione=1
c
        if( ifpot .eq. 1 .or. ifgrad .eq. 1 ) then
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,pot0,pre0,grad0)
        do j=1,nsource
        do i=1,nsource
c
        if (ifsingle .eq. 1 ) then
        if( i .eq. j ) then
        call supst3triadirectself
     $     (ione,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        else
        call supst3triadirecttarg
     $     (ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        endif
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        if( i .eq. j ) then 
        call stpst3triadirectself
     $     (ifdouble,ione,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        else
        call stpst3triadirecttarg
     $     (ifdouble,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),pot0,pre0,ifgrad,grad0)
        endif
        if (ifpot .eq. 1) then
        pot(1,j)=pot(1,j)+pot0(1)
        pot(2,j)=pot(2,j)+pot0(2)
        pot(3,j)=pot(3,j)+pot0(3)
        pre(j)=pre(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad(1,1,j)=grad(1,1,j)+grad0(1,1)
        grad(2,1,j)=grad(2,1,j)+grad0(2,1)
        grad(3,1,j)=grad(3,1,j)+grad0(3,1)
        grad(1,2,j)=grad(1,2,j)+grad0(1,2)
        grad(2,2,j)=grad(2,2,j)+grad0(2,2)
        grad(3,2,j)=grad(3,2,j)+grad0(3,2)
        grad(1,3,j)=grad(1,3,j)+grad0(1,3)
        grad(2,3,j)=grad(2,3,j)+grad0(2,3)
        grad(3,3,j)=grad(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
C$OMP END PARALLEL DO
c
        endif
c
c       ... targets
c
        if( ifpottarg .eq. 1 .or. ifgradtarg .eq. 1 ) then
c       
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,pot0,pre0,grad0)
        do j=1,ntarget
        do i=1,nsource
        if (ifsingle .eq. 1 ) then
        call supst3triadirecttarg
     $     (ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        if (ifdouble .ge. 1) then
        call stpst3triadirecttarg
     $     (ifdouble,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pottarg(1,j)=pottarg(1,j)+pot0(1)
        pottarg(2,j)=pottarg(2,j)+pot0(2)
        pottarg(3,j)=pottarg(3,j)+pot0(3)
        pretarg(j)=pretarg(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        gradtarg(1,1,j)=gradtarg(1,1,j)+grad0(1,1)
        gradtarg(2,1,j)=gradtarg(2,1,j)+grad0(2,1)
        gradtarg(3,1,j)=gradtarg(3,1,j)+grad0(3,1)
        gradtarg(1,2,j)=gradtarg(1,2,j)+grad0(1,2)
        gradtarg(2,2,j)=gradtarg(2,2,j)+grad0(2,2)
        gradtarg(3,2,j)=gradtarg(3,2,j)+grad0(3,2)
        gradtarg(1,3,j)=gradtarg(1,3,j)+grad0(1,3)
        gradtarg(2,3,j)=gradtarg(2,3,j)+grad0(2,3)
        gradtarg(3,3,j)=gradtarg(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
C$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
c***********************************************************************
c
c       Quadrature routines for Stokes single and double layers
c       Constant densities on flat triangles
c
c***********************************************************************
c
c
        subroutine supst3triadirecttarg_one
     $     (triangle,sigma_sl,ifself,target,
     $     pot,pre,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c
c     Direct evaluation of velocity and gradient due to constant
c     single layer Stokes kernel on a flat triangle.
c
c     Double layer Stokes kernel: constant-density on flat triangles
c
c     Computes velocity and gradient at arbitrary point TARGET
c     due to piecewise-constant single layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri                 number of triangles
c     sigma_sl(3)          SLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     pot(3)            velocity at TARGET
c     pre               pressure at TARGET
c     grad(3,3)         gradient at TARGET
c
c
c
        real *8 triangle(3,3),sigma_sl(3)
        real *8 pot(3),grad(3,3)
        real *8 w(20),vert1(3),vert2(3),vert3(3)
        real *8 vectout(3),vertout(3), ders(3,3)
        real *8 rtable(0:2,0:2,0:2),btable(0:4,0:4,0:4)
c
        do i=1,3
        pot(i)=0
        enddo
        pre=0
c
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
        call tri_for_vect(w,sigma_sl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0
c
        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif


        if( ifgrad .eq. 0 ) maxb=2
        if( ifgrad .eq. 0 ) maxr=0
        if( ifgrad .eq. 1 ) maxb=3
        if( ifgrad .eq. 1 ) maxr=1

c       ... needed for pressure calculation
        maxr=1

        do k=0,maxb
        do i=0,maxb
        do j=0,maxb
        if( i+j+k .le. maxb ) then
        call triabtable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        btable(i,j,k)=d
        endif
        enddo
        enddo
        enddo
        
        do k=0,maxr
        do i=0,maxr
        do j=0,maxr
        if( i+j+k .le. maxr ) then
        call triartable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        rtable(i,j,k)=d
        endif
        enddo
        enddo
        enddo


        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(2,0,0)
        valx = -d

        d=btable(1,1,0)
        valy = -d

        d=btable(1,0,1)
        valz = -d

        d=rtable(0,0,0)
        valx = valx+d*2

        if( ifgrad .eq. 1 ) then
        d=btable(3,0,0)
        valxx = -d

        d=btable(2,1,0)
        valxy = -d

        d=btable(2,0,1)
        valxz = -d

        d=btable(2,1,0)
        valyx = -d

        d=btable(1,2,0)
        valyy = -d

        d=btable(1,1,1)
        valyz = -d

        d=btable(2,0,1)
        valzx = -d

        d=btable(1,1,1)
        valzy = -d

        d=btable(1,0,2)
        valzz = -d

        d=rtable(1,0,0)
        valxx = valxx+d*2

        d=rtable(0,1,0)
        valxy = valxy+d*2

        d=rtable(0,0,1)
        valxz = valxz+d*2
        endif

        endif

        if ( k .eq. 2 ) then
        
        d=btable(1,1,0)
        valx = -d

        d=btable(0,2,0)
        valy = -d

        d=btable(0,1,1)
        valz = -d

        d=rtable(0,0,0)
        valy = valy+d*2

        if( ifgrad .eq. 1 ) then
        d=btable(2,1,0)
        valxx = -d

        d=btable(1,2,0)
        valxy = -d

        d=btable(1,1,1)
        valxz = -d

        d=btable(1,2,0)
        valyx = -d

        d=btable(0,3,0)
        valyy = -d

        d=btable(0,2,1)
        valyz = -d

        d=btable(1,1,1)
        valzx = -d

        d=btable(0,2,1)
        valzy = -d

        d=btable(0,1,2)
        valzz = -d

        d=rtable(1,0,0)
        valyx = valyx+d*2

        d=rtable(0,1,0)
        valyy = valyy+d*2

        d=rtable(0,0,1)
        valyz = valyz+d*2
        endif

        endif

        if ( k .eq. 3 ) then
        
        d=btable(1,0,1)
        valx = -d

        d=btable(0,1,1)
        valy = -d

        d=btable(0,0,2)
        valz = -d

        d=rtable(0,0,0)
        valz = valz+d*2

        if( ifgrad .eq. 1 ) then
        d=btable(2,0,1)
        valxx = -d

        d=btable(1,1,1)
        valxy = -d

        d=btable(1,0,2)
        valxz = -d

        d=btable(1,1,1)
        valyx = -d

        d=btable(0,2,1)
        valyy = -d

        d=btable(0,1,2)
        valyz = -d

        d=btable(1,0,2)
        valzx = -d

        d=btable(0,1,2)
        valzy = -d

        d=btable(0,0,3)
        valzz = -d

        d=rtable(1,0,0)
        valzx = valzx+d*2

        d=rtable(0,1,0)
        valzy = valzy+d*2

        d=rtable(0,0,1)
        valzz = valzz+d*2
        endif

        endif

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        pot(1)=pot(1)+vectout(k)*derx
        pot(2)=pot(2)+vectout(k)*dery
        pot(3)=pot(3)+vectout(k)*derz
c
c
        if( ifgrad .eq. 1 ) then
c
        if( 1 .eq. 2 ) then
c       ... symmetrize the derivative matrix, prepare to compute grad
c        
        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        grad(1,1)=grad(1,1)+vectout(k)*derxx
        grad(1,2)=grad(1,2)+vectout(k)*derxy
        grad(1,3)=grad(1,3)+vectout(k)*derxz
        grad(2,2)=grad(2,2)+vectout(k)*deryy
        grad(2,3)=grad(2,3)+vectout(k)*deryz
        grad(3,3)=grad(3,3)+vectout(k)*derzz
        
        grad(2,1)=grad(1,2)
        grad(3,1)=grad(1,3)
        grad(3,2)=grad(2,3)
        endif

        if( 2 .eq. 2 ) then
c       ... all derivatives
c        
        call rothess3d_arb(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,valyx,valzx,valzy,
     $      derxx,deryy,derzz,derxy,derxz,deryz,deryx,derzx,derzy)

        grad(1,1)=grad(1,1)+vectout(k)*derxx
        grad(1,2)=grad(1,2)+vectout(k)*derxy
        grad(1,3)=grad(1,3)+vectout(k)*derxz
        grad(2,2)=grad(2,2)+vectout(k)*deryy
        grad(2,3)=grad(2,3)+vectout(k)*deryz
        grad(3,3)=grad(3,3)+vectout(k)*derzz
        
        grad(2,1)=grad(2,1)+vectout(k)*deryx
        grad(3,1)=grad(3,1)+vectout(k)*derzx
        grad(3,2)=grad(3,2)+vectout(k)*derzy
        endif

        endif

        enddo
c
c
c       ... pre is scalar, no rotation 
        valx=rtable(1,0,0)
        valy=rtable(0,1,0)
        valz=rtable(0,0,1)
c
        pre=pre+valx*vectout(1)
        pre=pre+valy*vectout(2)
        pre=pre+valz*vectout(3)
c
c
        do i=1,3
        pot(i)=pot(i)/2
        enddo
c
c
        if( ifgrad .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        grad(i,j)=-grad(i,j)
        enddo
        enddo
c
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)/2
        enddo
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine supst3triadirecttarg
     $     (ntri,triangles,sigma_sl,
     1     target,pot,pre,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c     Single layer Stokes kernel: constant-densities on flat triangles
c
c     Computes velocity and gradient at arbitrary point TARGET not lying
c     on the surface due to piecewise-constant single layer density on
c     collection of triangles.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri                 number of triangles
c     sigma_sl(3,ntri)     array of SLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     target(3)            target location
c
c     OUTPUT:
c
c     pot(3)            velocity at TARGET
c     pre               pressure at TARGET
c     grad(3,3)         gradient at TARGET
c
c
        real *8 triangles(3,3,1),sigma_sl(3,1)
        real *8 pot0(3),grad0(3,3)
        real *8 pot(3),grad(3,3),target(3)


        do i=1,3
        pot(i)=0
        enddo
        pre=0
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo

        do k=1,ntri

        ifself=0
        call supst3triadirecttarg_one
     $     (triangles(1,1,k),sigma_sl(1,k),
     1     ifself,target,pot0,pre0,ifgrad,grad0)

        do i=1,3
        pot(i)=pot(i)+pot0(i)
        enddo
        pre=pre+pre0
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)+grad0(i,j)
        enddo
        enddo

        enddo

        return
        end
c
c
c
c
c
        subroutine supst3triadirectself
     $     (ipatch,ntri,triangles,sigma_sl,
     1     zparts,pot,pre,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c     Single layer Stokes kernel: constant-densities on flat triangles
c
c     Computes velocity and gradient at centroid zparts(*,ipatch) 
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles, numbered jpatch = 1,...,ntri. 
c
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c
c     INPUT:
c
c     ntri                 number of triangles
c     sigma_sl(3,ntri)     array of SLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     zparts(3,1)          array of triangle centroids
c
c     OUTPUT:
c
c     pot(3)            velocity at centroid zparts(*,ipatch)
c     pre               pressure at centroid zparts(*,ipatch)
c     grad(3,3)         gradient at centroid zparts(*,ipatch)
c
c
        real *8 triangles(3,3,1),sigma_sl(3,1),trinorm(3,1)
        real *8 pot0(3),grad0(3,3)
        real *8 pot(3),grad(3,3),zparts(3,1)
c
        do i=1,3
        pot(i)=0
        enddo
        pre=0
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo

        do k=1,ntri

        if( k .eq. ipatch ) ifself=1
        if( k .ne. ipatch ) ifself=0

        call supst3triadirecttarg_one
     $     (triangles(1,1,k),sigma_sl(1,k),
     1     ifself,zparts(1,ipatch),pot0,pre0,ifgrad,grad0)

        do i=1,3
        pot(i)=pot(i)+pot0(i)
        enddo
        pre=pre+pre0
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)+grad0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
        subroutine stpst3triadirecttarg_one
     $     (ifdouble,triangle,sigma_dl,trinorm,
     $     ifself,target,pot,pre,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c
c     Direct evaluation of velocity and gradient due to constant
c     double layer Stokes kernel on a flat triangle.
c
c     Double layer Stokes kernel: constant-density on flat triangles
c
c     Computes velocity and gradient at arbitrary point TARGET 
c     due to piecewise-constant double layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri                 number of triangles
c     sigma_dl(3)          DLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     trianorm(3)          triangle normal
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     pot(3)            velocity at TARGET
c     pre               pressure at TARGET
c     grad(3,3)         gradient at TARGET
c
c
c
        real *8 triangle(3,3),sigma_dl(3),trinorm(3)
        real *8 pot(3),grad(3,3)
        real *8 w(20),vert1(3),vert2(3),vert3(3)
        real *8 vectout(3),vertout(3), ders(3,3)
        real *8 tmatr(3,3),umatr(3,3,3)
        real *8 smatr(3,3,3),dmatr(3,3,3,3)
        real *8 rtable(0:2,0:2,0:2),btable(0:4,0:4,0:4)
c
        do i=1,3
        pot(i)=0
        enddo
        pre=0
c
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
ccc        call tri_for_vect(w,trinorm,vectout)
ccc        call prin2('inside stpst3triadirecttarg, trinorm=*',vectout,3)

        call tri_for_vect(w,sigma_dl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0

        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3
        do i=1,3
        do j=1,3
        umatr(i,j,k)=0
        enddo
        enddo
        enddo
        
        do m=1,3
        do k=1,3
        do i=1,3
        do j=1,3
        dmatr(i,j,k,m)=0
        enddo
        enddo
        enddo
        enddo


        if( ifgrad .eq. 0 ) maxb=3
        if( ifgrad .eq. 0 ) maxr=1
        if( ifgrad .eq. 1 ) maxb=4
        if( ifgrad .eq. 1 ) maxr=2

        do k=0,maxb
        do i=0,maxb
        do j=0,maxb
        if( i+j+k .le. maxb ) then
        call triabtable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        btable(i,j,k)=d
        endif
        enddo
        enddo
        enddo
        
        do k=0,maxr
        do i=0,maxr
        do j=0,maxr
        if( i+j+k .le. maxr ) then
        call triartable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        rtable(i,j,k)=d
        endif
        enddo
        enddo
        enddo


        
        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(3,0,0)
        valxx = -d

        d=btable(2,1,0)
        valxy = -d

        d=btable(2,0,1)
        valxz = -d

        d=btable(2,1,0)
        valyx = -d

        d=btable(1,2,0)
        valyy = -d

        d=btable(1,1,1)
        valyz = -d

        d=btable(2,0,1)
        valzx = -d

        d=btable(1,1,1)
        valzy = -d

        d=btable(1,0,2)
        valzz = -d

        d=rtable(1,0,0)
        valxx = valxx+d*2

        d=rtable(0,1,0)
        valxy = valxy+d*2

        d=rtable(0,0,1)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        d=btable(2,1,0)
        valxx = -d

        d=btable(1,2,0)
        valxy = -d

        d=btable(1,1,1)
        valxz = -d

        d=btable(1,2,0)
        valyx = -d

        d=btable(0,3,0)
        valyy = -d

        d=btable(0,2,1)
        valyz = -d

        d=btable(1,1,1)
        valzx = -d

        d=btable(0,2,1)
        valzy = -d

        d=btable(0,1,2)
        valzz = -d

        d=rtable(1,0,0)
        valyx = valyx+d*2

        d=rtable(0,1,0)
        valyy = valyy+d*2

        d=rtable(0,0,1)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        d=btable(2,0,1)
        valxx = -d

        d=btable(1,1,1)
        valxy = -d

        d=btable(1,0,2)
        valxz = -d

        d=btable(1,1,1)
        valyx = -d

        d=btable(0,2,1)
        valyy = -d

        d=btable(0,1,2)
        valyz = -d

        d=btable(1,0,2)
        valzx = -d

        d=btable(0,1,2)
        valzy = -d

        d=btable(0,0,3)
        valzz = -d

        d=rtable(1,0,0)
        valzx = valzx+d*2

        d=rtable(0,1,0)
        valzy = valzy+d*2

        d=rtable(0,0,1)
        valzz = valzz+d*2

        endif

        umatr(1,k,1)=valxx
        umatr(1,k,2)=valxy
        umatr(1,k,3)=valxz
        umatr(2,k,1)=valyx
        umatr(2,k,2)=valyy
        umatr(2,k,3)=valyz
        umatr(3,k,1)=valzx
        umatr(3,k,2)=valzy
        umatr(3,k,3)=valzz
c
        enddo
c
c
        if( ifgrad .eq. 1 ) then
c       ... and the derivatives 
c
        do m=1,3
c
        if( m .eq. 1 ) then
        ix=1
        iy=0
        iz=0
        endif
        if( m .eq. 2 ) then
        ix=0
        iy=1
        iz=0
        endif
        if( m .eq. 3 ) then
        ix=0
        iy=0
        iz=1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(3+ix,0+iy,0+iz)
        valxx = -d

        d=btable(2+ix,1+iy,0+iz)
        valxy = -d

        d=btable(2+ix,0+iy,1+iz)
        valxz = -d

        d=btable(2+ix,1+iy,0+iz)
        valyx = -d

        d=btable(1+ix,2+iy,0+iz)
        valyy = -d

        d=btable(1+ix,1+iy,1+iz)
        valyz = -d

        d=btable(2+ix,0+iy,1+iz)
        valzx = -d

        d=btable(1+ix,1+iy,1+iz)
        valzy = -d

        d=btable(1+ix,0+iy,2+iz)
        valzz = -d

        d=rtable(1+ix,0+iy,0+iz)
        valxx = valxx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valxy = valxy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        d=btable(2+ix,1+iy,0+iz)
        valxx = -d

        d=btable(1+ix,2+iy,0+iz)
        valxy = -d

        d=btable(1+ix,1+iy,1+iz)
        valxz = -d
        
        d=btable(1+ix,2+iy,0+iz)
        valyx = -d

        d=btable(0+ix,3+iy,0+iz)
        valyy = -d

        d=btable(0+ix,2+iy,1+iz)
        valyz = -d

        d=btable(1+ix,1+iy,1+iz)
        valzx = -d

        d=btable(0+ix,2+iy,1+iz)
        valzy = -d

        d=btable(0+ix,1+iy,2+iz)
        valzz = -d

        d=rtable(1+ix,0+iy,0+iz)
        valyx = valyx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valyy = valyy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        d=btable(2+ix,0+iy,1+iz)
        valxx = -d

        d=btable(1+ix,1+iy,1+iz)
        valxy = -d

        d=btable(1+ix,0+iy,2+iz)
        valxz = -d

        d=btable(1+ix,1+iy,1+iz)
        valyx = -d

        d=btable(0+ix,2+iy,1+iz)
        valyy = -d

        d=btable(0+ix,1+iy,2+iz)
        valyz = -d

        d=btable(1+ix,0+iy,2+iz)
        valzx = -d

        d=btable(0+ix,1+iy,2+iz)
        valzy = -d

        d=btable(0+ix,0+iy,3+iz)
        valzz = -d

        d=rtable(1+ix,0+iy,0+iz)
        valzx = valzx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valzy = valzy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valzz = valzz+d*2

        endif

        dmatr(1,k,1,m)=valxx
        dmatr(1,k,2,m)=valxy
        dmatr(1,k,3,m)=valxz
        dmatr(2,k,1,m)=valyx
        dmatr(2,k,2,m)=valyy
        dmatr(2,k,3,m)=valyz
        dmatr(3,k,1,m)=valzx
        dmatr(3,k,2,m)=valzy
        dmatr(3,k,3,m)=valzz
c
        enddo
        enddo
c
ccc        call prin2('umatr=*',umatr,3*3*3)
c
        endif
c
c
        do i=1,3
        do j=1,3
        tmatr(i,j)=0
        enddo
        enddo
c
        if( ifdouble .eq. 1 ) then
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+(umatr(j,i,3)+umatr(3,i,j))/2
        enddo
        enddo
c
        d=rtable(1,0,0)
        tmatr(1,3) = tmatr(1,3)+d
        d=rtable(0,1,0)
        tmatr(2,3) = tmatr(2,3)+d
        d=rtable(0,0,1)
        tmatr(3,3) = tmatr(3,3)+d
        endif

        if( ifdouble .eq. 2 ) then
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+(umatr(j,i,3)+umatr(3,i,j))/2
        enddo
        enddo
        endif
c
        if( ifdouble .eq. 3 ) then
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+(umatr(j,i,3)-umatr(3,i,j))/2
        enddo
        enddo
        endif
c
        if( ifdouble .eq. 4 ) then
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+umatr(j,i,3)
        enddo
        enddo
        endif
c
c
ccc        call prin2('tmatr=*',tmatr,3*3)
c
        if( ifgrad .eq. 1 ) then
c
        do m=1,3
c
        if( m .eq. 1 ) then
        ix=1
        iy=0
        iz=0
        endif
        if( m .eq. 2 ) then
        ix=0
        iy=1
        iz=0
        endif
        if( m .eq. 3 ) then
        ix=0
        iy=0
        iz=1
        endif
c
        do i=1,3
        do j=1,3
        smatr(i,j,m)=0
        enddo
        enddo
c
        if( ifdouble .eq. 1 ) then
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+(dmatr(j,i,3,m)+dmatr(3,i,j,m))/2
        enddo
        enddo
c
        d=rtable(1+ix,0+iy,0+iz)
        smatr(1,3,m) = smatr(1,3,m)+d
        d=rtable(0+ix,1+iy,0+iz)
        smatr(2,3,m) = smatr(2,3,m)+d
        d=rtable(0+ix,0+iy,1+iz)
        smatr(3,3,m) = smatr(3,3,m)+d
        endif

        if( ifdouble .eq. 2 ) then
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+(dmatr(j,i,3,m)+dmatr(3,i,j,m))/2
        enddo
        enddo
        endif

        if( ifdouble .eq. 3 ) then
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+(dmatr(j,i,3,m)-dmatr(3,i,j,m))/2
        enddo
        enddo
        endif

        if( ifdouble .eq. 4 ) then
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+dmatr(j,i,3,m)
        enddo
        enddo
        endif

        enddo
c
ccc        call prin2('smatr=*',tmatr,3*3*3)
c
        endif
c
        do k=1,3

        valx=tmatr(1,k)
        valy=tmatr(2,k)
        valz=tmatr(3,k)

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        pot(1)=pot(1)+vectout(k)*derx
        pot(2)=pot(2)+vectout(k)*dery
        pot(3)=pot(3)+vectout(k)*derz
c
c
        if( ifgrad .eq. 1 ) then
c
        if( 1 .eq. 2 ) then
c       ... symmetrize the derivative matrix, prepare to compute grad
c        
        valxx=smatr(1,k,1)
        valyx=smatr(2,k,1)
        valzx=smatr(3,k,1)
        valxy=smatr(1,k,2)
        valyy=smatr(2,k,2)
        valzy=smatr(3,k,2)
        valxz=smatr(1,k,3)
        valyz=smatr(2,k,3)
        valzz=smatr(3,k,3)

        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        grad(1,1)=grad(1,1)+vectout(k)*derxx
        grad(1,2)=grad(1,2)+vectout(k)*derxy
        grad(1,3)=grad(1,3)+vectout(k)*derxz
        grad(2,2)=grad(2,2)+vectout(k)*deryy
        grad(2,3)=grad(2,3)+vectout(k)*deryz
        grad(3,3)=grad(3,3)+vectout(k)*derzz
        
        grad(2,1)=grad(1,2) 
        grad(3,1)=grad(1,3)
        grad(3,2)=grad(2,3)
c
        endif

        if( 2 .eq. 2 ) then
c       ... all derivatives
c        
        valxx=smatr(1,k,1)
        valyx=smatr(2,k,1)
        valzx=smatr(3,k,1)
        valxy=smatr(1,k,2)
        valyy=smatr(2,k,2)
        valzy=smatr(3,k,2)
        valxz=smatr(1,k,3)
        valyz=smatr(2,k,3)
        valzz=smatr(3,k,3)

        call rothess3d_arb(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,valyx,valzx,valzy,
     $      derxx,deryy,derzz,derxy,derxz,deryz,deryx,derzx,derzy)

        grad(1,1)=grad(1,1)+vectout(k)*derxx
        grad(1,2)=grad(1,2)+vectout(k)*derxy
        grad(1,3)=grad(1,3)+vectout(k)*derxz
        grad(2,2)=grad(2,2)+vectout(k)*deryy
        grad(2,3)=grad(2,3)+vectout(k)*deryz
        grad(3,3)=grad(3,3)+vectout(k)*derzz
        
        grad(2,1)=grad(2,1)+vectout(k)*deryx
        grad(3,1)=grad(3,1)+vectout(k)*derzx
        grad(3,2)=grad(3,2)+vectout(k)*derzy
c
        endif

        endif
c
        enddo
c
c
        if( ifdouble .eq. 1 .or. ifdouble. eq. 2 
     $     .or. ifdouble .eq. 4 ) then
c
c       ... pre is scalar, no rotation 
        valxz=rtable(1,0,1)
        valyz=rtable(0,1,1)
        valzz=rtable(0,0,2)
c
        pre=pre+valxz*vectout(1)*2
        pre=pre+valyz*vectout(2)*2
        pre=pre+valzz*vectout(3)*2
c
        endif
c
c
        do i=1,3
        pot(i)=pot(i)
        enddo
c
        if( ifgrad .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        grad(i,j)=-grad(i,j)
        enddo
        enddo
c
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)
        enddo
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine stpst3triadirecttarg
     $     (ifdouble,ntri,triangles,sigma_dl,trinorm,
     1     target,pot,pre,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c     Double layer Stokes kernel: constant-densities on flat triangles
c
c     Computes velocity and gradient at arbitrary point TARGET not lying
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri                 number of triangles
c     sigma_dl(3,ntri)     array of DLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     trianorm(3,ntri)     array of triangle normals
c     target(3)            target location
c
c     OUTPUT:
c
c     pot(3)            velocity at TARGET
c     pre               pressure at TARGET
c     grad(3,3)         gradient at TARGET
c
c
        real *8 triangles(3,3,1),sigma_dl(3,1),trinorm(3,1)
        real *8 pot0(3),grad0(3,3)
        real *8 pot(3),grad(3,3),target(3)
c
        do i=1,3
        pot(i)=0
        enddo
        pre=0
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo

        do k=1,ntri

        ifself=0
        call stpst3triadirecttarg_one
     $     (ifdouble,triangles(1,1,k),sigma_dl(1,k),trinorm(1,k),
     1     ifself,target,pot0,pre0,ifgrad,grad0)

        do i=1,3
        pot(i)=pot(i)+pot0(i)
        enddo
        pre=pre+pre0
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)+grad0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
        subroutine stpst3triadirectself
     $     (ifdouble,ipatch,ntri,triangles,sigma_dl,trinorm,
     1     zparts,pot,pre,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c     Double layer Stokes kernel: constant-densities on flat triangles
c
c     Computes velocity and gradient at centroid zparts(*,ipatch) 
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles, numbered jpatch = 1,...,ntri. 
c
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c
c     INPUT:
c
c     ntri                 number of triangles
c     sigma_dl(3,ntri)     array of DLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     trianorm(3,ntri)     array of triangle normals
c     zparts(3,1)          array of triangle centroids
c
c     OUTPUT:
c
c     pot(3)            velocity at centroid zparts(*,ipatch)
c     pre               pressure at centroid zparts(*,ipatch)
c     grad(3,3)         gradient at centroid zparts(*,ipatch)
c
c
        real *8 triangles(3,3,1),sigma_dl(3,1),trinorm(3,1)
        real *8 pot0(3),grad0(3,3)
        real *8 pot(3),grad(3,3),zparts(3,1)
c
        do i=1,3
        pot(i)=0
        enddo
        pre=0
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo
        
        do k=1,ntri

        if( k .eq. ipatch ) ifself=1
        if( k .ne. ipatch ) ifself=0

        call stpst3triadirecttarg_one
     $     (ifdouble,triangles(1,1,k),sigma_dl(1,k),trinorm(1,k),
     1     ifself,zparts(1,ipatch),pot0,pre0,ifgrad,grad0)

        do i=1,3
        pot(i)=pot(i)+pot0(i)
        enddo
        pre=pre+pre0
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)+grad0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
