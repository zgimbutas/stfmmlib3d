cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-07-05 22:24:41 -0400 (Mon, 05 Jul 2010) $
c    $Revision: 1057 $
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The code is highly specialized, it can handle linear charge,
c     constant dipole strength and linear dipole vector densities on
c     flat triangles.  FAR FIELD CONTRIBUTION ONLY. Input densities are
c     defined on vertices of triangles, output is defined on centroids.
c
c
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Laplace potentials, fields, and
c     hessians on surfaces defined by a collection of positively
c     oriented flat triangles.  (FORTRAN 90 VERSION)
c
c     lfmm3dtrilhess - Laplace FMM in R^3: evaluate all pairwise triangle
c         interactions (including self-interaction)
c
c     lfmm3dtrilhesstarg - Laplace FMM in R^3: evaluate all pairwise
c         triangle interactions (including self-interaction) +
c         interactions with targets
c
c     The routines in this file permit the calculation of SECOND 
c     DERIVATIVES of the potentials. ARBITRARY ORIENTED dipole vectors 
c     are permitted.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine lfmm3dtrilhess(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the harmonic potential, field, and
c       hessian due to a collection of flat triangles with linear
c       single, and constant double layer densities. We use (1/r) for the
c       Green's function, without the (1/4 pi) scaling.
c
c       FAR FIELD CONTRIBUTION ONLY. 
c       Input densities are defined on vertices of triangles.
c       Output is defined on centroids of triangles.
c
c       The main FMM routine permits both evaluation on surface
c       and at a collection of off-surface targets. 
c
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c
        lused7=0
c
        ntarget=0
        ifpottarg=0
        iffldtarg=0
        ifhesstarg=0
c
        call lfmm3dtrilhesstarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtrilhesstarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c       
c
c       This is the principal subroutine for evaluating the harmonic
c       potential, field, and hessian due to a collection of flat
c       triangles with linear single, and constant double layer
c       densities. We use (1/r) for the Green's function, without the
c       (1/4 pi) scaling.  
c
c       FAR FIELD CONTRIBUTION ONLY. 
c       Input densities are defined on vertices of triangles.
c       Output is defined on centroids of triangles.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  
c   
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dtrilhesstargmain.
c
c       NOTE: In this routine, arbitrary oriented dipole vectors are permitted.
c
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c       charge: complex *16 (3,nsource): piecewise linear SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
c       dipstr: complex *16 (3,nsource): piecewise constant DLP strengths
c       dipvec: real *8 (3,3,nsource): piecewise constant dipole orientation 
c                                    vectors. 
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ifhess:  hessian flag (1=compute hessian, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c       ihesstarg:  target hessian flag 
c                   (1=compute hessian, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pot: complex *16 (nsource): potential at triangle centroids
c       fld: complex *16 (3,nsource): field (-gradient) at triangle centroids 
c       hess: complex *16 (6,nsource): hessian at triangle centroids 
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c       hesstarg: complex *16 (6,ntarget): hessian at target locations
c

        dimension triaflat(3,3,1)
        dimension trianorm(3,1)
        dimension source(3,1)
        dimension target(3,1)
        complex *16 charge(3,1)
        complex *16 dipstr(3,1)
        dimension dipvec(3,3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
        complex *16 hesstarg(6,1)        
c
        dimension timeinfo(10)
c       
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 evels of refinment, which is 
c     more than enough for any non-pathological case.
c
        dimension laddr(2,200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:) 
        real *8, allocatable :: wlists(:) 
        real *8, allocatable :: wrmlexp(:) 
        complex *16 ptemp,ftemp(3),htemp(6)
c
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
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
        if (ifprint.eq.1) call prin2('epsfmm=*',epsfmm,1)
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
        if (ifprint.eq.1) call prinf('nbox=*',nbox,1)
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
        if (ifprint.eq.1) call prin2('scale=*',scale,nlev+1)
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
        itriaflatsort = lused7 + 5
        ltriaflatsort = 3*3*nsource
        itrianormsort = itriaflatsort + ltriaflatsort
        ltrianormsort = 3*nsource
        isourcesort = itrianormsort + ltrianormsort 
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource *3
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 3*nsource *3
          ldipstr = 2*nsource *3
        else
          ldipvec = 3
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr       
c
c
c     allocate the potential and field arrays
c
        ipot = lused7
        lpot = 2*nsource
        lused7=lused7+lpot
c       
        ifld = lused7
        if( iffld .eq. 1) then
        lfld = 2*(3*nsource)
        else
        lfld=6
        endif
        lused7=lused7+lfld
c      
        ihess = lused7
        if( ifhess .eq. 1) then
        lhess = 2*(6*nsource)
        else
        lhess=12
        endif
        lused7=lused7+lhess
c      
        ipottarg = lused7
        lpottarg = 2*ntarget
        lused7=lused7+lpottarg
c       
        ifldtarg = lused7
        if( iffldtarg .eq. 1) then
        lfldtarg = 2*(3*ntarget)
        else
        lfldtarg=6
        endif
        lused7=lused7+lfldtarg
c      
        ihesstarg = lused7
        if( ifhesstarg .eq. 1) then
        lhesstarg = 2*(6*ntarget)
        else
        lhesstarg=12
        endif
        lused7=lused7+lhesstarg
c      
        if (ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c            
        nmax = 0
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call l3dterms(epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        nquad=2*nmax        
c       
c     ixnodes is pointer for quadrature nodes
c     iwhts is pointer for quadrature weights
c
        ixnodes = lused7 
        iwts = ixnodes + nquad
        lused7 = iwts + nquad
c
        if (ifprint.eq.1) call prinf('nterms=*',nterms,nlev+1)
        if (ifprint.eq.1) call prinf('nmax=*',nmax,1)
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c   
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7 
        imptemp = iiaddr + 2*nboxes
        lmptemp = (nmax+1)*(2*nmax+1)*2 
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                   lused7 is *',lused7,1)
           ier = 8
           return          
        endif
c
c     reorder triangles, centroids etc. so that each box holds
c     contiguous list of source numbers.
c
        call l3dreorder_linear
     $     (nsource,source,ifcharge,charge,wlists(iisource),
     $     ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 
c       
        call l3dreordertria(nsource,wlists(iisource),
     $     triaflat,w(itriaflatsort),trianorm,w(itrianormsort))
c
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
        if (ifprint.eq.1) then
          call prinf('finished reordering=*',ier,1)
          call prinf('ier=*',ier,1)
          call prinf('nboxes=*',nboxes,1)
          call prinf('nlev=*',nlev,1)
          call prinf('nboxes=*',nboxes,1)
          call prinf('lused7=*',lused7,1)
        endif
c
        ifinit=1
        call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
c
ccc        call prin2('xnodes=*',xnodes,nquad)
ccc        call prin2('wts=*',wts,nquad)
c
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call l3dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if (ifprint.eq.1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                   lused7 is *',lused7,1)
           ier = 16
           return          
        endif
c
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far field and local fields (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
c
c
        ifevalfar=1
        ifevalloc=0
c
        call lfmm3dtrilhesstargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,w(itriaflatsort),w(itrianormsort),
     $     w(isourcesort),w(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),iffld,w(ifld),ifhess,w(ihess),
     $     ntarget,w(itargetsort),w(iitarget),
     $     ifpottarg,w(ipottarg),iffldtarg,w(ifldtarg),
     $     ifhesstarg,w(ihesstarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     w(ixnodes),w(iwts),nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
        if (ifprint.eq.1) then
        call prinf('lwlists=*',lwlists,1)
        call prinf('lused total =*',lused7,1)       
        call prin2('memory / point = *',(lused7)/dble(nsource),1)
        endif
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
        if(iffld .eq. 1) 
     $     call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
        if(ifhess .eq. 1) 
     $     call l3dhsort(nsource,wlists(iisource),w(ihess),hess)
c
        if(ifpottarg .eq. 1 )
     $     call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
        if(iffldtarg .eq. 1) 
     $     call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
        if(ifhesstarg .eq. 1) 
     $     call l3dhsort(ntarget,wlists(iitarget),w(ihesstarg),hesstarg)
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dtrilhesstargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,triaflatsort,trianormsort,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ifhess,hess,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,xnodes,wts,nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension triaflatsort(3,3,1),trianormsort(3,1)
        dimension sourcesort(3,1), isource(1)
        complex *16 chargesort(3,1)
        complex *16 dipstrsort(3,1)
        dimension dipvecsort(3,3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        complex *16 hess(6,1)
        dimension targetsort(3,1), itarget(1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        complex *16 hesstarg(6,1)
        dimension wlists(1)
        dimension iaddr(2,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension xnodes(nquad),wts(nquad)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3),htemp(6)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(40 000)
        dimension nterms_eval(4,0:200)
c
        real *8, allocatable :: scarray_local(:)
        real *8, allocatable :: scarray_mpole(:)
c
        data ima/(0.0d0,1.0d0)/

c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
c
c       ... set the potential and field to zero
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
        endif
        if( ifhess .eq. 1) then
           hess(1,i)=0
           hess(2,i)=0
           hess(3,i)=0
           hess(4,i)=0
           hess(5,i)=0
           hess(6,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        if( ifhesstarg .eq. 1) then
           hesstarg(1,i)=0
           hesstarg(2,i)=0
           hesstarg(3,i)=0
           hesstarg(4,i)=0
           hesstarg(5,i)=0
           hesstarg(6,i)=0
        endif
        enddo
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        norder=1
        nqtri=1
c
        if( iprec .eq. -2 ) then
        norder=1
        nqtri=1
        endif
c
        if( iprec .eq. -1 ) then
        norder=2
        nqtri=2
        endif
c
        if( iprec .eq. 0 ) then
        norder=3
        nqtri=3
        endif
c
        if( iprec .ge. 1 ) then
        norder=6
        nqtri=6
        endif
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c       ... initialize Legendre function evaluation routines
c
        nlege=100
        lw7=40 000
        call ylgndrfwini(nlege,wlege,lw7,lused7)
c
        do i=0,nlev
        do itype=1,4
        call l3dterms_eval(itype,epsfmm,
     1       nterms_eval(itype,i),ier)
        enddo
        enddo
c
        if(ifprint .ge. 2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
        do ibox = 1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        enddo
c
c
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
        do 1300 ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        level=box(1)
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(14)
c        npts=box(15)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(15)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(14)),box(15))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(15) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	    radius = (corners0(1,1) - center0(1))**2
	    radius = radius + (corners0(2,1) - center0(2))**2
	    radius = radius + (corners0(3,1) - center0(3))**2
	    radius = sqrt(radius)
c
 	    if (ifcharge .eq. 1) then
               call l3dformmptris_linear_add(ier,scale(level),
     1  	  triaflatsort(1,1,box(14)),chargesort(1,box(14)),npts,
     2            center0,norder,nterms(level),rmlexp(iaddr(1,ibox)))        
            endif
c
            if (ifdipole .eq. 1 ) then
               call l3dformmptrid_linear_add(ier,scale(level),
     1           triaflatsort(1,1,box(14)),trianormsort(1,box(14)),
     2           dipstrsort(1,box(14)),dipvecsort(1,1,box(14)),
     $           npts,center0,norder,
     3           nterms(level),rmlexp(iaddr(1,ibox)))
            endif

         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
        if(ifprint .ge. 1) 
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and fields directly
c 
         do 3251 ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(level,npts,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            level1=box1(1)
c
            ifdirect3 = 0
c
c            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
c     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
c
               npts=box(15)
c
               if( ifcharge .eq. 1 ) then
                  call l3dformtatris_linear_add(ier,scale(level1),
     1  	    triaflatsort(1,1,box(14)),chargesort(1,box(14)),
     $              npts,center1,
     2              norder,nterms(level1),
     $              rmlexp(iaddr(2,jbox)))
	       endif
               if( ifdipole .eq. 1 ) then
                  call l3dformtatrid_linear_add(ier,scale(level1),
     1               triaflatsort(1,1,box(14)),trianormsort(1,box(14)),
     2               dipstrsort(1,box(14)),dipvecsort(1,1,box(14)),
     $               npts,center1,
     3               norder,nterms(level1),
     $               rmlexp(iaddr(2,jbox)))
               endif
c
            else

            call lfmm3dtrilhess_direct(nqtri,box,box1,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
            endif
         enddo
C$OMP END PARALLEL DO
c
 3251    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c
        if(ifprint .ge. 1) 
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0
        if (ifhess.eq.1) ifprune_list2 = 0
        call lfmm3d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2)
c
c
        allocate( scarray_mpole(0:100000) )
        call l3dmpevalhessdini(nterms(0),scarray_mpole)
c       
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and fields directly
c
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,level,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,htemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box(1)
c
            ifdirect4 = 0
c
c            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
c     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then
               do j=box1(14),box1(14)+box1(15)-1
                  if( ifhess .eq. 1 ) then
c                  call l3dmpevalhess(scale(level),center0,
c     $               rmlexp(iaddr(1,ibox)),nterms(level),
c     $               sourcesort(1,j),
c     $               ptemp,iffld,ftemp,ifhess,htemp,
c     $               ier)
                  call l3dmpevalhessd_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),
     $               sourcesort(1,j),
     $               ptemp,iffld,ftemp,ifhess,htemp,
     $               scarray_mpole,wlege,nlege)
                  else
                  call l3dmpeval_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),nterms(level),
     $               sourcesort(1,j),
     $               ptemp,iffld,ftemp,
     $               wlege,nlege,ier)
                  endif
                  if( ifpot .eq. 1 ) pot(j)=pot(j)+ptemp
                  if( iffld .eq. 1 ) then
                     fld(1,j)=fld(1,j)+ftemp(1)
                     fld(2,j)=fld(2,j)+ftemp(2)
                     fld(3,j)=fld(3,j)+ftemp(3)
                  endif
                  if (ifhess .eq. 1) then
                  hess(1,j)=hess(1,j)+htemp(1)
                  hess(2,j)=hess(2,j)+htemp(2)
                  hess(3,j)=hess(3,j)+htemp(3)
                  hess(4,j)=hess(4,j)+htemp(4)
                  hess(5,j)=hess(5,j)+htemp(5)
                  hess(6,j)=hess(6,j)+htemp(6)
                  endif
               enddo
               do j=box1(16),box1(16)+box1(17)-1
                  if( ifhesstarg .eq. 1 ) then
c                  call l3dmpevalhess(scale(level),center0,
c     $               rmlexp(iaddr(1,ibox)),nterms(level),
c     $               targetsort(1,j),
c     $               ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
c     $               ier)
                  call l3dmpevalhessd_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),
     $               targetsort(1,j),
     $            ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
     $            scarray_mpole,wlege,nlege)
                  else
                  call l3dmpeval_trunc(scale(level),center0,
     $               rmlexp(iaddr(1,ibox)),nterms(level),nterms(level),
     $               targetsort(1,j),
     $               ptemp,iffldtarg,ftemp,
     $               wlege,nlege,ier)
                  endif
                  if( ifpottarg .eq. 1 ) pottarg(j)=pottarg(j)+ptemp
                  if( iffldtarg .eq. 1 ) then
                     fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                     fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                     fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
                  endif
                  if (ifhesstarg .eq. 1) then
                  hesstarg(1,j)=hesstarg(1,j)+htemp(1)
                  hesstarg(2,j)=hesstarg(2,j)+htemp(2)
                  hesstarg(3,j)=hesstarg(3,j)+htemp(3)
                  hesstarg(4,j)=hesstarg(4,j)+htemp(4)
                  hesstarg(5,j)=hesstarg(5,j)+htemp(5)
                  hesstarg(6,j)=hesstarg(6,j)+htemp(6)
                  endif
               enddo
            else
            
            call lfmm3dtrilhess_direct(nqtri,box,box1,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
            endif
        enddo
C$OMP END PARALLEL DO
 3252   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c

        allocate( scarray_local(0:100000) )
        call l3dtaevalhessdini(nterms(0),scarray_local)

        if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$PRIVATE(i,j,ptemp,ftemp,htemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(15)
c       
        if (level .ge. 2) then
            do j=box(14),box(14)+box(15)-1
               if( ifhess .eq. 1 ) then
c                  call l3dtaevalhess(scale(level),center0,
c     $               rmlexp(iaddr(2,ibox)),nterms(level),
c     $               sourcesort(1,j),
c     $               ptemp,iffld,ftemp,ifhess,htemp,
c     $               ier)
               call l3dtaevalhessd_trunc(scale(level),center0,
     $              rmlexp(iaddr(2,ibox)),nterms(level),sourcesort(1,j),
     $              ptemp,iffld,ftemp,ifhess,htemp,
     $              scarray_local,wlege,nlege)
               else
               call l3dtaeval_trunc(scale(level),center0,
     $            rmlexp(iaddr(2,ibox)),nterms(level),nterms(level),
     $            sourcesort(1,j),
     $            ptemp,iffld,ftemp,
     $            wlege,nlege,ier)
               endif
               if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
               if (iffld .eq. 1) then
                  fld(1,j)=fld(1,j)+ftemp(1)
                  fld(2,j)=fld(2,j)+ftemp(2)
                  fld(3,j)=fld(3,j)+ftemp(3)
               endif
               if (ifhess .eq. 1) then
               hess(1,j)=hess(1,j)+htemp(1)
               hess(2,j)=hess(2,j)+htemp(2)
               hess(3,j)=hess(3,j)+htemp(3)
               hess(4,j)=hess(4,j)+htemp(4)
               hess(5,j)=hess(5,j)+htemp(5)
               hess(6,j)=hess(6,j)+htemp(6)
               endif
            enddo

            do j=box(16),box(16)+box(17)-1
               if( ifhesstarg .eq. 1 ) then
c                  call l3dtaevalhess(scale(level),center0,
c     $               rmlexp(iaddr(2,ibox)),nterms(level),
c     $               targetsort(1,j),
c     $               ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
c     $               ier)
               call l3dtaevalhessd_trunc(scale(level),center0,
     $              rmlexp(iaddr(2,ibox)),nterms(level),
     $              targetsort(1,j),
     $              ptemp,iffldtarg,ftemp,ifhesstarg,htemp,
     $              scarray_local,wlege,nlege)
               else
               call l3dtaeval_trunc(scale(level),center0,
     $            rmlexp(iaddr(2,ibox)),nterms(level),nterms(level),
     $            targetsort(1,j),
     $            ptemp,iffldtarg,ftemp,
     $            wlege,nlege,ier)
               endif
               if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
               if (iffldtarg .eq. 1) then
                  fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                  fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                  fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
               endif
               if (ifhesstarg .eq. 1) then
               hesstarg(1,j)=hesstarg(1,j)+htemp(1)
               hesstarg(2,j)=hesstarg(2,j)+htemp(2)
               hesstarg(3,j)=hesstarg(3,j)+htemp(3)
               hesstarg(4,j)=hesstarg(4,j)+htemp(4)
               hesstarg(5,j)=hesstarg(5,j)+htemp(5)
               hesstarg(6,j)=hesstarg(6,j)+htemp(6)
               endif
            enddo
        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
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
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        call lfmm3dtrilhess_direct_self(nqtri,box,
     $     triaflatsort,trianormsort,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
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
            call lfmm3dtrilhess_direct(nqtri,box1,box,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,ifhess,hess,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg,
     $         ifhesstarg,hesstarg)
c
 6203       continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if(ifprint .ge. 1) then
        call prin2('timeinfo=*',timeinfo,8)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        endif
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dtrilhess_direct(nqtri,box,box1,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),hess(6,1)
        complex *16 pottarg(1),fldtarg(3,1),hesstarg(6,1)
        complex *16 ptemp,ftemp(3),htemp(6)
c
c       ... do nothing
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtrilhess_direct_self(nqtri,box,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ifhess,hess,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     ifhesstarg,hesstarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),hess(6,1)
        complex *16 pottarg(1),fldtarg(3,1),hesstarg(6,1)
        complex *16 ptemp,ftemp(3),htemp(6)
c
c
c       ... do nothing
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       Low level FMM routines for triangles in Laplace regime.
c       (forming multipole expansions and forming taylor expansions).
c       
c       Linear densities on flat triangles
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     FORMMP routines: Individual triangle to multipole expansion
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine l3dformmptrisone_linear
     $     (ier, scale, triang,
     $     charge, x0y0z0, norder, nterms, mpole)
        implicit real *8 (a-h,o-z)
c
c
c     Form  multipole expansion due to SLP on  single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     charge       = linear density of SLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c     work(lw)     = workspace
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3), x0y0z0(3)
        real *8  scale
        real *8  rnodes(2,500), weights(500), zparts(3,500)
        real *8  vert1(2),vert2(2),vert3(2),vertout(3)
        real *8  rnodes0(2,500), weights0(500)
        real *8  vert10(2),vert20(2),vert30(2)
        real *8  w(20)
c
        complex *16 cval
c       
        complex *16 eye, charge(3), ctemp(500)
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
        data eye/(0.0d0,1.0d0)/
c
        vert10(1)=0
        vert10(2)=0
        vert20(1)=1
        vert20(2)=0
        vert30(1)=0
        vert30(2)=1
c
        call triasymq
     $     (norder,vert10,vert20,vert30,rnodes,weights,nnodes)
c       
        call triangle_area(triang,ds)
        ds=ds*2
c
        do i = 1,nnodes
        
        u=rnodes(1,i)
        v=rnodes(2,i)
        
        zparts(1,i)=triang(1,1)+
     $     u*(triang(1,2)-triang(1,1))+v*(triang(1,3)-triang(1,1))
        zparts(2,i)=triang(2,1)+
     $     u*(triang(2,2)-triang(2,1))+v*(triang(2,3)-triang(2,1))
        zparts(3,i)=triang(3,1)+
     $     u*(triang(3,2)-triang(3,1))+v*(triang(3,3)-triang(3,1))
        
        cval=charge(1)+
     $     u*(charge(2)-charge(1))+v*(charge(3)-charge(1))
        
        ctemp(i) = cval*weights(i) *ds

        enddo
        
        call l3dformmp(ier,scale,zparts,ctemp,nnodes,x0y0z0,
     1     nterms,mpole)
c       
        return
        end
c
c
c
c
c
        subroutine l3dformmptris_linear
     $     (ier, scale, triang, charge, 
     1     ntri, x0y0z0, norder, nterms, mpole, mptemp)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of linear SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), x0y0z0(3), scale
        complex *16  charge(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms), mptemp(1)
c
        do i = 1,ntri
        call l3dformmptrisone_linear(ier, scale,triang(1,1,i),
     $     charge(1,i), x0y0z0, norder, nterms, mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
        subroutine l3dformmptris_linear_add
     $     (ier, scale, triang, charge, 
     1     ntri, x0y0z0, norder, nterms, mpole)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of linear SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), x0y0z0(3), scale
        complex *16  charge(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mptemp(:,:)
c
        allocate( mptemp(0:nterms,-nterms:nterms) )
c
        do i = 1,ntri
        call l3dformmptrisone_linear(ier, scale,triang(1,1,i),
     $     charge(1,i), x0y0z0, norder, nterms, mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
c
c
c
c
        subroutine l3dformmptridone_linear
     $     (ier, scale, triang, trinorm, dipstr, dipvec,
     $     x0y0z0, norder, nterms, mpole)
        implicit real *8 (a-h,o-z)
c
c     Form multipole expansion due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal 
c                        (to be interpolated for curved triangles)
c     dipstr       = linear density of DLP strength
c     dipvec       = linear density of DLP vector
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c     work(lw)     = workspace
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
        real *8 triang(3,3), trinorm(3), x0y0z0(3)
        real *8 scale
        real *8 rnodes(2,500), weights(500), zparts(3,500)
        real *8 vert1(2),vert2(2),vert3(2),vertout(3)
        real *8 rnodes0(2,500), weights0(500)
        real *8 vert10(2),vert20(2),vert30(2)
        real *8 w(20)
c
        real *8 dipvec(3,3), dtemp(3,500)
        complex *16 eye, dipstr(3), ctemp(500)
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
        complex *16 cval
c       
        data eye/(0.0d0,1.0d0)/
c
        vert10(1)=0
        vert10(2)=0
        vert20(1)=1
        vert20(2)=0
        vert30(1)=0
        vert30(2)=1
c
        call triasymq
     $     (norder,vert10,vert20,vert30,rnodes,weights,nnodes)
c
        call triangle_area(triang,ds)
        ds=ds*2
c       
        do i = 1,nnodes

        u=rnodes(1,i)
        v=rnodes(2,i)

        zparts(1,i)=triang(1,1)+
     $     u*(triang(1,2)-triang(1,1))+v*(triang(1,3)-triang(1,1))
        zparts(2,i)=triang(2,1)+
     $     u*(triang(2,2)-triang(2,1))+v*(triang(2,3)-triang(2,1))
        zparts(3,i)=triang(3,1)+
     $     u*(triang(3,2)-triang(3,1))+v*(triang(3,3)-triang(3,1))
        
        cval=dipstr(1)+
     $     u*(dipstr(2)-dipstr(1))+v*(dipstr(3)-dipstr(1))
        
        ctemp(i)=cval*weights(i) *ds

        dtemp(1,i)=dipvec(1,1)+
     $     u*(dipvec(1,2)-dipvec(1,1))+v*(dipvec(1,3)-dipvec(1,1))
        dtemp(2,i)=dipvec(2,1)+
     $     u*(dipvec(2,2)-dipvec(2,1))+v*(dipvec(2,3)-dipvec(2,1))
        dtemp(3,i)=dipvec(3,1)+
     $     u*(dipvec(3,2)-dipvec(3,1))+v*(dipvec(3,3)-dipvec(3,1))
         
        enddo
c
        call l3dformmp_dp(ier,scale,zparts,ctemp,dtemp,nnodes,x0y0z0,
     1     nterms,mpole)
c
        return
        end
c
c
c
c
c
        subroutine l3dformmptrid_linear
     $     (ier,scale,triang,trinorm,dipstr,dipvec,
     1     ntri,x0y0z0,norder,nterms,mpole,mptemp)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of linear DLP densities
c     dipvec             = array of linear DLP vectors
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
        real *8  dipvec(3,3,1)
        complex *16  eye, dipstr(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms), mptemp(1)
c       
        do i = 1,ntri
        call l3dformmptridone_linear(ier,scale,triang(1,1,i),
     1     trinorm(1,i),dipstr(1,i),dipvec(1,1,i),
     $     x0y0z0, norder,nterms,mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
        subroutine l3dformmptrid_linear_add
     $     (ier,scale,triang,trinorm,dipstr,dipvec,
     1     ntri,x0y0z0,norder,nterms,mpole)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of linear DLP densities
c     dipvec             = array of linear DLP vectors
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
        real *8  dipvec(3,3,1)
        complex *16  eye, dipstr(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mptemp(:,:)
c
        allocate( mptemp(0:nterms,-nterms:nterms) )
c       
        do i = 1,ntri
        call l3dformmptridone_linear(ier,scale,triang(1,1,i),
     1     trinorm(1,i),dipstr(1,i),dipvec(1,1,i),
     $     x0y0z0, norder,nterms,mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     FORMTA routines: Individual triangle to taylor expansion
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine l3dformtatrisone_linear
     $     (ier, scale, triang,
     $     charge, x0y0z0, norder, nterms, mpole)
        implicit real *8 (a-h,o-z)
c
c
c     Form  local expansion due to SLP on  single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     charge       = linear density of SLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of local expansion
c     work(lw)     = workspace
c
c     OUTPUT:
c
c     mpole        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3), x0y0z0(3)
        real *8  scale
        real *8  rnodes(2,500), weights(500), zparts(3,500)
        real *8  vert1(2),vert2(2),vert3(2),vertout(3)
        real *8  rnodes0(2,500), weights0(500)
        real *8  vert10(2),vert20(2),vert30(2)
        real *8  w(20)
c
        complex *16 cval
c       
        complex *16 eye, charge(3), ctemp(500)
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
        data eye/(0.0d0,1.0d0)/
c
        vert10(1)=0
        vert10(2)=0
        vert20(1)=1
        vert20(2)=0
        vert30(1)=0
        vert30(2)=1
c       
        call triasymq
     $     (norder,vert10,vert20,vert30,rnodes,weights,nnodes)
c       
        call triangle_area(triang,ds)
        ds=ds*2
c
        do i = 1,nnodes
        
        u=rnodes(1,i)
        v=rnodes(2,i)
        
        zparts(1,i)=triang(1,1)+
     $     u*(triang(1,2)-triang(1,1))+v*(triang(1,3)-triang(1,1))
        zparts(2,i)=triang(2,1)+
     $     u*(triang(2,2)-triang(2,1))+v*(triang(2,3)-triang(2,1))
        zparts(3,i)=triang(3,1)+
     $     u*(triang(3,2)-triang(3,1))+v*(triang(3,3)-triang(3,1))
        
        cval=charge(1)+
     $     u*(charge(2)-charge(1))+v*(charge(3)-charge(1))
        
        ctemp(i) = cval*weights(i) *ds
        
        enddo
        
        call l3dformta(ier,scale,zparts,ctemp,nnodes,x0y0z0,
     1     nterms,mpole)
c       
        return
        end
c
c
c
c
c
        subroutine l3dformtatris_linear
     $     (ier, scale, triang, charge, 
     1     ntri, x0y0z0, norder, nterms, mpole, mptemp)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of linear SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of local expansion
c     mptemp             = work array to hold temp local expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), x0y0z0(3), scale
        complex *16  charge(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms), mptemp(1)
c
        do i = 1,ntri
        call l3dformtatrisone_linear(ier, scale,triang(1,1,i),
     $     charge(1,i), x0y0z0, norder, nterms, mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
        subroutine l3dformtatris_linear_add
     $     (ier, scale, triang, charge, 
     1     ntri, x0y0z0, norder, nterms, mpole)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of linear SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of local expansion
c     mptemp             = work array to hold temp local expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), x0y0z0(3), scale
        complex *16  charge(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mptemp(:,:)
c
        allocate( mptemp(0:nterms,-nterms:nterms) )
c
        do i = 1,ntri
        call l3dformtatrisone_linear(ier, scale,triang(1,1,i),
     $     charge(1,i), x0y0z0, norder, nterms, mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
c
c
c
c
        subroutine l3dformtatridone_linear
     $     (ier, scale, triang, trinorm, dipstr, dipvec,
     $     x0y0z0, norder, nterms, mpole)
        implicit real *8 (a-h,o-z)
c
c     Form local expansion due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal 
c                        (to be interpolated for curved triangles)
c     dipstr       = linear density of DLP strength
c     dipvec       = linear density of DLP vector
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of local expansion
c     work(lw)     = workspace
c
c     OUTPUT:
c
c     mpole        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
        real *8 triang(3,3), trinorm(3), x0y0z0(3)
        real *8 scale
        real *8 rnodes(2,500), weights(500), zparts(3,500)
        real *8 vert1(2),vert2(2),vert3(2),vertout(3)
        real *8 rnodes0(2,500), weights0(500)
        real *8 vert10(2),vert20(2),vert30(2)
        real *8 w(20)
c
        real *8 dipvec(3,3), dtemp(3,500)
        complex *16 eye, dipstr(3), ctemp(500)
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
        complex *16 cval
c       
        data eye/(0.0d0,1.0d0)/
c
        vert10(1)=0
        vert10(2)=0
        vert20(1)=1
        vert20(2)=0
        vert30(1)=0
        vert30(2)=1
c
        call triasymq
     $     (norder,vert10,vert20,vert30,rnodes,weights,nnodes)
c       
        call triangle_area(triang,ds)
        ds=ds*2
c
        do i = 1,nnodes

        u=rnodes(1,i)
        v=rnodes(2,i)

        zparts(1,i)=triang(1,1)+
     $     u*(triang(1,2)-triang(1,1))+v*(triang(1,3)-triang(1,1))
        zparts(2,i)=triang(2,1)+
     $     u*(triang(2,2)-triang(2,1))+v*(triang(2,3)-triang(2,1))
        zparts(3,i)=triang(3,1)+
     $     u*(triang(3,2)-triang(3,1))+v*(triang(3,3)-triang(3,1))
        
        cval=dipstr(1)+
     $     u*(dipstr(2)-dipstr(1))+v*(dipstr(3)-dipstr(1))
        
        ctemp(i)=cval*weights(i) *ds

        dtemp(1,i)=dipvec(1,1)+
     $     u*(dipvec(1,2)-dipvec(1,1))+v*(dipvec(1,3)-dipvec(1,1))
        dtemp(2,i)=dipvec(2,1)+
     $     u*(dipvec(2,2)-dipvec(2,1))+v*(dipvec(2,3)-dipvec(2,1))
        dtemp(3,i)=dipvec(3,1)+
     $     u*(dipvec(3,2)-dipvec(3,1))+v*(dipvec(3,3)-dipvec(3,1))
         
        enddo
c
        call l3dformta_dp(ier,scale,zparts,ctemp,dtemp,nnodes,x0y0z0,
     1     nterms,mpole)
c
        return
        end
c
c
c
c
c
        subroutine l3dformtatrid_linear
     $     (ier,scale,triang,trinorm,dipstr,dipvec,
     1     ntri,x0y0z0,norder,nterms,mpole,mptemp)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of linear DLP densities
c     dipvec             = array of linear DLP vectors
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of local expansion
c     mptemp             = work array to hold temp local expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
        real *8  dipvec(3,3,1)
        complex *16  eye, dipstr(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms), mptemp(1)
c       
        do i = 1,ntri
        call l3dformtatridone_linear(ier,scale,triang(1,1,i),
     1     trinorm(1,i),dipstr(1,i),dipvec(1,1,i),
     $     x0y0z0, norder,nterms,mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
        subroutine l3dformtatrid_linear_add
     $     (ier,scale,triang,trinorm,dipstr,dipvec,
     1     ntri,x0y0z0,norder,nterms,mpole)
        implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of linear DLP densities
c     dipvec             = array of linear DLP vectors
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of local expansion
c     mptemp             = work array to hold temp local expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
        real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
        real *8  dipvec(3,3,1)
        complex *16  eye, dipstr(3,1)
        complex *16  mpole(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mptemp(:,:)
c
        allocate( mptemp(0:nterms,-nterms:nterms) )
c
c       
        do i = 1,ntri
        call l3dformtatridone_linear(ier,scale,triang(1,1,i),
     1     trinorm(1,i),dipstr(1,i),dipvec(1,1,i),
     $     x0y0z0, norder,nterms,mptemp)
        call l3dadd(mptemp,mpole,nterms)
        enddo
        return
        end
c
c
c
c
c
        subroutine l3dreorder_linear(nsource,source,
     $     ifcharge,charge,isource,ifdipole,
     1     dipstr,dipvec,sourcesort,chargesort,dipvecsort,dipstrsort) 
        implicit real *8 (a-h,o-z)
        dimension source(3,1),sourcesort(3,1),isource(1)
        complex *16 charge(3,1),chargesort(3,1)
        complex *16 dipstr(3,1),dipstrsort(3,1)
        dimension dipvec(3,3,1),dipvecsort(3,3,1)
c       
ccc        call prinf('nsource=*',nsource,1)
        do i = 1,nsource
        sourcesort(1,i) = source(1,isource(i))
        sourcesort(2,i) = source(2,isource(i))
        sourcesort(3,i) = source(3,isource(i))
        if( ifcharge .eq. 1 ) then
        chargesort(1,i) = charge(1,isource(i))
        chargesort(2,i) = charge(2,isource(i))
        chargesort(3,i) = charge(3,isource(i))
        endif
        if (ifdipole .eq. 1) then
        dipstrsort(1,i) = dipstr(1,isource(i))
        dipstrsort(2,i) = dipstr(2,isource(i))
        dipstrsort(3,i) = dipstr(3,isource(i))
        dipvecsort(1,1,i) = dipvec(1,1,isource(i))
        dipvecsort(2,1,i) = dipvec(2,1,isource(i))
        dipvecsort(3,1,i) = dipvec(3,1,isource(i))
        dipvecsort(1,2,i) = dipvec(1,2,isource(i))
        dipvecsort(2,2,i) = dipvec(2,2,isource(i))
        dipvecsort(3,2,i) = dipvec(3,2,isource(i))
        dipvecsort(1,3,i) = dipvec(1,3,isource(i))
        dipvecsort(2,3,i) = dipvec(2,3,isource(i))
        dipvecsort(3,3,i) = dipvec(3,3,isource(i))
        endif
        enddo
        return
        end
c
c
c
c
c