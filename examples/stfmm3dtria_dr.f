c       
c       This is a simple driver to test the layer potential FMM routines
c       in R^3 using the free-space Stokes Green's functions.
c 
c       The geometry is assumed to consist of a collection of flat
c       triangles with piecewise constant source densities (tractions
c       for the single layer and/or jumps in displacement for the double
c       layer).
c
        implicit real *8 (a-h,o-z)
        parameter (lw=120 000 000)
        dimension w(lw)
c
        call stfmm3d_test(w,lw)
c
        stop
        end
c
c
c
c
c
        subroutine stfmm3d_test(w,lw)
c
c       Compute layer potentials by direct calculation and via FMM.
c
c
        implicit real *8 (a-h,o-z)
        parameter(nmax=100000)
c
        real *8 triangles(3,3,nmax),centroids(3,nmax)
        real *8 trinorm(3,nmax),triarea(nmax)
        real *8 verts(3,nmax),ifaces(3,nmax)       
        dimension itrivert(3,nmax)
c
        real *8 sigma_sl(3,nmax)
        real *8 sigma_dl(3,nmax)
c
        real *8 target(3,nmax)
c
        dimension vert1(3),vert2(3),vert3(3),vertout(3)
        dimension w(lw)
c
        dimension pot(3,nmax),grad(3,3,nmax)
        dimension pottarg(3,nmax),gradtarg(3,3,nmax)
c
        dimension pot0(3),grad0(3,3),stress0(3,3),tract0(3)
c
        dimension pot1(3,nmax),grad1(3,3,nmax)
        dimension pot2(3,nmax),grad2(3,3,nmax)
c
        dimension pre(nmax)
        dimension pretarg(nmax)
c
        dimension pre1(nmax)
        dimension pre2(nmax)
c
c
c       SET ALL PARAMETERS
c
c
c       PRINTING: prini determines whether output is printed or not
c       with subsequent calls to prin2 or prinf. prina(i1,i2) simply
c       prints to the Fortran unit numbers i1 and i2 assuming they are
c       nonzero. Setting i1 or i2 to zero suppresses printing. 
c
        call prini(6,13)
c
c       ... get scatterer geometry
c
        call getgeom(ntri,triangles,centroids,trinorm,triarea,nmax,
     1     verts,nmax,itrivert,nmax,iergeom)
        call prinf('geometry error flag is *',iergeom,1)
c
        call prinf('ntri=*',ntri,1)
c
ccc        call prin2('triangles=*',triangles,3*3*ntri)
ccc        call prin2('centroids=*',centroids,3*ntri)
ccc        call prin2('trinorm=*',trinorm,3*ntri)
c       
c
c       define (random) piecewise constant densities
c
        do i=1,ntri
        sigma_sl(1,i)=hkrand(0)
        sigma_sl(2,i)=hkrand(0)
        sigma_sl(3,i)=hkrand(0)
        sigma_dl(1,i)=hkrand(0)
        sigma_dl(2,i)=hkrand(0)
        sigma_dl(3,i)=hkrand(0)
        enddo        
c
c       define piecewise constant densities
c
        do i=1,ntri
c        sigma_sl(1,i)=1
c        sigma_sl(2,i)=1
c        sigma_sl(3,i)=1
c        sigma_dl(1,i)=1
c        sigma_dl(2,i)=1
c        sigma_dl(3,i)=1
        enddo        
c
c       set targets on one surface just outside the sphere and
c       one surface just inside the sphere.
c
        do i=1,ntri
        h=1d-4
        target(1,2*i-1)=centroids(1,i)+h*trinorm(1,i) 
        target(2,2*i-1)=centroids(2,i)+h*trinorm(2,i) 
        target(3,2*i-1)=centroids(3,i)+h*trinorm(3,i) 
        target(1,2*i  )=centroids(1,i)-h*trinorm(1,i)
        target(2,2*i  )=centroids(2,i)-h*trinorm(2,i)
        target(3,2*i  )=centroids(3,i)-h*trinorm(3,i)
        enddo
        ntargs=2*ntri
c        ntargs=0
c
        call prinf('ntargs=*',ntargs,1)
ccc        call prin2('target=*',target,3*ntargs)
c
        call prinf('==== targets ====*',i,0)
c       
c       since we are testing by direct calculation, only compute
c       direct calculation at m targets.
c
        m=min(ntargs,20)
        m=20
c
        m= min(m,ntri)
        m= min(m,ntargs)
c
c       turn on single and/or double layer potential
c
        ifsingle=1
        ifdouble=4
c
c       set whether displacement and/or gradient to be compute on surface
c       and at target locations.
c
        ifpot=1
        ifgrad=1
        ifpottarg=1
        ifgradtarg=1
c
        ifprint=0
c
        call prinf('m=*',m,1)
c
c       ... evaluate via direct quadrature on triangles
c
        do i=1,m
        do j=1,3
           pot1(j,i)=0
        enddo
        pre1(i)=0
        do j=1,3
        do k=1,3
           grad1(j,k,i)=0
        enddo
        enddo
        enddo
c
        do i=1,m
        do j=1,3
           pot2(j,i)=0
        enddo
        pre2(i)=0
        do j=1,3
        do k=1,3
           grad2(j,k,i)=0
        enddo
        enddo
        enddo
c
        t1=second()
C$        t1=omp_get_wtime()

        call st3dtriadirecttarg_test(M,
     $     TRIANGLES,TRINORM,NTRI,CENTROIDS,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifpot,pot2,pre2,ifgrad,grad2,NTARGS,
     $     target,ifpottarg,POT1,PRE1,
     $     ifgradtarg,GRAD1)

        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('after stfmm3triadirect, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (m)/(t2-t1),1)
        call prin2('after estimated time for direct=*',
     $     ntri/dble(m)*(t2-t1),1)
c
c       ... evaluate via FMM on triangles 
c
        iprec=1
c
ccc        call prini(0,13)
        t1=second()
C$        t1=omp_get_wtime()
        call stfmm3dtriatarg
     $     (ier,iprec,TRIANGLES,TRINORM,NTRI,centroids,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntargs,target,ifpottarg,POTtarg,PREtarg,
     $     ifgradtarg,GRADtarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prini(6,13)
c
        call prin2('after stfmm3triatarg, time=*',t2-t1,1)
        call prin2('speed, targets/sec=*',
     $     (ntargs)/(t2-t1),1)
c
c
        if( ifprint .eq. 1 ) then
c
        if( ifpot .eq. 1 ) then
        call prin2('after fmm, pot=*',pot,3*m)
        call prin2('after dir, pot2=*',pot2,3*m)
        call prin2('after fmm, pre=*',pre,m)
        call prin2('after dir, pre2=*',pre2,m)
        endif
        if( ifgrad .eq. 1 ) then
        call prin2('after fmm, grad=*',grad,3*3*m)
        call prin2('after dir, grad2=*',grad2,3*3*m)
        endif
c
        if( ifpottarg .eq. 1 ) then
        call prin2('after fmm, pottarg=*',pottarg,3*m)
        call prin2('after dir, pot1=*',pot1,3*m)
        call prin2('after fmm, pretarg=*',pretarg,m)
        call prin2('after dir, pre1=*',pre1,m)
        endif
        if( ifgradtarg .eq. 1 ) then
        call prin2('after fmm, gradtarg=*',gradtarg,3*3*m)
        call prin2('after dir, grad1=*',grad1,3*3*m)
        endif
c
        endif
c
c
        if( ifpot .eq. 1 ) then
        call d3error(pot2,pot,3*m,a,r)
ccc        call prin2('absolute error in pot=*',a,1)
        call prin2('relative error in pot=*',r,1)
        call d3error(pre2,pre,m,a,r)
        call prin2('absolute error in pre=*',a,1)
        call prin2('relative error in pre=*',r,1)
        endif
c
        if( ifgrad .eq. 1 ) then
        call d3error(grad2,grad,3*3*m,a,r)
ccc        call prin2('absolute error in grad=*',a,1)
        call prin2('relative error in grad=*',r,1)
        endif
c
        if( ifpottarg .eq. 1 ) then
        call d3error(pot1,pottarg,3*m,a,r)
ccc        call prin2('absolute error in target pot=*',a,1)
        call prin2('relative error in target pot=*',r,1)
        call d3error(pre1,pretarg,m,a,r)
        call prin2('absolute error in target pre=*',a,1)
        call prin2('relative error in target pre=*',r,1)
        endif
c
        if( ifgradtarg .eq. 1 ) then
        call d3error(grad1,gradtarg,3*3*m,a,r)
ccc        call prin2('absolute error in target grad=*',a,1)
        call prin2('relative error in target grad=*',r,1)
        endif
c
        stop
c
        if( ifpot .eq. 1 ) then
        call prin2('after stfmm3triatarg, self pot=*',
     $     pot,3)
        endif

        if( ifgrad .eq. 1 ) then
        call prin2('after stfmm3triatarg, self grad=*',
     $     grad,3*3)
        endif
c

        return
        end
c
c
c
c
c
        subroutine getgeom(ntri,triangles,centroids,trinorm,triarea,
     1     nmax,verts,lv,itrivert,li,iergeom)
        implicit real *8 (a-h,o-z)
        dimension triangles(3,3,1),centroids(3,1)
        dimension trinorm(3,1),triarea(1)
        dimension verts(3,lv),itrivert(3,li)
c       
c       INPUT 
c       
c       nmax = dimension of trisngles, centroids, trinorm, triarea arrays
c       verts =  work array to read in vertices 
c       lv = dimension of verts (work) array
c       itrivert = integer work array to read in vertices associated 
c          with each triangle
c       li = dimension of itrivert (work) array
c       
c       OUTPUT 
c
c       ntri = number of triangles
c       triangles = triangles array
c       centroids = centroids array
c       trinorm = triangle normals array
c       triarea = triangle areas array
c       iergeom = error return code
c      
c       iergeom = 0  normal execution 
c       iergeom = 1  error opening unit ir in atrireadchk
c       iergeom = 2  lv is too small to read all vertices
c       iergeom = 3  li is too small (must exceed ntri)
c       iergeom = 4  nmax is too small (must exceed ntri)
c
c-----------------------------------------------------------------------
c
        ir = 17
cc        open (unit = ir,file='sphere180.a.tri')
        open (unit = ir,file='sphere720.a.tri')
cc        open (unit = ir,file='sphere2880.a.tri')
c        
        call atrireadchk(ir,verts,lv,nverts,itrivert,li,ntri,iergeom)
c
        if (ntri.gt.nmax) then
        iergeom = 4
        return
        endif
c    
c       create triangles from tri format data
c
c
c       scale and offset if desired.
c
        scale=1
        dx=0
        dy=0
        dz=0
c       
        do itri = 1,ntri
        v11 = verts(1,itrivert(1,itri)) *scale + dx
        v21 = verts(2,itrivert(1,itri)) *scale + dy
        v31 = verts(3,itrivert(1,itri)) *scale + dz
        v12 = verts(1,itrivert(2,itri)) *scale + dx 
        v22 = verts(2,itrivert(2,itri)) *scale + dy
        v32 = verts(3,itrivert(2,itri)) *scale + dz
        v13 = verts(1,itrivert(3,itri)) *scale + dx
        v23 = verts(2,itrivert(3,itri)) *scale + dy
        v33 = verts(3,itrivert(3,itri)) *scale + dz
c        
        triangles(1,1,itri) = v11
        triangles(2,1,itri) = v21
        triangles(3,1,itri) = v31
        triangles(1,2,itri) = v12
        triangles(2,2,itri) = v22
        triangles(3,2,itri) = v32
        triangles(1,3,itri) = v13
        triangles(2,3,itri) = v23
        triangles(3,3,itri) = v33
        centroids(1,itri) = (v11+v12+v13)/3
        centroids(2,itri) = (v21+v22+v23)/3
        centroids(3,itri) = (v31+v32+v33)/3
        call triangle_norm(triangles(1,1,itri),trinorm(1,itri))
        call triangle_area(triangles(1,1,itri),triarea(itri))
c
        enddo
        return
        end
c
c
c
c
c
c      
        subroutine d3error(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        dimension pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
        d=d+abs(pot1(i)-pot2(i))**2
        a=a+abs(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
c
c
c
        subroutine st3dtriadirecttarg_test(M,
     $     TRIANGLE,TRINORM,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifpot,pot,pre,ifgrad,grad,NTARGET,
     $     target,ifpottarg,POTtarg,PREtarg,
     $     ifgradtarg,GRADtarg)
        implicit real *8 (a-h,o-z)
c
c
c       Stokes interactions in R^3: evaluate all pairwise triangle
c       interactions and interactions with targets using the direct
c       O(N^2) algorithm.
c
c       INPUT:
c
c       M - the number of triangles/target locations to be tested
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
c       grad(3,3,nsource) - gradient at source locations
c       pottarg(3,ntarget) - velocity at target locations
c       gradtarg(3,3,ntarget) - gradient at target locations
c
c
        dimension triangle(3,3,1),trinorm(3,1),source(3,1)
        dimension sigma_sl(3,1),sigma_dl(3,1),sigma_dv(3,1)
        dimension target(3,1)
c
        dimension pot(3,1),pre(1),grad(3,3,1)
        dimension pottarg(3,1),pretarg(1),gradtarg(3,3,1)
c
        dimension pot0(3),grad0(3,3)
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
        do j=1,min(m,nsource)
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
        do j=1,min(m,ntarget)
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
