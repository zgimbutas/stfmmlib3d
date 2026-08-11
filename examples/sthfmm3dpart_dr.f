c
c       This is a simple driver to test the particle FMM routines in R^3
c       using the half-space Stokes Green's functions, with the no-slip
c       (zero-velocity) boundary condition u = 0 at z=0.
c
        implicit real *8 (a-h,o-z)
c
        call sthfmm3d_test()
c
        stop
        end
c
c
c
c
c
        subroutine sthfmm3d_test()
        implicit real *8 (a-h,o-z)
c
c       Compute potentials by direct calculation and via FMM.
c
        parameter(nmax=100000)
c
        dimension sigma_sl(3,nmax)
        dimension sigma_dl(3,nmax)
        dimension sigma_dv(3,nmax)
c
c
        dimension source(3,nmax)
        dimension target(3,nmax)
c
        dimension pot(3,nmax),grad(3,3,nmax)
        dimension pottarg(3,nmax),gradtarg(3,3,nmax)
c
        dimension pot0(3),grad0(3,3)
c
        dimension pot1(3,nmax),grad1(3,3,nmax)
c
        dimension pre(nmax)
        dimension pretarg(nmax)
c
        dimension pre1(nmax)
c
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        print *, 'ENTER n'
        read *, nsource
c
c
        call prinf('nsource=*',nsource,1)
c
        done=1
        pi=4*atan(done)
c
        idist=3
c
        if( idist .eq. 1 ) then
c
c       ... construct randomly located charge distribution on a unit cube
c       
        do i=1,nsource
        source(1,i)=hkrand(0)
        source(2,i)=hkrand(0)
        source(3,i)=hkrand(0)
        source(1,i)=source(1,i)-0.5
        source(2,i)=source(2,i)-0.5
        source(3,i)=source(3,i)-0.5
        enddo
c
        endif
c
        if( idist .eq. 2 ) then
c
c       ... construct charge distribution on a curve in R^3
c       
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        source(1,i)=sin(1.1*a)
        source(2,i)=cos(2.1*a)
        source(3,i)=cos(3.1*a)
        enddo
c
        endif
c       
        if( idist .eq. 3 ) then
c
c       ... construct randomly located charge distribution on a unit sphere
c
        d=hkrand(0)
        do i=1,nsource
c
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        source(1,i)=.5d0*cos(phi)*sin(theta)
        source(2,i)=.5d0*sin(phi)*sin(theta)
        source(3,i)=.5d0*cos(theta)
c
        d=hkrand(0)
        d=hkrand(0)
        enddo
c
        endif
c
        if( idist .eq. 4 ) then
c
c       ... construct the grid of charges
c
        ngrid=nsource-1
        kk=0
        do i=1,ngrid+1
        do j=1,ngrid+1
        do k=1,ngrid+1
        kk=kk+1
        source(1,kk)=(i-1.0)/ngrid
        source(2,kk)=(j-1.0)/ngrid
        source(3,kk)=(k-1.0)/ngrid
        source(1,kk)=source(1,kk)+hkrand(0)*.0001
        source(2,kk)=source(2,kk)+hkrand(0)*.0001
        source(3,kk)=source(3,kk)+hkrand(0)*.0001
ccc        call prin2('source=*',source(1,kk),3)
        enddo
        enddo
        enddo
        nsource=kk
c       
        call prinf('after grid build, nsource=*',nsource,1)
c
        endif
c
c        do i=1,nsource
c        source(1,i)=source(1,i)*10
c        source(2,i)=source(2,i)*10
c        source(3,i)=source(3,i)*10
c        enddo


c       ... shift the sources below the wall z=0
c
        do i=1,nsource
        source(3,i)=source(3,i) - 3
        enddo

c       ... set up the targets
c
        nwall=0
c
        if( idist .eq. 1 .or. idist .eq. 4 ) then
        do i=1,nsource
        target(1,i)=source(1,i) + 0.1
        target(2,i)=source(2,i)
        target(3,i)=source(3,i)
        enddo
        ntarget=nsource
        do i=1,nsource
        target(1,i+nsource)=source(1,i) 
        target(2,i+nsource)=source(2,i) - 0.2
        target(3,i+nsource)=source(3,i)
        enddo
        ntarget=nsource*2
        endif
c
c
        if( idist .eq. 2 ) then
c
c       ... construct target distribution on a curve in R^3
c       
        ntarget=nsource*4
        do i=1,ntarget
        a=2*pi*dble(i)/ntarget
        target(1,i)=sin(1.1*a)/2
        target(2,i)=cos(2.1*a)/2
        target(3,i)=cos(3.1*a)/2
        enddo
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct target distribution on a unit sphere 
c       highly oversampled 
c
        ntarget=nsource*1
        do i=1,ntarget
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        target(1,i)=.5d0*cos(phi)*sin(theta) + 1
        target(2,i)=.5d0*sin(phi)*sin(theta)
        target(3,i)=.5d0*cos(theta) - 2
        d=hkrand(0)
        d=hkrand(0)
        enddo
c
c       ... append a batch of targets on the wall z=0,
c       the velocity must vanish there
c
        nwall=min(nsource,100)
        do i=1,nwall
        target(1,ntarget+i)=2*hkrand(0)-1
        target(2,ntarget+i)=2*hkrand(0)-1
        target(3,ntarget+i)=0
        enddo
        ntarget=ntarget+nwall
        endif
c
        call prinf('ntarget=*',ntarget,1)
c       
ccc        call prin2('source=*',source,3*nsource)
ccc        call prin2('target=*',target,3*ntarget)
        
c
c       define (random) piecewise constant densities
c
        do i=1,nsource
        sigma_sl(1,i)=hkrand(0)
        sigma_sl(2,i)=hkrand(0)
        sigma_sl(3,i)=hkrand(0)
        sigma_dl(1,i)=hkrand(0)
        sigma_dl(2,i)=hkrand(0)
        sigma_dl(3,i)=hkrand(0)
        enddo        
c
        do i=1,nsource
        sigma_dv(1,i)=hkrand(0)
        sigma_dv(2,i)=hkrand(0)
        sigma_dv(3,i)=hkrand(0)
        enddo        
c
c
c
c       turn on single and/or double layer potential
c
        ifsingle=1
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
c       ... test all double layer kernel types against the independent
c       half-space Green's function evaluators:
c       ifdouble = 1 stresslet, 2 symmetric stresslet, 3 rotlet, 4 doublet
c
        do 9000 ifdouble = 1,4
c
        call prinf('ifdouble=*',ifdouble,1)
c
c       ... evaluate via FMM
c
c       itype=1: include both direct arrival and image contribution
c
        itype=1
        iprec=3
c
ccc        call prini(0,13)
        t1=second()
C$        t1=omp_get_wtime()
        call sthfmm3dparttarg
     $     (ier,iprec,itype,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,SIGMA_DV,
     $     ifpot,pot,pre,ifgrad,grad,
     $     ntarget,target,ifpottarg,POTtarg,PREtarg,
     $     ifgradtarg,GRADtarg)
        t2=second()
C$        t2=omp_get_wtime()
        call prini(6,13)
c
        call prin2('after sthfmm3dparttarg, time=*',t2-t1,1)
        call prin2('speed, (sources+targets)/sec=*',
     $     (nsource+ntarget)/(t2-t1),1)
c
c


c       since we are testing by direct calculation, only compute
c       direct calculation at m sources.
c
        m=min(nsource,10)
c
c       ... direct evaluation
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
        if( ifprint .eq. 1 ) then
        call prin2('after sthfmm3dparttarg, pot=*',
     $     pot,3*m)
        call prin2('after sthfmm3dparttarg, grad=*',
     $     grad,3*3*m)
        endif
c
c       
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,itype0,pot0,pre0,grad0)
        do 6550 j=1,m
        do 6540 i=1,nsource
c
c       ... skip the direct arrival for the self-interaction,
c       the image contribution is included
c
        itype0=1
        if( i .eq. j ) itype0=2

        if (ifsingle .eq. 1 ) then
        call green3suph_eval(itype0,source(1,i),
     $     sigma_sl(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot1(1,j)=pot1(1,j)+pot0(1)
        pot1(2,j)=pot1(2,j)+pot0(2)
        pot1(3,j)=pot1(3,j)+pot0(3)
        pre1(j)=pre1(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad1(1,1,j)=grad1(1,1,j)+grad0(1,1)
        grad1(2,1,j)=grad1(2,1,j)+grad0(2,1)
        grad1(3,1,j)=grad1(3,1,j)+grad0(3,1)
        grad1(1,2,j)=grad1(1,2,j)+grad0(1,2)
        grad1(2,2,j)=grad1(2,2,j)+grad0(2,2)
        grad1(3,2,j)=grad1(3,2,j)+grad0(3,2)
        grad1(1,3,j)=grad1(1,3,j)+grad0(1,3)
        grad1(2,3,j)=grad1(2,3,j)+grad0(2,3)
        grad1(3,3,j)=grad1(3,3,j)+grad0(3,3)
        endif
        endif
c
        if (ifdouble .ge. 1) then
        call green3stph_arb_eval(itype0,ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     source(1,j),pot0,pre0,ifgrad,grad0)
        if (ifpot .eq. 1) then
        pot1(1,j)=pot1(1,j)+pot0(1)
        pot1(2,j)=pot1(2,j)+pot0(2)
        pot1(3,j)=pot1(3,j)+pot0(3)
        pre1(j)=pre1(j)+pre0
        endif
        if (ifgrad .eq. 1) then
        grad1(1,1,j)=grad1(1,1,j)+grad0(1,1)
        grad1(2,1,j)=grad1(2,1,j)+grad0(2,1)
        grad1(3,1,j)=grad1(3,1,j)+grad0(3,1)
        grad1(1,2,j)=grad1(1,2,j)+grad0(1,2)
        grad1(2,2,j)=grad1(2,2,j)+grad0(2,2)
        grad1(3,2,j)=grad1(3,2,j)+grad0(3,2)
        grad1(1,3,j)=grad1(1,3,j)+grad0(1,3)
        grad1(2,3,j)=grad1(2,3,j)+grad0(2,3)
        grad1(3,3,j)=grad1(3,3,j)+grad0(3,3)
        endif
        endif
 6540   continue
 6550   continue
C$OMP END PARALLEL DO

        if( ifprint .eq. 1 ) then
        call prin2('after directtarg, pot=*',pot1,3*m)
        call prin2('after directtarg, grad=*',grad1,3*3*m)
        endif
c
        t2=second()
C$        t2=omp_get_wtime()

        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (points/sec)=*',
     $     m/(t2-t1),1)

        if( ifpot .eq. 1 ) then
        call d3error(pot1,pot,3*m,a,r)
ccc        call prin2('absolute error in pot=*',a,1)
        call prin2('relative error in pot=*',r,1)
        call d3error(pre1,pre,m,a,r)
ccc        call prin2('absolute error in pre=*',a,1)
        call prin2('relative error in pre=*',r,1)
        endif
c
        if( ifgrad .eq. 1 ) then
        call d3error(grad1,grad,3*3*m,a,r)
ccc        call prin2('absolute error in grad=*',a,1)
        call prin2('relative error in grad=*',r,1)
        endif
c
c

c       since we are testing by direct calculation, only compute
c       direct calculation at m targets.
c
        m=min(ntarget,10)
c

c       ... direct evaluation
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
        if( ifprint .eq. 1 ) then
        call prin2('after sthfmm3dparttarg, pottarg=*',
     $     pottarg,3*m)
        call prin2('after sthfmm3dparttarg, gradtarg=*',
     $     gradtarg,3*3*m)
        endif
c

c       
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,itype0,pot0,pre0,grad0)
        do j=1,m
        do i=1,nsource
        itype0=1
        if (ifsingle .eq. 1 ) then
        call green3suph_eval(itype0,source(1,i),
     $     sigma_sl(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pot1(1,j)=pot1(1,j)+pot0(1)
        pot1(2,j)=pot1(2,j)+pot0(2)
        pot1(3,j)=pot1(3,j)+pot0(3)
        pre1(j)=pre1(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        grad1(1,1,j)=grad1(1,1,j)+grad0(1,1)
        grad1(2,1,j)=grad1(2,1,j)+grad0(2,1)
        grad1(3,1,j)=grad1(3,1,j)+grad0(3,1)
        grad1(1,2,j)=grad1(1,2,j)+grad0(1,2)
        grad1(2,2,j)=grad1(2,2,j)+grad0(2,2)
        grad1(3,2,j)=grad1(3,2,j)+grad0(3,2)
        grad1(1,3,j)=grad1(1,3,j)+grad0(1,3)
        grad1(2,3,j)=grad1(2,3,j)+grad0(2,3)
        grad1(3,3,j)=grad1(3,3,j)+grad0(3,3)
        endif
        endif
c
        if (ifdouble .ge. 1) then
        call green3stph_arb_eval(itype0,ifdouble,source(1,i),
     $     sigma_dl(1,i),sigma_dv(1,i),
     $     target(1,j),pot0,pre0,ifgradtarg,grad0)
        if (ifpottarg .eq. 1) then
        pot1(1,j)=pot1(1,j)+pot0(1)
        pot1(2,j)=pot1(2,j)+pot0(2)
        pot1(3,j)=pot1(3,j)+pot0(3)
        pre1(j)=pre1(j)+pre0
        endif
        if (ifgradtarg .eq. 1) then
        grad1(1,1,j)=grad1(1,1,j)+grad0(1,1)
        grad1(2,1,j)=grad1(2,1,j)+grad0(2,1)
        grad1(3,1,j)=grad1(3,1,j)+grad0(3,1)
        grad1(1,2,j)=grad1(1,2,j)+grad0(1,2)
        grad1(2,2,j)=grad1(2,2,j)+grad0(2,2)
        grad1(3,2,j)=grad1(3,2,j)+grad0(3,2)
        grad1(1,3,j)=grad1(1,3,j)+grad0(1,3)
        grad1(2,3,j)=grad1(2,3,j)+grad0(2,3)
        grad1(3,3,j)=grad1(3,3,j)+grad0(3,3)
        endif
        endif
        enddo
        enddo
C$OMP END PARALLEL DO

        t2=second()
C$        t2=omp_get_wtime()

        if( ifprint .eq. 1 ) then
        call prin2('after directtarg, pottarg=*',pot1,3*m)
        call prin2('after directtarg, gradtarg=*',grad1,3*3*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)

        if( ifpottarg .eq. 1 ) then
        call d3error(pot1,pottarg,3*m,a,r)
ccc        call prin2('absolute error in target pot=*',a,1)
        call prin2('relative error in target pot=*',r,1)
        call d3error(pre1,pretarg,m,a,r)
ccc        call prin2('absolute error in target pre=*',a,1)
        call prin2('relative error in target pre=*',r,1)
        endif
c
        if( ifgradtarg .eq. 1 ) then
        call d3error(grad1,gradtarg,3*3*m,a,r)
ccc        call prin2('absolute error in target grad=*',a,1)
        call prin2('relative error in target grad=*',r,1)
        endif
c
c       ... check the no-slip condition on the wall z=0
c
        if( nwall .gt. 0 .and. ifpottarg .eq. 1 ) then
        d=0
        do i=ntarget-nwall+1,ntarget
        do j=1,3
        if( abs(pottarg(j,i)) .gt. d ) d=abs(pottarg(j,i))
        enddo
        enddo
        call prin2('max wall velocity, should be zero=*',d,1)
        endif
c
 9000   continue
c
        return
        end
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
