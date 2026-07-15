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
c    $Date: 2012-04-10 19:48:49 -0400 (Tue, 10 Apr 2012) $
c    $Revision: 2891 $
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of half-space Stokes Green's functions in R^3,
c       that satisfy u=0 at the half space interface z=0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c**********************************************************************C
        subroutine green3suph_brute_v1(source,target,df,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Naive implementation of Stokes SLP with zero boundary condition
c       on lower half-space.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),uvec(3),df(3)
        dimension source(3),target(3),delta(3,3),rr(3),df_image(3)
c
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
c
c       === PART 1 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df(j)
        enddo
        enddo
c
        uout(1)=uvec(1) * (+0.5d0)
        uout(2)=uvec(2) * (+0.5d0)
        uout(3)=uvec(3) * (+0.5d0)
c
        pout=(dx*df(1)+dy*df(2)+dz*df(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df_image(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df_image(j)
        enddo
        enddo
c
        uout(1)=uout(1)-uvec(1) * (+0.5d0)
        uout(2)=uout(2)-uvec(2) * (+0.5d0)
        uout(3)=uout(3)-uvec(3) * (+0.5d0)
c
        pout=pout-(dx*df_image(1)+dy*df_image(2)+dz*df_image(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c
c
c       === PART 3 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
c
        uvec(i)=uvec(i)+(-delta(i,j)/cd**3) *df_image(j)
        uvec(i)=uvec(i)+(3*rr(i)*rr(j)/cd**5) *df_image(j)
c
        enddo
        enddo
c
        uout(1)=uout(1)+source(3)*target(3)*uvec(1)
        uout(2)=uout(2)+source(3)*target(3)*uvec(2)
        uout(3)=uout(3)+source(3)*target(3)*uvec(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(3)=uout(3)+source(3)*(rr(1)/cd**3) *df_image(1)
        uout(3)=uout(3)+source(3)*(rr(2)/cd**3) *df_image(2)
        uout(3)=uout(3)+source(3)*(rr(3)/cd**3) *df_image(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(1)=uout(1)-target(3)*(rr(1)/cd**3) *df(3)
        uout(2)=uout(2)-target(3)*(rr(2)/cd**3) *df(3)
        uout(3)=uout(3)-target(3)*(rr(3)/cd**3) *df(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(3)=uout(3)-1/cd *df(3)
c
c
        pout = pout + source(3)*uvec(3) *2
        pout = pout - rr(3)/cd**3*df(3) *2
c
ccc        call prin2('pout=*',pout,1)
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
c
        return
        end
c
c
c
c
c
c 
c**********************************************************************C
        subroutine green3suph_brute_v2(source,target,df,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Naive implementation of Stokes SLP with zero boundary condition
c       on lower half-space, yet another decomposition.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),uvec(3),df(3)
        dimension source(3),target(3),delta(3,3),rr(3),df_image(3)
c
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
        uout(1)=0
        uout(2)=0
        uout(3)=0
        pout=0

        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        
c
c       === PART 1 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df(j)
        enddo
        enddo
c
        uout(1)=uvec(1) * (+0.5d0)
        uout(2)=uvec(2) * (+0.5d0)
        uout(3)=uvec(3) * (+0.5d0)
c
        pout=+(dx*df(1)+dy*df(2)+dz*df(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
 5200   continue
c
c       === PART 2 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
        uvec(i)=uvec(i)+delta(i,j)/cd *df(j)
        uvec(i)=uvec(i)+rr(i)*rr(j)/cd**3 *df(j)
        enddo
        enddo
c
        uout(1)=uout(1)-uvec(1) * (+0.5d0)
        uout(2)=uout(2)-uvec(2) * (+0.5d0)
        uout(3)=uout(3)-uvec(3) * (+0.5d0)
c
        pout=pout-(dx*df(1)+dy*df(2)+dz*df(3))/cd**3
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c
 5300   continue
c
c       === PART 3 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
c
        do i=1,3
        do j=1,3
c
        uvec(i)=uvec(i)+(-delta(i,j)/cd**3) *df_image(j)
        uvec(i)=uvec(i)+(3*rr(i)*rr(j)/cd**5) *df_image(j)
c
        enddo
        enddo
c
        uout(1)=uout(1)+source(3)*rr(3)*uvec(1)
        uout(2)=uout(2)+source(3)*rr(3)*uvec(2)
        uout(3)=uout(3)+source(3)*rr(3)*uvec(3)
c
        uout(1)=uout(1)-source(3)*source(3)*uvec(1)
        uout(2)=uout(2)-source(3)*source(3)*uvec(2)
        uout(3)=uout(3)-source(3)*source(3)*uvec(3)
c
        uout(3)=uout(3)+source(3)*(rr(1)/cd**3) *df_image(1)
        uout(3)=uout(3)+source(3)*(rr(2)/cd**3) *df_image(2)
        uout(3)=uout(3)+source(3)*(rr(3)/cd**3) *df_image(3)
c
ccc        call prin2('uout=*',uout,3)
c
        uout(1)=uout(1)+source(3)*(rr(1)/cd**3) *df(3)
        uout(2)=uout(2)+source(3)*(rr(2)/cd**3) *df(3)
        uout(3)=uout(3)+source(3)*(rr(3)/cd**3) *df(3)
c
ccc        call prin2('uout=*',uout,3)
c
        pout = pout + source(3)*uvec(3) *2
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
c
        return
        end
c
c
c**********************************************************************C
        subroutine green3suph_brute(source,target,df,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes SLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),df(3)
        dimension source(3),target(3),df_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,pot,fld(3)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3sup(xyz,df,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... stokeslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        call green3sup(xyz,df_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
cc        call prin2('uout=*',uout,3)
cc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... dipole image
c
        dipstr = source(3)
        call lpotfld3d_dp(iffld,sourceim,dipstr,df_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... charge image
c
        charge = df_image(3)
        call lpotfld3d(iffld,sourceim,charge,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP, stresslet (type 1),
c       with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c       Same as green3stph_brute1.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... dlp stresslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
c       ... dipole image
c
        dipstr = rnorm_image(1)*du_image(1)+
     $           rnorm_image(2)*du_image(2)+
     $           rnorm_image(3)*du_image(3)
        dipstr = dipstr*2
        dipvec(1)=0
        dipvec(2)=0
        dipvec(3)=1
        call lpotfld3d_dp(iffld,sourceim,dipstr,dipvec,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_arb_brute(
     $     ifdouble,source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c       Dispatcher over the double layer kernel types:
c       ifdouble = 1: stresslet (green3stph_brute1)
c       ifdouble = 2: symmetric stresslet (green3stph_brute2)
c       ifdouble = 3: rotlet (green3stph_brute3)
c       ifdouble = 4: doublet (green3stph_brute4)
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       ifdouble        Double layer kernel type, see above.
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
        if( ifdouble .eq. 1 ) then
        call green3stph_brute1
     $     (source,target,du,rnorm,uout,pout)
        endif

        if( ifdouble .eq. 2 ) then
        call green3stph_brute2
     $     (source,target,du,rnorm,uout,pout)
        endif

        if( ifdouble .eq. 3 ) then
        call green3stph_brute3
     $     (source,target,du,rnorm,uout,pout)
        endif

        if( ifdouble .eq. 4 ) then
        call green3stph_brute4
     $     (source,target,du,rnorm,uout,pout)
        endif

c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute1(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP, stresslet (type 1),
c       with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... dlp stresslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
c       ... dipole image
c
        dipstr = rnorm_image(1)*du_image(1)+
     $           rnorm_image(2)*du_image(2)+
     $           rnorm_image(3)*du_image(3)
        dipstr = dipstr*2
        dipvec(1)=0
        dipvec(2)=0
        dipvec(3)=1
        call lpotfld3d_dp(iffld,sourceim,dipstr,dipvec,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute2(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP, symmetric stresslet (type 2),
c       with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c       The image system is quadrupole only: the Papkovich-Neuber
c       dipole corrections of type 1 and of the potential source
c       [-r_i/r^3](n.g) cancel exactly.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6),strain(3,3)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp_stresslet_sym(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... stresslet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_stresslet_sym(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute3(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP, rotlet (type 3),
c       with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp_rotlet(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... rotlet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_rotlet(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... dipole image
c
        dipstr = -rnorm_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,du_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... dipole image
c
        dipstr = du_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,rnorm_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_brute4(source,target,du,rnorm,uout,pout)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP, doublet (type 4),
c       with zero boundary condition on lower half-space,
c       check calls to free space Stokes and Laplace kernels.
c       Doublet = symmetric stresslet + rotlet, so this field
c       must equal green3stph_brute2 + green3stph_brute3.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c
        dimension xyz(3),rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3)
        complex *16 charge,dipstr,quadstr,pot,fld(3)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
        call green3stp_doublet(xyz,du,rnorm,uout,pout)
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
c       ... doublet image
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)+source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_doublet(xyz,du_image,rnorm_image,fvec,pvec)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
ccc        call prin2('uout=*',uout,3)
ccc        call prin2('pout=*',pout,1)
ccc        if (2.ne.3) return
c
c       === PART 3 ===   
c
        iffld = 1
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3d_quad(iffld,sourceim,quadvec,target,
     1                        pot,fld)
        fld(1)=fld(1)*quadstr
        fld(2)=fld(2)*quadstr
        fld(3)=fld(3)*quadstr
        pot=pot*quadstr
c
        uout(1)=uout(1)+target(3)*dreal(fld(1)) 
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... dipole image
c
        dipstr = -rnorm_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,du_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c       ... dipole image
c
        dipstr = du_image(3)*2
        call lpotfld3d_dp(iffld,sourceim,dipstr,rnorm_image,target,
     1                        pot,fld)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3suph_eval(itype,source,df,target,uout,pout,
     $     ifgrad,grad)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes SLP with zero boundary condition on lower half-space,
c       direct evaluation of velocity, pressure, and velocity gradient.
c       The image system matches green3suph_brute.
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       itype           Half space Green's function evaluation flag
c                       1 => include both direct arrival and
c                            image contribution
c                       2 => include image contribution only
c       source(3)       Source location in lower half-space (z<0).
c       df(3)           Strength of single force source
c       target(3)       Target point in lower half-space (z<0).
c       ifgrad          Velocity gradient computation flag
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c       grad (real *8) - the velocity gradient at the target,
c                        grad(i,j) = du_i/dx_j, if ifgrad=1
c
        dimension uout(3),fvec(3),df(3)
        dimension source(3),target(3),df_image(3)
        dimension sourceim(3),grad(3,3),dgrad(3,3)
        complex *16 charge,dipstr,pot,fld(3),hess(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        if( itype .eq. 1 ) then
        call green3sup_eval(source,df,target,uout,pout,ifgrad,grad)
        else
        uout(1)=0
        uout(2)=0
        uout(3)=0
        pout=0
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo
        endif
        endif
c
c
c       === PART 2 ===
c
c       ... stokeslet image
c
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
        df_image(1)=+df(1)
        df_image(2)=+df(2)
        df_image(3)=-df(3)
        call green3sup_eval(sourceim,df_image,target,fvec,pvec,
     $     ifgrad,dgrad)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)-dgrad(i,j)
        enddo
        enddo
        endif
c
c       === PART 3 ===   
c
        iffld = 1
        ifhess = ifgrad
c
c       ... dipole image
c
        dipstr = source(3)
        call lpotfld3dhess_dp(iffld,ifhess,sourceim,dipstr,df_image,
     1                        target,pot,fld,hess)
        call green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c
c       ... charge image
c
        charge = df_image(3)
        call lpotfld3dhess(iffld,ifhess,sourceim,charge,
     1                        target,pot,fld,hess)
        call green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3stph_arb_eval(itype,
     $     ifdouble,source,du,rnorm,target,uout,pout,ifgrad,grad)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Stokes DLP with zero boundary condition on lower half-space,
c       direct evaluation of velocity, pressure, and velocity gradient.
c       The image systems match green3stph_arb_brute:
c       ifdouble = 1: stresslet
c       ifdouble = 2: symmetric stresslet
c       ifdouble = 3: rotlet
c       ifdouble = 4: doublet
c
c       Half-space boundary condition is assumed,
c       u = 0 at z=0.
c
c       INPUT:
c
c       itype           Half space Green's function evaluation flag
c                       1 => include both direct arrival and
c                            image contribution
c                       2 => include image contribution only
c       ifdouble        Double layer kernel type, see above.
c       source(3)       Source location in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of double force source
c       target(3)       Target point in lower half-space (z<0).
c       ifgrad          Velocity gradient computation flag
c
c       OUTPUT:
c
c       uout (real *8) - the velocity field at the target
c       pout (real *8) - the pressure at the target
c       grad (real *8) - the velocity gradient at the target,
c                        grad(i,j) = du_i/dx_j, if ifgrad=1
c
        dimension rnorm(3),uout(3),fvec(3),du(3)
        dimension source(3),target(3),du_image(3),rnorm_image(3)
        dimension sourceim(3),grad(3,3),dgrad(3,3)
        complex *16 charge,dipstr,quadstr,pot,fld(3),hess(6)
        dimension dipvec(3),quadvec(6)
c
c
c       === PART 1 ===
c
c       ... direct arrival
c
        if( itype .eq. 1 ) then
        call green3stp_arb_eval(ifdouble,source,du,rnorm,
     $     target,uout,pout,ifgrad,grad)
        else
        uout(1)=0
        uout(2)=0
        uout(3)=0
        pout=0
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=0
        enddo
        enddo
        endif
        endif
c
c
c       === PART 2 ===
c
c       ... double layer image
c
        sourceim(1) = source(1)
        sourceim(2) = source(2)
        sourceim(3) = -source(3)
        du_image(1)=+du(1)
        du_image(2)=+du(2)
        du_image(3)=-du(3)
        rnorm_image(1)=+rnorm(1)
        rnorm_image(2)=+rnorm(2)
        rnorm_image(3)=-rnorm(3)
        call green3stp_arb_eval(ifdouble,sourceim,du_image,rnorm_image,
     $     target,fvec,pvec,ifgrad,dgrad)
c
        uout(1) = uout(1) - fvec(1)
        uout(2) = uout(2) - fvec(2)
        uout(3) = uout(3) - fvec(3)
        pout = pout - pvec
c
        if( ifgrad .eq. 1 ) then
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)-dgrad(i,j)
        enddo
        enddo
        endif
c
c       === PART 3 ===   
c
        iffld = 1
        ifhess = ifgrad
c
        if( ifdouble .eq. 1 .or. ifdouble .eq. 2
     $     .or. ifdouble .eq. 4 ) then
c
c       ... quadrupole image
c
        quadstr = source(3)*2
        quadvec(1)=du_image(1)*rnorm_image(1)
        quadvec(2)=du_image(2)*rnorm_image(2)
        quadvec(3)=du_image(3)*rnorm_image(3)
        quadvec(4)=du_image(1)*rnorm_image(2)
        quadvec(5)=du_image(1)*rnorm_image(3)
        quadvec(6)=du_image(2)*rnorm_image(3)
        quadvec(4)=quadvec(4)+du_image(2)*rnorm_image(1)
        quadvec(5)=quadvec(5)+du_image(3)*rnorm_image(1)
        quadvec(6)=quadvec(6)+du_image(3)*rnorm_image(2)
        call lpotfld3dhess_qp(iffld,ifhess,sourceim,quadstr,quadvec,
     1                        target,pot,fld,hess)
        call green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c
        endif
c
        if( ifdouble .eq. 1 ) then
c
c       ... dipole image
c
        dipstr = rnorm_image(1)*du_image(1)+
     $           rnorm_image(2)*du_image(2)+
     $           rnorm_image(3)*du_image(3)
        dipstr = dipstr*2
        dipvec(1)=0
        dipvec(2)=0
        dipvec(3)=1
        call lpotfld3dhess_dp(iffld,ifhess,sourceim,dipstr,dipvec,
     1                        target,pot,fld,hess)
        call green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c
        endif
c
        if( ifdouble .eq. 3 .or. ifdouble .eq. 4 ) then
c
c       ... dipole image
c
        dipstr = -rnorm_image(3)*2
        call lpotfld3dhess_dp(iffld,ifhess,sourceim,dipstr,du_image,
     1                        target,pot,fld,hess)
        call green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c
c       ... dipole image
c
        dipstr = du_image(3)*2
        call lpotfld3dhess_dp(iffld,ifhess,sourceim,dipstr,rnorm_image,
     1                        target,pot,fld,hess)
        call green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c
        endif
c
c
        return
        end
c
c
c
c
c**********************************************************************C
        subroutine green3pnh_eval_add(target,pot,fld,hess,
     $     uout,pout,ifgrad,grad)
c**********************************************************************C
        implicit real *8 (a-h,o-z)
c
c
c       Accumulate the Papkovich-Neuber correction
c
c       u_i = z fld_i + delta_i3 pot,  p = 2 fld_3
c
c       and, if ifgrad=1, its target gradient
c
c       du_i/dx_j = delta_j3 fld_i - z hess(i,j) - delta_i3 fld_j
c
c       for a single Laplace charge/dipole/quadrupole image source,
c       given pot, fld = -grad(pot) and 
c       hess = (potxx,potyy,potzz,potxy,potxz,potyz) at the target.
c
        dimension target(3),uout(3),grad(3,3)
        dimension hmatr(3,3)
        complex *16 pot,fld(3),hess(6)
c
        uout(1)=uout(1)+target(3)*dreal(fld(1))
        uout(2)=uout(2)+target(3)*dreal(fld(2))
        uout(3)=uout(3)+target(3)*dreal(fld(3))
        uout(3)=uout(3)+dreal(pot)
        pout = pout + 2*dreal(fld(3))
c
        if( ifgrad .eq. 1 ) then
c
        hmatr(1,1)=dreal(hess(1))
        hmatr(2,2)=dreal(hess(2))
        hmatr(3,3)=dreal(hess(3))
        hmatr(1,2)=dreal(hess(4))
        hmatr(1,3)=dreal(hess(5))
        hmatr(2,3)=dreal(hess(6))
        hmatr(2,1)=hmatr(1,2)
        hmatr(3,1)=hmatr(1,3)
        hmatr(3,2)=hmatr(2,3)
c
        do i=1,3
        do j=1,3
        grad(i,j)=grad(i,j)-target(3)*hmatr(i,j)
        enddo
        grad(i,3)=grad(i,3)+dreal(fld(i))
        grad(3,i)=grad(3,i)-dreal(fld(i))
        enddo
c
        endif
c
        return
        end
