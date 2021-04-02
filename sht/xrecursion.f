cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
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
c  This file contains a suite of evaluation codes for associated vector
c  Legendre functions with various scalings, arguent types, etc.
c  Following is a brief description of the user-callable subroutines.
c  (FORTRAN 77 VERSION).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the code for the evaluation of vector spherical harmonics
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This file contains a set of subroutines for the handling of
c       electromagnetic field multipoles. It contains 9 subroutines that
c       are user-callable. Following is a brief description of these
c       subroutines.
c
c
c       xlgndr - evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the first kind
c               in spherical coordinates
c
c       ulgndr - evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the second kind
c               in spherical coordinates
c
c       xulgndr - evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the first and second kind
c               in spherical coordinates
c
c
c       xlgndr_trunc - 
c               evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the first kind
c               in spherical coordinates, truncated 
c
c       ulgndr_trunc - 
c               evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the second kind
c               in spherical coordinates, truncated 
c
c       xulgndr_trunc - 
c               evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the first and second kind
c               in spherical coordinates, truncated 
c
c
c       zxlgndr - evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the first kind
c               in spherical coordinates. Complex argument
c
c       zulgndr - evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the second kind
c               in spherical coordinates. Complex argument
c
c       zxulgndr - evaluate normalized Legendre functions, their derivatives
c               and the spherical vector functions of the first and second kind
c               in spherical coordinates. Complex argument
c
c
c
c
        subroutine xlgndr(nmax,x,y,d,x2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        complex *16 x2(2,0:nmax,0:nmax)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call ylgndr2s(nmax,x,y,d)
c
        u=-sqrt(1-x*x)
c
        x2(1,0,0)=0
        x2(2,0,0)=0
c
        do n=1,nmax
        x2(1,n,0)=u*d(n,0)
        x2(2,n,0)=0
        enddo
c
        do n=0,nmax
        do m=1,n
        x2(1,n,m)=-d(n,m)
        x2(2,n,m)=-ima*m*y(n,m)
        enddo
        enddo
c
        do n=1,nmax
        scale=1/sqrt(dble(n)*dble(n+1))
        do m=0,n
        x2(1,n,m)=x2(1,n,m)*scale
        x2(2,n,m)=x2(2,n,m)*scale
        enddo
        enddo
c
        scale=sqrt(1-x*x)
        do n=0,nmax
        do m=1,n
        y(n,m)=y(n,m)*scale
        d(n,m)=d(n,m)/scale
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
        subroutine ulgndr(nmax,x,y,d,u2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        complex *16 u2(2,0:nmax,0:nmax)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call xlgndr(nmax,x,y,d,u2)
        call sphrot90i((nmax+1)**2,u2)
c
        return
        end
c
c
c
c
c
        subroutine xulgndr(nmax,x,y,d,x2,u2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        complex *16 x2(2,0:nmax,0:nmax)
        complex *16 u2(2,0:nmax,0:nmax)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call xlgndr(nmax,x,y,d,x2)
        call sphrot90u((nmax+1)**2,x2,u2)
c
        return
        end
c
c
c
c
c
        subroutine xlgndr_trunc(nmax,m2,x,y,d,x2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        complex *16 x2(2,0:nmax,0:nmax)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call ylgndr2s_trunc(nmax,m2,x,y,d)
c
        u=-sqrt(1-x*x)
c
        x2(1,0,0)=0
        x2(2,0,0)=0
c
        do n=1,nmax
        x2(1,n,0)=u*d(n,0)
        x2(2,n,0)=0
        enddo
c
        do n=0,nmax
        do m=1,min(n,m2)
        x2(1,n,m)=-d(n,m)
        x2(2,n,m)=-ima*m*y(n,m)
        enddo
        enddo
c
        do n=1,nmax
        scale=1/sqrt(dble(n)*dble(n+1))
        do m=0,min(n,m2)
        x2(1,n,m)=x2(1,n,m)*scale
        x2(2,n,m)=x2(2,n,m)*scale
        enddo
        enddo
c
        scale=sqrt(1-x*x)
        do n=0,nmax
        do m=1,min(n,m2)
        y(n,m)=y(n,m)*scale
        d(n,m)=d(n,m)/scale
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
        subroutine ulgndr_trunc(nmax,m2,x,y,d,u2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        complex *16 u2(2,0:nmax,0:nmax)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call xlgndr_trunc(nmax,m2,x,y,d,u2)
        call sphrot90i((nmax+1)**2,u2)
c
        return
        end
c
c
c
c
c
        subroutine xulgndr_trunc(nmax,m2,x,y,d,x2,u2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        complex *16 x2(2,0:nmax,0:nmax)
        complex *16 u2(2,0:nmax,0:nmax)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call xlgndr_trunc(nmax,m2,x,y,d,x2)
        call sphrot90u((nmax+1)**2,x2,u2)
c
        return
        end
c
c
c
c
c
        subroutine sphrot90i(nnout,unm2out)
        implicit real *8 (a-h,o-z)
        complex *16 unm2out(2,1)
        complex *16 cd
c
        do i=1,nnout
        cd=unm2out(1,i)
        unm2out(1,i)=-unm2out(2,i)
        unm2out(2,i)=cd
        enddo
c       
        return
        end
c
c
c
c
c
        subroutine sphrot90u(nnout,xnm2,unm2out)
        implicit real *8 (a-h,o-z)
        complex *16 xnm2(2,1),unm2out(2,1)
c
        do i=1,nnout
        unm2out(1,i)=-xnm2(2,i)
        unm2out(2,i)=+xnm2(1,i)
        enddo
c       
        return
        end
c
c
c
c
c
        subroutine zxlgndr(nmax,z,y,d,x2)
        implicit real *8 (a-h,o-z)
        complex *16 y(0:nmax,0:nmax)
        complex *16 d(0:nmax,0:nmax)
        complex *16 x2(2,0:nmax,0:nmax)
        complex *16 z,u,cscale
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call zylgndr2s(nmax,z,y,d)
c
        u=-sqrt(1-z*z)
c
        x2(1,0,0)=0
        x2(2,0,0)=0
c
        do n=1,nmax
        x2(1,n,0)=u*d(n,0)
        x2(2,n,0)=0
        enddo
c
        do n=0,nmax
        do m=1,n
        x2(1,n,m)=-d(n,m)
        x2(2,n,m)=-ima*m*y(n,m)
        enddo
        enddo
c
        do n=1,nmax
        scale=1/sqrt(dble(n)*dble(n+1))
        do m=0,n
        x2(1,n,m)=x2(1,n,m)*scale
        x2(2,n,m)=x2(2,n,m)*scale
        enddo
        enddo
c
        cscale=sqrt(1-z*z)
        do n=0,nmax
        do m=1,n
        y(n,m)=y(n,m)*cscale
        d(n,m)=d(n,m)/cscale
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
        subroutine zulgndr(nmax,z,y,d,u2)
        implicit real *8 (a-h,o-z)
        complex *16 y(0:nmax,0:nmax)
        complex *16 d(0:nmax,0:nmax)
        complex *16 u2(2,0:nmax,0:nmax)
        complex *16 z
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call zxlgndr(nmax,z,y,d,u2)
        call sphrot90i((nmax+1)**2,u2)
c
        return
        end
c
c
c
c
c
        subroutine zxulgndr(nmax,z,y,d,x2,u2)
        implicit real *8 (a-h,o-z)
        complex *16 y(0:nmax,0:nmax)
        complex *16 d(0:nmax,0:nmax)
        complex *16 x2(2,0:nmax,0:nmax)
        complex *16 u2(2,0:nmax,0:nmax)
        complex *16 z
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call zxlgndr(nmax,z,y,d,x2)
        call sphrot90u((nmax+1)**2,x2,u2)
c
        return
        end
c
c
c
c
c
