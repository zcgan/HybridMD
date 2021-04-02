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
c    $Date: 2010-01-05 12:57:52 -0500 (Tue, 05 Jan 2010) $
c    $Revision: 782 $
c
c
c     ROTATION VIA RECURSION  (FORTRAN 77 AND 90 VERSIONS).
c     good up to order 100 or so.
c
c     Multipole expansions for complex-valued functions.
c
c       rotviarecur3 is an f77 routine (needs a workspace).
c       rotviarecur3f90 is an f90 routine (no workspace needed).
c
c       rotviarecur3p_init  - precompute the coefficients of rotation matrix.
c       rotviarecur3p_apply - apply the rotation matrix.
c
c
C*****************************************************************
        subroutine rotviarecur3(theta,nterms,m1,m2,mpole,ld1,marray,
     1                         ld2,w,lw,lused)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c       The rotation matrices for each order (first index) are computed
c       from the lowest to the highest. As each one is generated, it
c       is applied to the input expansion "mpole" and overwritten.
c
c       As a result, it is sufficient to use two arrays rd1 and rd2 for
c       the two term recurrence, rather than storing them for all orders
c       as in the original code. There is some loss in speed
c       if the rotation operator is to be used multiple times, but the 
c       memory savings is often more critical.
c
c       Note that this is just a driver to rotviarecur3s and
c       rotviarecur3t routines, the main purpose of this routine is to
c       have the calling sequence compatible with rotproj.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multiple expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       w     :  work array 
c       lw    :  length of work array 
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c       lused    amount of workspace used.
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2,lw
        integer ldc,ird1,lrd,ird2,isqc,lused
        real *8 theta,w(lw)
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
c
        ldc = nterms
        ird1 = 1
        lrd  = (ldc+1)*(2*ldc+1) + 3
        ird2 = ird1 + lrd
        isqc = ird2 + lrd
        lused = isqc + 2*(2*ldc+1)
        if (lused.gt.lw) then
           write(6,*)' workspace failure in rotviarecur3'
           stop
        endif
c
c
        if( m1 .lt. nterms .or. m2 .lt. nterms ) then
        call rotviarecur3t(nterms,m1,m2,mpole,ld1,marray,
     1                    ld2,w(ird1),w(ird2),w(isqc),theta,ldc)
        else
        call rotviarecur3s(nterms,m1,m2,mpole,ld1,marray,
     1                    ld2,w(ird1),w(ird2),w(isqc),theta,ldc)
        endif
c
        return
        end
c        
c
C*****************************************************************
        subroutine rotviarecur3t(nterms,m1,m2,mpole,ld1,marray,
     1                         ld2,rd1,rd2,sqc,theta,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c       The rotation matrices for each order (first index) are computed
c       from the lowest to the highest. As each one is generated, it
c       is applied to the input expansion "mpole" and overwritten.
c
c       As a result, it is sufficient to use two arrays rd1 and rd2 for
c       the two term recurrence, rather than storing them for all orders
c       as in the original code. There is some loss in speed
c       if the rotation operator is to be used multiple times, but the 
c       memory savings is often more critical.
c
c       Use symmetry properties of rotation matrices
c
C---------------------------------------------------------------------
c       INPUT:
c
c       nterms: dimension parameter for d - the rotation matrix.
c       m1    : max m index for first expansion. 
c       m2    : max m index for second expansion.
C       mpole   coefficients of original multiple expansion
C       rd1     work space 
C       rd2     work space
c       sqc:    an array contains the square roots of the
c               binomial coefficients.
c       theta:  the rotate angle about the y-axis.
c       ldc     dimensions of sqc array
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2,ldc
        integer ij,im,imp,m,mp
        real *8 theta
	real *8 rd1(0:ldc,-ldc:ldc)
	real *8 rd2(0:ldc,-ldc:ldc)
	real *8 sqc(0:2*ldc,2)
        real *8 ww,done,ctheta,stheta,hsthta,cthtap,cthtan,d
        real *8 precis,scale
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
c
	data precis/1.0d-20/
c
        ww=1/sqrt(2.0d0)
        do m = 0, 2*ldc 
	   sqc(m,1) = dsqrt(m+0.0d0)
        enddo
	sqc(0,2) = 0.0d0
	if( ldc .gt. 0 ) sqc(1,2) = 0.0d0
        do m = 2, 2*ldc 
	   sqc(m,2) = dsqrt((m+0.0d0)*(m-1)/2.0d0)
        enddo
c
	done=1
	ctheta=dcos(theta)
	if (dabs(ctheta).le.precis) ctheta=0.0d0
	stheta=dsin(-theta)
	if (dabs(stheta).le.precis) stheta=0.0d0
	hsthta=ww*stheta
	cthtap=+2.0d0*ww*dcos(theta/2.0d0)**2
	cthtan=-2.0d0*ww*dsin(theta/2.0d0)**2
c
c     Compute the (0,0,0) term.
c
      rd1(0,0)=done
      marray(0,0)=mpole(0,0)*rd1(0,0)
c
c ... Loop over first index ij=1,nterms, constructing
c     rotation matrices recursively.
c
      do ij=1,nterms
c
c     For mprime=0, use formula (1).
c
         do im=-ij,-1
	    rd2(0,im)=-sqc(ij-im,2)*rd1(0,im+1)
	    if (im.gt.(1-ij)) then
	       rd2(0,im)=rd2(0,im)+sqc(ij+im,2)*rd1(0,im-1)
	    endif
	    rd2(0,im)=rd2(0,im)*hsthta
	    if (im.gt.-ij) then
	       rd2(0,im)=rd2(0,im)+
     1         rd1(0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
	    endif
	    rd2(0,im)=rd2(0,im)/ij
         enddo
c
	 rd2(0,0)=rd1(0,0)*ctheta
	 if (ij.gt.1) then
	    rd2(0,0)=rd2(0,0)+hsthta*sqc(ij,2)*(2*rd1(0,-1))/ij
	 endif
c         
	 do im=1,ij
	    rd2(0,im)=rd2(0,-im)
            if( mod(im,2) .eq. 0 ) then
            rd2(im,0)=+rd2(0,im)
            else
            rd2(im,0)=-rd2(0,im)
            endif
         enddo
c
c        For 0<mprime<=j (2nd index) case, use formula (2).
c
	 do imp=1,min(ij,min(m1,m2))         
c
            scale=(ww/sqc(ij+imp,2))
c
	    do im=imp,ij
c
	       rd2(imp,+im)=rd1(imp-1,+im-1)*(cthtap*sqc(ij+im,2))
	       rd2(imp,-im)=rd1(imp-1,-im+1)*(cthtan*sqc(ij+im,2))
c
	       if (im.lt.(ij-1)) then
	          rd2(imp,+im)=rd2(imp,+im)-rd1(imp-1,+im+1)*
     $               (cthtan*sqc(ij-im,2))
	          rd2(imp,-im)=rd2(imp,-im)-rd1(imp-1,-im-1)*
     $               (cthtap*sqc(ij-im,2))
	       endif
c
	       if (im.lt.ij) then
                  d=(stheta*sqc(ij+im,1)*sqc(ij-im,1))
	          rd2(imp,+im)=rd2(imp,+im)+rd1(imp-1,+im)*d
	          rd2(imp,-im)=rd2(imp,-im)+rd1(imp-1,-im)*d
	       endif
c
	       rd2(imp,+im)=rd2(imp,+im)*scale
	       rd2(imp,-im)=rd2(imp,-im)*scale
c
	       if (im.gt.imp) then
               if( mod(imp+im,2) .eq. 0 ) then
               rd2(im,+imp)=+rd2(imp,+im)
               rd2(im,-imp)=+rd2(imp,-im)
               else
               rd2(im,+imp)=-rd2(imp,+im)
               rd2(im,-imp)=-rd2(imp,-im)
               endif
               endif
c
            enddo
c
         enddo
c
         do m=-ij,ij
            marray(ij,m)=0
         enddo
c
         do m=-min(ij,m2),min(ij,m2)
            marray(ij,m)=mpole(ij,0)*rd2(0,m) 
            do mp=1,min(ij,m1)
               marray(ij,m)=marray(ij,m)+
     1	       mpole(ij,mp)*rd2(mp,m)+
     1         mpole(ij,-mp)*rd2(mp,-m)
            enddo
         enddo
c
         do m=-ij,ij
            do mp=0,min(ij,min(m1,m2))
	       rd1(mp,m) = rd2(mp,m)
            enddo
         enddo
      enddo
      return
      end
c
c
C*****************************************************************
        subroutine rotviarecur3s(nterms,m1,m2,mpole,ld1,marray,
     1                         ld2,rd1,rd2,sqc,theta,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c       The rotation matrices for each order (first index) are computed
c       from the lowest to the highest. As each one is generated, it
c       is applied to the input expansion "mpole" and overwritten.
c
c       As a result, it is sufficient to use two arrays rd1 and rd2 for
c       the two term recurrence, rather than storing them for all orders
c       as in the original code. There is some loss in speed
c       if the rotation operator is to be used multiple times, but the 
c       memory savings is often more critical.
c
c       Use symmetry properties of rotation matrices
c
c       Parameters m1 and m2 are NOT USED by this routine, use
c       rotviarecur3t instead
c
C---------------------------------------------------------------------
c       INPUT:
c
c       nterms: dimension parameter for d - the rotation matrix.
c       m1    : max m index for first expansion. NOT USED in this routine
c       m2    : max m index for second expansion. NOT USED in this routine
C       mpole   coefficients of original multiple expansion
C       rd1     work space 
C       rd2     work space
c       sqc:    an array contains the square roots of the
c               binomial coefficients.
c       theta:  the rotate angle about the y-axis.
c       ldc     dimensions of sqc array
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2,ldc
        integer ij,im,imp,m,mp
        real *8 theta
        real *8 ww,done,ctheta,stheta,hsthta,cthtap,cthtan,d
        real *8 precis,scale
	real *8 rd1(0:ldc,-ldc:ldc)
	real *8 rd2(0:ldc,-ldc:ldc)
	real *8 sqc(0:2*ldc,2)
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
c
	data precis/1.0d-20/
c
        ww=1/sqrt(2.0d0)
        do m = 0, 2*ldc 
	   sqc(m,1) = dsqrt(m+0.0d0)
        enddo
	sqc(0,2) = 0.0d0
	if( ldc .gt. 0 ) sqc(1,2) = 0.0d0
        do m = 2, 2*ldc 
	   sqc(m,2) = dsqrt((m+0.0d0)*(m-1)/2.0d0)
        enddo
	done=1
	ctheta=dcos(theta)
	if (dabs(ctheta).le.precis) ctheta=0.0d0
	stheta=dsin(-theta)
	if (dabs(stheta).le.precis) stheta=0.0d0
	hsthta=ww*stheta
	cthtap=+2.0d0*ww*dcos(theta/2.0d0)**2
	cthtan=-2.0d0*ww*dsin(theta/2.0d0)**2
c
c     Compute the (0,0,0) term.
c
      rd1(0,0)=done
      marray(0,0)=mpole(0,0)*rd1(0,0)
c
c ... Loop over first index ij=1,nterms, constructing
c     rotation matrices recursively.
c
      do ij=1,nterms
c
c     For mprime=0, use formula (1).
c
         do im=-ij,-1
	    rd2(0,im)=-sqc(ij-im,2)*rd1(0,im+1)
	    if (im.gt.(1-ij)) then
	       rd2(0,im)=rd2(0,im)+sqc(ij+im,2)*rd1(0,im-1)
	    endif
	    rd2(0,im)=rd2(0,im)*hsthta
	    if (im.gt.-ij) then
	       rd2(0,im)=rd2(0,im)+
     1         rd1(0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
	    endif
	    rd2(0,im)=rd2(0,im)/ij
         enddo
c
	 rd2(0,0)=rd1(0,0)*ctheta
	 if (ij.gt.1) then
	    rd2(0,0)=rd2(0,0)+hsthta*sqc(ij,2)*(2*rd1(0,-1))/ij
	 endif
c         
	 do im=1,ij
	    rd2(0,im)=rd2(0,-im)
            if( mod(im,2) .eq. 0 ) then
            rd2(im,0)=+rd2(0,im)
            else
            rd2(im,0)=-rd2(0,im)
            endif
         enddo
c
c        For 0<mprime<=j (2nd index) case, use formula (2).
c
	 do imp=1,ij
c
            scale=(ww/sqc(ij+imp,2))
c
	    do im=imp,ij
c
	       rd2(imp,+im)=rd1(imp-1,+im-1)*(cthtap*sqc(ij+im,2))
	       rd2(imp,-im)=rd1(imp-1,-im+1)*(cthtan*sqc(ij+im,2))
c
	       if (im.lt.(ij-1)) then
	          rd2(imp,+im)=rd2(imp,+im)-rd1(imp-1,+im+1)*
     $               (cthtan*sqc(ij-im,2))
	          rd2(imp,-im)=rd2(imp,-im)-rd1(imp-1,-im-1)*
     $               (cthtap*sqc(ij-im,2))
	       endif
c
	       if (im.lt.ij) then
                  d=(stheta*sqc(ij+im,1)*sqc(ij-im,1))
	          rd2(imp,+im)=rd2(imp,+im)+rd1(imp-1,+im)*d
	          rd2(imp,-im)=rd2(imp,-im)+rd1(imp-1,-im)*d
	       endif
c
	       rd2(imp,+im)=rd2(imp,+im)*scale
	       rd2(imp,-im)=rd2(imp,-im)*scale
c
	       if (im.gt.imp) then
               if( mod(imp+im,2) .eq. 0 ) then
               rd2(im,+imp)=+rd2(imp,+im)
               rd2(im,-imp)=+rd2(imp,-im)
               else
               rd2(im,+imp)=-rd2(imp,+im)
               rd2(im,-imp)=-rd2(imp,-im)
               endif
               endif
c
            enddo
c
         enddo
c
         do m=-ij,ij
            marray(ij,m)=mpole(ij,0)*rd2(0,m)
            do mp=1,ij
               marray(ij,m)=marray(ij,m)+
     1	       mpole(ij,mp)*rd2(mp,m)+
     1         mpole(ij,-mp)*rd2(mp,-m)
            enddo
         enddo
c
         do m=-ij,ij
            do mp=0,ij
	       rd1(mp,m) = rd2(mp,m)
            enddo
         enddo
c
c         if( ij.eq.nterms ) then
c         call prinf('rd1=*',ij,1)
c         do m=-ij,ij
c         call prin2('=*',rd1(0,m),ij+1)
c         enddo
c         do m=0,ij
c         call prin2('_*',rd1(0,m),m+1)
c         enddo
c         endif
c
      enddo
      return
      end
c
c
C*****************************************************************
        subroutine rotviarecur3x(nterms,m1,m2,mpole,ld1,marray,
     1                         ld2,rd1,rd2,sqc,theta,ldc)
C*****************************************************************
c
c       WARNING: this code fails for values m2 .ne. nterms. DO NOT USE!
c
c       Purpose:
c
c       Fast, recursive method for applying rotation matrix about
c       the y-axis determined by angle theta.
c       The rotation matrices for each order (first index) are computed
c       from the lowest to the highest. As each one is generated, it
c       is applied to the input expansion "mpole" and overwritten.
c
c       As a result, it is sufficient to use two arrays rd1 and rd2 for
c       the two term recurrence, rather than storing them for all orders
c       as in the original code. There is some loss in speed
c       if the rotation operator is to be used multiple times, but the 
c       memory savings is often more critical.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       nterms: dimension parameter for d - the rotation matrix.
c       m1    : max m index for first expansion.
c       m2    : max m index for second expansion.
C       mpole   coefficients of original multiple expansion
C       rd1     work space 
C       rd2     work space
c       sqc:    an array contains the square roots of the
c               binomial coefficients.
c       theta:  the rotate angle about the y-axis.
c       ldc     dimensions of sqc array
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2,ldc
        integer ij,im,imp,m,mp
        real *8 theta
        real *8 ww,done,ctheta,stheta,hsthta,cthtap,cthtan,d
        real *8 precis,scale
        real *8 rd1(0:ldc,-ldc:ldc)
        real *8 rd2(0:ldc,-ldc:ldc)
        real *8 sqc(0:2*ldc,2)
        complex *16 mpole(0:ld1,-ld1:ld1)
        complex *16 marray(0:ld2,-ld2:ld2)
        data precis/1.0d-20/
c
        ww=1/sqrt(2.0d0)
        do m = 0, 2*ldc 
           sqc(m,1) = dsqrt(m+0.0d0)
        enddo
        sqc(0,2) = 0.0d0
        if( ldc .gt. 0 ) sqc(1,2) = 0.0d0
        do m = 2, 2*ldc 
           sqc(m,2) = dsqrt((m+0.0d0)*(m-1)/2.0d0)
        enddo
        do m = 0, ldc 
           do mp = -ldc, ldc
              rd1(m,mp) = 0d0
              rd2(m,mp) = 0d0
           enddo
        enddo
c
        done=1
        ctheta=dcos(theta)
        if (dabs(ctheta).le.precis) ctheta=0.0d0
        stheta=dsin(-theta)
        if (dabs(stheta).le.precis) stheta=0.0d0
        hsthta=ww*stheta
        cthtap=2.0d0*ww*dcos(theta/2.0d0)**2
        cthtan=-2.0d0*ww*dsin(theta/2.0d0)**2
c
c     Compute the (0,0,0) term.
c
      rd1(0,0)=done
      marray(0,0)=mpole(0,0)*rd1(0,0)
c
c ... Loop over first index ij=1,nterms, constructing
c     rotation matrices recursively.
c
      do ij=1,nterms
c
c     For mprime=0, use formula (1).
c
         do im=-min(ij,m2),-1
            rd2(0,im)=-sqc(ij-im,2)*rd1(0,im+1)
            if (im.gt.(1-ij)) then
               rd2(0,im)=rd2(0,im)+sqc(ij+im,2)*rd1(0,im-1)
            endif
            rd2(0,im)=rd2(0,im)*hsthta
            if (im.gt.-ij) then
               rd2(0,im)=rd2(0,im)+
     1         rd1(0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
            endif
            rd2(0,im)=rd2(0,im)/ij
         enddo
c
         rd2(0,0)=rd1(0,0)*ctheta
         if (ij.gt.1) then
            rd2(0,0)=rd2(0,0)+hsthta*sqc(ij,2)*(rd1(0,-1)+
     1                rd1(0,1))/ij
         endif
c
         do im=1,min(ij,m2)
            rd2(0,im)=-sqc(ij+im,2)*rd1(0,im-1)
            if (im.lt.(ij-1)) then
               rd2(0,im)=rd2(0,im)+sqc(ij-im,2)*rd1(0,im+1)
            endif
            rd2(0,im)=rd2(0,im)*hsthta
            if (im.lt.ij) then
               rd2(0,im)=rd2(0,im)+
     1         rd1(0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
            endif
            rd2(0,im)=rd2(0,im)/ij
         enddo
c
c        For 0<mprime<=j (2nd index) case, use formula (2).
c
         do imp=1,min(ij,m1)
c
            do im=-min(ij,m2),-1
               rd2(imp,im)=rd1(imp-1,im+1)*cthtan*sqc(ij-im,2)
               if (im.gt.(1-ij)) then
                  rd2(imp,im)=rd2(imp,im)-
     1            rd1(imp-1,im-1)*cthtap*sqc(ij+im,2)
               endif
               if (im.gt.-ij) then
                  rd2(imp,im)=rd2(imp,im)+
     1            rd1(imp-1,im)*stheta*sqc(ij+im,1)*sqc(ij-im,1)
               endif
               rd2(imp,im)=rd2(imp,im)*ww/sqc(ij+imp,2)
            enddo
c
            rd2(imp,0)=ij*stheta*rd1(imp-1,0)
            if (ij.gt.1) then
               rd2(imp,0)=rd2(imp,0)-sqc(ij,2)*(
     1         rd1(imp-1,-1)*cthtap+rd1(imp-1,1)*cthtan)
            endif
            rd2(imp,0)=rd2(imp,0)*ww/sqc(ij+imp,2)
c
            do im=1,min(ij,m2)
               rd2(imp,im)=rd1(imp-1,im-1)*cthtap*sqc(ij+im,2)
               if (im.lt.(ij-1)) then
                  rd2(imp,im)=rd2(imp,im)-
     1            rd1(imp-1,im+1)*cthtan*sqc(ij-im,2)
               endif
               if (im.lt.ij) then
                  rd2(imp,im)=rd2(imp,im)+
     1           rd1(imp-1,im)*stheta*sqc(ij+im,1)*sqc(ij-im,1)
               endif
               rd2(imp,im)=rd2(imp,im)*ww/sqc(ij+imp,2)
            enddo
c
         enddo
c
         do m=-ij,ij
            marray(ij,m)=mpole(ij,0)*rd2(0,m)
            do mp=1,ij
               marray(ij,m)=marray(ij,m)+
     1         mpole(ij,mp)*rd2(mp,m)+
     1         mpole(ij,-mp)*rd2(mp,-m)
            enddo
         enddo
c
         do mp=0,ij
            do m=-ij,ij
               rd1(mp,m) = rd2(mp,m)
            enddo
         enddo
      enddo
      return
      end
c
c
C*****************************************************************
        subroutine rotviarecur3f90(theta,nterms,m1,m2,mpole,ld1,
     1                         marray,ld2)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c       The rotation matrices for each order (first index) are computed
c       from the lowest to the highest. As each one is generated, it
c       is applied to the input expansion "mpole" and overwritten.
c
c       As a result, it is sufficient to use two arrays rd1 and rd2 for
c       the two term recurrence, rather than storing them for all orders
c       as in the original code. There is some loss in speed
c       if the rotation operator is to be used multiple times, but the 
c       memory savings is often more critical.
c
c       Note that this is just a driver to rotviarecur3s and
c       rotviarecur3t routines, the main purpose of this routine is to
c       have the calling sequence compatible with rotproj.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multiple expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c       lused    amount of workspace used.
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2
        integer ldc,ird1,lrd,ird2,isqc,lused,ier
        real *8 theta
        real *8, allocatable :: w(:)
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
c
        ldc = nterms
        ird1 = 1
        lrd  = (ldc+1)*(2*ldc+1) + 3
        ird2 = ird1 + lrd
        isqc = ird2 + lrd
        lused = isqc + 2*(2*ldc+1) 
        allocate (w(lused), stat=ier)
        if( ier.ne.0 ) return
c
        if( m1 .lt. nterms .or. m2 .lt. nterms ) then
        call rotviarecur3t(nterms,m1,m2,mpole,ld1,marray,
     1                    ld2,w(ird1),w(ird2),w(isqc),theta,ldc)
        else
        call rotviarecur3s(nterms,m1,m2,mpole,ld1,marray,
     1                    ld2,w(ird1),w(ird2),w(isqc),theta,ldc)
        endif
c
        return
        end
c        
c
c
c
C*****************************************************************
        subroutine rotviarecur3p_init(ier,rotmat,ldc,theta)
C*****************************************************************
c
c       Purpose:
c
c       Precompute the coefficients of rotation matrix. 
c
c       Use symmetry properties of rotation matrices. Fortran 90 code.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotate angle about the y-axis.
c       ldc     dimensions of rotation matrix array
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       ier: error return code
c       rotmat: coefficients of rotation matrix
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2,ldc
        integer ier,ij,im,imp,m,mp
        real *8 theta
        real *8 ww,done,ctheta,stheta,hsthta,cthtap,cthtan,d
        real *8 precis,scale
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
c
	real *8, allocatable :: rd1(:,:)
	real *8, allocatable :: rd2(:,:)
	real *8, allocatable :: sqc(:,:)
c
	data precis/1.0d-20/
c
	allocate (rd1(0:ldc,-ldc:ldc),stat=ier)
        if(ier .ne. 0 ) return
	allocate (rd2(0:ldc,-ldc:ldc),stat=ier)
        if(ier .ne. 0 ) return
	allocate (sqc(0:2*ldc,2),stat=ier)
        if(ier .ne. 0 ) return
c
        ww=1/sqrt(2.0d0)
        do m = 0, 2*ldc 
	   sqc(m,1) = dsqrt(m+0.0d0)
        enddo
	sqc(0,2) = 0.0d0
	if( ldc .gt. 0 ) sqc(1,2) = 0.0d0
        do m = 2, 2*ldc 
	   sqc(m,2) = dsqrt((m+0.0d0)*(m-1)/2.0d0)
        enddo
c
	done=1
	ctheta=dcos(theta)
	if (dabs(ctheta).le.precis) ctheta=0.0d0
	stheta=dsin(-theta)
	if (dabs(stheta).le.precis) stheta=0.0d0
	hsthta=ww*stheta
	cthtap=+2.0d0*ww*dcos(theta/2.0d0)**2
	cthtan=-2.0d0*ww*dsin(theta/2.0d0)**2
c
c     Compute the (0,0,0) term.
c
      rd1(0,0)=done
      rotmat(0,0,0)=rd1(0,0)
c
c ... Loop over first index ij=1,nterms, constructing
c     rotation matrices recursively.
c
      nterms=ldc
      do ij=1,nterms
c
c     For mprime=0, use formula (1).
c
         do im=-ij,-1
	    rd2(0,im)=-sqc(ij-im,2)*rd1(0,im+1)
	    if (im.gt.(1-ij)) then
	       rd2(0,im)=rd2(0,im)+sqc(ij+im,2)*rd1(0,im-1)
	    endif
	    rd2(0,im)=rd2(0,im)*hsthta
	    if (im.gt.-ij) then
	       rd2(0,im)=rd2(0,im)+
     1         rd1(0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
	    endif
	    rd2(0,im)=rd2(0,im)/ij
         enddo
c
	 rd2(0,0)=rd1(0,0)*ctheta
	 if (ij.gt.1) then
	    rd2(0,0)=rd2(0,0)+hsthta*sqc(ij,2)*(2*rd1(0,-1))/ij
	 endif
c         
	 do im=1,ij
	    rd2(0,im)=rd2(0,-im)
            if( mod(im,2) .eq. 0 ) then
            rd2(im,0)=+rd2(0,im)
            else
            rd2(im,0)=-rd2(0,im)
            endif
         enddo
c
c        For 0<mprime<=j (2nd index) case, use formula (2).
c
	 do imp=1,ij
c
            scale=(ww/sqc(ij+imp,2))
c
	    do im=imp,ij
c
	       rd2(imp,+im)=rd1(imp-1,+im-1)*(cthtap*sqc(ij+im,2))
	       rd2(imp,-im)=rd1(imp-1,-im+1)*(cthtan*sqc(ij+im,2))
c
	       if (im.lt.(ij-1)) then
	          rd2(imp,+im)=rd2(imp,+im)-rd1(imp-1,+im+1)*
     $               (cthtan*sqc(ij-im,2))
	          rd2(imp,-im)=rd2(imp,-im)-rd1(imp-1,-im-1)*
     $               (cthtap*sqc(ij-im,2))
	       endif
c
	       if (im.lt.ij) then
                  d=(stheta*sqc(ij+im,1)*sqc(ij-im,1))
	          rd2(imp,+im)=rd2(imp,+im)+rd1(imp-1,+im)*d
	          rd2(imp,-im)=rd2(imp,-im)+rd1(imp-1,-im)*d
	       endif
c
	       rd2(imp,+im)=rd2(imp,+im)*scale
	       rd2(imp,-im)=rd2(imp,-im)*scale
c
	       if (im.gt.imp) then
               if( mod(imp+im,2) .eq. 0 ) then
               rd2(im,+imp)=+rd2(imp,+im)
               rd2(im,-imp)=+rd2(imp,-im)
               else
               rd2(im,+imp)=-rd2(imp,+im)
               rd2(im,-imp)=-rd2(imp,-im)
               endif
               endif
c
            enddo
c
         enddo
c
         do m=-ij,ij
            do mp=0,ij
	       rd1(mp,m) = rd2(mp,m)
               rotmat(mp,ij,m) = rd2(mp,m)
            enddo
         enddo
c
      enddo
      return
      end
c
c
C*****************************************************************
        subroutine rotviarecur3p_apply
     $     (theta,nterms,m1,m2,mpole,ld1,marray,ld2,rotmat,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multiple expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       rotmat:  real *8 (0:ldc,0:ldc,-ldc:ldc): rotation matrix 
c       ldc   :  leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
	integer nterms,m1,m2,ld1,ld2,ldc
	integer ij,m,mp
	real *8 theta
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
c
        if (m1.lt.nterms .or. m2.lt.nterms) then
c
c       ... truncate multipole expansions
c
           do ij=0,nterms

              do m=-ij,ij
                  marray(ij,m)=0
              enddo
c
              do m=-min(ij,m2),min(ij,m2)
                 marray(ij,m)=mpole(ij,0)*rotmat(0,ij,m) 
                 do mp=1,min(ij,m1)
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*rotmat(mp,ij,m)+
     1              mpole(ij,-mp)*rotmat(mp,ij,-m)
                 enddo
              enddo
           enddo
        else
c
c       ... apply rotation matrix
c
           do ij=0,nterms
c
c              do m=-ij,ij
c                 marray(ij,m)=0
c              enddo
c
              do m=-ij,ij
                 marray(ij,m)=mpole(ij,0)*rotmat(0,ij,m) 
                 do mp=1,ij
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*rotmat(mp,ij,m)+
     1              mpole(ij,-mp)*rotmat(mp,ij,-m)
                 enddo
              enddo
           enddo
        endif
        return
        end
c        
c
c
C*****************************************************************
        subroutine rotviarecur3p_init_old(ier,rotmat,ldc,theta)
C*****************************************************************
c
c       Purpose:
c
c       Precompute the coefficients of rotation matrix. 
c
c       Use symmetry properties of rotation matrices. Fortran 90 code.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotate angle about the y-axis.
c       ldc     dimensions of rotation matrix array
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       ier: error return code
c       rotmat: coefficients of rotation matrix
c
C---------------------------------------------------------------------
	implicit none
        integer ld1,ld2,nterms,m1,m2,ldc
        integer ier,ij,im,imp,m,mp
        real *8 theta
        real *8 ww,done,ctheta,stheta,hsthta,cthtap,cthtan,d
        real *8 precis,scale
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
c
	real *8, allocatable :: rd1(:,:)
	real *8, allocatable :: rd2(:,:)
	real *8, allocatable :: sqc(:,:)
c
	data precis/1.0d-20/
c
	allocate (rd1(0:ldc,-ldc:ldc),stat=ier)
        if(ier .ne. 0 ) return
	allocate (rd2(0:ldc,-ldc:ldc),stat=ier)
        if(ier .ne. 0 ) return
	allocate (sqc(0:2*ldc,2),stat=ier)
        if(ier .ne. 0 ) return
c
        ww=1/sqrt(2.0d0)
        do m = 0, 2*ldc 
	   sqc(m,1) = dsqrt(m+0.0d0)
        enddo
	sqc(0,2) = 0.0d0
	if( ldc .gt. 0 ) sqc(1,2) = 0.0d0
        do m = 2, 2*ldc 
	   sqc(m,2) = dsqrt((m+0.0d0)*(m-1)/2.0d0)
        enddo
c
	done=1
	ctheta=dcos(theta)
	if (dabs(ctheta).le.precis) ctheta=0.0d0
	stheta=dsin(-theta)
	if (dabs(stheta).le.precis) stheta=0.0d0
	hsthta=ww*stheta
	cthtap=+2.0d0*ww*dcos(theta/2.0d0)**2
	cthtan=-2.0d0*ww*dsin(theta/2.0d0)**2
c
c     Compute the (0,0,0) term.
c
      rd1(0,0)=done
      rotmat(0,0,0)=rd1(0,0)
c
c ... Loop over first index ij=1,nterms, constructing
c     rotation matrices recursively.
c
      nterms=ldc
      do ij=1,nterms
c
c     For mprime=0, use formula (1).
c
         do im=-ij,-1
	    rd2(0,im)=-sqc(ij-im,2)*rd1(0,im+1)
	    if (im.gt.(1-ij)) then
	       rd2(0,im)=rd2(0,im)+sqc(ij+im,2)*rd1(0,im-1)
	    endif
	    rd2(0,im)=rd2(0,im)*hsthta
	    if (im.gt.-ij) then
	       rd2(0,im)=rd2(0,im)+
     1         rd1(0,im)*ctheta*sqc(ij+im,1)*sqc(ij-im,1)
	    endif
	    rd2(0,im)=rd2(0,im)/ij
         enddo
c
	 rd2(0,0)=rd1(0,0)*ctheta
	 if (ij.gt.1) then
	    rd2(0,0)=rd2(0,0)+hsthta*sqc(ij,2)*(2*rd1(0,-1))/ij
	 endif
c         
	 do im=1,ij
	    rd2(0,im)=rd2(0,-im)
            if( mod(im,2) .eq. 0 ) then
            rd2(im,0)=+rd2(0,im)
            else
            rd2(im,0)=-rd2(0,im)
            endif
         enddo
c
c        For 0<mprime<=j (2nd index) case, use formula (2).
c
	 do imp=1,ij
c
            scale=(ww/sqc(ij+imp,2))
c
	    do im=imp,ij
c
	       rd2(imp,+im)=rd1(imp-1,+im-1)*(cthtap*sqc(ij+im,2))
	       rd2(imp,-im)=rd1(imp-1,-im+1)*(cthtan*sqc(ij+im,2))
c
	       if (im.lt.(ij-1)) then
	          rd2(imp,+im)=rd2(imp,+im)-rd1(imp-1,+im+1)*
     $               (cthtan*sqc(ij-im,2))
	          rd2(imp,-im)=rd2(imp,-im)-rd1(imp-1,-im-1)*
     $               (cthtap*sqc(ij-im,2))
	       endif
c
	       if (im.lt.ij) then
                  d=(stheta*sqc(ij+im,1)*sqc(ij-im,1))
	          rd2(imp,+im)=rd2(imp,+im)+rd1(imp-1,+im)*d
	          rd2(imp,-im)=rd2(imp,-im)+rd1(imp-1,-im)*d
	       endif
c
	       rd2(imp,+im)=rd2(imp,+im)*scale
	       rd2(imp,-im)=rd2(imp,-im)*scale
c
	       if (im.gt.imp) then
               if( mod(imp+im,2) .eq. 0 ) then
               rd2(im,+imp)=+rd2(imp,+im)
               rd2(im,-imp)=+rd2(imp,-im)
               else
               rd2(im,+imp)=-rd2(imp,+im)
               rd2(im,-imp)=-rd2(imp,-im)
               endif
               endif
c
            enddo
c
         enddo
c
         do m=-ij,ij
            do mp=0,ij
	       rd1(mp,m) = rd2(mp,m)
               rotmat(ij,mp,m) = rd2(mp,m)
            enddo
         enddo
c
      enddo
      return
      end
c
c
C*****************************************************************
        subroutine rotviarecur3p_apply_old
     $     (theta,nterms,m1,m2,mpole,ld1,marray,ld2,rotmat,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       theta:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
c
c       m1    :  max m index for first expansion  
c       m2    :  max m index for second expansion 
C       mpole :  coefficients of original multiple expansion
C       ld1   :  leading dim for mpole (must exceed nterms)
C       ld2   :  leading dim for marray (must exceed nterms)
c       rotmat:  real *8 (0:ldc,0:ldc,-ldc:ldc): rotation matrix 
c       ldc   :  leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
	integer nterms,m1,m2,ld1,ld2,ldc
	integer ij,m,mp
	real *8 theta
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
c
        if (m1.lt.nterms .or. m2.lt.nterms) then
c
c       ... truncate multipole expansions
c
           do ij=0,nterms

              do m=-ij,ij
                  marray(ij,m)=0
              enddo
c
              do m=-min(ij,m2),min(ij,m2)
                 marray(ij,m)=mpole(ij,0)*rotmat(ij,0,m) 
                 do mp=1,min(ij,m1)
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*rotmat(ij,mp,m)+
     1              mpole(ij,-mp)*rotmat(ij,mp,-m)
                 enddo
              enddo
           enddo
        else
c
c       ... apply rotation matrix
c
           do ij=0,nterms
c
c              do m=-ij,ij
c                 marray(ij,m)=0
c              enddo
c
              do m=-ij,ij
                 marray(ij,m)=mpole(ij,0)*rotmat(ij,0,m) 
                 do mp=1,ij
                    marray(ij,m)=marray(ij,m)+
     1	            mpole(ij,mp)*rotmat(ij,mp,m)+
     1              mpole(ij,-mp)*rotmat(ij,mp,-m)
                 enddo
              enddo
           enddo
        endif
        return
        end
c        
c
c
