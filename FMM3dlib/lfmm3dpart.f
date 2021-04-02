cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-07-05 09:56:35 -0400 (Mon, 05 Jul 2010) $
c    $Revision: 1048 $
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Laplace potentials and fields due to
c     point charges and dipoles.  (FORTRAN 90 VERSION)
c
c     lfmm3dpart - Laplace FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm3dpartself - Laplace FMM in R^3: evaluate all pairwise particle
c         interactions (ignoring self-interaction)
c
c     lfmm3dparttarg - Laplace FMM in R^3: evaluate all pairwise
c         particle interactions (ignoring self-interaction) +
c         interactions with targets
c
c     l3dpartdirect - Laplace interactions in R^3: evaluate all
c         pairwise particle interactions (ignoring self-interaction) +
c       interactions with targets via direct O(N^2) algorithm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        subroutine lfmm3dpart(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c              
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use (1/r) for the Green's function, without the 
c       (1/4 pi) scaling. Self-interactions are not included.
c   
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
c
        implicit real *8 (a-h,o-z)
        dimension source(3,1)
        complex *16 charge(1)
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
c
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        iffldtarg=0
c
        call lfmm3dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
        return
        end
c
c
c
c
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        subroutine lfmm3dpartself(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c              
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction). 
c       We use (1/r) for the Green's function, without the 
c       (1/4 pi) scaling. Self-interactions are not included.
c   
c       The main FMM routine permits both evaluation at sources
c       and at a collection of targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
c
        implicit real *8 (a-h,o-z)
        dimension source(3,1)
        complex *16 charge(1)
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
c
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
        ifpottarg=0
        iffldtarg=0
c
        call lfmm3dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
        return
        end
c
c
c
c
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        subroutine lfmm3dparttarg(ier,iprec,nsource,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c       
c       Laplace FMM in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dparttargmain.
c
c       INPUT PARAMETERS:
c
c       iprec:     FMM precision flag
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
c       nsource          : number of sources                           (integer)
c       source(3,nsource): source locations 			       (real *8)  
c       ifcharge         : charge computation flag                     (integer)
c                          ifcharge = 1 => include charge contribution
c                                     otherwise do not
c       charge(nsource)  : charge strengths                        (complex *16)
c       ifdipole         : dipole computation flag                     (integer)
c                          ifdipole = 1 =>  include dipole contribution
c                                     otherwise do not
c       dipstr(nsource)  : dipole strengths                        (complex *16) 
c       dipvec(3,nsource): dipole orientation vectors                  (real *8) 
c
c       ifpot            : potential flag                              (integer)
c                          (1=compute potential, otherwise no)
c       iffld            : field flag                                  (integer) 
c                          (1=compute field, otherwise no)
c       ntarget          : number of targets                           (integer)  
c       target(3,ntarget): target locations                            (real *8) 
c       ifpottarg        : target potential flag                       (integer)
c                          (1=compute potential, otherwise no)
c       iffldtarg        : target field flag                           (integer)
c                          (1=compute field, otherwise no)
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
c       pot(nsource)       : potential at source locations         (complex *16)
c       fld(3,nsource)     : field (-gradient) at source locations (complex *16)
c       pottarg(ntarget)   : potential at target locations         (complex *16)
c       fldtarg(3,ntarget) : field (-gradient) at target locations (complex *16) 
c-----------------------------------------------------------------------
c
cf2py   intent(out) ier
cf2py   intent(in) iprec
cf2py   intent(in) nsource, source
cf2py   intent(in) ifcharge,charge
cf2py   check(!ifcharge || (shape(charge,0) == nsource))  charge
cf2py   depend(nsource)  charge
cf2py   intent(in) ifdipole,dipvec,dipstr
cf2py   check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
cf2py   depend(nsource)  dipstr
cf2py   intent(in) ifpot,iffld
cf2py   intent(out) pot,fld
cf2py   intent(in) ifpottarg, iffldtarg
cf2py   intent(in) target
cf2py   intent(in) ntarget
cf2py   check((!ifpottarg && !iffldtarg) || (shape(target,0)==3 && shape(target,1) == ntarget))  target
cf2py   check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
cf2py   check((!iffldtarg) || (shape(fldtarg,0)==3 && shape(fldtarg,1) == ntarget))  fldtarg
c
c       (F2PY workaround: pottarg, fldtarg must be input because f2py
c       refuses to allocate zero-size output arrays.)
c
cf2py   intent(in,out) pottarg,fldtarg
c
        implicit real *8 (a-h,o-z)
        dimension source(3,nsource)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        dimension dipvec(3,nsource)
        complex *16 ima
        complex *16 pot(nsource)
        complex *16 fld(3,nsource)
        dimension target(3,nsource)
        complex *16 pottarg(ntarget)
        complex *16 fldtarg(3,ntarget)
c
        dimension timeinfo(10)
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 evels of refinment, which is 
c     more than enough for any non-pathological case.
c
 
        dimension laddr(2,200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
        real *8, allocatable :: wrmlexp(:)
        complex *16 ptemp,ftemp(3)
c       
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
        lused7 = 0
c       
        done=1
        pi=4*atan(done)
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
        ifprint=1
c
c     set fmm tolerance based on iprec flag.
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
        if (ifprint .ge. 1) call prin2('epsfmm=*',epsfmm,1)
c
c
c     set criterion for box subdivision (number of sources per box)
c
        if( iprec .eq. -2 ) nbox=40*1.0
        if( iprec .eq. -1 ) nbox=50*1.0
        if( iprec .eq. 0 ) nbox=80*1.0
        if( iprec .eq. 1 ) nbox=160*1.0
        if( iprec .eq. 2 ) nbox=400*1.0
        if( iprec .eq. 3 ) nbox=800*1.0
        if( iprec .eq. 4 ) nbox=1200*1.0
        if( iprec .eq. 5 ) nbox=1400*1.0
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        if (ifprint .ge. 1) call prinf('nbox=*',nbox,1)
c
c
c     create oct-tree data structure
c
        t1=second()
C$        t1=omp_get_wtime()
        ntot = 100*(nsource+ntarget)+10000
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
        t2=second()
C$        t2=omp_get_wtime()
        if( ifprint .eq. 1 ) call prin2('time in d3tstrcr=*',t2-t1,1)
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
c       
        if (ifprint .ge. 1) call prin2('scale=*',scale,nlev+1)
c       
c
c       carve up workspace further
c
c     isourcesort is pointer for sorted source coordinates
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
        isourcesort = lused7 + 5
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 3*nsource
          ldipstr = 2*nsource
        else
          ldipvec = 3
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr
c
c       ... allocate the potential and field arrays
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
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c      
        nmax = 0
        call l3dterms(epsfmm, nterms_lap, ier)
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           nterms(i)=nterms_lap
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        if (ifprint .ge. 1) call prinf('nterms=*',nterms,nlev+1)
        if (ifprint .ge. 1) call prinf('nmax=*',nmax,1)
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
     1                  lused7 is *',lused7,1)
           ier = 8
           return
        endif
c
c     reorder sources, targets so that each box holds
c     contiguous list of source/target numbers.
c
        call l3dreorder(nsource,source,ifcharge,charge,wlists(iisource),
     $     ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 
c       
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     1       w(itargetsort))
c
        if (ifprint .ge. 1) call prinf('finished reordering=*',ier,1)
        if (ifprint .ge. 1) call prinf('ier=*',ier,1)
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('nlev=*',nlev,1)
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('lused7=*',lused7,1)
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call l3dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if (ifprint .ge. 1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                  lused7 is *',lused7,1)
           ier = 16
           return
        endif
c
c       
ccc        do i=lused7+1,lused7+1+100
ccc        w(i)=777
ccc        enddo
c
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far field and local fields (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
        ifevalfar=1
        ifevalloc=1
c
        t1=second()
C$        t1=omp_get_wtime()
        call lfmm3dparttargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,w(isourcesort),wlists(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),iffld,w(ifld),
     $     ntarget,w(itargetsort),wlists(iitarget),
     $     ifpottarg,w(ipottarg),iffldtarg,w(ifldtarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
        t2=second()
C$        t2=omp_get_wtime()
        if( ifprint .eq. 1 ) call prin2('time in fmm main=*',t2-t1,1)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if (ifprint .ge. 1) call prinf('lwlists=*',lused,1)
        if (ifprint .ge. 1) call prinf('lused total =*',lused7,1)
c       
        if (ifprint .ge. 1) 
     $      call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
        if(iffld .eq. 1) 
     $     call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
c
        if(ifpottarg .eq. 1 )
     $     call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
        if(iffldtarg .eq. 1) 
     $     call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dparttargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension sourcesort(3,1), isource(1)
        complex *16 chargesort(1)
        complex *16 dipstrsort(1)
        dimension dipvecsort(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        dimension targetsort(3,1), itarget(1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        dimension wlists(1)
        dimension iaddr(2,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(40 000)
        dimension nterms_eval(4,0:200)
c
        data ima/(0.0d0,1.0d0)/

c
c
c     INPUT PARAMETERS:
c
c     iprec        precision flag (see above)   ELIMINATE???
c     ifevalfar    far field flag (1 means compute far field, 
c                                  else dont)
c     ifevalloc    local field flag (1 means compute local field, 
c                                    else dont)
c     nsource      number of sources
c     sourcesort   sorted source coordinates
c     isource      sorting index for sources
c     ifcharge     flag indicating potential includes contribution
c                  from charges
c     chargesort   sorted charge values
c     ifdipole     flag indicating potential includes contribution
c                  from dipoles
c     dipstrsort   sorted dipole strengths
c     dipvecsort   sorted dipole orientation vectors
c     ifpot        potential flag (1 => compute, else do not)
c     iffld        field flag (1 => compute, else do not)
c     ntarget      number of targets
c     targetsort   sorted array of target locations
c     itarget      sorting index for targets
c     ifpottarg    target potential flag (1 => compute, else do not)
c     iffldtarg    target field flag (1 => compute, else do not)
c     epsfmm       FMM tolerance
c     iaddr        iaddr(2,nboxes) array points to mpole/local
c                     expansions for each box
c     rmlexp       workspace to contain mpole/local expansions.
c     nboxes       number of boxes in FMM hierarchy
c     laddr        indexing array for FMM data structure
c     nlev         number of levels in FMM hierarchy
c     scale        array of scaling parameters
c     bsize        box dimension for FMM
c     nterms       array of nterms needed at each level
c     wlists       FMM data structure (real array)
c     lw           length of wlists
c
c
c     OUTPUT PARAMETERS:
c
c     pot          surface potential (if ifpot=1)
c     fld          surface field=-gradient(potential) (if iffld=1)
c     pottarg      target potential (if ifpot=1)
c     fldtarg      target field=-gradient(potential) (if iffld=1)
c     ier          error return code
c                  ier = 0    =>   normal execution
c                  ier = 4    =>   cannot allocate tree workspace
c                  ier = 8    =>   cannot alocate bulk FMM workspace
c                  ier = 16   =>   cannot allocate mpole exp workspace
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
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        enddo
c
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
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
        if (ifprint .ge. 2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,box,center0,corner0,level)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(8) 
        do ibox = 1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        enddo
C$OMP END PARALLEL DO
c
c
        if (ifprint .ge. 1) call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
ccc        do 1300 ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(8) 
        do 1200 ibox=1,nboxes
ccc        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
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
            call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))

            if( ifcharge .eq. 1 ) then
c
            call l3dformmp_add_trunc(ier,scale(level),
     1         sourcesort(1,box(14)),chargesort(box(14)),
     $         npts,center0,
     $         nterms(level),nterms_eval(1,level),
     2         rmlexp(iaddr(1,ibox)),wlege,nlege)        
c
            endif
c 
            if (ifdipole .eq. 1 ) then

            call l3dformmp_dp_add_trunc(ier,scale(level),
     $         sourcesort(1,box(14)),
     1         dipstrsort(box(14)),dipvecsort(1,box(14)),
     $         npts,center0,nterms(level),nterms_eval(1,level),
     2         rmlexp(iaddr(1,ibox)),wlege,nlege)
            
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
        if (ifprint .ge. 1) call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c-----------------------------------------------------------------------
c       Step 2: In the adaptive FMM, a large leaf node may need to interact
c       with separated boxes at finer levels. This is called <list 3> in the
c       FMM. One takes individual sources in the large leaf node and 
c       maps them to loca expansions in the target boxes.
c-----------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist,npts) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(8) 
         do 3251 ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
ccc         if( box(15) .eq. 0 ) goto 3251
c
         itype=4
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
ccc         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
ccc         call prinf('nlist3=*', nlist,1)
         do 3250 ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            npts=box1(15)            
            if( npts .eq. 0 ) goto 3250
c
            level0=box(1)
            level1=box1(1)
c
c            ifdirect3 = 0
c            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
c     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
c
               if( ifcharge .eq. 1 ) then
c
               call l3dformta_add_trunc(ier,scale(level0),
     1            sourcesort(1,box1(14)),chargesort(box1(14)),
     $            npts,center0,
     $            nterms(level0),nterms_eval(1,level0),
     2            rmlexp(iaddr(2,ibox)),wlege,nlege)
c
               endif
c
               if( ifdipole .eq. 1 ) then

               call l3dformta_dp_add_trunc(ier,scale(level0),
     1            sourcesort(1,box1(14)),dipstrsort(box1(14)),
     2            dipvecsort(1,box1(14)),npts,center0,
     3            nterms(level0),nterms_eval(1,level0),
     $            rmlexp(iaddr(2,ibox)),wlege,nlege)

               endif
c
            else

            call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

            endif
 3250    continue
c
 3251    continue
C$OMP END PARALLEL DO
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c-----------------------------------------------------------------------
c       Steps 3,4,5 are carried out by the routine lfmm3d_list2.
c       Step 3 is the merging of multipole expansions at every level.
c       Step 4 is the mapping of multipole expansions to local expansions
c              using <list 2>.
c       Step 5 is the recursive mapping of local expansions from parent to 
c              child.
c-----------------------------------------------------------------------
c
        if (ifprint .ge. 1) call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0
        call lfmm3d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,ifprune_list2)
c
c
c
        if (ifprint .ge. 1) call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c-----------------------------------------------------------------------
c       Step 6: In the adaptive FMM, a small leaf node may need to interact
c       with separated boxes at coarser levels. This is the dual of 
c       Step 2 and is called <list 4> in the FMM. 
c       The multipole expansion for the small leaf node is evaluated directly
c       at targets in the large target node.
c-----------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,itype,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist,level) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(8) 
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
ccc         if( box(15) .eq. 0 ) goto 3252
c
         itype=3
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
ccc         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
ccc         call prinf('nlist4=*', nlist,1)
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box1(1)
c
c            ifdirect4 = 0
c
c            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
c     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
c           for future optimization - here, we just evaluate the 
c           multipole expansion, regardless of the number of sources
c           in the source box.
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then

            if( box(15) .gt. 0 ) 
     $         call l3dmpevalall_trunc(scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),
     $         nterms(level),nterms_eval(1,level),
     $         sourcesort(1,box(14)),box(15),
     $         ifpot,pot(box(14)),
     $         iffld,fld(1,box(14)),
     $         wlege,nlege,ier)

            if( box(17) .gt. 0 ) 
     $         call l3dmpevalall_trunc(scale(level),center1,
     $         rmlexp(iaddr(1,jbox)),
     $         nterms(level),nterms_eval(1,level),
     $         targetsort(1,box(16)),box(17),
     $         ifpottarg,pottarg(box(16)),
     $         iffldtarg,fldtarg(1,box(16)),
     $         wlege,nlege,ier)

            else
            
            call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

            endif
        enddo
 3252   continue
C$OMP END PARALLEL DO
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c
c-----------------------------------------------------------------------
c       Step 7: Evaluate the local expansions for all relevant sources/targets.
c-----------------------------------------------------------------------
c
        if (ifprint .ge. 1) call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(8) 
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

        call l3dtaevalall_trunc(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),
     $     nterms(level),nterms_eval(1,level),
     $     sourcesort(1,box(14)),box(15),
     $     ifpot,pot(box(14)),
     $     iffld,fld(1,box(14)),
     $     wlege,nlege,ier)

        call l3dtaevalall_trunc(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),
     $     nterms(level),nterms_eval(1,level),
     $     targetsort(1,box(16)),box(17),
     $     ifpottarg,pottarg(box(16)),
     $     iffldtarg,fldtarg(1,box(16)),
     $     wlege,nlege,ier)
        
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
        if (ifprint .ge. 1) call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c-----------------------------------------------------------------------
c       Step 8: Evaluate direct interactions locally
c-----------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(8) 
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
        if( box(15) .gt. 0 ) 
     $     call lfmm3dpart_direct_self(box,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,
     $     targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
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
            call lfmm3dpart_direct_targ(box1,box,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
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
        if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
c       
        d=0
        do i=1,8
        d=d+timeinfo(i)
        enddo
c       
        if (ifprint .ge. 1) call prin2('sum(timeinfo)=*',d,1)
c
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('nsource=*',nsource,1)
        if (ifprint .ge. 1) call prinf('ntarget=*',ntarget,1)
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dpart_direct_self(box,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
        do 6160 j=box(14),box(14)+box(15)-1
        do 6150 i=box(14),box(14)+box(15)-1
            if (i .eq. j) goto 6150
            if (ifcharge .eq. 1 ) then
            call lpotfld3d(iffld,source(1,i),charge(i),
     1           source(1,j),ptemp,ftemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (iffld .eq. 1) then
               fld(1,j)=fld(1,j)+ftemp(1)
               fld(2,j)=fld(2,j)+ftemp(2)
               fld(3,j)=fld(3,j)+ftemp(3)
            endif
            endif
            if (ifdipole .eq. 1) then
               call lpotfld3d_dp(iffld,source(1,i),
     $              dipstr(i),dipvec(1,i),
     $              source(1,j),ptemp,ftemp)
               if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
               if (iffld .eq. 1) then
                  fld(1,j)=fld(1,j)+ftemp(1)
                  fld(2,j)=fld(2,j)+ftemp(2)
                  fld(3,j)=fld(3,j)+ftemp(3)
               endif
            endif
 6150   continue
 6160   continue
        endif

        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
        do j=box(16),box(16)+box(17)-1
        if (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
        call lpotfld3dall_targ
     $     (iffldtarg,source(1,box(14)),charge(box(14)),
     1     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        if (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
        call lpotfld3dall_dp_targ(iffldtarg,source(1,box(14)),
     $     dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
        call lpotfld3dall_sdp_targ(iffldtarg,source(1,box(14)),
     $     charge(box(14)),dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        enddo
        endif
c
        return
c       
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
        do j=box(16),box(16)+box(17)-1
        if (ifcharge .eq. 1 ) then
        call lpotfld3dall_targ
     $     (iffldtarg,source(1,box(14)),charge(box(14)),
     1     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        if (ifdipole .eq. 1) then
        call lpotfld3dall_dp_targ(iffldtarg,source(1,box(14)),
     $     dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        enddo
        endif
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dpart_direct_targ(box,box1,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
        do j=box1(14),box1(14)+box1(15)-1
        if (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
        call lpotfld3dall_targ
     $     (iffld,source(1,box(14)),charge(box(14)),
     1     box(15),source(1,j),ptemp,ftemp)
        if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
        if (iffld .eq. 1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
        if (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
        call lpotfld3dall_dp_targ(iffld,source(1,box(14)),
     $     dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),source(1,j),ptemp,ftemp)
        if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
        if (iffld .eq. 1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
        if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
        call lpotfld3dall_sdp_targ(iffld,source(1,box(14)),
     $     charge(box(14)),dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),source(1,j),ptemp,ftemp)
        if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
        if (iffld .eq. 1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
        enddo
        endif
c
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
        do j=box1(16),box1(16)+box1(17)-1
        if (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
        call lpotfld3dall_targ
     $     (iffldtarg,source(1,box(14)),charge(box(14)),
     1     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        if (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
        call lpotfld3dall_dp_targ(iffldtarg,source(1,box(14)),
     $     dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
        call lpotfld3dall_sdp_targ(iffldtarg,source(1,box(14)),
     $     charge(box(14)),dipstr(box(14)),dipvec(1,box(14)),
     $     box(15),target(1,j),ptemp,ftemp)
        if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg .eq. 1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
        enddo
        endif
c       
        return
        end
c
c
c
c
c
        subroutine l3dpartdirect(nsource,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarget,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
c       Laplace interactions in R^3: evaluate all pairwise particle
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not-included.
c   
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of sources
c       source: real *8 (3,nsource):  source locations
c       ifcharge:  charge computation flag
c                  ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): charge strengths
c       ifdipole:  dipole computation flag
c                  ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): dipole strengths
c       dipvec: real *8 (3,nsource): dipole orientation vectors. 
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at source locations
c       fld: complex *16 (3,nsource): field (-gradient) at source locations
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
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
        enddo
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp) 
        do 6160 j=1,nsource
        do 6150 i=1,nsource
            if (i .eq. j) goto 6150
            if (ifcharge .eq. 1 ) then
            call lpotfld3d(iffld,source(1,i),charge(i),
     1           source(1,j),ptemp,ftemp)
            if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
            if (iffld .eq. 1) then
               fld(1,j)=fld(1,j)+ftemp(1)
               fld(2,j)=fld(2,j)+ftemp(2)
               fld(3,j)=fld(3,j)+ftemp(3)
            endif
            endif
            if (ifdipole .eq. 1) then
               call lpotfld3d_dp(iffld,source(1,i),
     $              dipstr(i),dipvec(1,i),
     $              source(1,j),ptemp,ftemp)
               if (ifpot .eq. 1) pot(j)=pot(j)+ptemp
               if (iffld .eq. 1) then
                  fld(1,j)=fld(1,j)+ftemp(1)
                  fld(2,j)=fld(2,j)+ftemp(2)
                  fld(3,j)=fld(3,j)+ftemp(3)
               endif
            endif
 6150   continue
 6160   continue
C$OMP END PARALLEL DO
        endif

        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp) 
        do j=1,ntarget
        do i=1,nsource
            if (ifcharge .eq. 1 ) then
            call lpotfld3d(iffldtarg,source(1,i),charge(i),
     1           target(1,j),ptemp,ftemp)
            if (ifpottarg .eq. 1) pottarg(j)=pottarg(j)+ptemp
            if (iffldtarg .eq. 1) then
               fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
               fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
               fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
            endif
            endif
            if (ifdipole .eq. 1) then
               call lpotfld3d_dp(iffldtarg,source(1,i),
     $              dipstr(i),dipvec(1,i),
     $              target(1,j),ptemp,ftemp)
               if (ifpottarg .eq. 1 ) pottarg(j)=pottarg(j)+ptemp
               if (iffldtarg .eq. 1) then
                  fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
                  fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
                  fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
               endif
            endif
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
