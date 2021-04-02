
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Vector spherical harmonic transforms
c
c       Fortran F90 version, allocatable arrays
c
c       NOTE #1: heavy use of allocatable arrays to simplify programming 
c       We use automatic deallocation feature of Fortran F95
c
c       vshgrid - form the nodes and weights for quadrature on the unit sphere
c                 legendre nodes along theta, and equispaced point along phi
c
c       vshgfun - precompute the vsh functions needed in vsh transforms
c                 memory intensive O(p^3) storage, suitable for repeated
c                 calls to vsh transform.
c
c       vshcoef - compute the coefficients of vector spherical harmonics
c                 given a vector field on the unit sphere grid
c                 needs to call vshgfun first to obtain vsh function values.
c       
c       vshevalg - evaluate the vector field on the unit sphere grid given 
c                 the coefficients of its vector spherical harmonic expansion
c                 needs to call vshgfun first to obtain vsh function values.
c
c       vshcoef2 - compute the coefficients of vector spherical harmonics
c                 given a vector field on the unit sphere grid
c                 compute vsh function on the fly.
c                 adapted from em3ehformex_fast
c       
c       vshevalg2 - evaluate the vector field on the unit sphere grid given 
c                 the coefficients of its vector spherical harmonic expansion
c                 compute vsh function on the fly
c                 adapted from em3exevalehz_fast
c
c       vsheval - evaluate the vector field at a given point in R^3 
c                 given the coefficients of its vector spherical 
c                 harmonic expansion
c
c
c
        subroutine vshgrid(nphi,ntheta,rnodes,weights,nnodes)
c       This subroutine returns the xyz coordinates of the grid 
c       on the unit sphere,
c       the grid is used in (vector) spherical harmonic transform.
c
c       Input: nphi, ntheta 
c       Output: rnodes, weights, nnodes
c
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nphi,ntheta),weights(nphi,ntheta)
ccc        dimension ts(10 000), ws(10 000)
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )

        done=1
        pi=4*atan(done)
c
c       ... construct the Gaussian nodes and weights on the interval [-1,1]
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        do 1400 i=1,ntheta
c
        z=ts(i)
        r=sqrt(1-z*z)
c
c       ... construct FFT compatible angular quadratures
c
        do 1200 j=1,nphi
c
        phi=(j-1)*(2*pi/nphi)
        x=r*cos(phi)
        y=r*sin(phi)
        rnodes(1,j,i)=x
        rnodes(2,j,i)=y
        rnodes(3,j,i)=z
        weights(j,i)=ws(i)*(2*pi/nphi)
 1200   continue
 1400   continue
c
        nnodes=nphi*ntheta
c       
        return
        end
c
c
c
        subroutine vshgfun
     $     (nterms,ntheta,ynm,xnm2,unm2)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the vsh function values needed
c       in the vsh transforms.
c
c       O(p^3) work and storage routine
c
c          Input parameters:
c
c       nterms - the number of terms in vsh expansion
c       ntheta - the number of grid points in theta direction
c
c          Output parameters:
c
c       ynm (complex*16)(0:nterms,-nterms:nterms,ntheta)
c       xnm2 (complex*16)(2,0:nterms,-nterms:nterms,ntheta) 
c       unm2 (complex*16)(2,0:nterms,-nterms:nterms,ntheta)
c           - the values of vector spherical harmonic functions
c        
        complex *16 ynm(0:nterms,-nterms:nterms,ntheta)
        complex *16 xnm2(2,0:nterms,-nterms:nterms,ntheta)
        complex *16 unm2(2,0:nterms,-nterms:nterms,ntheta)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)

        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c       
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )

        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 

        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        do 1400 itheta=1,ntheta
c
        costheta=ts(itheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c

        do n=0,nterms
        do m=0,n
        ynm(n,m,itheta)=pnm(n,m)
        xnm2(1,n,m,itheta)=pxnm2(1,n,m)
        xnm2(2,n,m,itheta)=pxnm2(2,n,m)
        unm2(1,n,m,itheta)=punm2(1,n,m)
        unm2(2,n,m,itheta)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m,itheta)=conjg(ynm(n,-m,itheta))
        xnm2(1,n,m,itheta)=conjg(pxnm2(1,n,-m))
        xnm2(2,n,m,itheta)=conjg(pxnm2(2,n,-m))
        unm2(1,n,m,itheta)=conjg(punm2(1,n,-m))
        unm2(2,n,m,itheta)=conjg(punm2(2,n,-m))
        enddo
        enddo
c
 1400   continue
c
        return
        end
c
c
c
c
        subroutine vshcoef
     $     (fgrid,nterms,nphi,ntheta,ynm,xnm2,unm2,
     $      ycoef,xcoef,ucoef)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the coefficients of vector spherical
c       harmonic expansion given the values of the vector field on the grid
c
c       Fast O(p^3) routine
c
c       Fast algorithm, perform FFT along all parallels
c
c          Input parameters:
c
c       fgrid (complex *16)(3,nphi,ntheta) - H field grid
c       nterms - the number of terms in vsh expansion
c       nphi, ntheta - the number of grid points in each direction
c       ynm,xnm2,unm2 - precomputed vsh function values
c
c          Output parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c       xcoef (complex*16)(0:nterms,-nterms:nterms) 
c       ucoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of vector spherical harmonic expansion
c        
        complex *16 ycoef(0:nterms,-nterms:nterms)
        complex *16 xcoef(0:nterms,-nterms:nterms)
        complex *16 ucoef(0:nterms,-nterms:nterms)
c
        complex *16 fgrid(3,nphi,ntheta)
c
        complex *16 ynm(0:nterms,-nterms:nterms,ntheta)
        complex *16 xnm2(2,0:nterms,-nterms:nterms,ntheta)
        complex *16 unm2(2,0:nterms,-nterms:nterms,ntheta)

        complex *16, allocatable :: fcgrid(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c       
        allocate( fcgrid(3,nphi,ntheta) )

        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do n=0,nterms
        do m=-nterms,nterms
           ycoef(n,m)=0
           xcoef(n,m)=0
           ucoef(n,m)=0
        enddo
        enddo
c
c       ... first, convert vector field grids to spherical coordinates
c       and perform FFT along each parallel
c
        done=1
        pi4=16*atan(done)

        do i=1,ntheta
           do j=1,nphi
              do k=1,3
                 fcgrid(k,j,i)=fgrid(k,j,i)/pi4
              enddo
           enddo
        enddo
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        call vshsphfft(fcgrid,nphi,ntheta,ts,ws)
c
c
        do 1400 itheta=1,ntheta
c
c       ... compute coefficients via inner product
c
c

        do n=0,nterms
        do m=-n,n
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ycoef(n,m)=ycoef(n,m)+ws(itheta)*
     $     fcgrid(3,mf,itheta)*conjg(ynm(n,m,itheta))

        xcoef(n,m)=xcoef(n,m)+ws(itheta)*
     $     (
     $     fcgrid(1,mf,itheta)*conjg(xnm2(1,n,m,itheta))+
     $     fcgrid(2,mf,itheta)*conjg(xnm2(2,n,m,itheta))
     $     )

        ucoef(n,m)=ucoef(n,m)+ws(itheta)*
     $     (
     $     fcgrid(1,mf,itheta)*conjg(unm2(1,n,m,itheta))+
     $     fcgrid(2,mf,itheta)*conjg(unm2(2,n,m,itheta))
     $     )
        enddo
        enddo
c
c
 1400   continue
c
c
        return
        end
c
c
c
c
        subroutine vshcoef2
     $     (fgrid,nterms,nphi,ntheta,ycoef,xcoef,ucoef)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the coefficients of vector spherical
c       harmonic expansion given the values of the vector field on the grid
c
c       Fast O(p^3) routine
c
c       Fast algorithm, perform FFT along all parallels
c
c          Input parameters:
c
c       fgrid (complex *16)(3,nphi,ntheta) - H field grid
c       nterms - the number of terms in vsh expansion
c       nphi, ntheta - the number of grid points in each direction
c
c          Output parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c       xcoef (complex*16)(0:nterms,-nterms:nterms) 
c       ucoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of vector spherical harmonic expansion
c        
        complex *16 ycoef(0:nterms,-nterms:nterms)
        complex *16 xcoef(0:nterms,-nterms:nterms)
        complex *16 ucoef(0:nterms,-nterms:nterms)
c
        complex *16 fgrid(3,nphi,ntheta)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: fcgrid(:,:,:)

        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c       
        allocate( fcgrid(3,nphi,ntheta) )

        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 

        allocate( ynm(0:nterms,-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do n=0,nterms
        do m=-nterms,nterms
           ycoef(n,m)=0
           xcoef(n,m)=0
           ucoef(n,m)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        nnodes=nphi*ntheta
c
c       ... first, convert vector field grids to spherical coordinates
c       and perform FFT along each parallel
c
        done=1
        pi4=16*atan(done)

        do i=1,ntheta
           do j=1,nphi
              do k=1,3
                 fcgrid(k,j,i)=fgrid(k,j,i)/pi4
              enddo
           enddo
        enddo
c


        call vshsphfft(fcgrid,nphi,ntheta,ts,ws)
c
c
        do 1400 itheta=1,ntheta
c
        costheta=ts(itheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m))
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m))
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m))
        unm2(1,n,m)=conjg(punm2(1,n,-m))
        unm2(2,n,m)=conjg(punm2(2,n,-m))
        enddo
        enddo
c
c
c       ... compute coefficients via inner product
c
c
        do n=0,nterms
        do m=-n,n
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ycoef(n,m)=ycoef(n,m)+ws(itheta)*
     $     fcgrid(3,mf,itheta)*conjg(ynm(n,m))

        xcoef(n,m)=xcoef(n,m)+ws(itheta)*
     $     (
     $     fcgrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     fcgrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )

        ucoef(n,m)=ucoef(n,m)+ws(itheta)*
     $     (
     $     fcgrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     fcgrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        enddo
        enddo
c
c
 1400   continue
c
c
ccc        call prin2('bmpole=*',bmpole,(nterms+1)*(2*nterms+1)*2)
c

        return
        end
c
c
c
c
        subroutine vshevalg
     $     (ycoef,xcoef,ucoef,nterms,nphi,ntheta,
     $     ynm,xnm2,unm2,fgrid)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates the vector field at the grids 
c       given its coefficients of vector spherical harmonic expansion
c
c       Fast O(p^3) routine
c
c          Input parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c       xcoef (complex*16)(0:nterms,-nterms:nterms) 
c       ucoef (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of vector spherical harmonic expansion
c       nterms - the number of terms in vsh expansion
c
c       nphi,ntheta - the number of points along each direction
c       ynm,xnm2,unm2 - precomputed vsh function values
c
c          Output parameters:
c
c       fgrid (complex*16) - the values of the vector field at the grids
c
c
        complex *16 ycoef(0:nterms,-nterms:nterms)
        complex *16 xcoef(0:nterms,-nterms:nterms)
        complex *16 ucoef(0:nterms,-nterms:nterms)
c
        complex *16 ynm(0:nterms,-nterms:nterms,ntheta)
        complex *16 xnm2(2,0:nterms,-nterms:nterms,ntheta)
        complex *16 unm2(2,0:nterms,-nterms:nterms,ntheta)
c
        complex *16 fgrid(3,nphi,ntheta)
c
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
           fgrid(1,iphi,itheta)=0
           fgrid(2,iphi,itheta)=0
           fgrid(3,iphi,itheta)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
c
        do 2000 k=1,ntheta
c
c
c       ... evaluate vector fields via vector spherical harmonic expansions
c
c
        do n=0,nterms
        do m=-n,n
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        fgrid(1,mf,k)=fgrid(1,mf,k)+xcoef(n,m)*xnm2(1,n,m,k)
        fgrid(2,mf,k)=fgrid(2,mf,k)+xcoef(n,m)*xnm2(2,n,m,k)
c
        fgrid(1,mf,k)=fgrid(1,mf,k)+ucoef(n,m)*unm2(1,n,m,k)
        fgrid(2,mf,k)=fgrid(2,mf,k)+ucoef(n,m)*unm2(2,n,m,k)
c
        fgrid(3,mf,k)=fgrid(3,mf,k)+ycoef(n,m)*ynm(n,m,k)
c
        endif
c
        enddo
        enddo
c
c
 2000   continue
c       
c
        zshift=0
        radius=1
        call vshfftcar(zshift,radius,fgrid,nphi,ntheta,ts,ws)
c

        return
        end
c
c
c
c
        subroutine vshevalg2
     $     (ycoef,xcoef,ucoef,nterms,nphi,ntheta,
     $     fgrid)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates the vector field at the grids 
c       given its coefficients of vector spherical harmonic expansion
c
c       Fast O(p^3) routine
c
c          Input parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c       xcoef (complex*16)(0:nterms,-nterms:nterms) 
c       ucoef (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of vector spherical harmonic expansion
c       nterms - the number of terms in vsh expansion
c
c       nphi,ntheta - the number of points along each direction
c
c
c          Output parameters:
c
c       fgrid (complex*16) - the values of the vector field at the grids
c
c
        complex *16 ycoef(0:nterms,-nterms:nterms)
        complex *16 xcoef(0:nterms,-nterms:nterms)
        complex *16 ucoef(0:nterms,-nterms:nterms)

        complex *16 fgrid(3,nphi,ntheta)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )

        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 

        allocate( ynm(0:nterms,-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       

        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
           fgrid(1,iphi,itheta)=0
           fgrid(2,iphi,itheta)=0
           fgrid(3,iphi,itheta)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
c
        do 2000 k=1,ntheta
c
c       ... evaluate xnm3, unm3, ynm
c
        costheta=ts(k)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m)) 
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m)) 
        unm2(1,n,m)=conjg(punm2(1,n,-m)) 
        unm2(2,n,m)=conjg(punm2(2,n,-m)) 
        enddo
        enddo
c
c
c
c       ... evaluate vector fields via vector spherical harmonic expansions
c
c
        do n=0,nterms
        do m=-n,n
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        fgrid(1,mf,k)=fgrid(1,mf,k)+xcoef(n,m)*xnm2(1,n,m)
        fgrid(2,mf,k)=fgrid(2,mf,k)+xcoef(n,m)*xnm2(2,n,m)
c
        fgrid(1,mf,k)=fgrid(1,mf,k)+ucoef(n,m)*unm2(1,n,m)
        fgrid(2,mf,k)=fgrid(2,mf,k)+ucoef(n,m)*unm2(2,n,m)
c
        fgrid(3,mf,k)=fgrid(3,mf,k)+ycoef(n,m)*ynm(n,m)
c
        endif
c
        enddo
        enddo
c
c
 2000   continue
c       
c
        zshift=0
        radius=1
        call vshfftcar(zshift,radius,fgrid,nphi,ntheta,ts,ws)
c

        return
        end
c
c
c
c        
        subroutine vsheval
     $     (ycoef,xcoef,ucoef,nterms,target,
     $     fval)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates the vector field at the target point
c       given its coefficients of vector spherical harmonic expansion
c
c          Input parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c       xcoef (complex*16)(0:nterms,-nterms:nterms) 
c       ucoef (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of vector spherical harmonic expansion
c       nterms - the number of terms in vsh expansion
c       target - the xyz coordinates of the target point
c
c
c          Output parameters:
c
c       fval (complex*16) - the value of the vector field at the target
c
c
        complex *16 ycoef(0:nterms,-nterms:nterms)
        complex *16 xcoef(0:nterms,-nterms:nterms)
        complex *16 ucoef(0:nterms,-nterms:nterms)

        dimension target(3)

        complex *16 fval(3),gval(3),ephi,ephim,ima
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)

        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        data ima/(0.0d0,1.0d0)/
c
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )

        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 

        allocate( ynm(0:nterms,-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
c
        x=target(1)
        y=target(2)
        z=target(3)

        r=sqrt(x**2+y**2+z**2)
        costheta = z/r
        
        rho=sqrt(x**2+y**2)
        if (rho .gt. 1d-14) then
           ephi=x/rho+ima*y/rho
        else
           ephi=1
        endif
c
c       ... evaluate xnm3, unm3, ynm
c
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c

        do n=0,nterms
        ephim=1

        do m=0,n
        ynm(n,m)=pnm(n,m)*ephim
        xnm2(1,n,m)=pxnm2(1,n,m)*ephim
        xnm2(2,n,m)=pxnm2(2,n,m)*ephim
        unm2(1,n,m)=punm2(1,n,m)*ephim
        unm2(2,n,m)=punm2(2,n,m)*ephim
        
        ephim=ephim*ephi
        enddo
        
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m))    
        xnm2(1,n,m)=conjg(xnm2(1,n,-m)) 
        xnm2(2,n,m)=conjg(xnm2(2,n,-m)) 
        unm2(1,n,m)=conjg(unm2(1,n,-m)) 
        unm2(2,n,m)=conjg(unm2(2,n,-m)) 
        enddo
        enddo
c
c
c
c       ... evaluate vector fields via vector spherical harmonic expansions
c
c

        do i=1,3
           gval(i)=0
        enddo

        do n=0,nterms
        do m=-n,n
c
        gval(1)=gval(1)+xcoef(n,m)*xnm2(1,n,m)
        gval(2)=gval(2)+xcoef(n,m)*xnm2(2,n,m)
c
        gval(1)=gval(1)+ucoef(n,m)*unm2(1,n,m)
        gval(2)=gval(2)+ucoef(n,m)*unm2(2,n,m)
c
        gval(3)=gval(3)+ycoef(n,m)*ynm(n,m)
        enddo
        enddo
c
        itype=2
        call vshsphcar2(itype,target,gval,fval)
c
        return
        end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Auxilary files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine vshsphfft(cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 work(nphi)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        real *8, allocatable :: wsave(:)
c
        allocate( wsave(10*nphi+15) )
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        call zffti(nphi,wsave)
c
        done=1
        pi=4*atan(done)
c
c
c       ... project and convert to spherical coordinates
c
        zshift=0
        zradius=1
        itype=1
        call vshsphcar(itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
c
c       ... perform the Fourier transform along each parallel
c
        do 1200 i=1,ntheta
c
        do j=1,nphi
        work(j)=cvals(1,j,i)
        enddo
        call zfftf(nphi,work,wsave)
        do j=1,nphi
        cvals(1,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(2,j,i)
        enddo
        call zfftf(nphi,work,wsave)
        do j=1,nphi
        cvals(2,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(3,j,i)
        enddo
        call zfftf(nphi,work,wsave)
        do j=1,nphi
        cvals(3,j,i)=work(j)
        enddo
c
 1200   continue
c
c
c       ... multiply by quadrature weights
c
        scale=(2*pi)/dble(nphi)
c
        do 1500 i=1,ntheta
        do 1400 j=1,nphi
        cvals(1,j,i)=cvals(1,j,i)*scale
        cvals(2,j,i)=cvals(2,j,i)*scale
        cvals(3,j,i)=cvals(3,j,i)*scale
 1400   continue
 1500   continue
c
c
        return
        end
c
c
c
c
c
        subroutine vshctheta(z,r,ctheta)
        implicit real *8 (a-h,o-z)
c
        if( abs(r) .gt. 0 ) then
        ctheta = z/r
        else
        ctheta = 0.0d0
        endif
c
        return
        end
c
c
c
c
c
        subroutine vshfftcar(zshift,zradius,cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        real *8, allocatable :: wsave(:)
        complex *16, allocatable :: work(:)
c
        allocate( wsave(10*nphi+15) )
        allocate( work(nphi) )
c
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        call zffti(nphi,wsave)
c
        done=1
        pi=4*atan(done)
c
c
c       ... perform the Fourier transform along each parallel
c
        do 1200 i=1,ntheta
c
        do j=1,nphi
        work(j)=cvals(1,j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(1,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(2,j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(2,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(3,j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(3,j,i)=work(j)
        enddo
c
 1200   continue
c
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do 1500 i=1,ntheta
        do 1400 j=1,nphi
        cvals(1,j,i)=cvals(1,j,i)*scale
        cvals(2,j,i)=cvals(2,j,i)*scale
        cvals(3,j,i)=cvals(3,j,i)*scale
 1400   continue
 1500   continue
c
c
c
c       ... convert the spherical coordinates to cartesian
c
        itype=2
        call vshsphcar(itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine vshsphcar
     $     (itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
        dimension ts(1),ws(1)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        done=1
        pi=4*atan(done)
c
        do i=1,ntheta
c
        z=ts(i)
c
        z=zshift+ts(i)*zradius
        x=sqrt(1-ts(i)*ts(i))*zradius
        r=sqrt(x*x+z*z)
        z=z/r
c
        do j=1,nphi
c
           phi=(2*pi)*(j-1)/dble(nphi)
           x=cos(phi)*sqrt(1-z**2)
           y=sin(phi)*sqrt(1-z**2)
c
           if( itype .eq. 1 ) then
c
c       ... convert the cartesian coordinates to spherical
c
           c(1,1)=-sin(phi)
           c(2,1)=cos(phi)
           c(3,1)=0
c
           c(1,2)=+cos(phi)*z
           c(2,2)=+sin(phi)*z
           c(3,2)=-sqrt(1-z**2) 
c
           c(1,3)=x
           c(2,3)=y
           c(3,3)=z
c
           endif
c
           if( itype .eq. 2 ) then
c
c       ... convert the spherical coordinates to cartesian
c
c       ... this is a transpose (inverse) matrix
c
           c(1,1)=-sin(phi)
           c(1,2)=cos(phi)
           c(1,3)=0
c
           c(2,1)=+cos(phi)*z
           c(2,2)=+sin(phi)*z
           c(2,3)=-sqrt(1-z**2) 
c
           c(3,1)=x
           c(3,2)=y
           c(3,3)=z
c
           endif
c
           c1=
     $        cvals(1,j,i)*c(1,1)+
     $        cvals(2,j,i)*c(2,1)+
     $        cvals(3,j,i)*c(3,1)
c
           c2=
     $        cvals(1,j,i)*c(1,2)+
     $        cvals(2,j,i)*c(2,2)+
     $        cvals(3,j,i)*c(3,2)
c
           c3=
     $        cvals(1,j,i)*c(1,3)+
     $        cvals(2,j,i)*c(2,3)+
     $        cvals(3,j,i)*c(3,3)
c
           cvals(1,j,i)=c1
           cvals(2,j,i)=c2
           cvals(3,j,i)=c3
c
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
        subroutine vshsphcar2
     $     (itype,xyz,vin,vout)
        implicit real *8 (a-h,o-z)
        complex *16 vin(3),vout(3)
        dimension c(3,3),xyz(3)
c
        done=1
        pi=4*atan(done)
c
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
        r=sqrt(x*x+y*y+z*z)
c
        rho=sqrt(x*x+y*y)

        if (rho .lt. 1d-14) then
           if (z.gt.0) then
           vout(1)=vin(2)
           vout(2)=vin(1)
           vout(3)=vin(3)
           else
           vout(1)=-vin(2)
           vout(2)=vin(1)
           vout(3)=-vin(3)
           endif
        else

           sinphi=y/rho
           cosphi=x/rho

           x=x/r
           y=y/r
           z=z/r
c
           if( itype .eq. 1 ) then
c
c       ... convert the cartesian coordinates to spherical
c
           c(1,1)=-sinphi
           c(2,1)=cosphi
           c(3,1)=0
c
           c(1,2)=+cosphi*z
           c(2,2)=+sinphi*z
           c(3,2)=-sqrt(1-z**2) 
c
           c(1,3)=x
           c(2,3)=y
           c(3,3)=z
c
           endif
c
           if( itype .eq. 2 ) then
c
c       ... convert the spherical coordinates to cartesian
c
c       ... this is a transpose (inverse) matrix
c
           c(1,1)=-sinphi
           c(1,2)=cosphi
           c(1,3)=0
c
           c(2,1)=+cosphi*z
           c(2,2)=+sinphi*z
           c(2,3)=-sqrt(1-z**2) 
c
           c(3,1)=x
           c(3,2)=y
           c(3,3)=z
c
           endif
c
           vout(1)=
     $        vin(1)*c(1,1)+
     $        vin(2)*c(2,1)+
     $        vin(3)*c(3,1)
c
           vout(2)=
     $        vin(1)*c(1,2)+
     $        vin(2)*c(2,2)+
     $        vin(3)*c(3,2)
c
           vout(3)=
     $        vin(1)*c(1,3)+
     $        vin(2)*c(2,3)+
     $        vin(3)*c(3,3)
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
c
