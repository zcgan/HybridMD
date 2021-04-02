
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Scalar spherical harmonic transforms
c
c       Fortran F90 version, allocatable arrays
c
c       NOTE #1: heavy use of allocatable arrays to simplify programming 
c       We use automatic deallocation feature of Fortran F95
c
c       sshgrid - form the nodes and weights for quadrature on the unit sphere
c                 legendre nodes along theta, and equispaced point along phi
c
c       sshgfun - precompute the ssh functions needed in ssh transforms
c                 memory intensive O(p^3) storage, suitable for repeated
c                 calls to ssh transform.
c
c       sshcoef - compute the coefficients of scalar spherical harmonics
c                 given a scalar function on the unit sphere grid
c                 needs to call sshgfun first to obtain ssh function values.
c       
c       sshevalg - evaluate the scalar function on the unit sphere grid given 
c                 the coefficients of its scalar spherical harmonic expansion
c                 needs to call sshgfun first to obtain ssh function values.
c
c       sshcoef2 - compute the coefficients of scalar spherical harmonics
c                 given a scalar function on the unit sphere grid
c                 compute ssh function on the fly.
c                 adapted from em3ehformex_fast
c       
c       sshevalg2 - evaluate the scalar function on the unit sphere grid given 
c                 the coefficients of its scalar spherical harmonic expansion
c                 compute ssh function on the fly
c                 adapted from em3exevalehz_fast
c
c       ssheval - evaluate the scalar function at a given point in R^3 
c                 given the coefficients of its scalar spherical 
c                 harmonic expansion
c
c
c
        subroutine sshgrid(nphi,ntheta,rnodes,weights,nnodes)
c       This subroutine returns the xyz coordinates of the grid 
c       on the unit sphere,
c       the grid is used in (scalar) spherical harmonic transform.
c
c       Input: nphi, ntheta 
c       Output: rnodes, weights, nnodes
c
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nphi,ntheta),weights(nphi,ntheta)
c
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
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
        subroutine sshgfun(nterms,ntheta,ynm)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the ssh function values needed
c       in the ssh transforms.
c
c       O(p^3) work and storage routine
c
c          Input parameters:
c
c       nterms - the number of terms in ssh expansion
c       ntheta - the number of grid points in theta direction
c
c          Output parameters:
c
c       ynm (complex*16)(0:nterms,-nterms:nterms,ntheta)
c           - the values of scalar spherical harmonic functions
c        
        complex *16 ynm(0:nterms,-nterms:nterms,ntheta)
c
        real *8, allocatable :: pnm(:,:)

c
        real *8, allocatable :: ts(:), ws(:)
c       
        allocate( pnm(0:nterms,0:nterms) )

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
        call ylgndr(nterms,costheta,pnm)
c

        do n=0,nterms
        do m=0,n
        ynm(n,m,itheta)=pnm(n,m)
        enddo
        do m=-n,-1
        ynm(n,m,itheta)=ynm(n,-m,itheta)
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
        subroutine sshcoef
     $     (fgrid,nterms,nphi,ntheta,ynm,ycoef)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the coefficients of scalar spherical
c       harmonic expansion given the values of the scalar function on the grid
c
c       Fast O(p^3) routine
c
c       Fast algorithm, perform FFT along all parallels
c
c          Input parameters:
c
c       fgrid (complex *16)(nphi,ntheta) - scalar function grid
c       nterms - the number of terms in ssh expansion
c       nphi, ntheta - the number of grid points in each direction
c       ynm, - precomputed ssh function values
c
c          Output parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of scalar spherical harmonic expansion
c        
        complex *16 ycoef(0:nterms,-nterms:nterms)
c
        complex *16 fgrid(nphi,ntheta)
c
        complex *16 ynm(0:nterms,-nterms:nterms,ntheta)

        complex *16, allocatable :: fcgrid(:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c       
        allocate( fcgrid(nphi,ntheta) )

        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do n=0,nterms
        do m=-nterms,nterms
           ycoef(n,m)=0
        enddo
        enddo
c
c       ... first, convert scalar function grids to spherical coordinates
c       and perform FFT along each parallel
c
        done=1
        pi4=16*atan(done)

        do i=1,ntheta
           do j=1,nphi
                 fcgrid(j,i)=fgrid(j,i)/pi4
           enddo
        enddo
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        call sshsphfft(fcgrid,nphi,ntheta,ts,ws)
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
     $     fcgrid(mf,itheta)*conjg(ynm(n,m,itheta))
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
        subroutine sshcoef2
     $     (fgrid,nterms,nphi,ntheta,ycoef)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the coefficients of scalar spherical
c       harmonic expansion given the values of the scalar function on the grid
c
c       Fast O(p^3) routine
c
c       Fast algorithm, perform FFT along all parallels
c
c          Input parameters:
c
c       fgrid (complex *16)(nphi,ntheta) - scalar function grid
c       nterms - the number of terms in ssh expansion
c       nphi, ntheta - the number of grid points in each direction
c
c          Output parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of scalar spherical harmonic expansion
c        
        complex *16 ycoef(0:nterms,-nterms:nterms)
c
        complex *16 fgrid(nphi,ntheta)
c
        real *8, allocatable :: pnm(:,:)
        complex *16, allocatable :: fcgrid(:,:)

        complex *16, allocatable :: ynm(:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c       
        allocate( fcgrid(nphi,ntheta) )

        allocate( pnm(0:nterms,0:nterms) )

        allocate( ynm(0:nterms,-nterms:nterms) )
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do n=0,nterms
        do m=-nterms,nterms
           ycoef(n,m)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        nnodes=nphi*ntheta
c
c       ... first, convert scalar function grids to spherical coordinates
c       and perform FFT along each parallel
c
        done=1
        pi4=16*atan(done)

        do i=1,ntheta
           do j=1,nphi
                 fcgrid(j,i)=fgrid(j,i)/pi4
           enddo
        enddo
c

        call sshsphfft(fcgrid,nphi,ntheta,ts,ws)
c
c
        do 1400 itheta=1,ntheta
c
        costheta=ts(itheta)
        call ylgndr(nterms,costheta,pnm)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=ynm(n,-m)
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
     $     fcgrid(mf,itheta)*conjg(ynm(n,m))
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
        subroutine sshevalg
     $     (ycoef,nterms,nphi,ntheta,ynm,fgrid)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates the scalar function at the grids 
c       given its coefficients of scalar spherical harmonic expansion
c
c       Fast O(p^3) routine
c
c          Input parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of scalar spherical harmonic expansion
c       nterms - the number of terms in ssh expansion
c
c       nphi,ntheta - the number of points along each direction
c       ynm - precomputed ssh function values
c
c          Output parameters:
c
c       fgrid (complex*16) - the values of the scalar function at the grids
c
c
        complex *16 ycoef(0:nterms,-nterms:nterms)
c
        complex *16 ynm(0:nterms,-nterms:nterms,ntheta)
c
        complex *16 fgrid(nphi,ntheta)
c
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
           fgrid(iphi,itheta)=0
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
c       ... evaluate scalar functions via scalar spherical harmonic expansions
c
c
        do n=0,nterms
        do m=-n,n
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        fgrid(mf,k)=fgrid(mf,k)+ycoef(n,m)*ynm(n,m,k)
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
        call sshfftcar(zshift,radius,fgrid,nphi,ntheta,ts,ws)
c

        return
        end
c
c
c
c
        subroutine sshevalg2
     $     (ycoef,nterms,nphi,ntheta,
     $     fgrid)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates the scalar function at the grids 
c       given its coefficients of scalar spherical harmonic expansion
c
c       Fast O(p^3) routine
c
c          Input parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of scalar spherical harmonic expansion
c       nterms - the number of terms in ssh expansion
c
c       nphi,ntheta - the number of points along each direction
c
c
c          Output parameters:
c
c       fgrid (complex*16) - the values of the scalar function at the grids
c
c
        complex *16 ycoef(0:nterms,-nterms:nterms)

        complex *16 fgrid(nphi,ntheta)
c
        real *8, allocatable :: pnm(:,:)
c
        complex *16, allocatable :: ynm(:,:)
c
c
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( pnm(0:nterms,0:nterms) )

        allocate( ynm(0:nterms,-nterms:nterms) )
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
           fgrid(iphi,itheta)=0
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
        call ylgndr(nterms,costheta,pnm)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=ynm(n,-m) 
        enddo
        enddo
c
c
c
c       ... evaluate scalar functions via scalar spherical harmonic expansions
c
c
        do n=0,nterms
        do m=-n,n
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
c
        fgrid(mf,k)=fgrid(mf,k)+ycoef(n,m)*ynm(n,m)
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
        call sshfftcar(zshift,radius,fgrid,nphi,ntheta,ts,ws)
c

        return
        end
c
c
c
c        
        subroutine ssheval(ycoef,nterms,target,fval)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates the scalar function at the target point
c       given its coefficients of scalar spherical harmonic expansion
c
c          Input parameters:
c
c       ycoef (complex*16)(0:nterms,-nterms:nterms)
c           - the coefficients of scalar spherical harmonic expansion
c       nterms - the number of terms in ssh expansion
c       target - the xyz coordinates of the target point
c
c
c          Output parameters:
c
c       fval (complex*16) - the value of the scalar function at the target
c
c
        complex *16 ycoef(0:nterms,-nterms:nterms)

        dimension target(3)

        complex *16 fval,gval,ephi,ephim,ima
c
        real *8, allocatable :: pnm(:,:)

        complex *16, allocatable :: ynm(:,:)
c
        data ima/(0.0d0,1.0d0)/
c
        allocate( pnm(0:nterms,0:nterms) )

        allocate( ynm(0:nterms,-nterms:nterms) )
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
        call ylgndr(nterms,costheta,pnm)
c

        do n=0,nterms
        ephim=1

        do m=0,n
        ynm(n,m)=pnm(n,m)*ephim
        ephim=ephim*ephi
        enddo
        
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m))    
        enddo
        enddo
c
c
c
c       ... evaluate scalar functions via scalar spherical harmonic expansions
c
c
        fval=0
        do n=0,nterms
        do m=-n,n
c
        fval=fval+ycoef(n,m)*ynm(n,m)
        enddo
        enddo
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
        subroutine sshsphfft(cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 work(nphi)
        complex *16 cvals(nphi,ntheta)
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
c       ... perform the Fourier transform along each parallel
c
        do 1200 i=1,ntheta
c
        do j=1,nphi
        work(j)=cvals(j,i)
        enddo

        call zfftf(nphi,work,wsave)

        do j=1,nphi
        cvals(j,i)=work(j)
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
        cvals(j,i)=cvals(j,i)*scale
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
        subroutine sshfftcar(zshift,zradius,cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 cvals(nphi,ntheta)
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
        work(j)=cvals(j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(j,i)=work(j)
        enddo
c
 1200   continue
c
        return
        end
c
c
c
c
c
