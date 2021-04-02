        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),center(3)
        
        complex *16 rk,fval(3)

        complex *16, allocatable :: ycoef(:,:)
        complex *16, allocatable :: xcoef(:,:)
        complex *16, allocatable :: ucoef(:,:)

        complex *16, allocatable :: ynm(:,:,:)
        complex *16, allocatable ::  xnm2(:,:,:,:)
        complex *16, allocatable ::  unm2(:,:,:,:)

        complex *16, allocatable :: fgrid(:,:,:)
        complex *16, allocatable :: fcgrid(:,:,:)

        real *8, allocatable :: rnodes(:,:,:)
        real *8, allocatable :: weights(:,:)

        complex *16 cjvec(3),evec(3),hvec(3),d,df1,df2,df3,dn
        complex *16 ima
        
        data ima/(0.0d0,1.0d0)/
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
c
        print *, 'ENTER nterms'
        read *, nterms
c
        call prinf('nterms=*',nterms,1)
c
        rk=5
c
        center(1)=0
        center(2)=0
        center(3)=0
c
        source(1)=center(1)+0.2
        source(2)=center(2)+0.67
        source(3)=center(3)+0.1
c       
        cjvec(1)=3+ima*2
        cjvec(2)=2+ima*1.33
        cjvec(3)=1+ima*1
c

        nphi=2*nterms
        ntheta=2*nterms

        allocate( ycoef(0:nterms,-nterms:nterms) )
        allocate( xcoef(0:nterms,-nterms:nterms) )
        allocate( ucoef(0:nterms,-nterms:nterms) )

        allocate ( ynm(0:nterms,-nterms:nterms,ntheta) )
        allocate ( xnm2(2,0:nterms,-nterms:nterms,ntheta) )
        allocate ( unm2(2,0:nterms,-nterms:nterms,ntheta) )

        allocate( fgrid(3,nphi,ntheta) )
        allocate( fcgrid(3,nphi,ntheta) )

        allocate( rnodes(3,nphi,ntheta) )
        allocate( weights(nphi,ntheta) )
c
c
c
        call vshgrid(nphi,ntheta,rnodes,weights,nnodes)

        do i=1,ntheta
           do j=1,nphi
              target(1)=rnodes(1,j,i)
              target(2)=rnodes(2,j,i)
              target(3)=rnodes(3,j,i)

              call dipole3et(rk,source,target,cjvec,evec)

              fgrid(1,j,i)=evec(1)
              fgrid(2,j,i)=evec(2)
              fgrid(3,j,i)=evec(3)
           enddo
        enddo

        call vshgfun(nterms,ntheta,ynm,xnm2,unm2)
        call vshcoef
     $     (fgrid,nterms,nphi,ntheta,ynm,xnm2,unm2,
     $      ycoef,xcoef,ucoef)

c        call vshcoef2(fgrid,nterms,nphi,ntheta,ycoef,xcoef,ucoef)

        call vshevalg(ycoef,xcoef,ucoef,nterms,nphi,ntheta,
     $     ynm,xnm2,unm2,fcgrid)

c        call vshevalg2(ycoef,xcoef,ucoef,nterms,nphi,ntheta,
c     $     fcgrid)

        d=0
        dn=0
        do i=1,ntheta
           do j=1,nphi
              df1=fgrid(1,j,i)-fcgrid(1,j,i)
              df2=fgrid(2,j,i)-fcgrid(2,j,i)
              df3=fgrid(3,j,i)-fcgrid(3,j,i)

              d=d+df1*conjg(df1)+df2*conjg(df2)+df3*conjg(df3)
              
              dn=dn+fgrid(1,j,i)*conjg(fgrid(1,j,i))
     $           +fgrid(2,j,i)*conjg(fgrid(2,j,i))
     $           +fgrid(3,j,i)*conjg(fgrid(3,j,i))

          enddo
        enddo

        d=dble(d)
        d=sqrt(d)

        dn=dble(dn)
        dn=sqrt(dn)

        call prin2('relative l2 error=*',d/dn,1)

        done=1
        pi=4*atan(done)
c
        d=hkrand(second())
        do i=1,5
c        theta=hkrand(second())*pi
c        phi=hkrand(second())*2*pi

        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        r=1

        target(1)=cos(phi)*sin(theta)*r
        target(2)=sin(phi)*sin(theta)*r
        target(3)=cos(theta)*r
c
        call prin2('target=*',target,3)

c        target(1)=0
c        target(2)=1
c        target(3)=0

        call dipole3et(rk,source,target,cjvec,evec)
        
        call vsheval
     $     (ycoef,xcoef,ucoef,nterms,target,
     $     fval)


        df1=fval(1)-evec(1)
        df2=fval(2)-evec(2)
        df3=fval(3)-evec(3)
        d=df1*conjg(df1)+df2*conjg(df2)+df3*conjg(df3)
        d=sqrt(dble(d))
        call prin2('relative l2 error=*',d,1)
        enddo
c        
        return
        end
c
c
c
c
c
        subroutine dipole3et(rk,source,target,cjvec,evec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the arbitrary 
c       source location
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),evec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3e(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine green3e(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the electric field 
c       (free space, Lorenz gauge)
c
c                 [        \grad  \grad  ]   exp(I*rk*R)
c       green3e = [ I_3 +  ------------  ]   -----------
c                 [            rk^2      ]        R
c
c       where I_3 is a 3x3 identity matrix
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,qmat(3,3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=exp(ima*rk*cd)/cd
c
        fvec(1)=fout*cjvec(1)
        fvec(2)=fout*cjvec(2)
        fvec(3)=fout*cjvec(3)
c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)*(1-ima*rk*cd)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)*(1-ima*rk*cd)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)*(1-ima*rk*cd)
c
        qmat(1,1)=qmat(1,1)+(-rk**2*dx**2*cd**2)
        qmat(2,2)=qmat(2,2)+(-rk**2*dy**2*cd**2)
        qmat(3,3)=qmat(3,3)+(-rk**2*dz**2*cd**2)
c
        qmat(1,2)=(3-rk**2*cd**2-3*ima*rk*cd)*(dx*dy)
        qmat(2,3)=(3-rk**2*cd**2-3*ima*rk*cd)*(dy*dz)
        qmat(3,1)=(3-rk**2*cd**2-3*ima*rk*cd)*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=exp(ima*rk*cd)/cd**5/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
        return
        end
c
c
c
c
c
