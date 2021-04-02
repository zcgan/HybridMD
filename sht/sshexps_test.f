        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),center(3)
        
        complex *16 rk,fval,f1(0:10 000),f2(0:10 000),f3(0:10 000)

        complex *16, allocatable :: ycoef(:,:)

        complex *16, allocatable :: ynm(:,:,:)

        complex *16, allocatable :: fgrid(:,:)
        complex *16, allocatable :: fcgrid(:,:)

        real *8, allocatable :: rnodes(:,:,:)
        real *8, allocatable :: weights(:,:)

        complex *16 evec,d,df1,df2,df3,dn,dnp,work(1 000 000)
        complex *16 ima,beta(2 000 000),z

        external besseljcoef,besseljfun01,intexppncoef,intexppnfun01

        data ima/(0.0d0,1.0d0)/
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
c
        d=hkrand(second())

        n=100
        ifder=0
        lwork=1 000 000

        do i=1,5
           z=hkrand(0)+ima*hkrand(0)
           z=-z*8**(i-1)
           scale=sqrt(dble(z*conjg(z)))
           if (scale .gt. 1) scale=1.0d0

           call cjfuns3d(ier,n,z,scale,f2,ifder,f3,work,lwork)
           call prin2('f2=*',f2,20)
           call recurrence(ier,n,z,scale,f1,besseljcoef,besseljfun01)
           call prin2('f1=*',f1,20)

           d=0
           dn=0
           do j=0,n
              d=d+(f1(j)-f2(j))*conjg(f1(j)-f2(j))
              dn=dn+f2(j)*conjg(f2(j))
           enddo
           d=sqrt(dble(d))
           dn=sqrt(dble(dn))
           call prin2('z=*',z,2)
           call prin2('relative l2 error=*',d/dn,1)
        enddo

c        return

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

        nphi=2*nterms
        ntheta=2*nterms

        allocate( ycoef(0:nterms,-nterms:nterms) )

        allocate ( ynm(0:nterms,-nterms:nterms,ntheta) )

        allocate( fgrid(nphi,ntheta) )
        allocate( fcgrid(nphi,ntheta) )

        allocate( rnodes(3,nphi,ntheta) )
        allocate( weights(nphi,ntheta) )
c
c
c
        call sshgrid(nphi,ntheta,rnodes,weights,nnodes)

        do i=1,ntheta
           do j=1,nphi
              target(1)=rnodes(1,j,i)
              target(2)=rnodes(2,j,i)
              target(3)=rnodes(3,j,i)

              x=target(1)-source(1)
              y=target(2)-source(2)
              z=target(3)-source(3)

              r=sqrt(x*x+y*y+z*z)

              fgrid(j,i)=exp(ima*rk*r)/r
           enddo
        enddo

c        call sshgfun(nterms,ntheta,ynm)
c        call sshcoef(fgrid,nterms,nphi,ntheta,ynm,ycoef)

        call sshcoef2(fgrid,nterms,nphi,ntheta,ycoef)

c        call sshevalg(ycoef,nterms,nphi,ntheta,ynm,fcgrid)

        call sshevalg2(ycoef,nterms,nphi,ntheta,fcgrid)

        d=0
        dn=0
        do i=1,ntheta
           do j=1,nphi
              df1=fgrid(j,i)-fcgrid(j,i)

              d=d+df1*conjg(df1)
              
              dn=dn+fgrid(j,i)*conjg(fgrid(j,i))

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

        x=target(1)-source(1)
        y=target(2)-source(2)
        z=target(3)-source(3)
        
        r=sqrt(x*x+y*y+z*z)
        
        evec=exp(ima*rk*r)/r
        
        call ssheval(ycoef,nterms,target,fval)


        df1=fval-evec
        d=df1*conjg(df1)
        d=sqrt(dble(d))
        call prin2('relative l2 error=*',d,1)
        enddo
c        

        ifder=1
        open(unit=26,file="Dzeros.txt")        

        kk=1
        do i=1,100
           do j=1,i+1
              read(26,*) beta(kk)
              if (i.lt. 10) then
                 call besselD(i,beta(kk),dn,ifder,dnp)
                 call prin2('dn=*',dn,2)
              endif
              kk=kk+1
           enddo
        enddo

        close(26)

        print *, 'ENTER n'
        read *, n
c
        call prinf('n=*',n,1)

        print *,beta(100)
        z=beta(100)*1.0d1

        print *, 'z=',z
        scale=sqrt(dble(z*conjg(z)))
        if (scale .gt. 1) scale=1.0d0

        call recurrence(ier,n,z,scale,f1,intexppncoef,intexppnfun01)
        print *,'ier=',ier
        print *, 'f1(n)=', f1(n)

c        call prin2('beta=*',beta,100)
        call prin2('f1=*',f1,100)

        return
        end
c
c
c
c
c
