c
c Driver for the code ROWMAP with example NILIDI
c
c Version of January 10, 1997.
c
c link driver rowmap dblas
c
c-----!----------------------------------------------------------------- 

      program driver
      implicit none
      integer n,liwork,lwork,nx,ipar,idid
      parameter(liwork=100, lwork=500000, nx=5000)
      real*8 u(nx),t,tend,atol,rtol,work(lwork),hs,rpar
      external fex,jacv,solout,fdt
      integer iwork(liwork),itol,ifcn,ijacv,ifdt,iout,i

c  set the initial values
      call init(n,u,t,tend,ifcn,rpar,ipar)
c  tolerances
      atol=1d-4
      rtol=1d-4
c  atol and rtol are scalars.
      itol=0
c  use default initial step size
      hs=0d0
c  use finite differences for jacobian-vector products
      ijacv=0
c  problem is autonomous, don't call a subroutine computing the
c  derivative f_t:
      ifdt=0
c  don't call solout
      iout=0 
c  clear work and iwork, i.e. use the default values for optional input.
      do 100 i=1,20
        iwork(i)=0
        work(i)=0d0
 100  continue
      rpar=0d0
      ipar=0

c
c Call to ROWMAP.
c 
      call rowmap(n,fex,ifcn,t,u,tend,hs,rtol,atol,itol,
     1                  jacv,ijacv,fdt,ifdt,solout,iout,work,
     2                  lwork,iwork,liwork,rpar,ipar,idid)

      write (*,9990) iwork(5),iwork(6),iwork(7),iwork(8),
     1                float(iwork(8))/float(iwork(5))
      write (*,9992) iwork(9),iwork(10)
 9990 format (/ 'Statistics: '/,
     1 'number of steps     =',i5,5x,'number of rejected steps  =',i5/,
     2 'number of f evals.  =',i5,5x,'number of jac-vector op.  =',i5/,
     3 'average krylov dim. =',f7.2)
 9992 format(
     1 'minimum for lwork   =',i7,3x,'minimum for liwork        =',i5)

      end


c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   F D T 
c-----!-----------------------------------------------------------------
c A dummy only, the problem is autonomous.
c
      subroutine fdt(n,t,u,ft,rpar,ipar) 
      real*8 t,u(n),ft(n),rpar(*)
      integer ipar(*)
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   S O L O U T
c-----!-----------------------------------------------------------------
c This subroutine is for output of the numerical solution. It applies
c the dense output subroutine "ROWCON" at times t=0.1,0.2,..,1.0;
c ROWMAP calls this subroutine solout only if "iout" is set equal to 1 above.
c
      subroutine solout(n,told,tnew,uold,unew,fold,fnew,ucon,
     1                  intr,rpar,ipar)
      implicit none
      integer n,intr,ipar(*)
      real*8 told,tnew,uold(n),unew(n),fold(n),fnew(n),ucon(n)
      real*8 rpar(*),s
      s=rpar(1)
      if(s.eq.0d0) write (6,9990)  
      if (s.ge.tnew) return
c Call for dense output.
1     call rowcon(n,s,told,tnew,uold,unew,fold,fnew,ucon)
      write (6,9991) s,ucon(1)
      s=s+1d-1
      if (s.lt.tnew) goto 1
c Save last evaluation point:
      rpar(1)=s
 9990 format ('Solution u(t,x,y) at grid point (1,1) at time ' )
 9991 format ('t=',f5.1,2x,'is',2x,d12.6) 
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   F E X 
c-----!-----------------------------------------------------------------
c An ODE system is generated from the following two-dimensional 
c nonlinear diffusion equation
c
c     du           [  d^2 u     d^2 u  ]
c     --  = exp(u) [  -----  +  -----  ] + u*(18 exp(u) - 1),
c     dt           [  dx^2      dy^2   ]
c
c     t in [0,1], (x,y) in [0,pi/3]^2
c
c with Dirichlet boundary conditions and initial values chosen
c to yield the exact solution
c
c     u(t,x,y) = exp(-t) * sin(3x) * sin(3y).
c
c The PDE system is treated by central differences on a uniform
c mesh with 69 x 69 inner grid points.
c

      subroutine fex(n,t,u,fw,rpar,ipar)
      implicit none
      integer n,m, i,j,ipar(*)
      parameter (m=69)
      real*8 t,u(m,m),fw(m,m), pi,x1,x2,fk2,eu,dx,d2,rpar(*)
      data pi/ 3.141592654D0/
      fk2 = (m+1)*(m+1)*9D0/(pi*pi)
      dx = 1D0/(m+1)
      do 100 i=1,m
       x1 = i*dx
       do 100 j=1,m
       x2 = j*dx
       d2 = -4D0*u(i,j)
       if (i.gt.1) d2 = d2+u(i-1,j)
       if (i.lt.m) d2 = d2+u(i+1,j)
       if (j.gt.1) d2 = d2+u(i,j-1)
       if (j.lt.m) d2 = d2+u(i,j+1)
       eu = dexp(u(i,j))
       fw(i,j) = eu*fk2*d2+u(i,j)*(18D0*eu-1D0)
  100 continue
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   J A C V 
c-----!-----------------------------------------------------------------
c This subroutine computes z:=A*y with A=Jacobian of fex.
c It is only called if "ijacv" is set equal 1 above.
c
c
      subroutine jacv(n,t,u,y,z,rpar,ipar)
      implicit none
      integer n,m, i,j,ipar(*)
      parameter (m=69)
      real*8 t,u(m,m),y(m,m),z(m,m), pi,fk2,eu,d2,du,rpar(*)
      data pi/3.14159265358979D0/
      fk2 = (m+1)*(m+1)*9D0/(pi*pi)
      do 100 i=1,m
       do 100 j=1,m
       eu = dexp(u(i,j))
       d2 = -4D0*y(i,j)
       du = -4D0*u(i,j)
       if (i.gt.1) then
        d2 = d2+y(i-1,j)
        du = du+u(i-1,j)
       endif
       if (i.lt.m) then
        d2 = d2+y(i+1,j)
        du = du+u(i+1,j)
       endif
       if (j.gt.1) then
        d2 = d2+y(i,j-1)
        du = du+u(i,j-1)
       endif
       if (j.lt.m) then
        d2 = d2+y(i,j+1)
        du = du+u(i,j+1)
       endif
       z(i,j) = eu*fk2*(d2+du*y(i,j))+(18D0*eu*(1D0+u(i,j))-1D0)*y(i,j)
  100 continue
      return
      end

c-----!-----------------------------------------------------------------
c --- S U B R O U T I N E   I N I T 
c-----!-----------------------------------------------------------------
c
c
      subroutine init(n,u,t0,te,ifcn,rpar,ipar)
      implicit none
      integer n,m, i,j,ifcn,ipar(*)
      parameter (m=69)
      real*8 u(m,m),t0,te, x1,x2,dx,pi,rpar(*)
      data pi/3.14159265358979D0/
c The problem is autonomous.
      ifcn=0
c Number of equations.
      n = m*m
c Initial values.
      t0 = 0D0
      dx = pi/(m+1)
      do 100 i=1,m
       x1 = i*dx
       do 100 j=1,m
       x2 = j*dx
  100  u(i,j) = dsin(x1)*dsin(x2)
c Endpoint of integration.
      te = 1d0
      return
      end

