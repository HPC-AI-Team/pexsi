      real*8 function seekeig(cdf, n, x, y, mu)
      implicit none
!
!     --- use Newton's method or bisection to find a mu
!         with a prescribed cdf value ---
!
      integer :: n, cdf
      real*8 :: mpci
      real*8 :: x(n), y(n)
!
!     initial guess of mu
!
      real*8 :: mu 

      real*8 :: f, r, d, lb, ub
      integer, parameter :: maxiter = 100
      real*8, parameter :: h = 0.1, ftol = 1.5;
      integer :: iter
!
      do iter = 1, maxiter
         f = mpci(mu, n, x, y)
         r = f - dble(cdf)
         if (abs(r) < ftol) then
            write(6,*) 'converged!'
            goto 100 
         endif
         !
         ! derivative
         !
         d = (mpci(mu+h, n, x, y) - f)/h
         if (abs(d) .gt. 1e-6) then
            !
            ! use Newton's method
            !
            mu = mu - r/d
         else
            print*, 'derivative too small, d = ', d
            !
            ! use bisection
            !
            if ( abs(f-dble(cdf)) < ftol) then
               print*, 'terminating...'
               goto 100 
            endif

            if (mu .gt. x(n)) mu = x(n)
            if (mu .lt. x(1)) mu = x(1)
            f = mpci(mu,n, x, y)
            print*, 'f = ', f, 'cdf = ', cdf
            if (mu .le. x(1) .and. f .gt. dble(cdf)) then
               print*, 'solution is out of bound, increase nbins!'
               goto 100
            endif
            if (mu .ge. x(n) .and. f .le. dble(cdf)) then
               print*, 'solution is out of bound, increase nbins!'
               goto 100
            endif

            if (f .gt. dble(cdf)) then
              lb = x(1)
              ub = mu  
            else
              lb = mu
              ub = x(n) 
            endif
            if (ub-lb .gt. lb*1e-6) then
               print*, 'call bisection...'
               call bisect(n, x, y, dble(cdf), lb, ub, mu)
            endif
            goto 100
         endif
         write(6,200) iter, mu, r
 200     format(1x, ' iter = ', I2, ' mu = ', 1pe11.3, ' r = ',
     &          1pe11.3)
      end do
 100  continue
      seekeig = mu
      return
      end function seekeig
!
      subroutine bisect(n, x, y, f, lb, ub, mu)
      implicit none
      integer :: n
      real*8 :: x(n), y(n)
      real*8 :: lb, ub, mu, f, r
      integer, parameter :: maxitr = 100
      integer :: iter
      real*8 :: mpci
!
      mu = (lb+ub)/2.0
      r = mpci(mu, n, x, y) - f
      iter = 1
      print*, 'r = ', r, ' lb = ', lb, ' ub = ', ub 
      do while ( abs(r) > 2 .and. iter .lt. maxitr .and. 
     &           ub-lb .gt. lb*1.0e-6) 
         if (r .gt. 0) then
            ub = mu
         else
            lb = mu
         endif
         mu = (lb + ub)/2.0
         r = mpci(mu, n, x, y) - f
         write(6,111) iter, mu, r
 111     format('bisect iter = ', I5, ' mu = ', 1pe15.6, 
     &          ' r = ', 1pe11.3e) 
         iter = iter + 1
      end do
      end subroutine bisect

