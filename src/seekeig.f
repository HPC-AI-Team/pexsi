      real*8 function seekeig(cdf, n, x, y, mu)
      implicit none
!
!     --- use bisection to find a mu
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
      real*8, parameter :: ftol = 0.01;
      integer :: iter
!
      if (mu .gt. x(n)) mu = x(n)
      if (mu .lt. x(1)) mu = x(1)
      f = mpci(mu, n, x, y)
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
      real*8, parameter :: ftol = 0.01
      integer :: iter
      real*8 :: mpci
!
      mu = (lb+ub)/2.0
      r = mpci(mu, n, x, y) - f
      iter = 1
      print*, 'r = ', r, ' lb = ', lb, ' ub = ', ub 
      do while ( abs(r) > ftol .and. iter .lt. maxitr .and. 
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
     &          ' r = ', 1pe11.3) 
         iter = iter + 1
      end do
      end subroutine bisect

