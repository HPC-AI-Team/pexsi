      ! NOTE: This function is superceded by the MonotoneRootFinding
      ! routine in interface.cpp.

      real*8 function mpci(t, n, x, y)
      implicit none
      real*8, intent(in) :: t
      integer, intent(in) :: n
      real*8, intent(in) :: x(n), y(n)

      real*8 :: d(n), h(n), delta(n), xi(n) 
      integer :: i, j
      real*8 :: s

      s = 0.0
      d(1:n) = 0.0
      h(1:n) = 0.0
      delta(1:n) = 0.0
      xi(1:n) = 0.0

      do i = 1, n-1
         h(i) = x(i+1)-x(i) 
         if (h(i) .gt. 0.0) then
            delta(i) = (y(i+1)-y(i))/h(i)
         else
            write(6,*) 'i = ', i, 'x(i) not less than x(i+1)',
     &                 'x(i) = ', x(i), ' x(i+1) = ', x(i+1)
         endif
      end do 
      do i = 1, n-2
         xi(i) = (h(i)+2.0*h(i+1))/(3.0*(h(i)+h(i+1)))
      end do

      ! derivative at x_j
      
      do j = 2, n-1
         if (delta(j-1)*delta(j) .gt. 0) then
            d(j) = (delta(j-1)*delta(j))/(xi(j-1)*delta(j)
     &           + (1-xi(j-1))*delta(j-1))
         endif
      end do
!
!     not-a-knot condition
!
      d(1) = h(1)*(3.0*delta(2)-(2.0*d(2)+d(3)))/h(2)
     &     + 3.0*delta(1)-2*d(2)
      d(n) = h(n-1)*(3.0*delta(n-1)-(2.0*d(n-1)+d(n-2)))/h(n-2)
     &     + 3.0*delta(n-1) - 2.0*d(n-1)
      if (d(1) .lt. 0.0)           d(1) = 0.0
      if (d(1) .gt. 3*delta(1))    d(1) = 3.0*delta(1)
      if (d(n) .lt. 0)             d(n) = 0.0
      if (d(n) .gt. 3*delta(n-1))  d(n) = 3*delta(n-1)
!
!     evaluate the spline
!
      if (t < x(1)) then
         s = y(1)
         goto 100
      endif 

      do j = 1, n-1
         if (x(j) .le. t .and. t .le. x(j+1)) then
            s = ((d(j)+d(j+1)-2.0*delta(j))/(h(j)*h(j)))*(t-x(j))**3.0
     &        + ((-2.0*d(j)-d(j+1)+3.0*delta(j))/h(j))*(t-x(j))**2.0
     &        + d(j)*(t-x(j))+y(j) 
            goto 100
         else
            cycle
         endif
      end do
      
      if (t .ge. x(n)) s = y(n)
 100  continue
      mpci = s
      return
      end function mpci
