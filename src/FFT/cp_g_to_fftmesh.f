      subroutine cp_g_to_fftmesh(g,k1d,k2d,k3d,f,nfft)

!     copies the array g(-k1d:k1d,-k2d:k2d,-k3d:k3d) to the
!     fft mesh defined by nfft(3).

      use mod_comp

      implicit none

      integer, intent(in) :: k1d,k2d,k3d,nfft(3)
!--->    note: lower bounds of g are defined here
      complex(DP), intent(in)    :: g(-k1d:,-k2d:,-k3d:)
      complex(DP), intent(inout) :: f(:,:,:)

      integer                         :: i1,i2,i3,m1,m2,m3,n1,n2,n3

!--->    get the limits and shifts
      i1 = min( k1d, (nfft(1) - 1)/2 )
      i2 = min( k2d, (nfft(2) - 1)/2 )
      i3 = min( k3d, (nfft(3) - 1)/2 )

      m1 = nfft(1) + 1 - i1
      n1 = nfft(1)
      m2 = nfft(2) + 1 - i2
      n2 = nfft(2)
      m3 = nfft(3) + 1 - i3
      n3 = nfft(3)

      f = 0.0_DP

      f( 1:i1+1,  1:i2+1,  1:i3+1)  =  g(  0:i1,   0:i2,   0:i3 )

      f(m1:n1  ,  1:i2+1,  1:i3+1)  =  g(-i1:-1,   0:i2,   0:i3 )
      f( 1:i1+1, m2:n2  ,  1:i3+1)  =  g(  0:i1, -i2:-1,   0:i3 )
      f( 1:i1+1,  1:i2+1, m3:n3  )  =  g(  0:i1,   0:i2, -i3:-1 )

      f( 1:i1+1, m2:n2  , m3:n3  )  =  g(  0:i1, -i2:-1, -i3:-1 )
      f(m1:n1  ,  1:i2+1, m3:n3  )  =  g(-i1:-1,   0:i2, -i3:-1 )
      f(m1:n1  , m2:n2  ,  1:i3+1)  =  g(-i1:-1, -i2:-1,   0:i3 )

      f(m1:n1  , m2:n2  , m3:n3  )  =  g(-i1:-1, -i2:-1, -i3:-1 )

      return
      end subroutine cp_g_to_fftmesh

      subroutine cp_fftmesh_to_g(g,k1d,k2d,k3d,f,nfft)

!     copies the values from the fft mesh defined by nfft(3) to an array
!     g(-k1d:k1d,-k2d:k2d,-k3d:k3d)

      use mod_comp

      implicit none

      integer, intent(in) :: k1d,k2d,k3d,nfft(3)
!--->    note: lower bounds of g are defined here
      complex(DP), intent(inout) :: g(-k1d:,-k2d:,-k3d:)
      complex(DP), intent(in)    :: f(:,:,:)

      integer                         :: i1,i2,i3,m1,m2,m3,n1,n2,n3

!--->    get the limits and shifts
      i1 = min( k1d, ( nfft(1) - 1)/2 )
      i2 = min( k2d, ( nfft(2) - 1)/2 )
      i3 = min( k3d, ( nfft(3) - 1)/2 )

      m1 = nfft(1) + 1 - i1
      n1 = nfft(1)
      m2 = nfft(2) + 1 - i2
      n2 = nfft(2)
      m3 = nfft(3) + 1 - i3
      n3 = nfft(3)

      g = 0.0_DP

      g(  0:i1,   0:i2,   0:i3 ) = f( 1:i1+1,  1:i2+1,  1:i3+1)

      g(-i1:-1,   0:i2,   0:i3 ) = f(m1:n1  ,  1:i2+1,  1:i3+1)
      g(  0:i1, -i2:-1,   0:i3 ) = f( 1:i1+1, m2:n2  ,  1:i3+1)
      g(  0:i1,   0:i2, -i3:-1 ) = f( 1:i1+1,  1:i2+1, m3:n3  )

      g(  0:i1, -i2:-1, -i3:-1 ) = f( 1:i1+1, m2:n2  , m3:n3  )
      g(-i1:-1,   0:i2, -i3:-1 ) = f(m1:n1  ,  1:i2+1, m3:n3  )
      g(-i1:-1, -i2:-1,   0:i3 ) = f(m1:n1  , m2:n2  ,  1:i3+1)

      g(-i1:-1, -i2:-1, -i3:-1 ) = f(m1:n1  , m2:n2  , m3:n3  )

      return
      end subroutine cp_fftmesh_to_g
