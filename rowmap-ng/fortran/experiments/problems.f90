! Test problems for the ROWMAP experiments. Definitions: docs/problems.md.
!
! The module holds the problem state; the plain external subroutines
! pfcn / pjacv / pfdt at the end have the F77 interfaces expected by rowmap.F.
module problems
  implicit none
  private
  public :: setup, pinit, prhs, pjv, pft, prob_name, prob_auton, nprob
  public :: id_of

  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

  integer, parameter :: HEAT = 1, NILIDI = 2, BRUSS = 3, FOOD = 4, FOOD3 = 5, &
                        VDPL = 6, DIFKON = 7, DIFFU2 = 8
  integer, parameter :: nprob = 8
  character(len=8), parameter :: names(nprob) = [character(len=8) :: &
       'heat', 'nilidi', 'bruss', 'food', 'food3', 'vdpl', 'difkon', 'diffu2']

  integer, save :: prob = 0, m = 0, np = 0   ! problem id, grid parameter, points per component
  real(dp), save :: h = 0, x0 = 0
  real(dp), allocatable, save :: xs(:)       ! coordinates per direction (length m)
  real(dp), allocatable, save :: bfood3(:)   ! FOOD3 b(x,y,z)
  real(dp), allocatable, save :: wrk1(:), wrk2(:)

contains

  integer function id_of(name)
    character(len=*), intent(in) :: name
    integer :: k
    id_of = 0
    do k = 1, nprob
      if (trim(names(k)) == trim(name)) id_of = k
    end do
  end function id_of

  character(len=8) function prob_name()
    prob_name = names(prob)
  end function prob_name

  logical function prob_auton()
    prob_auton = .not. (prob == DIFKON .or. prob == DIFFU2)
  end function prob_auton

  ! Set problem and grid; returns the ODE dimension n.
  subroutine setup(id, mm, n)
    integer, intent(in) :: id, mm
    integer, intent(out) :: n
    integer :: i, j, l
    prob = id
    m = mm
    if (allocated(xs)) deallocate(xs)
    allocate(xs(m))
    select case (prob)
    case (HEAT, DIFKON, DIFFU2)
      x0 = 0; h = 1.0_dp / (m + 1)
      do i = 1, m; xs(i) = x0 + i * h; end do
      np = m * m; n = np
    case (NILIDI)
      x0 = 0; h = (pi / 3) / (m + 1)
      do i = 1, m; xs(i) = x0 + i * h; end do
      np = m * m; n = np
    case (BRUSS, FOOD, VDPL)
      h = 1.0_dp / (m - 1)
      do i = 1, m; xs(i) = (i - 1) * h; end do
      np = m * m; n = 2 * np
    case (FOOD3)
      h = 1.0_dp / (m - 1)
      do i = 1, m; xs(i) = (i - 1) * h; end do
      np = m * m * m; n = 2 * np
      if (allocated(bfood3)) deallocate(bfood3)
      allocate(bfood3(np))
      do l = 1, m
        do j = 1, m
          do i = 1, m
            bfood3(i + (j - 1) * m + (l - 1) * m * m) = &
                 (1 + 0.2_dp * xs(i) * xs(j) * xs(l)) * (1d6 - 1 + 1d-6)
          end do
        end do
      end do
    case default
      stop 'setup: unknown problem'
    end select
    if (allocated(wrk1)) deallocate(wrk1, wrk2)
    allocate(wrk1(np), wrk2(np))
  end subroutine setup

  ! ---------------------------------------------------------------- Laplacians
  ! Dirichlet, zero boundary values, m x m inner points.
  subroutine lap_dir(u, out)
    real(dp), intent(in) :: u(m, m)
    real(dp), intent(out) :: out(m, m)
    integer :: i, j
    real(dp) :: s
    do j = 1, m
      do i = 1, m
        s = -4 * u(i, j)
        if (i > 1) s = s + u(i - 1, j)
        if (i < m) s = s + u(i + 1, j)
        if (j > 1) s = s + u(i, j - 1)
        if (j < m) s = s + u(i, j + 1)
        out(i, j) = s / (h * h)
      end do
    end do
  end subroutine lap_dir

  ! Central differences D_x, D_y, Dirichlet zero boundary.
  subroutine dxy_dir(u, dx, dy)
    real(dp), intent(in) :: u(m, m)
    real(dp), intent(out) :: dx(m, m), dy(m, m)
    integer :: i, j
    real(dp) :: a, b
    do j = 1, m
      do i = 1, m
        a = 0; b = 0
        if (i < m) a = u(i + 1, j)
        if (i > 1) b = u(i - 1, j)
        dx(i, j) = (a - b) / (2 * h)
        a = 0; b = 0
        if (j < m) a = u(i, j + 1)
        if (j > 1) b = u(i, j - 1)
        dy(i, j) = (a - b) / (2 * h)
      end do
    end do
  end subroutine dxy_dir

  integer pure function mlo(i)
    integer, intent(in) :: i
    if (i > 1) then; mlo = i - 1; else; mlo = 2; end if
  end function mlo

  integer pure function mhi(i, mm)
    integer, intent(in) :: i, mm
    if (i < mm) then; mhi = i + 1; else; mhi = mm - 1; end if
  end function mhi

  ! Neumann 2D, vertex grid with ghost mirror.
  subroutine lap_neu2(u, out)
    real(dp), intent(in) :: u(m, m)
    real(dp), intent(out) :: out(m, m)
    integer :: i, j
    do j = 1, m
      do i = 1, m
        out(i, j) = (u(mlo(i), j) + u(mhi(i, m), j) + u(i, mlo(j)) + u(i, mhi(j, m)) &
                     - 4 * u(i, j)) / (h * h)
      end do
    end do
  end subroutine lap_neu2

  ! Neumann 3D, vertex grid with ghost mirror.
  subroutine lap_neu3(u, out)
    real(dp), intent(in) :: u(m, m, m)
    real(dp), intent(out) :: out(m, m, m)
    integer :: i, j, l
    do l = 1, m
      do j = 1, m
        do i = 1, m
          out(i, j, l) = (u(mlo(i), j, l) + u(mhi(i, m), j, l) + u(i, mlo(j), l) &
                          + u(i, mhi(j, m), l) + u(i, j, mlo(l)) + u(i, j, mhi(l, m)) &
                          - 6 * u(i, j, l)) / (h * h)
        end do
      end do
    end do
  end subroutine lap_neu3

  subroutine lap_neu(u, out)
    real(dp), intent(in) :: u(*)
    real(dp), intent(out) :: out(*)
    if (prob == FOOD3) then
      call lap_neu3(u, out)
    else
      call lap_neu2(u, out)
    end if
  end subroutine lap_neu

  ! ------------------------------------------------------ DIFKON / DIFFU2 sources
  subroutine sources(t, g, gt)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: g(m, m), gt(m, m)
    integer :: i, j
    real(dp) :: x, y, s, cx, cy, e, g0, g1, eps
    eps = 0.1_dp
    e = exp(-2 * eps * pi * pi * t)
    do j = 1, m
      do i = 1, m
        x = xs(i); y = xs(j)
        s = sin(pi * x) * sin(pi * y)
        cx = cos(pi * x) * sin(pi * y)
        cy = sin(pi * x) * cos(pi * y)
        if (prob == DIFKON) then
          g0 = 2 * eps * pi * pi * s + pi * (y * cx - x * cy)
          g1 = -8 * eps * pi * (y * cx + x * cy) + 4 * pi * x * y * (y * cx - x * cy) &
               + 4 * s * (y * y - x * x)
          g(i, j) = g0 + e * g1
          gt(i, j) = -2 * eps * pi * pi * e * g1
        else
          g(i, j) = 4 * x * y * s * cos(t) + 2 * pi * pi * s * (1 + 4 * x * y * sin(t)) &
                    - 8 * pi * sin(t) * (y * cx + x * cy)
          gt(i, j) = -4 * x * y * s * sin(t) + 8 * pi * pi * x * y * s * cos(t) &
                     - 8 * pi * cos(t) * (y * cx + x * cy)
        end if
      end do
    end do
  end subroutine sources

  ! ------------------------------------------------------------ initial values
  subroutine pinit(n, y, t0, tend)
    integer, intent(in) :: n
    real(dp), intent(out) :: y(n), t0, tend
    integer :: i, j, l, k
    real(dp) :: x, yy, z
    t0 = 0
    select case (prob)
    case (HEAT)
      tend = 10.24_dp
      do j = 1, m; do i = 1, m
        x = xs(i); yy = xs(j)
        y(i + (j - 1) * m) = 16 * x * (1 - x) * yy * (1 - yy)
      end do; end do
    case (NILIDI)
      tend = 1
      do j = 1, m; do i = 1, m
        y(i + (j - 1) * m) = sin(3 * xs(i)) * sin(3 * xs(j))
      end do; end do
    case (BRUSS)
      tend = 1
      do j = 1, m; do i = 1, m
        k = i + (j - 1) * m
        y(k) = 0.5_dp + xs(j)
        y(np + k) = 1 + 5 * xs(i)
      end do; end do
    case (FOOD)
      tend = 3
      do j = 1, m; do i = 1, m
        k = i + (j - 1) * m
        x = xs(i); yy = xs(j)
        y(k) = 10 - 5 * cos(pi * x) * cos(10 * pi * yy)
        y(np + k) = 17 + 5 * cos(10 * pi * x) * cos(pi * yy)
      end do; end do
    case (FOOD3)
      tend = 10
      do l = 1, m; do j = 1, m; do i = 1, m
        k = i + (j - 1) * m + (l - 1) * m * m
        x = xs(i); yy = xs(j); z = xs(l)
        y(k) = 500 + 250 * cos(pi * x) * cos(3 * pi * yy) * cos(10 * pi * z)
        y(np + k) = 200 + 150 * cos(10 * pi * x) * cos(pi * yy) * cos(3 * pi * z)
      end do; end do; end do
    case (VDPL)
      tend = 1000
      do j = 1, m; do i = 1, m
        k = i + (j - 1) * m
        x = xs(i); yy = xs(j)
        y(k) = 1 - cos(pi * x) * cos(2 * pi * yy)
        y(np + k) = 1 + 2 * cos(2 * pi * x) * cos(pi * yy)
      end do; end do
    case (DIFKON)
      tend = 1
      do j = 1, m; do i = 1, m
        x = xs(i); yy = xs(j)
        y(i + (j - 1) * m) = sin(pi * x) * sin(pi * yy) * (1 + 4 * x * yy)
      end do; end do
    case (DIFFU2)
      tend = 1
      do j = 1, m; do i = 1, m
        y(i + (j - 1) * m) = sin(pi * xs(i)) * sin(pi * xs(j))
      end do; end do
    end select
  end subroutine pinit

  ! ------------------------------------------------------------- rhs f(t,y)
  subroutine prhs(n, t, y, f)
    integer, intent(in) :: n
    real(dp), intent(in) :: t, y(n)
    real(dp), intent(out) :: f(n)
    integer :: k, i, j
    real(dp) :: u, v, eu
    select case (prob)
    case (HEAT)
      call lap_dir(y, f)
    case (NILIDI)
      call lap_dir(y, f)
      do k = 1, n
        eu = exp(y(k))
        f(k) = eu * f(k) + y(k) * (18 * eu - 1)
      end do
    case (BRUSS)
      call lap_neu(y(1), f(1))
      call lap_neu(y(np + 1), f(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        f(k) = 0.02_dp * f(k) + 1 + u * u * v - 4 * u
        f(np + k) = 0.02_dp * f(np + k) + 3 * u - u * u * v
      end do
    case (FOOD)
      call lap_neu(y(1), f(1))
      call lap_neu(y(np + 1), f(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        f(k) = 0.05_dp * f(k) + u * (1 - 0.1_dp * v)
        f(np + k) = f(np + k) + v * (-1000 + 100 * u)
      end do
    case (FOOD3)
      call lap_neu(y(1), f(1))
      call lap_neu(y(np + 1), f(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        f(k) = 0.05_dp * f(k) + u * (bfood3(k) - 1d6 * u - v)
        f(np + k) = f(np + k) + v * (bfood3(k) - (1d6 - 1) * u - 1d6 * v)
      end do
    case (VDPL)
      call lap_neu(y(1), f(1))
      call lap_neu(y(np + 1), f(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        f(k) = 0.05_dp * f(k) + v
        f(np + k) = f(np + k) + 100 * (1 - u * u) * v - u
      end do
    case (DIFKON)
      call conv_diff(y, f)
      call sources(t, wrk1, wrk2)
      f(1:np) = f(1:np) + wrk1
    case (DIFFU2)
      call lap_dir(y, f)
      call sources(t, wrk1, wrk2)
      f(1:np) = f(1:np) + wrk1
    end select
  contains
    ! eps*Lap(u) - y D_x u + x D_y u   (DIFKON, eps = 0.1)
    subroutine conv_diff(uu, out)
      real(dp), intent(in) :: uu(*)
      real(dp), intent(out) :: out(*)
      real(dp), parameter :: nu = 0.1_dp
      call lap_dir(uu, out)
      call dxy_dir(uu, wrk1, wrk2)
      do j = 1, m
        do i = 1, m
          k = i + (j - 1) * m
          out(k) = nu * out(k) - xs(j) * wrk1(k) + xs(i) * wrk2(k)
        end do
      end do
    end subroutine conv_diff
  end subroutine prhs

  ! ------------------------------------------------------- z = J(t,y) v
  subroutine pjv(n, t, y, vv, z)
    integer, intent(in) :: n
    real(dp), intent(in) :: t, y(n), vv(n)
    real(dp), intent(out) :: z(n)
    integer :: k, i, j
    real(dp) :: u, v, eu, nu
    select case (prob)
    case (HEAT, DIFFU2)
      call lap_dir(vv, z)
    case (NILIDI)
      call lap_dir(y, wrk1)
      call lap_dir(vv, z)
      do k = 1, n
        eu = exp(y(k))
        z(k) = eu * (z(k) + wrk1(k) * vv(k)) + (18 * eu * (1 + y(k)) - 1) * vv(k)
      end do
    case (BRUSS)
      call lap_neu(vv(1), z(1))
      call lap_neu(vv(np + 1), z(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        z(k) = 0.02_dp * z(k) + (2 * u * v - 4) * vv(k) + u * u * vv(np + k)
        z(np + k) = 0.02_dp * z(np + k) + (3 - 2 * u * v) * vv(k) - u * u * vv(np + k)
      end do
    case (FOOD)
      call lap_neu(vv(1), z(1))
      call lap_neu(vv(np + 1), z(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        z(k) = 0.05_dp * z(k) + (1 - 0.1_dp * v) * vv(k) - 0.1_dp * u * vv(np + k)
        z(np + k) = z(np + k) + 100 * v * vv(k) + (-1000 + 100 * u) * vv(np + k)
      end do
    case (FOOD3)
      call lap_neu(vv(1), z(1))
      call lap_neu(vv(np + 1), z(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        z(k) = 0.05_dp * z(k) + (bfood3(k) - 2d6 * u - v) * vv(k) - u * vv(np + k)
        z(np + k) = z(np + k) - (1d6 - 1) * v * vv(k) &
                    + (bfood3(k) - (1d6 - 1) * u - 2d6 * v) * vv(np + k)
      end do
    case (VDPL)
      call lap_neu(vv(1), z(1))
      call lap_neu(vv(np + 1), z(np + 1))
      do k = 1, np
        u = y(k); v = y(np + k)
        z(k) = 0.05_dp * z(k) + vv(np + k)
        z(np + k) = z(np + k) + (-200 * u * v - 1) * vv(k) + 100 * (1 - u * u) * vv(np + k)
      end do
    case (DIFKON)
      nu = 0.1_dp
      call lap_dir(vv, z)
      call dxy_dir(vv, wrk1, wrk2)
      do j = 1, m
        do i = 1, m
          k = i + (j - 1) * m
          z(k) = nu * z(k) - xs(j) * wrk1(k) + xs(i) * wrk2(k)
        end do
      end do
    end select
  end subroutine pjv

  ! ---------------------------------------------------- ft = df/dt (t,y)
  subroutine pft(n, t, y, ft)
    integer, intent(in) :: n
    real(dp), intent(in) :: t, y(n)
    real(dp), intent(out) :: ft(n)
    select case (prob)
    case (DIFKON, DIFFU2)
      call sources(t, wrk1, wrk2)
      ft(1:np) = wrk2(1:np)
    case default
      ft = 0
    end select
  end subroutine pft

end module problems
