! eps-sweep near an exact Krylov breakdown, one ROWMAP step per eps
! (AGENTS.md §5.3 item 4, Fortran comparison).
!
! Usage: sweep <input> <output>
!   input  : results/sweep_input.txt (written by `cargo test sweep_export -- --ignored`)
!   output : one block of n values u(h) per eps, then a line "# jv = ..."
!
! Problem y' = (D + eps G) y, one step of size h from u0. Loose tolerances
! (rtol = atol = 1) make the step accepted and the stage-1 residual test stop
! at the minimal dimension 4, comparable to the Rust run with kappa1 = 4.
module sweepdata
  implicit none
  integer, parameter :: dp = kind(1d0)
  integer, save :: n
  real(dp), allocatable, save :: d(:), g(:, :)
  real(dp), save :: eps
end module sweepdata

subroutine sfcn(nn, t, u, f, rpar, ipar)
  use sweepdata
  implicit none
  integer :: nn, ipar(*)
  real(dp) :: t, u(nn), f(nn), rpar(*)
  f = d * u + eps * matmul(g, u)
end subroutine sfcn

subroutine sjacv(nn, t, u, v, z, rpar, ipar)
  use sweepdata
  implicit none
  integer :: nn, ipar(*)
  real(dp) :: t, u(nn), v(nn), z(nn), rpar(*)
  z = d * v + eps * matmul(g, v)
end subroutine sjacv

subroutine sdummy(nn, t, u, ft, rpar, ipar)
  implicit none
  integer :: nn, ipar(*)
  double precision :: t, u(nn), ft(nn), rpar(*)
end subroutine sdummy

subroutine ssolout(nn, told, tnew, uold, unew, fold, fnew, ucon, intr, rpar, ipar)
  implicit none
  integer :: nn, intr, ipar(*)
  double precision :: told, tnew, uold(nn), unew(nn), fold(nn), fnew(nn), ucon(nn), rpar(*)
end subroutine ssolout

program sweep
  use sweepdata
  implicit none
  external sfcn, sjacv, sdummy, ssolout
  character(len=512) :: fin, fout, line
  integer :: neps, i, j, k, ui, uo, lwork, liwork, idid, ipar(1), mx
  real(dp) :: h, t, hs, rtol(1), atol(1), rpar(1)
  real(dp), allocatable :: epsv(:), u0(:), u(:), work(:)
  integer, allocatable :: iwork(:)

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  open (newunit=ui, file=trim(fin), status='old', action='read')
  read (ui, '(a)') line            ! comment
  read (ui, *) n, neps, h
  allocate(epsv(neps), d(n), g(n, n), u0(n), u(n))
  do k = 1, neps
    read (ui, *) epsv(k)
  end do
  do i = 1, n
    read (ui, *) d(i)
  end do
  do i = 1, n                      ! G stored row-major
    do j = 1, n
      read (ui, *) g(i, j)
    end do
  end do
  do i = 1, n
    read (ui, *) u0(i)
  end do
  close (ui)

  mx = 70
  lwork = 10 + n * (mx + 11) + mx * (mx + 4) + 100
  liwork = mx + 40
  allocate(work(lwork), iwork(liwork))
  open (newunit=uo, file=trim(fout), status='replace')
  write (uo, '(a)') '# rowmap.F, one step per eps: n values u(h), then "# jv idid"'
  do k = 1, neps
    eps = epsv(k)
    u = u0
    t = 0
    hs = h
    work = 0
    iwork = 0
    iwork(3) = mx
    rtol(1) = 1
    atol(1) = 1
    call rowmap(n, sfcn, 0, t, u, h, hs, rtol, atol, 0, sjacv, 1, sdummy, 0, &
                ssolout, 0, work, lwork, iwork, liwork, rpar, ipar, idid)
    do i = 1, n
      write (uo, '(es25.16e3)') u(i)
    end do
    write (uo, '(a,i0,a,i0,a,i0)') '# jv ', iwork(8), ' idid ', idid, ' steps ', iwork(5)
  end do
  close (uo)
end program sweep
