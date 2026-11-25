! Experiment driver for rowmap.F (paper experiments, docs/problems.md).
!
! Usage:
!   driver run  <problem> <m> <tol> <method> <maxkry> <sft> <reffile|-> <out.toml>
!   driver dump <problem> <m> <out.txt>
!
! sft (paper §4): 2 = nonautonomous form (ifcn=1, exact f_t)
!                 3 = augmented autonomous system y_{n+1}' = 1
!                 4 = like 3, plus y_{n+1} := t after each step   (patched build only)
!                 5 = like 4, plus k_{i,n+1} := 1                  (patched build only)
! For autonomous problems sft is ignored (recorded as 0).
module drvstate
  implicit none
  integer, parameter :: dp = kind(1d0)
  logical, save :: aug = .false.
  integer, save :: nn = 0          ! dimension of the original ODE
  real(dp), allocatable, save :: tmp(:)
end module drvstate

subroutine fcn(n, t, u, udot, rpar, ipar)
  use drvstate
  use problems
  implicit none
  integer :: n, ipar(*)
  real(dp) :: t, u(n), udot(n), rpar(*)
  if (aug) then
    call prhs(nn, u(n), u, udot)
    udot(n) = 1
  else
    call prhs(n, t, u, udot)
  end if
end subroutine fcn

subroutine jacv(n, t, u, v, z, rpar, ipar)
  use drvstate
  use problems
  implicit none
  integer :: n, ipar(*)
  real(dp) :: t, u(n), v(n), z(n), rpar(*)
  if (aug) then
    call pjv(nn, u(n), u, v, z)
    call pft(nn, u(n), u, tmp)
    z(1:nn) = z(1:nn) + v(n) * tmp
    z(n) = 0
  else
    call pjv(n, t, u, v, z)
  end if
end subroutine jacv

subroutine fdt(n, t, u, ft, rpar, ipar)
  use drvstate
  use problems
  implicit none
  integer :: n, ipar(*)
  real(dp) :: t, u(n), ft(n), rpar(*)
  call pft(n, t, u, ft)
end subroutine fdt

subroutine solout(n, told, tnew, uold, unew, fold, fnew, ucon, intr, rpar, ipar)
  implicit none
  integer :: n, intr, ipar(*)
  double precision :: told, tnew, uold(n), unew(n), fold(n), fnew(n), ucon(n), rpar(*)
end subroutine solout

program driver
  use drvstate
  use problems
  implicit none
  external fcn, jacv, fdt, solout
  character(len=512) :: mode, pname, arg, reffile, outfile
  integer :: id, m, method, maxkry, sft, n, naug, ifcn, idid, lwork, liwork, i, u
  integer :: ipar(1), isft, nrep
  real(dp), parameter :: tmin = 0.5_dp   ! minimal accumulated CPU time for timing
  real(dp) :: tol, t, t0, tend, hs, rpar(1), rtol(1), atol(1), c0, c1, err
  real(dp), allocatable :: y(:), work(:), f(:), z(:), v(:), ft(:)
  integer, allocatable :: iwork(:)
  character(len=6), parameter :: mnames(6) = [character(len=6) :: &
       'shamp', 'grk4t', 'velds', 'veldd', 'lstab', 'grk4a']
  common /sftpatch/ isft

  call get_command_argument(1, mode)
  call get_command_argument(2, pname)
  id = id_of(trim(pname))
  if (id == 0) stop 'unknown problem'
  call get_command_argument(3, arg); read (arg, *) m
  call setup(id, m, n)
  nn = n

  if (trim(mode) == 'dump') then
    call get_command_argument(4, outfile)
    allocate(y(n), f(n), z(n), v(n), ft(n))
    call pinit(n, y, t0, tend)
    do i = 1, n
      v(i) = sin(dble(i))
    end do
    call prhs(n, t0, y, f)
    call pjv(n, t0, y, v, z)
    call pft(n, t0, y, ft)
    open (newunit=u, file=trim(outfile), status='replace')
    write (u, '(a,a)') '# problem = ', trim(pname)
    write (u, '(a,i0,a,i0)') '# m = ', m, '  n = ', n
    write (u, '(a,es24.16e3,a,es24.16e3)') '# t0 = ', t0, '  tend = ', tend
    write (u, '(a)') '# columns: y0  f(t0,y0)  J(t0,y0)*v (v_k = sin k)  f_t(t0,y0)'
    do i = 1, n
      write (u, '(4es25.16e3)') y(i), f(i), z(i), ft(i)
    end do
    close (u)
    stop
  end if

  if (trim(mode) /= 'run') stop 'mode must be run or dump'
  call get_command_argument(4, arg); read (arg, *) tol
  call get_command_argument(5, arg); read (arg, *) method
  call get_command_argument(6, arg); read (arg, *) maxkry
  call get_command_argument(7, arg); read (arg, *) sft
  call get_command_argument(8, reffile)
  call get_command_argument(9, outfile)

  if (prob_auton()) sft = 0
  aug = sft >= 3
  isft = 0
  if (sft >= 4) isft = sft
  naug = n
  if (aug) naug = n + 1
  allocate(y(naug), tmp(n))
  call pinit(n, y, t0, tend)
  if (aug) y(naug) = t0
  t = t0

  ifcn = 1
  if (prob_auton() .or. aug) ifcn = 0
  lwork = 10 + naug * (maxkry + 11) + maxkry * (maxkry + 4) + 100
  liwork = maxkry + 40
  allocate(work(lwork), iwork(liwork))
  work = 0
  iwork = 0
  iwork(1) = 1000000       ! max steps (paper default 10000 too small for VDPL)
  iwork(2) = method
  iwork(3) = maxkry
  rtol(1) = tol
  atol(1) = tol
  hs = 0
  rpar = 0
  ipar = 0

  ! Repeat the whole integration until >= tmin CPU seconds have accumulated
  ! and report the mean: single runs take only milliseconds and timer noise
  ! would dominate the work-precision diagrams. Every repetition is
  ! identical (fresh initial values, work arrays and options).
  nrep = 0
  call cpu_time(c0)
  do
    call pinit(n, y, t0, tend)
    if (aug) y(naug) = t0
    t = t0
    hs = 0
    work = 0
    iwork = 0
    iwork(1) = 1000000       ! max steps (paper default 10000 too small for VDPL)
    iwork(2) = method
    iwork(3) = maxkry
    call rowmap(naug, fcn, ifcn, t, y, tend, hs, rtol, atol, 0, &
                jacv, 1, fdt, 1, solout, 0, work, lwork, iwork, liwork, rpar, ipar, idid)
    nrep = nrep + 1
    call cpu_time(c1)
    if (c1 - c0 >= tmin .or. idid < 0) exit
  end do
  c1 = c0 + (c1 - c0) / nrep

  err = ieee_nan()
  if (idid > 0 .and. trim(reffile) /= '-') err = errref(trim(reffile), y(1:n))

  open (newunit=u, file=trim(outfile), status='replace')
  write (u, '(a)') '[run]'
  write (u, '(a)') 'producer = "rowmap-f77"'
  write (u, '(a,a,a)') 'problem = "', trim(pname), '"'
  write (u, '(a,i0)') 'm = ', m
  write (u, '(a,i0)') 'n = ', n
  write (u, '(a,a,a)') 'method = "', trim(mnames(method)), '"'
  write (u, '(a,es24.16e3)') 'rtol = ', tol
  write (u, '(a,es24.16e3)') 'atol = ', tol
  write (u, '(a,i0)') 'max_krylov = ', maxkry
  write (u, '(a,i0)') 'sft = ', sft
  write (u, '(a,es24.16e3)') 't_end = ', tend
  write (u, '(a,i0)') 'idid = ', idid
  write (u, '(a)') ''
  write (u, '(a)') '[stats]'
  if (err /= err) then
    write (u, '(a)') 'err = nan'
  else
    write (u, '(a,es24.16e3)') 'err = ', err
  end if
  write (u, '(a,i0)') 'nstep = ', iwork(5)
  write (u, '(a,i0)') 'nfeval = ', iwork(7)
  write (u, '(a,i0)') 'njv = ', iwork(8)        ! diagnostic: Jacobian-vector products
  write (u, '(a,es24.16e3)') 'cpu_time = ', c1 - c0   ! mean over nrep repetitions
  close (u)

contains

  real(dp) function ieee_nan()
    real(dp) :: zero
    zero = 0
    ieee_nan = zero / zero
  end function ieee_nan

  ! Paper ERR against a reference file (lines starting with '#' are skipped).
  real(dp) function errref(fname, yy)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: yy(:)
    character(len=256) :: line
    integer :: uu, k, ios
    real(dp) :: r, s
    open (newunit=uu, file=fname, status='old', action='read')
    k = 0
    s = 0
    do
      read (uu, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      k = k + 1
      if (k > size(yy)) stop 'reference file too long'
      read (line, *) r
      s = s + ((yy(k) - r) / (1 + abs(r)))**2
    end do
    close (uu)
    if (k /= size(yy)) stop 'reference file length mismatch'
    errref = sqrt(s / size(yy))
  end function errref

end program driver
