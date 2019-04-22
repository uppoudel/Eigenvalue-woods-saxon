program energy
implicit none
  integer(8) n, i
  real(8) beta, lambda, V, a, En
  real(8), parameter :: hbar=1.05457180d-34
  real(8), parameter :: mu=1.62747d-27
  real(8), parameter :: Ese=1.0d+6
  real(8), parameter :: Lse=1.0d-15
  open(1,file='in.dat',action='read')
   open(2,file='eigval.dat',action='write')
  read(1,*) V
  read(1,*) a
   close(1)
  V=V*Ese
  a=a*Lse
  beta=(mu*a**2)/hbar**2
  lambda=beta/2
  !n=floor(1.0d0/kappa-0.50d0)     ! determine the number of the energy states
  do i=0, n
    En=-lambda*[(beta*V/(i+1.0d0))**2+((i+1.0d0)/2)**2+beta*V]  ! determine the eigenvalue
    En=En/Ese;  ! energy unit conversion, back to MeV
  write(2,*) En
  end do
close(2)
end program energy
