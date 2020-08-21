subroutine PoroPerFoamMech
  use VarModFoamMech
  implicit none
  integer(kind=4)              :: irand,ir,clock,ic
  integer(kind=4), allocatable :: RandSeed(:)
  real(kind=8),    allocatable :: PoroAux(:)
  real(kind=8)                 :: RanNumGen,U1,U2,val1

      Pi = 4*atan(1.d0)
      call random_seed()
      call random_seed(size = irand)
      allocate(RandSeed(irand),PoroAux(iMax) , stat=status)
      call system_clock(count=clock)
      RandSeed = clock + 37 * (/ (ic-1, ic = 1, irand) /)
      call random_seed(put = RandSeed)
      !call random_seed(get = RandSeed(1:irand))

      do i=1,iMax
            call random_number(RanNumGen);    U1=RanNumGen
            call random_number(RanNumGen);    U2=RanNumGen
            val1 = dsqrt(-2*log(U1)) * dcos(2*pi*U2)
            PoroAux(i) = abs(val1*DesStanPoro + MediaPoro)/MediaPoro
      end do
      Poro(1) = (PoroAux(1) + PoroAux(2)+PoroAux(3) ) /3.d0
      Kp(1) = Poro(1)**(nPer+1)*(1 - MediaPoro)**nPEr/(1 - MediaPoro*Poro(1))**nPer
      do i=2,iMax-1
         Poro(i) = (PoroAux(i-1) + PoroAux(i) + PoroAux(i+1) ) /3.d0
         Kp(i) = Poro(i)**(nPer+1)*(1 - MediaPoro)**nPEr/(1 - MediaPoro*Poro(i))**nPer
      end do
      Poro(iMax) = (PoroAux(iMax-2) + PoroAux(iMax-1)+PoroAux(iMax) ) /3.d0
      Kp(iMax) = Poro(iMax)**(nPer+1)*(1 - MediaPoro)**nPEr/(1 - MediaPoro*Poro(iMax))**nPer


      open(unit=88, file='poro.dat')
        do i=1,iMax
           write(88,*)x(i),Poro(i),Kp(i)
        end do
      close(unit=88)


end subroutine PoroPerFoamMech
