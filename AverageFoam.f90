subroutine AverageFoam
  use VarModFoamMech
  implicit none

!$OMP PARALLEL

   AveSw = 0.d0
   !$OMP DO
   do i=2,iMax-1
     AveSw = AveSw + Sw(i,k+1)*x(i)
   end do
   AveSw = AveSw + (Sw(1,k+1)*x(1)+Sw(iMax,k+1)*x(iMax))/2.d0
   AveSw = AveSw*2*dx/(1.d0-rw**2)
   !$OMP END DO

   AveNf = 0.d0
   !$OMP DO
   do i=2,iMax-1
     AveNf = AveNf + Nf(i,k+1)*x(i)
   end do
   AveNf = AveNf + (Nf(1,k+1)*x(1)+Nf(iMax,k+1)*x(iMax))/2.d0
   AveNf = AveNf*2*dx/(1.d0-rw**2)
   !$OMP END DO

   AveUg = 0.d0
   !$OMP DO
   do i=2,iMax-1
     AveUg = AveUg + Ug(i,1)*x(i)
   end do
   AveUg = AveUg +( Ug(1,1)*x(1)+Ug(iMax,1)*x(iMax))/2.d0
   AveUg = AveUg*2*dx/(1.d0-rw**2)
   !$OMP END DO

   AveUw = 0.d0
   !$OMP DO
   do i=2,iMax-1
     AveUw = AveUw + Uw(i,1)*x(i)
   end do
   AveUw = AveUw +( Uw(1,1)*x(1)+Uw(iMax,1)*x(iMax))/2.d0
   AveUw = AveUw*2*dx/(1.d0-rw**2)
   !$OMP END DO

   Averg = 0.d0
   !$OMP DO
   do i=2,iMax-1
     Averg = Averg + rg(i,1)*x(i)
   end do
   Averg = Averg +( rg(1,1)*x(1)+rg(iMax,1)*x(iMax))/2.d0
   Averg = Averg*2*dx/(1.d0-rw**2)
   !$OMP END DO

   Averc = 0.d0
   !$OMP DO
   do i=1,iMax
     Averc = Averc + rc(i,1)*x(i)
   end do
   Averc = Averc +( rc(1,1)*x(1)+rc(iMax,1)*x(iMax))/2.d0
   Averc = Averc*2*dx/(1.d0-rw**2)
   !$OMP END DO

!$OMP END PARALLEL
end subroutine AverageFoam
