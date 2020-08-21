subroutine KrsFoamMech
  use VarModFoamMech
  implicit none
!$OMP PARALLEL

    !$OMP DO
    do i=1,iMax
      if(Sw(i,k+1)>1) Sw(i,k+1)=0.999999d0
      Krg0(i,1) = Krg00*((1-Sw(i,k+1)-Sgr)/(1-Swc-Sgr))**nkrg  ! gas permeability whitout foam
      if ( 1-Sw(i,k+1)<Sgr ) then
        write(*,*)'en Krg0',i,Krg0(i,1),Sw(i,k+1)
        stop
      end if
    end do
    !$OMP END DO

    !$OMP DO
    do i=1,iMax
      if ( Sw(i,k+1)<=Swc ) Sw(i,k+1)=Swc+eps
      Krw(i,1) = Krw0*((Sw(i,k+1)-Swc)/(1-Swc-Sgr))**nkrw  ! water permeability
    end do
    !$OMP END DO
!$OMP END PARALLEL
end subroutine KrsFoamMech
