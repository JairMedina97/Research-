subroutine PcFoamMech
  use VarModFoamMech
  implicit none
!$OMP PARALLEL
    !$OMP DO
    do i=1,iMax
      Pc(i,1) = rw*(Kp(i)*Ca)**(-1)*(Sw(i,k+1)/(1-Swc-Sgr))**(-n2Pc)
    end do
    !$OMP END DO
!$OMP END PARALLEL
end subroutine PcFoamMech
