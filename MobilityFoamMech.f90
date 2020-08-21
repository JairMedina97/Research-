subroutine MobilityFoamMech
  use VarModFoamMech
  implicit none
!$OMP PARALLEL
    !$OMP DO
    do i=1,iMax
!      lamdaG(i,k+1) = Krg0(i,1)/NuG(i,k+1)
!      lamdaW(i,k+1) = Krw(i,1)/Nuw0
      lamdaG(i,k+1) = relax*lamdaG(i,k+1) + (1-relax)*Krg0(i,1)/NuG(i,k+1)
      lamdaW(i,k+1) = relax*lamdaW(i,k+1) + (1-relax)*Krw(i,1)/Nuw0
    end do
    !$OMP END DO
!$OMP END PARALLEL
end subroutine MobilityFoamMech
