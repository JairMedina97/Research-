subroutine PromMov
    use VarModFoamMech
    implicit none

    !$OMP PARALLEL
    !$OMP DO
    do i=2,iMax-1
       !nfaux(i)=nf(i-1,k+1)*0.25+ nf(i,k+1)*0.5+nf(i+1,k+1)*0.25
!       Swaux(i)=Sw(i-1,k+1)*0.25+Sw(i,k+1)*0.5+Sw(i+1,k+1)*0.25
!       rcaux(i)=rc(i-1,1)*0.25+rc(i,1)*0.5+rc(i+1,1)*0.25
!       Ugaux(i)=Ug(i-1,1)*0.25+Ug(i,1)*0.5+Ug(i+1,1)*0.25
!       rgaux(i)=rg(i-1,1)*0.25+rg(i,1)*0.5+rg(i+1,1)*0.25
        Paux(i)=P(i-1,k+1)*0.25+P(i,k+1)*0.5+P(i+1,k+1)*0.25
       !NuGaux(i)=NuG(i-1,k+1)*0.25+NuG(i,k+1)*0.5+NuG(i+1,k+1)*0.25
!       Uwaux(i)=Uw(i-1,1)*0.25+Uw(i,1)*0.5+Uw(i+1,1)*0.25
    end do
   !$OMP END DO
    !$OMP DO
    do i=2,iMax-1
       !nf(i,k+1)=nfaux(i)
!       Sw(i,k+1)=Swaux(i)
        P(i,k+1)=Paux(i)
!       rg(i,1)=rgaux(i)
       !NuG(i,k+1)=NuGaux(i)
!       rc(i,1)=rcaux(i)
!       Ug(i,1)=Ugaux(i)
!       Uw(i,1)=Uwaux(i)
    end do
    !$OMP END DO
!$OMP END PARALLEL

end subroutine PromMov
