subroutine RheologyFoamMech
  use VarModFoamMech
  implicit none

 !$OMP PARALLEL
    !$OMP DO
    do i=1,iMax
      NuG(i,k+1) = NuG(i,k+1)*relax + (1-relax)*(Br + nf(i,k+1)/(Ug(i,1)**(1.d0/3.d0) + Br2) )
      if (NuG(i,k+1)<Br) NuG(i,k+1)=Br
      if (NuG(i,k+1)>Br+1.d0/Br2) NuG(i,k+1)=Br+1.d0/Br2
      if(NuG(i,k+1)/=NuG(i,k+1)) then
      !NuG(i,k+1)=1.d0
      write(*,*)'Nug aqui',iter,i,NuG(i,k+1),nf(i,k+1),Ug(i,1),Uw(i,1)
      stop
     end if
    end do
    !$OMP END DO
 !$OMP END PARALLEL

end subroutine RheologyFoamMech
