subroutine TextureFoamMech
    use VarModFoamMech
    implicit none
    real(kind=8)    :: step
    step=2000.d0

      Da = Dam*Poro(1)
      Ga = Gam*Kp(1)
      nf(1,k+1)=  Sg(1,k)*nf(1,k)/Sg(1,k+1) - &
                  dt*2*nf(1,k)*Ug(1,0)/(rw*dx*Sg(1,k+1)) + &
                  dt*Da*Sg(1,k)*(Ga*rg(1,0)-rc(1,0))/Sg(1,k+1)
      nf(1,k+1)= nf(1,k+1) /( 1+exp(-step*(Sg(1,k+1) - Sgcr)) )
      !if(Sg(1,k+1) <= 0.02d0) nf(1,k+1) = 0.d0

    do i=2,iMax-1
      Da = Dam*Poro(i)
      Ga = Gam*Kp(i)
      nf(i,k+1)=  Sg(i,k)*nf(i,k)/Sg(i,k+1) + &
                  dt*x(i-1)*nf(i-1,k)*Ug(i-1,0)/(x(i)*rw*dx*Sg(i,k+1)) - &
                  dt*nf(i,k)*Ug(i,0)/(rw*dx*Sg(i,k+1)) + &
                  dt*Da*Sg(i,k)*(Ga*rg(i,0)-rc(i,0))/Sg(i,k+1)
      nf(i,k+1)= nf(i,k+1) /( 1+exp(-step*(Sg(i,k+1) - Sgcr)) )
      !if(Sg(i,k+1) <= 0.02d0) nf(i,k+1) = 0.d0

      if(nf(i,k+1)/=nf(i,k+1)) then
        write(*,*)'nf sin Crank aqui',i,nf(i,     k+1),Ug(i,0),rg(i,0),rc(i,0),Sg(i,k)
        stop
      end if
    end do

    if(Sg(iMax,k+1) <= 0.02d0) nf(iMax-1,k+1) = 0.d0
    nf(iMax,k+1)= nf(iMax-1,k+1)

end subroutine TextureFoamMech

