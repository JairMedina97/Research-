subroutine TextureFoamMechIter

    use VarModFoamMech
    implicit none
    real(kind=8)    :: step
    step=2000.d0
      Da = Dam*Poro(1)
      Ga = Gam*Kp(1)
    div=Sg(1,k+1)+dt*(1-gama)*Ug(1,1)/(rw*dx)
     nf(1,k+1)=  (Sg(1,k)*nf(1,k))/div - &
                  dt*2*gama**nf(1,k)*Ug(1,0)/(div*rw*dx) + &
                  dt*gama*Da*Sg(1,k)*(Ga*rg(1,0)-rc(1,0))/div +&
                  dt*(1-gama)*Da*Sg(1,k+1)*(Ga*rg(1,1)-rc(1,1))/div
     nf(1,k+1)= nf(1,k+1) /( 1+exp(-step*(Sg(1,k+1) - Sgcr)) )

     !if(Sg(1,k+1) <= 0.02d0) nf(1,k+1) = 0.d0
     if (nf(1,k+1)>1.d0) nf(1,k+1)=1.d0
     if (nf(1,k+1)<0.d0) nf(1,k+1)=0.d0

   do i=2,iMax-1
         Da = Dam*Poro(i)
         Ga = Gam*Kp(i)
        div=Sg(i,k+1)+dt*(1-gama)*Ug(i,1)/(rw*dx)

        nf(i,k+1)=  Sg(i,k)*nf(i,k)/div + &
                    dt*gama*x(i-1)*nf(i-1,k)*Ug(i-1,0)/(div*rw*dx*x(i)) - &
                    dt*gama*nf(i,k)*Ug(i,0)/(div*rw*dx) + &
                    dt*(1-gama)*x(i-1)*nf(i-1,k+1)*Ug(i-1,1)/(div*rw*dx*x(i)) + &
                    dt*gama*Da*Sg(i,k)*(Ga*rg(i,0)-rc(i,0))/div + &
                    dt*(1-gama)*Da*Sg(i,k+1)*(Ga*rg(i,1)-rc(i,1))/div
         nf(i,k+1)= nf(i,k+1) /( 1+exp(-step*(Sg(i,k+1) - Sgcr)) )
         !if(Sg(i,k+1) <= 0.02d0) nf(i,k+1) = 0.d0
         if (nf(i,k+1)>1.d0) nf(i,k+1)=1.d0
         if (nf(i,k+1)<0.d0) nf(i,k+1)=0.d0
         if(nf(i,k+1)/=nf(i,k+1)) then !Revisa que si es NaN
          write(*,*)'nfiter',i,1.d0/( 1+exp(-step*(Sg(i-1,k+1) - Sgcr)) ),nf(i,k+1)
          stop
         end if
     end do

     nf(iMax,k+1)= nf(iMax-1,k+1)

end subroutine TextureFoamMechIter
