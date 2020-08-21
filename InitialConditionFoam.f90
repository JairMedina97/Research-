subroutine InitialConditionFoam
    use VarModFoamMech
    implicit none

   !$OMP PARALLEL
   !$OMP DO
    do i=1,iMax
      Sw(i,1)=Swini
      Sg(i,1) = 1.d0 - Sw(i,1)
    end do
   !$OMP END DO
   !$OMP DO
    do i=1,iMax
      nf(i,1)=nfini
    end do
   !$OMP END DO
   !$OMP DO
    do i=1,iMax
      pc(i,0) = rw*(Kp(i)*Ca)**(-1)*(Sw(i,1)/(1-Swc-Sgr))**(-n2Pc)
    end do
  !$OMP END DO
  !$OMP DO
    do i=1,iMax
      Krg0(i,0) = Krg00*((1-Sw(i,1)-Sgr)/(1-Swc-Sgr))**nkrg

      if(Krg0(i,0) /= Krg0(i,0)) then
         write(*,*)'Krg0.ciniciales aqui',i,Krg0(i,0),Sw(i,1)
         stop
      end if
    end do
  !$OMP END DO
  !$OMP DO
    do i=1,iMax
      if ( Sw(i,1)<Swc ) Sw(i,1)=Swc+eps
      Krw(i,0) = Krw0*((Sw(i,1)-Swc)/(1-Swc-Sgr))**nkrw  ! water permeability
    end do
  !$OMP END DO
  !$OMP DO
    do i=1,iMax
     lamdaW(i,1) = Krw(i,0)/Nuw0
    end do
  !$OMP END DO
  !$OMP END PARALLEL

    NuG(1:iMax,1) = Br !0.05d0  !Pa.s Programar una mejor elecci�n, din�mica
    NuGAverage=sum(Nug(1:iMax,1))/iMax
    NuGAverageOld=NuGAverage

     Uw(1,0) =  FwJ*rw/x(1)
     Ug(1,0) =  rw/x(1)-Uw(1,0)
     itera: do iter=1,iterMax
        !$OMP PARALLEL
        !$OMP DO
        do i=1,iMax
          lamdaG(i,1) = Krg0(i,0)/NuG(i,1)
            !write(*,*)'lamdaG',Krg0(i,0),NuG(i,1)

          if(NuG(i,1)/=NuG(i,1) .or. lamdaG(i,1)/=lamdaG(i,1) ) then
            write(*,*)'lamdaG.ciniciales aqui',i,Krg0(i,0),NuG(i,1)
            stop
          end if
        end do
        !$OMP END DO
        !$OMP DO
        do i=1,iMax
          Fw(i,0) = lamdaW(i,1) / (lamdaW(i,1) + lamdaG(i,1))
        end do
        !$OMP END DO


       !$OMP DO
       do i=2,iMax-1
         if(iter==1)then
            Uw(i,0) =  Fw(i,0)*rw/x(i) + &
                       lamdaW(i,1)*(pc(i,0)-x(i-1)*pc(i-1,0)/x(i))/dx- &
                       lamdaW(i,1)*Fw(i,0)*(pc(i,0)-x(i-1)*pc(i-1,0)/x(i))/dx
         else
            Uw(i,0) =  Uw(i,0)*relax+(1-relax)* (Fw(i,0)*rw/x(i) + &
                       lamdaW(i,1)*(pc(i,0)-x(i-1)*pc(i-1,0)/x(i))/dx- &
                       lamdaW(i,1)*Fw(i,0)*(pc(i,0)-x(i-1)*pc(i-1,0)/x(i))/dx)
              if(Fw(i,0)/=Fw(i,0) .or. lamdaG(i,1)/=lamdaG(i,1) ) then
                 write(*,*)'velocity aqui',i,Fw(i,0),lamdaG(i,1),pc(i,0)
                 stop
              end if
         end if
       end do
       !$OMP END DO

       !$OMP DO
       do i=2,iMax-1
         if(iter==1)then
           Ug(i,0) = rw/x(i)-Uw(i,0)
         else
           Ug(i,0) = Ug(i,0)*relax + (1-relax)*(rw/x(i)-Uw(i,0))
         end if
       end do
       !$OMP END DO

       if(iter==1)then    !The pressure-capilar in x=L is asumed zero
          Uw(iMax,0) = Fw(iMax,0)*rw + &
                       lamdaW(iMax,1)*(-x(iMax-1)*pc(iMax-1,0))/dx- &
                       lamdaW(iMax,1)*Fw(iMax,0)*(-x(iMax-1)*pc(iMax-1,0))/dx
          if ( Uw(iMax,0) < 0.d0 ) Uw(iMax,0) = 0.d0
          Ug(iMax,0) = rw/x(iMax)-Uw(iMax,0)

       else
          Uw(iMax,0) =Uw(iMax,0)*relax+(1-relax)*( Fw(iMax,0)*rw + &
                      lamdaW(iMax,1)*(-x(iMax-1)*pc(iMax-1,0))/dx- &
                      lamdaW(iMax,1)*Fw(iMax,0)*(-x(iMax-1)*pc(iMax-1,0))/dx)
          if ( Uw(iMax,0) < 0.d0 ) Uw(iMax,0) = 0.d0
          Ug(iMax,0) =Ug(iMax,0)*relax+(1-relax)*( rw/x(iMax)-Uw(iMax,0))
       end if

       P(1,1) = PwJ - (dx*UwJ)/(2*lamdaW(1,1))
       !$OMP DO
       do i=1,iMax-1
          P(i+1,1) = P(i,1) - (dx*Uw(i,0))/lamdaW(i,1)
       end do
       !$OMP END

       !$OMP DO
       do i=1,iMax
         NuG(i,1) = NuG(i,1)*relax + (1-relax)*(Br + nf(i,1)/(Ug(i,0)**(1.d0/3.d0) + Br2 ) )
         if (NuG(i,1)<Br) NuG(i,1)=Br

         if(NuG(i,1)/=NuG(i,1) .or. nf(i,1)/=nf(i,1) ) then
         write(*,*)'NuG.ciniciales aqui',i,nf(i,1),NuG(i,1),Ug(i,0),Uw(i,0)
         stop
         end if
       end do
       !$OMP END
       !$OMP END PARALLEL
       NuGAverage=sum(Nug(1:iMax,1))/iMax
       Conver=abs((NuGAverage-NuGAverageOld)/NuGAverageOld)
       NuGAverageOld=NuGAverage
       if(Conver < ConverMax) exit itera

       if(iter==iterMax) then
         write(*,*)'No Convergence in initialCondition'
         stop
       end if

    end do itera

    a(1)=rw
    iAvaOld=0

    rg(1,0) = Sw(1,1)*( ( PwJ-0.5d0*(P(1,1)+P(2,1)) ) /dx )**m
    !$OMP PARALLEL
    !$OMP DO
    do i=2,iMax-1
      rg(i,0) = Sw(i,1)*((P(i-1,1)-P(i+1,1))/(2*dx))**m
    end do
    !$OMP END DO
    rg(iMax,0) = rg(iMax-1,1)
    rg(1:iMax,1) = rg(1:iMax,0)

    !$OMP DO
    do i=1,iMax
       rc(i,0) = 0.d0
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    mark=0

    nfAverageOld = 0.d0

end subroutine InitialConditionFoam
