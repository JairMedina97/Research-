program MainModFoamMech
  use VarModFoamMech
  implicit none
  !$OMP  call OMP_SET_DYNAMIC(.true.)
     open (unit=11,file='Lamella.dat')
     call InputModFoamMech
     call GridModFoamMech
     call PoroPerFoamMech
     call initialConditionFoam
     write(*,*)'aqui voy'


 open(unit=9,file='slope.dat')

     icont=1
     time: do k=1,kMax-1
        dt=t(k+1)-t(k)
        call SatModFoamMech
        call TextureFoamMech
        call KrsFoamMech
        call PcFoamMech
        SwAverage=sum(Sw(1:iMax,k+1))/iMax
        SwAverageOld=SwAverage
        nfAverage=sum(nf(1:iMax,k+1))/iMax
        nfAverageOld=nfAverage
        NuG(1:iMax,k+1)=NuG(1:iMax,k)
        NuGAverage=sum(Nug(1:iMax,k+1))/iMax
        NuGAverageOld=NuGAverage

        itera: do iter=1,iterMax

           call MobilityFoamMech
           call FractionalFoamMech
           call VelocityFoamMech
           call PressureFoamMech
           call LamellaCreationCoalescenceFoamMech
           call RheologyFoamMech
           call SatModFoamMechIter
           call TextureFoamMechIter
           call ConverFoamMech

           if(Conver < ConverMax) exit itera
           if(iter==iterMax) then
              write(*,*)'Not Convergence in time =',t(k),convernug, conversw, convernf,iter
              call OutputFoam
              stop
           end if

        end do itera
        call PromMov
        call AverageFoam
        call Slope

         write(9,*)t(k+1),a(k+1)

        if(a(k+1)/=a(k))icont=0

        if(a(k+1)==a(k))then
           icont=icont+1
           if(icont>icontMax .and. mark==0)then
           write(*,*)'Calidad max en r',a(k+1),'al tiempo',t(k+1)
           mark = 1
           !exit time
          end if

        end if
        !write(11,*)t(k+1),AveSw,AveNf,AveUg,AveUw,Averg,Averc,DeltaP,NuGAverage
        iAvaOld=iAva
        kcount = kcount + 1
        if(kcount == 20 .or. k==1) then
             write(* ,*)t(k+1),AveSw,AveNf,AveUg,AveUw,Averg,Averc,DeltaP,NuGAverage
             write(11,*)t(k+1),AveSw,AveNf,AveUg,AveUw,Averg,Averc,DeltaP,NuGAverage
            kcount = 0
        end if
        do i=1,imax
          rg(i,0) = rg(i,1)
          rc(i,0) = rc(i,1)
          pc(i,0) = pc(i,1)
          Uw(i,0) = Uw(i,1)
          Ug(i,0) = Ug(i,1)
          Krw(i,0) = Krw(i,1)
          Krg0(i,0) = Krg0(i,1)
          Fw(i,0) = Fw(i,1)
        end do
      end do time

      close(unit=9)
      call OutputFoam

end program MainModFoamMech
