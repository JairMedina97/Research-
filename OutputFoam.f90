subroutine OutputFoam
    use VarModFoamMech
    implicit none

       open(unit=8,file='sw.dat')
       do i=1,iMax
          write(8,102)x(i), (Sw(i,k), k=1,kMax,kMax/100)
       end do
       close(unit=8)
      !nfMax=Maxval(nf(1:iMax,1:kMax))

      open(unit=8,file='nf.dat')
       do i=1,iMax
          write(8,102)x(i), (nf(i,k), k=1,kMax,kMax/100)
       end do
       close(unit=8)

     open(unit=7,file='Pres.dat')
       do i=1,iMax
          write(7,102)x(i), (P(i,k), k=1,kMax,kMax/100)
       end do
      close(unit=7)

     open(unit=9,file='NuG.dat')
       do i=1,iMax
         write(9,102)x(i), (NuG(i,k), k=1,kMax,kMax/100)
       end do
     close(unit=9)

      open(unit=9,file='MovG.dat')
       do i=1,iMax
         write(9,102)x(i), (lamdaG(i,k), k=1,kMax,kMax/100)
       end do
     close(unit=9)

     open(unit=9,file='MovW.dat')
       do i=1,iMax
         write(9,102)x(i), (lamdaW(i,k), k=1,kMax,kMax/100)
       end do
     close(unit=9)

     102 format(900e18.8)

end subroutine OutputFoam
