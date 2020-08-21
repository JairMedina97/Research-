subroutine InputModFoamMech
  use VarModFoamMech
  implicit none

    open(unit=8,file='InputModFoamMech.dat')


      read(8,*)m
      read(8,*)n
      read(8,*)SwAster
      read(8,*)PermAbs
      read(8,*)phy
      read(8,*)Nuw0
      read(8,*)Nug0
      read(8,*)Swc
      read(8,*)Sgr
      read(8,*)iMax
      read(8,*)kMax
      read(8,*)tMax
      read(8,*)Swini
      read(8,*)nfini
      read(8,*)FwJ
      read(8,*)Krg00
      read(8,*)Krw0
      read(8,*)PwJ
      read(8,*)Swc
      read(8,*)Sgr
      read(8,*)nkrg
      read(8,*)nkrw
      read(8,*)n2Pc
      read(8,*)iterMax
      read(8,*)ConverMax
      read(8,*)relax
      read(8,*)gama
      read(8,*)Dam
      read(8,*)Gam
      read(8,*)Ca
      read(8,*)Br
      read(8,*)rw
      read(8,*)eps
      read(8,*)icontMax
      read(8,*)Br2
      read(8,*)Sgcr
      read(8,*)MediaPoro
      read(8,*)DesStanPoro
      read(8,*)nPer
      UwJ=FwJ

  write(*,*)'UwJ',UwJ

  kcount = 0

    close(unit=8)

end subroutine InputModFoamMech
