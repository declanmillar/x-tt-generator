subroutine ZZwidthMulti(nzp,mass,rlambdaQCD4,nloop,width,neumass,ineu)
! NOTE 23 Jan 2014 found a factor two wrong in the convention for gv and ga
! in the E6 models coming from the factor 2 in the feynman rule for Elenas paper
! All widths have been rescaled by 1/4
! calculates Zp width contributions from light fermions
  implicit real*8 (a-h,o-z)
  real*8 :: mass(10),width(10)
! Lorenzos extra vars
  real*8 :: NeuMass
  common/Z/rmZ,gamZ
  common/t/rmt,gamt
  common/b/rmb
  common/LR/gZZu(2,10),gZZd(2,10),gZZe(2,10),gZZn(2,10), &
  gZZc(2,10),gZZs(2,10),gZZt(2,10),gZZb(2,10), &
  gZZta(2,10),gZZnt(2,10)
! ccccccccccccccccccc
! initialise
  do n=1,10
      width(n)=0.d0
  enddo
! cccccccccccccccccccc
! couplings.
  pi=dacos(-1.d0)
  do n = 1,nzp
  ! ZZ width.
      ZZwidth=0.d0
      temptt=0.d0
      tempqq=0.d0
      if (mass(n) > 0.d0) then
          a_s=alfas(mass(n),rlambdaQCD4,nloop)
      else
          goto 124
      endif
      do i=1,6
          tempZZ=0.d0
          if((i == 2) .OR. (i == 4))then
          ! u-quark.
              rmq=0.d0
              gv = (gzzu(1,n)+gzzu(2,n))
              ga = (gzzu(1,n)-gzzu(2,n))
          else if((i == 1) .OR. (i == 3))then
          ! d-quark.
              rmq=0.d0
              gv = (gzzd(1,n)+gzzd(2,n))
              ga = (gzzd(1,n)-gzzd(2,n))
          else if(i == 5)then
              rmq=rmb
              gv = (gzzb(1,n)+gzzb(2,n))
              ga = (gzzb(1,n)-gzzb(2,n))
          else if(i == 6)then
              rmq=rmt
              gv = (gzzt(1,n)+gzzt(2,n))
              ga = (gzzt(1,n)-gzzt(2,n))
          end if
          if(mass(n) <= 2.d0*rmq)goto 123
      ! with QCD kfactor
          tempZZ=3.d0/48.d0/pi*mass(n) &
          *sqrt(1.d0-4.d0*rmq**2/mass(n)**2) &
          *(gv**2*(1.d0+2.d0*rmq**2/mass(n)**2) &
          +ga**2*(1.d0-4.d0*rmq**2/mass(n)**2)) &
          *(1.d0+1.045d0*a_s/pi)
          ZZwidth=ZZwidth+tempZZ
      ! without QCD kfactor
      !        tempZZ = 3.d0/12.d0/pi*mass(n)
      !     &               *sqrt(1.d0-4.d0*rmq**2/mass(n)**2)
      !     &               *(gv**2*(1.d0+2.d0*rmq**2/mass(n)**2)
      !     &                +ga**2*(1.d0-4.d0*rmq**2/mass(n)**2))
      !        ZZwidth=ZZwidth+tempZZ
      !             print*,i,tempZZ
      ! cccccccccccccccccccc
          if(i == 6) then
              temptt=tempZZ
          endif
          if((i == 5) .OR. (i == 4) .OR. (i == 3) .OR. (i == 2) .OR. (i == 1)) then
              tempqq=tempqq+tempZZ
          endif
                
      ! cccccccccccccccccccc
          123 continue
      end do

  !      print *, 'ineu = ', ineu
  !       print*,'Z'' width due to quarks:',ZZwidth,' GeV'
  ! ccccccccccccccccccc
      temp=0.d0
      temp1=0.d0
      temp2=0.d0
      temp3=0.d0
  ! cccccccccccccccccccc
      do i=1,9
          tempZZ=0.d0
          rmq=0.d0
          if((i == 2) .OR. (i == 4))then
          ! neutrino.
              if(iNeu == 0) then
                  gV = (gzzn(1,n)+gzzn(2,n))
                  gA = (gzzn(1,n)-gzzn(2,n))
                  dof=1.d0
              else if(iNeu == 1) then
                  gV=0.d0
                  gA = gzzn(1,n)
                  dof=0.5d0
              endif
          else if(i == 6)then ! Tau neutrino
              rmq=NeuMass
              if(iNeu == 0) then
                  gV = (gzznt(1,n)+gzznt(2,n))
                  gA = (gzznt(1,n)-gzznt(2,n))
                  dof=1.d0
              else if(iNeu == 1) then
                  gV=0.d0
                  gA = gzznt(1,n)
                  dof=0.5d0
              endif
          else if((i == 1) .OR. (i == 3))then
          ! lepton.
              gv = (gzze(1,n)+gzze(2,n))
              ga = (gzze(1,n)-gzze(2,n))
              dof=1.d0
          else if(i == 5)then ! Tau
              gv = (gzzta(1,n)+gzzta(2,n))
              ga = (gzzta(1,n)-gzzta(2,n))
              dof=1.d0
          else if((i == 7) .OR. (i == 8) .OR. (i == 9))then
          ! heavy neutrinos
              rmq=NeuMass
              if(iNeu == 0) then
                  gV = (gzzn(1,n)+gzzn(2,n))
                  gA = (gzzn(1,n)-gzzn(2,n))
                  dof=1.d0
              else if(iNeu == 1) then
                  gV=0.d0
                  gA = gzzn(1,n)
                  dof=0.5d0
              endif
          end if
          if(mass(n) <= 2.d0*rmq)goto 456
          tempZZ = dof/48.d0/pi*mass(n) &
          *sqrt(1.d0-4.d0*rmq**2/mass(n)**2) &
          *(gv**2*(1.d0+2.d0*rmq**2/mass(n)**2) &
          +ga**2*(1.d0-4.d0*rmq**2/mass(n)**2))
          ZZwidth=ZZwidth+tempZZ
                  
      ! cccccccccccccccccccc
          temp=temp+tempZZ
          if((i == 1) .OR. (i == 3) .OR. (i == 5))then
              temp2=temp2+tempZZ
          else if((i == 2) .OR. (i == 4) .OR. (i == 6))then
              temp1=temp1+tempZZ
          else if((i == 7) .OR. (i == 8) .OR. (i == 9))then
              temp3=temp3+tempZZ
          endif
          tempZZ=0.d0
      ! cccccccccccccccccccc
      end do
      456 continue
      width(n)=ZZwidth
      124 continue
  ! cccccccccccccccccccc
      print*,'########### Width and BRs ###########'
      print *,'Gamma(Zp(',n,')->ff)=',width(n),' GeV'
      print *,'Gamma(Zp(',n,')->tt)=',temptt,' GeV'
      print *,'Gamma(Zp(',n,')->qq)=',tempqq,' GeV'
      print *,'Gamma(Zp(',n,')->ll)=',temp,' GeV'
      print *,'Gamma(Zp(',n,')->e/mu/tau)=',temp2,' GeV'
      print *,'Gamma(Zp(',n,')->nu)=',temp1,' GeV'
      print *,'Gamma(Zp(',n,')->nuH)=',temp3,' GeV'
  ! cccccccccccccccccccc
  end do
  return
end subroutine ZZwidthMulti
