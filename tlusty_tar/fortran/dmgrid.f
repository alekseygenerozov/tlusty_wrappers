      program dmgrid
      IMPLICIT double precision (A-H,O-Z), LOGICAL*1 (L)
      PARAMETER(MDEPTH=200)
      common/ints/nd,imod,nmod
      common/reals/dm1,dmtot
      COMMON/MODPAR/DM(MDEPTH),DMM(MDEPTH)

      call input
c      ND=70
c      DM1=0.001
c      DMTOT=1.0D4
      DM(ND)=DMTOT
      DM(ND-1)=0.99D0*DMTOT
      DML=LOG(DM(ND-1)/DM1)/(ND-2)
      DML1=LOG(DM1)
      DO ID=1,ND-1
         DM(ID)=EXP(DML1+(ID-1)*DML)
      END DO
      if(imod .eq. 1) then
         do id=41,ND
            dmm(id+nmod)=dm(id)
         enddo
         dml=LOG(DM(40)/DM1)/(40.0+nmod-1.0)
         do id=1,40+nmod
            DMM(ID)=EXP(DML1+(ID-1)*DML)
         enddo
         ND=ND+nmod
         do id=1,ND
            dm(id)=dmm(id)
         enddo
      endif

      call output
 
      end

c ***********************************************************************
c ***********************************************************************
 
      subroutine output
 
      IMPLICIT double precision (A-H,O-Z), LOGICAL*1 (L)
      PARAMETER(MDEPTH=200)
      common/ints/nd,imod,nmod
      COMMON/MODPAR/DM(MDEPTH),DMM(MDEPTH)
 
c      WRITE(8,501) ND
      OPEN(8,FILE='dmgrid.out',STATUS='NEW')
      write(8,501) 1
      WRITE(8,502) (DM(ID),ID=1,ND)
      CLOSE(8)
 
  501 FORMAT(I5)
  502 FORMAT(1P6D13.6)
  503 FORMAT(1P5D15.6)
 
      end

c ***********************************************************************
c ***********************************************************************
 
      subroutine input
 
      IMPLICIT double precision (A-H,O-Z), LOGICAL*1 (L)
      PARAMETER(MDEPTH=200)
      common/ints/nd,imod,nmod
      common/reals/dm1,dmtot
 
c      WRITE(8,501) ND
      OPEN(8,FILE='dmgrid.in',STATUS='OLD')
      READ(8,*) ND,dm1,dmtot,imod,nmod
      CLOSE(8)
 

      end

c ***********************************************************************
c ***********************************************************************
 
