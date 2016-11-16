C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit none
C  Input variables
      integer NatomsA, NatomsB, NelmtsA, NelmtsB
      integer NA1, NA2, NA3, NB1, NB2, NB3, homoA, homoB
      real*8 xa0, ya0, za0, xa1, ya1, za1, xa2, ya2, za2, xa3, ya3, za3
      real*8 xb0, yb0, zb0, xb1, yb1, zb1, xb2, yb2, zb2, xb3, yb3, zb3
      real*8 cena(3), cenb(3), cenaang(3), cenbang(3), center(3)
      real*8 expA, expB, abssuma, abssumb
      character*80 titleA, titleB, dummy
	real*8 sumint

C  Dynamically allocated array for transition density cube elements and
C    for atomic number and coordinates
      real*8, allocatable, dimension(:,:,:) :: qA, qB      
      real*8, allocatable, dimension(:) :: pxA, pyA, pzA, pxB, pyB, pzB,
     & atnoA, atnoB
      integer, allocatable, dimension(:) :: IA, IB

C  Distance vectors
      real*8 R0(3),R,RAB(3)
      real*8 riA(3), rjB(3), rij, Rcutoff

C  Work
      real*8 sumA, sumB, volA, volB, correctA, correctB
      real*8 exp_correctA, exp_correctB, exp_correct
C      real*8 Vdd, Vdq, Vqd, Vqq, Vdo, Vod
C      real*8 deltaab, deltaaa, deltabb, Ra, Rb
      real*8 dotAB, dotAR0, dotBR0
 
C  Output
      real*8 dipoleA(3), dipoleB(3), dipoleAexp(3), dipoleBexp(3) 
      real*8 dipoleAD(3), dipoleBD(3)
      real*8 quadA(6), quadB(6), quadAexp(6), quadBexp(6)
      real*8 dipoleADexp(3), dipoleBDexp(3)
      real*8 kappa, muA, muAD, muB, muBD
      real*8 muAexp, muBexp, muBDexp, muADexp
      real*8 Texact, Texactexp
C      real*8 Tdd, Tdq, Tqd, Tqq, Tdo, Tod
C      real*8 TabTddexp, Tdqexp, Tqdexp, Tqqexp
C      real*8 Tdoexp, Todexp, Tabexp
      real*8 Vdipdip, Vdipdipexp

C  Timer variables
C      real*4 timearray(2), calctime            !This line for Unix
C      real*8 calctime, TIMEF                  !This line for Win32

C  Constants
      real*8 :: BohrtoAng=0.529177249
      real*8 :: echarge=1.602188E-19
      real*8 :: ep0=8.8542e-12
      real*8 :: pi=3.1415926535
      real*8 :: HtoeV=27.2116
      real*8 :: JtoHartree=4.3597482e-18
      real*8 :: eVtocm=8065.
      real*8 :: eAngtoDebye=4.803

C  Loop variables, flags
      integer i,i1,i2,i3,j1,j2,j3
      integer irflag


C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C      expA = 13.      !RG1
C      expB = 6.13     !Bchl Qy
C      expB = 3.29     !Bchl Qx
      expA = 1.
      expB = 1.
      sumA = 0.
      abssumA = 0.
      center(1) = 0.
      center(2) = 0.
      center(3) = 0.
      dipoleA(1) = 0.
      dipoleA(2) = 0.
      dipoleA(3) = 0.
      muA = 0.
      muB = 0.
      sumB = 0.
      abssumB = 0.
      dipoleB(1) = 0.
      dipoleB(2) = 0.
      dipoleB(3) = 0.
      irflag = 0
C      quadB(1) = 0.
C      quadB(2) = 0.
C      quadB(3) = 0.
C      quadB(4) = 0.
C      quadB(5) = 0.
C      quadB(6) = 0.
      quadA(1) = 0.
      quadA(2) = 0.
      quadA(3) = 0.
      quadA(4) = 0.
      quadA(5) = 0.
      quadA(6) = 0.

C  Open and process cubeA
  3   open (unit=16, file='cubeA', err=999)      
      read (16,5) titleA
      read (16,5) dummy
      write (6,5) titleA
  5   format(A80)
  
C  Number of atoms and origin
      read (16,*) NatomsA, xa0, ya0, za0
C  Number of increments and size in each direction
      read (16,*) NA1, xa1, ya1, za1
      read (16,*) NA2, xa2, ya2, za2
      read (16,*) NA3, xa3, ya3, za3
      write (6,*) NatomsA, NA1, NA2, NA3
C  Allocate arrays for atomic coords, etc.
      allocate( IA(NatomsA), atnoA(NatomsA), pxA(NatomsA), pyA(NatomsA),
     &  pzA(NatomsA) )
C  Atomic number, charge and atomic positions
      read (16,'(I5,4F12.6)') (IA(i),atnoA(i),pxA(i),pyA(i)
     $   ,pzA(i), i=1,NatomsA)

C  Read in the density cube 
C  (Note that units are es and Bohr.)
      allocate( qA(NA3, NA2, NA1) )     
      qA = 0.            

      do i1=1,NA1
        do i2=1,NA2
          read (16,'(6E13.5)') (qA(i3,i2,i1), i3=1,NA3)          
        end do
      end do

      close (unit=16)
      

C  ----------------------------------------------------------------
C  Before we do anything with the cubes, we must correct for the small 
C    residual charge present.  This charge arises because the cube
C    files only store five significant figures.
      sumA = 0.
      sumB = 0.

      do i1 = 1, NA1
      do i2 = 1, NA2
      do i3 = 1, NA3
        sumA = sumA + qA(i3,i2,i1)
      end do
      end do
      end do
      NelmtsA = NA1 * NA2 * NA3
c      if (homoA .eq. 0) then
         correctA = sumA / NelmtsA
c      else
c         correctA = (sumA-2.) / NelmtsA
c      endif
C  This is the element volume correction term.  The sign is
C    negative to account for positive electron density
C    being actually negative charge density.
      VolA = -DSqrt((xa1*xa1+ya1*ya1+za1*za1)*
     $              (xa2*xa2+ya2*ya2+za2*za2)*
     $              (xa3*xa3+ya3*ya3+za3*za3))

      qA = qA - correctA
c      sumA = 0.

C  Define a volume dependent Rcutoff.  This is a bit arbitrary, but
C    is intended to allow some significant amount of overlap, but to not
C    let R get too small.  If VolA and VolB are the same and are cube-shaped
C    with sides of length s, then Rcutoff becomes (1/2Sqrt(2)) * s ~= .35s.

      Rcutoff = 0.25*DSqrt(abs(VolA)**(2.0/3.0)+abs(VolB)**(2.0/3.0))
C  
C  Now let's determine the center of 'charge' of the transition
C    densities rather than just assuming that the first point of 
C    the file is the center as we did before.
C
      do i1=1,NA1
      do i2=1,NA2
      do i3=1,NA3
        abssumA = abssumA + abs(qA(i3,i2,i1))
        riA(1) = xa0 + (i1-1)*xa1 + (i2-1)*xa2 + (i3-1)*xa3 - pxa(1)
        riA(2) = ya0 + (i1-1)*ya1 + (i2-1)*ya2 + (i3-1)*ya3 - pya(1)
        riA(3) = za0 + (i1-1)*za1 + (i2-1)*za2 + (i3-1)*za3 - pza(1)
        
        center =  center + (abs(qA(i3,i2,i1))*riA)
      end do
      end do
      end do

      center = center / abssumA

      cena(1) = pxa(1) + center(1)
      cena(2) = pya(1) + center(2)
      cena(3) = pza(1) + center(3)
C  Convert to angstroms
      cenaang = cena*BohrtoAng

C
C  Now determine charge and dipole moments of cubeA and cubeB
C

C  Charge and transition dipole moment from cubeA

	sumint=0.d0

      do i1=1,NA1
      do i2=1,NA2
      do i3=1,NA3

        sumA = sumA + qA(i3,i2,i1)

	  sumint=sumint+qA(i3,i2,i1)**2
c        riA(1) = xa0 + (i1-1)*xa1 + (i2-1)*xa2 + (i3-1)*xa3 - cena(1)
c        riA(2) = ya0 + (i1-1)*ya1 + (i2-1)*ya2 + (i3-1)*ya3 - cena(2)
c        riA(3) = za0 + (i1-1)*za1 + (i2-1)*za2 + (i3-1)*za3 - cena(3)
        
c        dipoleA = dipoleA + qA(i3,i2,i1)*riA
C  Quadrupole moment relative to centre A

c       quadA(1) = quadA(1) + qA(i3,i2,i1)*riA(1)*riA(1)
c       quadA(2) = quadA(2) + qA(i3,i2,i1)*riA(2)*riA(2)
c       quadA(3) = quadA(3) + qA(i3,i2,i1)*riA(3)*riA(3)
c       quadA(4) = quadA(4) + qA(i3,i2,i1)*riA(1)*riA(2)
c       quadA(5) = quadA(5) + qA(i3,i2,i1)*riA(1)*riA(3)
c       quadA(6) = quadA(6) + qA(i3,i2,i1)*riA(2)*riA(3)

      end do
      end do
      end do
C  Now convert from e*Bohr to e*Angstroms and normalize.      
c     dipoleA = dipoleA * BohrtoAng * volA
c     quadA = quadA * BohrtoAng*BohrtoAng * volA
	sumint=sumint*volA
	
      write (6,*) sumA
c      write (6,*) dipoleA(1), dipoleA(2), dipoleA(3)
	 write(6,*) 'volA=', volA
	 write(6,*) 'sumint=', sumint
	
c      dipoleAD = dipoleA*eAngtoDebye
c      write (6,*) quadA(1), quadA(2), quadA(3)
c      write (6,*) quadA(4), quadA(5), quadA(6)


C  Dipole-dipole coupling from point dipole model,
C  with separation, R, taken to be between the calculated centers of
C  charge of each transition density.

C  Unit vectors and magnitude defining the intermolecular separation
C  vector about which the Taylor expansion of the interaction potential
C  is taken.  This is assumed to be between the centers of charge of 
C  each transition density.

c      muA = DSQRT(DOT_PRODUCT(dipoleA,dipoleA))
c      muAD = DSQRT(DOT_PRODUCT(dipoleAD,dipoleAD))
 
c      write (6,*) ' '
c      write (6,*) 'Monomer A:'
c      write (6,100) titleA
c      write (6,*) ' SumA = ', sumA,' e in',NelmtsA,' elements'
c      write (6,96) cena(1), cena(2), cena(3)
c      write (6,97) dipoleA(1),dipoleA(2),dipoleA(3)
c      write (6,*) '         |Dipole| = ', muA, ' e.Ang'
c      write (6,98) dipoleAD(1),dipoleAD(2),dipoleAD(3)
c      write (6,*) '         |Dipole| = ', muAD, ' Debye'
c      write (6,99) (quadA(i), i=1,6)

        goto 1000

C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Error traps
 990  write (6,*) 'Error opening expAB file.'
      write (6,*) 'Continuing without use of experimental transition
     & moments.'
      goto 3
 999  write (6,*) ' Ahhhhhhh!!!'

 1000 continue
  96  format (' Center of Charge  = [',3(1X,E12.5),'  ]  bohr')
  97  format (' Calculated Dipole = [',3(1X,E12.5),'  ]  e.Ang')
  98  format (' Calculated Dipole = [',3(1X,E12.5),'  ]  Debye')
  99  format (' Calculated Quadrupole transition moments (e.Ang^2)'/
     &        '   XX=',F10.3,'       YY=',F10.3,'       ZZ=',F10.3/
     &        '   XY=',F10.3,'       XZ=',F10.3,'       YZ=',F10.3/)
 100  format (A80)
 101  format (A80)
 102  format (' Scaled Dipole = [',3(1X,E12.5),'  ]  e.Ang')
 103  format (' Scaled Dipole = [',3(1X,E12.5),'  ]  Debye')
 104  format (' Scaled Quadrupole transition moments (e.Ang^2)'/
     &        '   XX=',F10.3,'       YY=',F10.3,'       ZZ=',F10.3/
     &        '   XY=',F10.3,'       XZ=',F10.3,'       YZ=',F10.3/)
       end


