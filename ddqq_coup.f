c modified by Cherri Hsu to do dipole-quadrupole, dip-octapole, and quad-quad 
c couplings.  12/2000
c results are not quite the same as those from coupling4.f mainly due to the
c proper acount (I think) of Taylor expansion factors.  
c There is 
c      a factor of 2 in dip-quad and quad-dip
c      a factor of -4 in quad-quad and quad-quad
c      a factor of 6 in dip-oct and oct-dip
c the agreement between Tab (sum of those multipole interaction) and Texact 
c is improved with this version.
C
C                      DIPOLE-DIPOLE-ONLY version of 
C                      COUPLING for fort90
C  Multipolar coupling between charge transition density cubes (TDCs)
C     Original work by Greg Scholes, May 1997 (University of Chicago)
C
C                      Version 5   January 1999
C 
C  Transition density cubes output by Gaussian are used
C  to calculate the "Coulombic" (1/rij) coupling between large
C  molecules.  The results are compared to the dipole approximation
C  and a multipole expansion of the interaction potential which
C  Includes dipole-dipole, dipole-quadrupole, quadrupole-quadrupole
C  and dipole-octapole couplings between the transition densities.
C
C  Note that the density cube must be generated using cube=full
C  so that the integrated density = 0 and must be formatted.
C
C  The first atomic coordinates in each file are taken to be the
C  molecular centres.  It is important that these coordinates
C  (which may be of dummy atoms) are for the centre of the leading
C  transition moment.  (This is now obsolete... see version 4.)
C
C  Version 2 is now fortran90.  All large arrays are dynamically
C  allocated to make it flexible with use of different sized cubes.
C  Also, outputs are now also given in 'normal' units like Angstroms,
C  Debye and wavenumbers.                (Brent Krueger -- June 1997)
C
C  Version 3 now properly normalizes cube properties (and coupling)
C  by taking element volume into account (as in cubman.F from Gaussian).
C  Note that this is a key difference between TDC couplings and the
C  monopole method.  Monopolers can't normalize properly.
C  Also, this version allows user to input experimental transition
C  dipole moment and outputs normalized coupling along with 
C  coupling scaled to the experimental transition dipole magnitude.
C  The experimental transition dipoles must be in a file called
C  coupling.in which contains only two numbers, the experimental transition
C  dipole moments (in Debye) for cubeA and cubeB in that order (i.e.
C  a two line file with one number per line).  Note that there is also
C  an bit of untested code that is intended for use if one of the 
C  cubes is to be an orbital (HOMO, HOMO-1 etc.) then the experimental 
C  transition moment should be entered as zero signalling program to 
C  properly scale total charge to 2 rather than zero.  
C                                       (Brent Krueger -- October 1997)
C
C  Version 4 now calculates the center of 'charge' of the transition
C  density and uses this point as the center to expand the multipolar
C  transition moments about rather than just the first atom in the cubefile.
C  Exact couplings are, of course, unchanged.
C                                       (Brent Krueger -- June 1998)
C
C  Version 5 just calculates exact (TDC) coupling and dipole-dipole
C  coupling.  The higher-order multipole stuff has been removed in favor
C  of speed.  Also, Rcutoff which determines the distance at which cubes
C  are considered to be overlapping is now determined by the element 
C  volumes of cubeA and cubeB rather than aribtrarily set as a constant.
C  Also, some vector manipulation is now handled by intrinsic commands
C  (e.g. DOT-PRODUCT) rather than with clumsy loops.
C                                       (Brent Krueger -- January 1999)
C
C  Dipole-Dipole-Only version Does not calculate TDC coupling, just the 
C  individual transition dipoles and their interactions, for quick work.
C  This does properly center the dipoles on the center of charge, etc.
C                                       (Brent Krueger -- February 1999)
C
C 
C  Multipolar couplings are described in:
C   G.D.Scholes and D.L.Andrews, J.Chem.Phys. (1997) submitted.
C  For a description of this TDC coupling technique:
C   B.P.Krueger, G.D.Scholes, and G.R.Fleming J.Phys.Chem.B. (1998)
C   also  B.P.Krueger Doctoral Dissertation, U of Chicago (1999)
C
C  Input files:    cubeA, cubeB  (which should be formatted)
C                  coupling.in (contains the experimental dipole magnitudes)
C  Output file:    coupling.out
C
C  This code has been used with SGI, HP, Cray and Windows.  Some modification
C  is necessary for use in Windows.  These lines are commented in the text.
C
C  The compile line on SGI:  f90 -O3 -o coupling5 coupling5.f
C  The compile line on HP:   f90 +O4 +U77 -o coupling5 coupling5.f
C
C  *Note on HP-UX use +U77 flag to get timer to work
C
C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit none
C  Input variables
      integer NatomsA, NatomsB, NelmtsA, NelmtsB
      real*8  rNelmtsA, rNelmtsB
      integer NA1, NA2, NA3, NB1, NB2, NB3, homoA, homoB
      real*8 xa0, ya0, za0, xa1, ya1, za1, xa2, ya2, za2, xa3, ya3, za3
      real*8 xb0, yb0, zb0, xb1, yb1, zb1, xb2, yb2, zb2, xb3, yb3, zb3
      real*8 cena(3), cenb(3), cenaang(3), cenbang(3), center(3)
      real*8 expA, expB, abssuma, abssumb
      character*80 titleA, titleB, dummy

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
      real*8 Vdd, Vdq, Vqd, Vqq, Vdo, Vod
      real*8 deltaab, deltaaa, deltabb, Ra, Rb
      real*8 dotAB, dotAR0, dotBR0
 
C  Output
      real*8 dipoleA(3), dipoleB(3), dipoleAexp(3), dipoleBexp(3) 
      real*8 dipoleAD(3), dipoleBD(3)
      real*8 quadA(6), quadB(6), quadAexp(6), quadBexp(6)
      real*8 octaA(10), octaB(10)
      real*8 dipoleADexp(3), dipoleBDexp(3)
	real*8 quA(3,3), quB(3,3), oA(3,3,3), oB(3,3,3)
	real*8 k2(3,3), k3(3,3,3), k4(3,3,3,3),testdd
      real*8 kappa, muA, muAD, muB, muBD
      real*8 muAexp, muBexp, muBDexp, muADexp
      real*8 Texact, Texactexp
      real*8 Tdd, Tdq, Tqd, Tqq, Tdo, Tod, Tab
      real*8 TabTddexp, Tdqexp, Tqdexp, Tqqexp
      real*8 Tdoexp, Todexp, Tabexp
      real*8 Vdipdip, Vdipdipexp
	real*8 Vquaddip, Vdipquad, Vquadquad, Voctdip, Vdipoct

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
      integer i,i1,i2,i3,j1,j2,j3,j,k,l
C      integer irflag


C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C      expA = 13.      !RG1
C      expB = 6.13     !Bchl Qy
C      expB = 3.29     !Bchl Qx
      expA = 1.
      expB = 1.
      sumA = 0.
      center(1) = 0.
      center(2) = 0.
      center(3) = 0.
      dipoleA(1) = 0.
      dipoleA(2) = 0.
      dipoleA(3) = 0.
      muA = 0.
      muB = 0.
      sumB = 0.
      dipoleB(1) = 0.
      dipoleB(2) = 0.
      dipoleB(3) = 0.
C      irflag = 0

      quadB(1) = 0.
      quadB(2) = 0.
      quadB(3) = 0.
      quadB(4) = 0.
      quadB(5) = 0.
      quadB(6) = 0.
      quadA(1) = 0.
      quadA(2) = 0.
      quadA(3) = 0.
      quadA(4) = 0.
      quadA(5) = 0.
      quadA(6) = 0.
	DO i=1,10
	  octaA(i) = 0.
	  octaB(i) = 0.
	ENDDO

C      calctime = TIMEF()      !This line for Win32

C  The transition density cube files (from Gaussian)
C  are called cubeA and cubeB.  The file with experimental 
C  transition dipole moments is called coupling.in.
C
C  Open coupling.in if it exists.
      open (unit=18, file='coupling.in', err=990)
      read (18,*) expA
      read (18,*) expB
      close (unit=18)
      if ((expA < 1e-10).and.(expA > -1e-10)) then
         homoA = 1
      else
         homoA = 0
      endif
      if ((expB < 1e-10).and.(expB > -1e-10)) then
         homoB = 1
      else
         homoB = 0
      endif
C
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

      do i1=1,NA1
        do i2=1,NA2
          read (16,'(6E13.5)') (qA(i3,i2,i1), i3=1,NA3)          
        end do
      end do

      close (unit=16)
      
C  --------------------------------------------------------------
C  Now do the same for cubeB

      open (unit=14, file='cubeB', err=999)
      
      read (14,5) titleB
      read (14,5) dummy
      write (6,5) titleB
  
      read (14,*) NatomsB, xb0, yb0, zb0
      read (14,*) NB1, xb1, yb1, zb1
      read (14,*) NB2, xb2, yb2, zb2
      read (14,*) NB3, xb3, yb3, zb3
      write (6,*) NatomsB, NB1, NB2, NB3

      allocate( IB(NatomsB), atnoB(NatomsB), pxB(NatomsB), pyB(NatomsB),
     &  pzB(NatomsB) )            
 
      read (14,'(I5,4F12.6)') (IB(i),atnoB(i),pxB(i),pyB(i)
     $   ,pzB(i), i=1,NatomsB)

      allocate( qB(NB3, NB2, NB1) )

      do i1=1,NB1
        do i2=1,NB2
          read (14,'(6E13.5)') (qB(i3,i2,i1), i3=1,NB3)          
        end do
      end do

      close (unit=14)


C  ----------------------------------------------------------------
C  Before we do anything with the cubes, we must correct for the small 
C    residual charge present.  This charge arises because the cube
C    files only store five significant figures.
      do i1 = 1, NA1
      do i2 = 1, NA2
      do i3 = 1, NA3
        sumA = sumA + qA(i3,i2,i1)
      end do
      end do
      end do
      NelmtsA = NA1 * NA2 * NA3
      rNelmtsA = NelmtsA
      if (homoA .eq. 0) then
         correctA = sumA / rNelmtsA
      else
         correctA = (sumA-2.) / rNelmtsA
      endif
C  This is the element volume correction term.  The sign is
C    negative to account for the positive electron density
C    actually negative charge density.
      VolA = -DSqrt((xa1*xa1+ya1*ya1+za1*za1)*
     $              (xa2*xa2+ya2*ya2+za2*za2)*
     $              (xa3*xa3+ya3*ya3+za3*za3))

      sumA = 0.
      qA = qA - correctA

C  Correct residual charge for cubeB

      do i1 = 1, NB1
      do i2 = 1, NB2
      do i3 = 1, NB3
        sumB = sumB + qB(i3,i2,i1)
      end do
      end do
      end do
      NelmtsB = NB1 * NB2 * NB3
      rNelmtsB = NelmtsB
      if (homoB .eq. 0) then
         correctB = sumB / rNelmtsB
      else
         correctB = (sumB-2.) / rNelmtsB
      endif
      VolB = -DSqrt((xb1*xb1+yb1*yb1+zb1*zb1)*
     $              (xb2*xb2+yb2*yb2+zb2*zb2)*
     $              (xb3*xb3+yb3*yb3+zb3*zb3))

      sumB = 0.
      qB = qB - correctB


C  Define a volume dependent Rcutoff.  This is a bit arbitrary, but
C    is intended to allow some significant amount of overlap, but to not
C    let R get too small.  If VolA and VolB are the same and are cube-shaped
C    with sides of length s, then Rcutoff becomes Sqrt(2)/10 * s.

C      Rcutoff = 0.1*DSqrt(abs(VolA)**(2.0/3.0)+abs(VolB)**(2.0/3.0))
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
        
        center =  center + abs(qA(i3,i2,i1))*riA
      end do
      end do
      end do

      center = center / abssumA

      cena(1) = pxa(1) + center(1)
      cena(2) = pya(1) + center(2)
      cena(3) = pza(1) + center(3)
C  Convert to angstroms
      cenaang = cena*BohrtoAng

C  Reset center variable
      center = 0.
C   again for cubeB
      do i1=1,NB1
      do i2=1,NB2
      do i3=1,NB3
        abssumB = abssumB + abs(qB(i3,i2,i1))
        rjB(1) = xb0 + (i1-1)*xb1 + (i2-1)*xb2 + (i3-1)*xb3 - pxb(1)
        rjB(2) = yb0 + (i1-1)*yb1 + (i2-1)*yb2 + (i3-1)*yb3 - pyb(1)
        rjB(3) = zb0 + (i1-1)*zb1 + (i2-1)*zb2 + (i3-1)*zb3 - pzb(1)
        
        center =  center + abs(qB(i3,i2,i1))*rjB
      end do
      end do
      end do

      center = center / abssumB

      cenb(1) = pxb(1) + center(1)
      cenb(2) = pyb(1) + center(2)
      cenb(3) = pzb(1) + center(3)
C  Convert to angstroms
      cenbang = cenb*BohrtoAng

C
C  Now determine charge and dipole moments of cubeA and cubeB
C

C  Charge and transition dipole moment from cubeA


      do i1=1,NA1
      do i2=1,NA2
      do i3=1,NA3

        sumA = sumA + qA(i3,i2,i1)

        riA(1) = xa0 + (i1-1)*xa1 + (i2-1)*xa2 + (i3-1)*xa3 - cena(1)
        riA(2) = ya0 + (i1-1)*ya1 + (i2-1)*ya2 + (i3-1)*ya3 - cena(2)
        riA(3) = za0 + (i1-1)*za1 + (i2-1)*za2 + (i3-1)*za3 - cena(3)
        
        dipoleA = dipoleA + qA(i3,i2,i1)*riA

C  Quadrupole moment relative to centre A

        quadA(1) = quadA(1) + qA(i3,i2,i1)*riA(1)*riA(1)
        quadA(2) = quadA(2) + qA(i3,i2,i1)*riA(1)*riA(2)
        quadA(3) = quadA(3) + qA(i3,i2,i1)*riA(1)*riA(3)
        quadA(4) = quadA(4) + qA(i3,i2,i1)*riA(2)*riA(2)
        quadA(5) = quadA(5) + qA(i3,i2,i1)*riA(2)*riA(3)
        quadA(6) = quadA(6) + qA(i3,i2,i1)*riA(3)*riA(3)

C  octapole moment relative to centre A

	  octaA(1) = octaA(1) + qA(i3,i2,i1)*riA(1)*riA(1)*riA(1)
	  octaA(2) = octaA(2) + qA(i3,i2,i1)*riA(1)*riA(1)*riA(2)
	  octaA(3) = octaA(3) + qA(i3,i2,i1)*riA(1)*riA(1)*riA(3)
	  octaA(4) = octaA(4) + qA(i3,i2,i1)*riA(1)*riA(2)*riA(2)
	  octaA(5) = octaA(5) + qA(i3,i2,i1)*riA(1)*riA(2)*riA(3)
	  octaA(6) = octaA(6) + qA(i3,i2,i1)*riA(1)*riA(3)*riA(3)
	  octaA(7) = octaA(7) + qA(i3,i2,i1)*riA(2)*riA(2)*riA(2)
	  octaA(8) = octaA(8) + qA(i3,i2,i1)*riA(2)*riA(2)*riA(3)
	  octaA(9) = octaA(9) + qA(i3,i2,i1)*riA(2)*riA(3)*riA(3)
	  octaA(10)=octaA(10) + qA(i3,i2,i1)*riA(3)*riA(3)*riA(3)

      end do
      end do
      end do
C  Now convert from e*Bohr to e*Angstroms and normalize.      
      dipoleA = dipoleA * BohrtoAng * volA
      quadA = quadA * BohrtoAng*BohrtoAng * volA
      octaA = octaA * BohrtoAng*BohrtoAng*BohrtoAng * volA

      write (6,*) sumA
      write (6,*) dipoleA(1), dipoleA(2), dipoleA(3)
      write (6,*) quadA(1), quadA(2), quadA(3)
      write (6,*) quadA(4), quadA(5), quadA(6)
	write (6,*) octaA(1), octaA(2), octaA(3), octaA(4), octaA(5)
	write (6,*) octaA(6), octaA(7), octaA(8), octaA(9), octaA(10)

	quA(1,1) = quadA(1)
	quA(1,2) = quadA(2)
	quA(1,3) = quadA(3)
	quA(2,1) = quA(1,2)
	quA(2,2) = quadA(4)
	quA(2,3) = quadA(5)
	quA(3,1) = quA(1,3)
	quA(3,2) = quA(2,3)
	quA(3,3) = quadA(6)

	oA(1,1,1) = octaA(1)
	oA(1,1,2) = octaA(2)
	oA(1,1,3) = octaA(3)
	oA(1,2,1) = oA(1,1,2)
	oA(1,2,2) = octaA(4)
	oA(1,2,3) = octaA(5)
	oA(1,3,1) = oA(1,1,3)
	oA(1,3,2) = oA(1,2,3)
	oA(1,3,3) = octaA(6)
	oA(2,1,1) = oA(1,1,2)
	oA(2,1,2) = oA(1,2,2)
	oA(2,1,3) = oA(1,2,3)
	oA(2,2,1) = oA(1,2,2)
	oA(2,2,2) = octaA(7)
	oA(2,2,3) = octaA(8)
	oA(2,3,1) = oA(1,2,3)
	oA(2,3,2) = oA(2,2,3)
	oA(2,3,3) = octaA(9)
	oA(3,1,1) = oA(1,1,3)
	oA(3,1,2) = oA(1,2,3)
	oA(3,1,3) = oA(1,3,3)
	oA(3,2,1) = oA(1,2,3)
	oA(3,2,2) = oA(2,2,3)
	oA(3,2,3) = oA(2,3,3)
	oA(3,3,1) = oA(1,3,3)
	oA(3,3,2) = oA(2,3,3)
	oA(3,3,3) = octaA(10)

C  Charge and transition dipole moment from cubeB

      do i1=1,NB1
      do i2=1,NB2
      do i3=1,NB3
        sumB = sumB + qB(i3,i2,i1)
          
        rjB(1) = xb0 + (i1-1)*xb1 + (i2-1)*xb2 + (i3-1)*xb3 - cenb(1)
        rjB(2) = yb0 + (i1-1)*yb1 + (i2-1)*yb2 + (i3-1)*yb3 - cenb(2)
        rjB(3) = zb0 + (i1-1)*zb1 + (i2-1)*zb2 + (i3-1)*zb3 - cenb(3)
        
        dipoleB = dipoleB + qB(i3,i2,i1)*rjB

C  Quadrupole moment relative to centre B

        quadB(1) = quadB(1) + qB(i3,i2,i1)*rjB(1)*rjB(1)
        quadB(2) = quadB(2) + qB(i3,i2,i1)*rjB(1)*rjB(2)
        quadB(3) = quadB(3) + qB(i3,i2,i1)*rjB(1)*rjB(3)
        quadB(4) = quadB(4) + qB(i3,i2,i1)*rjB(2)*rjB(2)
        quadB(5) = quadB(5) + qB(i3,i2,i1)*rjB(2)*rjB(3)
        quadB(6) = quadB(6) + qB(i3,i2,i1)*rjB(3)*rjB(3)

C  octapole moment relative to centre B

	  octaB(1) = octaB(1) + qB(i3,i2,i1)*rjB(1)*rjB(1)*rjB(1)
	  octaB(2) = octaB(2) + qB(i3,i2,i1)*rjB(1)*rjB(1)*rjB(2)
	  octaB(3) = octaB(3) + qB(i3,i2,i1)*rjB(1)*rjB(1)*rjB(3)
	  octaB(4) = octaB(4) + qB(i3,i2,i1)*rjB(1)*rjB(2)*rjB(2)
	  octaB(5) = octaB(5) + qB(i3,i2,i1)*rjB(1)*rjB(2)*rjB(3)
	  octaB(6) = octaB(6) + qB(i3,i2,i1)*rjB(1)*rjB(3)*rjB(3)
	  octaB(7) = octaB(7) + qB(i3,i2,i1)*rjB(2)*rjB(2)*rjB(2)
	  octaB(8) = octaB(8) + qB(i3,i2,i1)*rjB(2)*rjB(2)*rjB(3)
	  octaB(9) = octaB(9) + qB(i3,i2,i1)*rjB(2)*rjB(3)*rjB(3)
	  octaB(10)=octaB(10) + qB(i3,i2,i1)*rjB(3)*rjB(3)*rjB(3)


      end do
      end do
      end do
C  Now convert from e*Bohr to e*Angstroms and normalize.      
      dipoleB = dipoleB * BohrtoAng * volB
      quadB = quadB * BohrtoAng*BohrtoAng * volB
      octaB = octaB * BohrtoAng*BohrtoAng*BohrtoAng * volB

c      write (6,*) sumB
c      write (6,*) dipoleB(1), dipoleB(2), dipoleB(3)
c      write (6,*) quadB(1), quadB(2), quadB(3)
c      write (6,*) quadB(4), quadB(5), quadB(6)
c	write (6,*) octaB(1), octaB(2), octaB(3), octaB(4), octaB(5)
c	write (6,*) octaB(6), octaB(7), octaB(8), octaB(9), octaB(10)
      
	quB(1,1) = quadB(1)
	quB(1,2) = quadB(2)
	quB(1,3) = quadB(3)
	quB(2,1) = quB(1,2)
	quB(2,2) = quadB(4)
	quB(2,3) = quadB(5)
	quB(3,1) = quB(1,3)
	quB(3,2) = quB(2,3)
	quB(3,3) = quadB(6)

	oB(1,1,1) = octaB(1)
	oB(1,1,2) = octaB(2)
	oB(1,1,3) = octaB(3)
	oB(1,2,1) = oB(1,1,2)
	oB(1,2,2) = octaB(4)
	oB(1,2,3) = octaB(5)
	oB(1,3,1) = oB(1,1,3)
	oB(1,3,2) = oB(1,2,3)
	oB(1,3,3) = octaB(6)
	oB(2,1,1) = oB(1,1,2)
	oB(2,1,2) = oB(1,2,2)
	oB(2,1,3) = oB(1,2,3)
	oB(2,2,1) = oB(1,2,2)
	oB(2,2,2) = octaB(7)
	oB(2,2,3) = octaB(8)
	oB(2,3,1) = oB(1,2,3)
	oB(2,3,2) = oB(2,2,3)
	oB(2,3,3) = octaB(9)
	oB(3,1,1) = oB(1,1,3)
	oB(3,1,2) = oB(1,2,3)
	oB(3,1,3) = oB(1,3,3)
	oB(3,2,1) = oB(1,2,3)
	oB(3,2,2) = oB(2,2,3)
	oB(3,2,3) = oB(2,3,3)
	oB(3,3,1) = oB(1,3,3)
	oB(3,3,2) = oB(2,3,3)
	oB(3,3,3) = octaB(10)

C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Calculate dipole-dipole coupling
C
      RAB(1) = cenb(1)-cena(1)
      RAB(2) = cenb(2)-cena(2)
      RAB(3) = cenb(3)-cena(3)
      R = DSQRT(DOT_PRODUCT(RAB,RAB))
      R0 = RAB/R

      R = R*BohrtoAng

C  Convert transition moments from e.Ang into Debye

      dipoleAD = dipoleA*eAngtoDebye
      dipoleBD = dipoleB*eAngtoDebye


C  Dipole-dipole coupling from point dipole model,
C  with separation, R, taken to be between the calculated centers of
C  charge of each transition density.

C  Unit vectors and magnitude defining the intermolecular separation
C  vector about which the Taylor expansion of the interaction potential
C  is taken.  This is assumed to be between the centers of charge of 
C  each transition density.

	muA = DSQRT(DOT_PRODUCT(dipoleA,dipoleA))
      muB = DSQRT(DOT_PRODUCT(dipoleB,dipoleB))
      muAD = DSQRT(DOT_PRODUCT(dipoleAD,dipoleAD))
      muBD = DSQRT(DOT_PRODUCT(dipoleBD,dipoleBD))
            

C  Orientation factor
      dotAB = DOT_PRODUCT(dipoleA,dipoleB)
      dotAR0 = DOT_PRODUCT(dipoleA,R0)
      dotBR0 = DOT_PRODUCT(dipoleB,R0)

      kappa = (dotAB - 3.* dotAR0 * dotBR0)/(muA*muB)
     
      Vdipdip = 1.E10*muA*muB*echarge**2.*kappa*HtoeV/
     $            (4.*JtoHartree*pi*ep0*R**3.)

c	write(*,*) 'vdipdip before converting units:', muA*muB*kappa/R**3
c	write(*,*) 'vdipdip=', (dotAB-3*DotAR0*dotBR0)/R**3
C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
C  Convert coupling to eV and correct for element volume
C      Texact = Texact*volA*volB*1.e10*HtoeV*echarge**2./(4.*pi*ep0*
C     $             BohrtoAng*JtoHartree)
      
      
C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Scale couplings to experimental transition moments
      exp_correctA = expA / muAD
      exp_correctB = expB / muBD
      exp_correct = exp_correctA * exp_correctB
      Vdipdipexp = Vdipdip * exp_correct
C      Texactexp = Texact * exp_correct
C      Tddexp = Tdd * exp_correct
C      Tdqexp = Tdq * exp_correct
C      Tqdexp = Tqd * exp_correct
C      Tdoexp = Tdo * exp_correct
C      Todexp = Tod * exp_correct
C      Tqqexp = Tqq * exp_correct
C      Tabexp = Tab * exp_correct

      dipoleAexp = dipoleA * exp_correctA
      dipoleBexp = dipoleB * exp_correctB
      dipoleADexp = dipoleAexp * eAngtoDebye
      dipoleBDexp = dipoleBexp * eAngtoDebye
      muAexp = muA * exp_correctA
      muBexp = muB * exp_correctB
      muADexp = muAD * exp_correctA
      muBDexp = muBD * exp_correctB

C      quadAexp = quadA * exp_correctA
C      quadBexp = quadB * exp_correctB

c Cherri: calculate d-d, d-q and q-q+d-o couplings.
c First test the d-d result.
	k2= 0.*k2
	k3 = 0.*k3
	k4=0.*k4

c kernal matrix for 2nd order: d-d coupling
	do i=1,3
	  do j=1,3
	    k2(i,j) = 3.*R0(i)*R0(j)
	    IF(i.eq.j) k2(i,j) = k2(i,j)-1.
	    write(*,*) i,j,k2(i,j)
	    k2(i,j) = -k2(i,j)/R**3
	  enddo
	enddo

	testdd=0.d0
	Do i=1,3
	   DO j=1,3
	     testdd = testdd+dipoleA(i)*k2(i,j)*dipoleB(j)
	   enddo
	enddo
	testdd = testdd*1.E10*echarge**2.*HtoeV/(4.*JtoHartree*pi*ep0)	
	write(*,*) 'dipdip coupling:', testdd

c kernal matrix for 3rd order: q-d coupling

	do i=1,3
	  do j=1,3
	    do k=1,3
	      k3(i,j,k) = 15*R0(i)*R0(j)*R0(k)
		if (i.eq.j) k3(i,j,k) = k3(i,j,k) - 3.* R0(k)
		if (i.eq.k) k3(i,j,k) = k3(i,j,k) - 3.* R0(j)
		if (k.eq.j) k3(i,j,k) = k3(i,j,k) - 3.* R0(i)
		k3(i,j,k) = k3(i,j,k)/R**4/6  ! taylor expan factor 
		k3(i,j,k) = k3(i,j,k)*3  ! for dipole-quadrupole there are 3 terms
 	    enddo
	  enddo
	enddo

c calculate Quadrupole(A)-dipole(B) (Vquaddip)
c calculate Quadrupole(B)-dipole(A) (Vdipquad):k3 must be used with a minus sign
	Vquaddip = 0.d0
	Vdipquad = 0.d0
	do i=1,3
	  do j=1,3
	    do k=1,3
		Vquaddip = Vquaddip + QuA(i,j)*dipoleB(k)*k3(i,j,k)
		Vdipquad = Vdipquad - QuB(i,j)*dipoleA(k)*k3(i,j,k)
	    enddo
	  enddo
	enddo
	Vquaddip = Vquaddip*1.E10*echarge**2.*HtoeV/(4.*JtoHartree*pi*ep0)	
	Vdipquad = Vdipquad*1.E10*echarge**2.*HtoeV/(4.*JtoHartree*pi*ep0)	
	write(*,*) 'A quad-B dip coupling:', Vquaddip
	write(*,*) 'B quad-A dip coupling:', Vdipquad

c setup the kernal matrix for the 4th order: k4 (q-q, o-d couplings)
	do i=1,3
	  do j=1,3
	    do k=1,3
	      do l=1,3
	      k4(i,j,k,l) = 105*R0(i)*R0(j)*R0(k)*R0(l)
		if (i.eq.j) k4(i,j,k,l) = k4(i,j,k,l) - 15.* R0(k)* R0(l)
		if (i.eq.k) k4(i,j,k,l) = k4(i,j,k,l) - 15.* R0(j)* R0(l)
		if (k.eq.j) k4(i,j,k,l) = k4(i,j,k,l) - 15.* R0(i)* R0(l)
		if (i.eq.l) k4(i,j,k,l) = k4(i,j,k,l) - 15.* R0(j)* R0(k)
		if (j.eq.l) k4(i,j,k,l) = k4(i,j,k,l) - 15.* R0(i)* R0(k)
		if (k.eq.l) k4(i,j,k,l) = k4(i,j,k,l) - 15.* R0(j)* R0(i)
		if((i.eq.j).and.(k.eq.l)) k4(i,j,k,l) = k4(i,j,k,l) + 3.
		if((i.eq.k).and.(j.eq.l)) k4(i,j,k,l) = k4(i,j,k,l) + 3.
		if((i.eq.l).and.(k.eq.j)) k4(i,j,k,l) = k4(i,j,k,l) + 3.
		k4(i,j,k,l) = k4(i,j,k,l)/R**5/24  ! taylor expan factor 4!
		enddo
 	    enddo
	  enddo
	enddo

c calculate Quadrupole(A)-quadpole(B) (Vquadquad) k4*6
c calculate Octapole(A)-dipole(B) (Voctdip): k4 *(-4)
c calculate Octpole(B)-dipole(A) (Vdipoct):k4 *(-4)
	Vquadquad = 0.d0
	Voctdip = 0.d0
	Vdipoct = 0.d0
	do i=1,3
	  do j=1,3
	    do k=1,3
		do l=1,3
		Vquadquad = Vquadquad + QuA(i,j)*QuB(k,l)*k4(i,j,k,l)*6.
		Voctdip = Voctdip + oA(i,j,k)*dipoleB(l)*k4(i,j,k,l)*(-4.)
		Vdipoct = Vdipoct + dipoleA(i)*oB(j,k,l)*k4(i,j,k,l)*(-4.)
		enddo
	    enddo
	  enddo
	enddo
	Vquadquad = Vquadquad*1.E10*echarge**2.*HtoeV/(4.*JtoHartree*pi*ep0)	
	Voctdip = Voctdip*1.E10*echarge**2.*HtoeV/(4.*JtoHartree*pi*ep0)	
	Vdipoct = Vdipoct*1.E10*echarge**2.*HtoeV/(4.*JtoHartree*pi*ep0)	
	write(*,*) 'A quad-B quad coupling:', Vquadquad
	write(*,*) 'A oct-B dip coupling:', Voctdip
	write(*,*) 'B oct-A dip coupling:', Vdipoct
C
C  Stop the clock and convert from seconds to minutes

C  The two lines with ETIME are for UNIX
C      CALL ETIME(timearray)
C      calctime = ( timearray(1) + timearray(2) ) / 60.

C  The two lines with TIMEF are for Win95 MS Fortran 
C      calctime = TIMEF()
C      calctime = calctime / 60.
C  Deallocate the transition density arrays
      deallocate (qA, qB)
	
	Tdd = Vdipdip
	Tqd = Vquaddip
	Tdq = Vdipquad
	Tod = Voctdip
	Tdo = Vdipoct
	Tqq = Vquadquad
      Tab = Tdd + Tdq + Tqd + Tdo + Tod + Tqq

C
C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Output the results (couplings and cartesian coordinates)
      open (unit=12, file='coupling.out', err=999)

C  Monomer details (to double-check input)
      write (12,*) ' '
      write (12,*) 'Monomer A:'
      write (12,100) titleA
      write (12,*) ' SumA = ', sumA,' e in',NelmtsA,' elements'
      write (12,96) cena(1), cena(2), cena(3)
      write (12,97) dipoleA(1),dipoleA(2),dipoleA(3)
      write (12,*) '         |Dipole| = ', muA, ' e.Ang'
      write (12,98) dipoleAD(1),dipoleAD(2),dipoleAD(3)
      write (12,*) '         |Dipole| = ', muAD, ' Debye'
      write (12,99) (quadA(i), i=1,6)
      write (12,*) 'Scaled to experimental dipole moment of ',expA,
     &               ' Debye:'
      write (12,102) dipoleAexp(1),dipoleAexp(2),dipoleAexp(3)
      write (12,*) '         |Dipole| = ', muAexp, ' e.Ang'
      write (12,103) dipoleADexp(1),dipoleADexp(2),dipoleADexp(3)
      write (12,*) '         |Dipole| = ', muADexp, ' Debye'
C      write (12,104) (quadAexp(i), i=1,6)
      write (12,*) ' '
      write (12,*) 'Monomer B:'
      write (12,101) titleB
      write (12,*) ' SumB = ', sumB,' e in',NelmtsB,' elements'
      write (12,96) cenb(1), cenb(2), cenb(3)
      write (12,97) dipoleB(1),dipoleB(2),dipoleB(3)
      write (12,*) '         |Dipole| = ', muB, ' e.Ang'
      write (12,98) dipoleBD(1),dipoleBD(2),dipoleBD(3)
      write (12,*) '         |Dipole| = ', muBD, ' Debye'
      write (12,99) (quadB(i), i=1,6)
      write (12,*) 'Scaled to experimental dipole moment of ',expB,
     &               ' Debye:'
      write (12,102) dipoleBexp(1),dipoleBexp(2),dipoleBexp(3)
      write (12,*) '         |Dipole| = ', muBexp, ' e.Ang'
      write (12,103) dipoleBDexp(1),dipoleBDexp(2),dipoleBDexp(3)
      write (12,*) '         |Dipole| = ', muBDexp, ' Debye'
C      write (12,104) (quadBexp(i), i=1,6)
      write (12,*) ' '
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

C  Output couplings
      write (12,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
C      write (12,139) calctime
      write (12,*) ' Coulombic couplings'
      write (12,141) R, kappa
      write (12,150) Vdipdip, Vdipdip*eVtocm
C      write (12,155) Texact, Texact*eVtocm
      write (12,*) ' Tab = Tdd + Tdq + Tqd + Tdo + Tod + Tqq:'
      write (12,160) Tab, Tdd, Tdq, Tqd, Tdo, Tod, Tqq
      write (12,*) ' '
      write (12,*) ' Couplings after scaling to experimental transition
     & moments'
      write (12,170) Vdipdipexp, Vdipdipexp*eVtocm
C      write (12,175) Texactexp, Texactexp*eVtocm
C      write (12,180) Tabexp,Tddexp,Tdqexp,Tqdexp,Tdoexp,Todexp,Tqqexp
C      write (12,*) ' Overlapping charge? ', irflag
      write (12,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
C      write (12,*) ' Diverging terms? ', idflag
 139  format ('  Elapsed time of calculation ', F9.2, ' minutes')
 141  format ('  R = ',F6.2, ' Angstrom'/'  Orientation factor = ',F6.3)
 150  format (/'  Ideal dipole-dipole: Vdipdip = ',E12.5,' eV'/
     $            31X,'= ',E12.5,' cm-1'/)
 155  format ('  Complete interaction: Texact = ',E12.5,' eV'/
     $            31X,'= ',E12.5,' cm-1'/)
 160  format ('  Multipole expansion (scaled to exp): Tab = '
     $           , E12.5,' eV'/
     $           30X,'   { Tdip-dip =',E12.5,' eV}'/
     $           30X,'  { Tdip-quad =',E12.5,' eV}'/
     $           30X,'  { Tquad-dip =',E12.5,' eV}'/
     $           30X,'   { Tdip-oct =',E12.5,' eV}'/
     $           30X,'   { Toct-dip =',E12.5,' eV}'/
     $           30X,' { Tquad-quad =',E12.5,' eV}'/)
 170  format ('  Ideal dip-dip scaled to experiment:  ',E12.5,' eV'/
     &            37x,'= ',E12.5,' cm-1'/)
 175  format ('  Complete interaction scaled to exp:  ',E12.5,' eV'/
     &            37x,'= ',E12.5,' cm-1'/)
 180  format ('  Multipole expansion scaled to exp: Tabexp = '
     $           , E12.5,' eV'/
     $           28X,'   { Tdip-dipexp =',E12.5,' eV}'/
     $           28X,'  { Tdip-quadexp =',E12.5,' eV}'/
     $           28X,'  { Tquad-dipexp =',E12.5,' eV}'/
     $           28X,'   { Tdip-octexp =',E12.5,' eV}'/
     $           28X,'   { Toct-dipexp =',E12.5,' eV}'/
     $           28X,' { Tquad-quadexp =',E12.5,' eV}'/)



C  Cartesian coordinates (convert to PDB file using "newzmat -ixyz -opdb")
c      write (12,*) ' Cartesian coordinates for dimer'
c      write (12,200) ( IA(i),pxA(i),pyA(i),pzA(i), i=1,NatomsA )
c      write (12,200) ( IB(i),pxB(i),pyB(i),pzB(i), i=1,NatomsB )
c 200  format(2X,I5,3F12.6)

      close (unit=12)
C  Deallocate atomic coordinates and such
      deallocate (pxA, pyA, pzA, pxB, pyB, pzB, atnoA, atnoB, IA, IB)

      goto 1000

C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Error traps
 990  write (6,*) 'Error opening expAB file.'
      write (6,*) 'Continuing without use of experimental transition
     & moments.'
      goto 3
 999  write (6,*) ' Ahhhhhhh!!!'

 1000 continue
      end

