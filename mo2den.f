c Copied from coupling5.f, this program is intend to convert MO's to transitio
c densities.
c
C
C  The compile line on SGI:  f90 -OPT:fold_arith_limit=1500 -O3 -o coupling5 coupling5.f
C  The compile line on HP:   f90 +O4 -o coupling5 coupling5.f
C
C  *Note on HP-UX use +U77 flag to get timer functions to work
C
C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit none
C  Input variables
      integer NatomsA, NatomsB, NelmtsA, NelmtsB
      integer NA1, NA2, NA3, NB1, NB2, NB3, homoA, homoB
      real*8 xa0, ya0, za0, xa1, ya1, za1, xa2, ya2, za2, xa3, ya3, za3
      real*8 xb0, yb0, zb0, xb1, yb1, zb1, xb2, yb2, zb2, xb3, yb3, zb3
      real*8 cena(3), cenb(3), cenaang(3), cenbang(3), center(3)
      real*8 expA, expB, abssuma, abssumb, density(200)
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
C      real*8 Vdd, Vdq, Vqd, Vqq, Vdo, Vod
C      real*8 deltaab, deltaaa, deltabb, Ra, Rb
      real*8 dotAB, dotAR0, dotBR0
 
C  Output
      real*8 dipoleA(3), dipoleB(3), dipoleAexp(3), dipoleBexp(3) 
      real*8 dipoleAD(3), dipoleBD(3)
C      real*8 quadA(6), quadB(6), quadAexp(6), quadBexp(6)
      real*8 dipoleADexp(3), dipoleBDexp(3)
      real*8 kappa, muA, muAD, muB, muBD
      real*8 muAexp, muBexp, muBDexp, muADexp
      real*8 Texact, Texactexp
C      real*8 Tdd, Tdq, Tqd, Tqq, Tdo, Tod
C      real*8 TabTddexp, Tdqexp, Tqdexp, Tqqexp
C      real*8 Tdoexp, Todexp, Tabexp
      real*8 Vdipdip, Vdipdipexp
	real*8 sqr2

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

C      calctime = TIMEF()      !This line for Win32

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
      qA = 0.            

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
      qB = 0.

      do i1=1,NB1
        do i2=1,NB2
          read (14,'(6E13.5)') (qB(i3,i2,i1), i3=1,NB3)          
        end do
      end do

      close (unit=14)

	IF(((NA1-NB1)**2+(NA2-NB2)**2+(NA3-NB3)**2).GT.0) THEN
	  write(*,*) 'files are not the same!'
	  STOP
	endif
	
      OPEN(12, file='cubeAB', status='unknown')
	write(12,'(A)') title A
	write(12,*)
      write(12,'(I5,4F12.6)') NatomsA, xa0, ya0, za0
C  Number of increments and size in each direction
      write(12,'(I5,4F12.6)') NA1, xa1, ya1, za1
      write(12,'(I5,4F12.6)') NA2, xa2, ya2, za2
      write(12,'(I5,4F12.6)') NA3, xa3, ya3, za3
C  Atomic number, charge and atomic positions
      write (12,'(I5,4F12.6)') (IA(i),atnoA(i),pxA(i),pyA(i)
     $   ,pzA(i), i=1,NatomsA)

      VolA = DSqrt((xa1*xa1+ya1*ya1+za1*za1)*
     $              (xa2*xa2+ya2*ya2+za2*za2)*
     $              (xa3*xa3+ya3*ya3+za3*za3))

      VolB = DSqrt((xb1*xb1+yb1*yb1+zb1*zb1)*
     $              (xb2*xb2+yb2*yb2+zb2*zb2)*
     $              (xb3*xb3+yb3*yb3+zb3*zb3))
	sqr2=DSQRT(2.d0)
      do i1=1,NA1
        do i2=1,NA2
	    do i3 = 1,NA3
	      density(i3) = qA(i3,i2,i1)*qB(i3,i2,i1)*sqr2 !*VolA
c	      write(*,*) density(i3)
	    enddo
          write (12,'(6E13.5)') (density(i), i=1,NA3)          
        end do
      end do

	CLOSE(12)
	go to 1000
C  Error traps
 990  write (6,*) 'Error opening expAB file.'
c      write (6,*) 'Continuing without use of experimental transition
c     & moments.'
c      goto 3
 999  write (6,*) ' Ahhhhhhh!!!'

 1000 continue
      end

