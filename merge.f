c A rather incorrect way to reduce the number of cubes for testing
c 12/21/00 Cherri Hsu
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

C  Dynamically allocated array for transition density cube elements and
C    for atomic number and coordinates (fort90)
      real*8, allocatable, dimension(:,:,:) :: qA, qB      
      real*8, allocatable, dimension(:) :: pxA, pyA, pzA, pxB, pyB, pzB,
     & atnoA, atnoB
      integer, allocatable, dimension(:) :: IA, IB

C  Statically allocated array for transition density cube elements and
C    for atomic number and coordinates (fort77)
C      real*8 qA(51,62,171), qB(94,65,92)      
C      real*8 pxA(123), pyA(123), pzA(123), pxB(81), pyB(81), pzB(81),
C     & atnoA(123), atnoB(81)
C      integer IA(123), IB(81)

C  Distance vectors
      real*8 R0(3),R,RAB(3)
      real*8 riA(3), rjB(3), rij, riAmag, rjBmag

C  Work
      real*8 sumA, sumB, volA, volB, correctA, correctB
      real*8 exp_correctA, exp_correctB, exp_correct
      real*8 Vdd, Vdq, Vqd, Vqq, Vdo, Vod
      real*8 deltaab, deltaaa, deltabb, Ra, Rb
      real*8 dotAB, dotAR0, dotBR0
 
C  Output
      real*8 dipoleA(3), dipoleB(3), quadA(6), quadB(6)
      real*8 dipoleAD(3), dipoleBD(3)
      real*8 dipoleAexp(3), dipoleBexp(3), quadAexp(6), quadBexp(6)
      real*8 dipoleADexp(3), dipoleBDexp(3)
      real*8 kappa, muA, muAD, muB, muBD
      real*8 muAexp, muBexp, muBDexp, muADexp
      real*8 Texact, Tdd, Tdq, Tqd, Tqq, Tdo, Tod, Tab
      real*8 Texactexp, Tddexp, Tdqexp, Tqdexp, Tqqexp
      real*8 Tdoexp, Todexp, Tabexp
      real*8 Vdipdip, Vdipdipexp
C  Timer variables
      real*4 timearray(2), calctime            !This line for Unix
C      real*8 calctime, TIMEF                  !This line for Win32

C  Constants
      real*8 :: BohrtoAng=0.529177249
      real*8 :: echarge=1.602188E-19
      real*8 :: ep0=8.8542e-12
      real*8 :: pi=3.1415926535
      real*8 :: HtoeV=27.2116
      real*8 :: JtoHartree=4.3597482e-18
      real*8 :: eVtocm=8065.
      real*8 :: Rcutoff=0.1
      real*8 :: eAngtoDebye=4.803

C  Loop variables, flags
      integer i,i1,i2,i3,j1,j2,j3
      integer :: nmerge=7
	integer NA1m, NA2m, NA3m
      real*8, allocatable, dimension(:) :: qmerge
	real*8 sum


C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C      expA = 13.      !RG1
C      expB = 3.29     !Bchl Qx
      expA = 1.
      expB = 1.
      sumA = 0.
      center(1) = 0.
      center(2) = 0.
      center(3) = 0.
      muA = 0.
      muB = 0.
      sumB = 0.

C      calctime = TIMEF()      !This line for Win32

C  The transition density cube files (from Gaussian)
C  are called cubeA and cubeB.  The file with experimental 
C  transition dipole moments is called expAB.
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

C  Read in the density cube and calculate charge and dipole moment
C  (Note that units are es and Bohr.)
      allocate( qA(NA3, NA2, NA1) )     

      do i1=1,NA1
        do i2=1,NA2
          read (16,'(6E13.5)') (qA(i3,i2,i1), i3=1,NA3)          
        end do
      end do

      close (unit=16)
      
	write(*,*) 'nmerge=', nmerge
	NA1m = NA1/nmerge
	write(*,*) 'NA1=', NA1, 'NA1m=', NA1m
	NA2m = NA2/nmerge
	write(*,*) 'NA2=', NA2, 'NA2m=', NA2m
	NA3m = NA3/nmerge
	write(*,*) 'NA3=', NA3, 'NA3m=', NA3m
      allocate( qmerge(NA3m) )     

	open (17, file='cubeA.merge')
	write(17,5) titleA
	write(17,5) dummy
100	format(i5,3f12.6)
      write(17,100) NatomsA, xa0, ya0, za0
      write(17,100) NA1m, xa1*nmerge, ya1*nmerge, za1*nmerge
      write(17,100) NA2m, xa2*nmerge, ya2*nmerge, za2*nmerge
      write(17,100) NA3m, xa3*nmerge, ya3*nmerge, za3*nmerge
      write(17,'(I5,4F12.6)') (IA(i),atnoA(i),pxA(i),pyA(i)
     $   ,pzA(i), i=1,NatomsA)

	  do i1 = 1,NA1m
	    Do i2 = 1, NA2m
		  do i3=1, NA3m
		    sum=0.d0
		    do j1=nmerge*(i1-1)+1,nmerge*i1
			  do j2=nmerge*(i2-1)+1,nmerge*i2
			    do j3=nmerge*(i3-1)+1,nmerge*i3
				  sum = sum + qA(j3,j2,j1)
				enddo !j3
			  enddo !j2
			enddo  !j1
			qmerge(i3) = sum/nmerge**3
		  enddo !i3
          write(17,'(6E13.5)') (qmerge(i3), i3=1,NA3m)          
	    enddo
	  enddo

	  write(17,*)
	  close(17)

	  goto 1000

 999  write (6,*) ' Ahhhhhhh!!!'

 1000 continue
      end
