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
C    for atomic number and coordinates
      real*8, allocatable, dimension(:,:,:) :: qA, qB      
      real*8, allocatable, dimension(:) :: pxA, pyA, pzA, pxB, pyB, pzB,
     & atnoA, atnoB
      integer, allocatable, dimension(:) :: IA, IB
	CHARACTER*6 DUMMY6
	CHARACTER*4 DUMMY4
	CHARACTER*2, allocatable, dimension(:) :: ATOM
	character*2 GROUP
	INTEGER iDUMMY, IMOL
	real one

C  Distance vectors
      real*8 R0(3),R,RAB(3)
      real*8 riA(3), rjB(3), rij, Rcutoff
 
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

C  Open and process cubeA
  3   open (unit=16, file='cubeA', err=999)      
      read (16,5) titleA
      read (16,5) dummy
c      write (6,5) titleA
  5   format(A80)
  
C  Number of atoms and origin
      read (16,*) NatomsA, xa0, ya0, za0
C  Number of increments and size in each direction
      read (16,*) NA1, xa1, ya1, za1
      read (16,*) NA2, xa2, ya2, za2
      read (16,*) NA3, xa3, ya3, za3
c      write (6,*) NatomsA, NA1, NA2, NA3
C  Allocate arrays for atomic coords, etc.
      allocate( IA(NatomsA), atnoA(NatomsA), pxA(NatomsA), pyA(NatomsA),
     &  pzA(NatomsA),ATOM(NATOMSA))
C  Atomic number, charge and atomic positions
      read (16,'(I5,4F12.6)') (IA(i),atnoA(i),pxA(i),pyA(i)
     $   ,pzA(i), i=1,NatomsA)

C  Read in the density cube 
C  (Note that units are es and Bohr.)
c      allocate( qA(NA3, NA2, NA1) )     
c      qA = 0.            

c      do i1=1,NA1
c        do i2=1,NA2
c          read (16,'(6E13.5)') (qA(i3,i2,i1), i3=1,NA3)          
c        end do
c      end do
      close (unit=16)

	do i=1,natomsa
	  IF(ia(i).eq.1) then
	    atom(i)=' H'
	  ELSEIF(ia(i).eq.6) then
	    atom(i)=' C'
	  ELSEIF(ia(i).eq.7) then
	    atom(i)=' N'
	  ELSEIF(ia(i).eq.8) then
	    atom(i)=' O'
	  ELSEIF(ia(i).eq.12) then
	    atom(i)='MG'
	  ELSE
	    atom(i)=' X'
	  ENDIF
	enddo
9	format(A6,I5,1x,A2,I2,A4,A2,I4,A4,3F8.3,2f6.2,10X,A2,2X)  
	dummy6='HETATM'
	GROUP=' P'
	IMOL=11
	dummy4='    '
	one=1.0
	do i=1,natomsA
	  write(*,9) DUMMY6,I,ATOM(I),i,DUMMY4,GROUP,IMOL,DUMMY4,
     .     PXA(i)*bohrtoang,PYA(I)*bohrtoang,PZA(I)*bohrtoang,	
     .     one,one,ATOM(I)
	enddo
			
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


