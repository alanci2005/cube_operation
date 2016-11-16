C  transform.f   (for fortran 90)
C  Brent Krueger Nov 1997 The University of Chicago -- Fleming Group
C
C
C  This program is for the re-orientation (transformation) of a Gaussian
C    transition density cube.  A cube can be transformed from Gaussian
C    internal coordinates to some molecular coordinates and/or transformed
C    by any arbitrary translation or rotation.  The program works by 
C    finding three atoms in the cube's coordinates and transforming 
C    onto three reference atomic positions specified by the user.  Atomic
C    coordinates for many of the chromophores from the Rps. acidophila 
C    structure are built in, so use of these coordinates is simple.
C
C  The program always provides a default response to any query.  The 
C    default is chosen simply by pressing enter.  So, if a transition
C    density cube for BChl is to be read in and transformed onto the
C    location of BChl #5 from chain P of the LH2 crystal structure, 
C    the user simply needs to name the input cube 'cube.transform.start'
C    run the program, and enter 2 for BChl.  The program automatically
C    reads the correct atomic coordinates from the input file and uses
C    them as the default, which the user accepts by striking enter for 
C    queries concerning characteristic atoms of the input cube.  When
C    asked for the destination, the user types 5 to indicate BChl #5P
C    and then more enter keys to accept the default responses the rest
C    of the way through the program.  Not using the default entries
C    allows the user to translate or rotate the cube and to specify
C    any possible coordinates for the three reference atoms.
C	
C
C  This version of the code runs on both the HP and SGI.
C
C  Compile line on the SGI: f90 -o transform transform.f
C  Compile line on the HP:  f90 -o transform transform.f
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C  Declarations:
C
      implicit none
C
C  Constants
      real*8 :: BohrtoAng = 0.529177249
      real*8 xhat(3), zhat(3)
      data xhat/1.,0.,0./
      data zhat/0.,0.,1./
C    The following molecule coordinates are taken from lh2_non.pdb
C      complete crystal structure of acidophila B800-B850 LH2
C    Chain P chromophores
      real*8 bchlp1_1(3), bchlp1_2(3), bchlp1_3(3)
      real*8 bchlp2_1(3), bchlp2_2(3), bchlp2_3(3)
      real*8 bchlp3_1(3), bchlp3_2(3), bchlp3_3(3)
      real*8 bchlp4_1(3), bchlp4_2(3), bchlp4_3(3)
      real*8 bchlp5_1(3), bchlp5_2(3), bchlp5_3(3)
      real*8 bchlp6_1(3), bchlp6_2(3), bchlp6_3(3)
      real*8 bchlp7_1(3), bchlp7_2(3), bchlp7_3(3)
      real*8 bchlp8_1(3), bchlp8_2(3), bchlp8_3(3)
      real*8 bchlp9_1(3), bchlp9_2(3), bchlp9_3(3)
      real*8 rgp10_1(3), rgp10_2(3), rgp10_3(3), rgp10_c(3)
      real*8 rgp11_1(3), rgp11_2(3), rgp11_3(3), rgp11_c(3)
      real*8 rgp12_1(3), rgp12_2(3), rgp12_3(3), rgp12_c(3)
      data bchlp1_1/25.468,-5.207,18.174/		!Mg
      data bchlp1_2/25.430,-5.481,20.211/		!NA
      data bchlp1_3/24.485,-6.934,17.893/		!NB
      data bchlp2_1/26.785,4.227,18.132/
      data bchlp2_2/26.476,4.597,20.134/
      data bchlp2_3/27.253,6.140,17.825/
      data bchlp3_1/22.857,12.380,18.174/
      data bchlp3_2/23.004,12.146,20.211/
      data bchlp3_3/23.214,10.426,17.893/
      data bchlp4_1/17.802,20.454,18.132/
      data bchlp4_2/17.327,20.538,20.134/  
      data bchlp4_3/16.930,22.219,17.825/    
      data bchlp5_1/9.551,24.175,18.174/
      data bchlp5_2/9.814,24.09,20.211/
      data bchlp5_3/11.081,22.907,17.893/
      data bchlp6_1/.488,27.11,18.132/
      data bchlp6_2/.07,26.87,20.134/
      data bchlp6_3/-1.315,27.902,17.825/
      data bchlp7_1/28.764,11.644,34.748/
      data bchlp7_2/27.393,10.205,34.350/
      data bchlp7_3/30.263,10.381,34.212/
      data bchlp8_1/14.549,27.407,34.748/
      data bchlp8_2/14.424,25.423,34.350/
      data bchlp8_3/16.509,27.403,34.212/
      data bchlp9_1/-6.474,30.346,34.748/
      data bchlp9_2/-5.294,28.746,34.350/
      data bchlp9_3/-4.969,31.602,34.212/
      data rgp10_1/25.583,3.086,31.491/       !C16
      data rgp10_2/26.079,1.620,33.337/       !C14
      data rgp10_3/25.860,1.801,31.971/	  !C15
      data rgp10_c/24.871,2.578,32.069/	  !avg of C5-C26
      data rgp11_1/17.614,18.807,31.491/
      data rgp11_2/18.937,18.003,33.337/
      data rgp11_3/18.652,18.000,31.971/
      data rgp11_c/17.395,17.960,32.069/
      data rgp12_1/1.403,25.728,31.491/
      data rgp12_2/2.933,25.962,33.337/
      data rgp12_3/2.717,25.778,31.971/
      data rgp12_c/1.779,24.938,32.069/
C    Chain R chromophores
      real*8 bchlr1_1(3), bchlr1_2(3), bchlr1_3(3)
      real*8 bchlr2_1(3), bchlr2_2(3), bchlr2_3(3)
      real*8 bchlr3_1(3), bchlr3_2(3), bchlr3_3(3)
      real*8 bchlr4_1(3), bchlr4_2(3), bchlr4_3(3)
      real*8 bchlr5_1(3), bchlr5_2(3), bchlr5_3(3)
      real*8 bchlr6_1(3), bchlr6_2(3), bchlr6_3(3)
      real*8 bchlr9_1(3), bchlr9_2(3), bchlr9_3(3)
      real*8 rgr12_1(3), rgr12_2(3), rgr12_3(3), rgr12_c(3)
      data bchlr1_1/-17.243,-19.452,18.174/
      data bchlr1_2/-17.462,-19.282,20.211/
      data bchlr1_3/-18.248,-17.738,17.893/
      data bchlr2_1/-9.732,-25.310,18.132/
      data bchlr2_2/-9.257,-25.227,20.134/
      data bchlr2_3/-8.309,-26.672,17.825/
      data bchlr3_1/-0.707,-25.985,18.174/
      data bchlr3_2/-0.983,-25.995,20.211/
      data bchlr3_3/-2.578,-25.317,17.893/
      data bchlr4_1/8.813,-25.644,18.132/
      data bchlr4_2/9.123,-25.275,20.134/
      data bchlr4_3/10.777,-25.771,17.825/
      data bchlr5_1/16.161,-20.359,18.174/
      data bchlr5_2/15.956,-20.544,20.211/
      data bchlr5_3/14.298,-21.050,17.893/
      data bchlr6_1/23.134,-13.829,18.107/
      data bchlr6_2/23.235,-13.495,20.134/
      data bchlr6_3/24.821,-12.812,17.825/
      data bchlr9_1/29.536,-9.532,34.522/
      data bchlr9_2/27.542,-9.788,34.350/
      data bchlr9_3/29.853,-11.497,34.212/
      data rgr12_1/21.580,-14.079,31.491/
      data rgr12_2/21.017,-15.521,33.337/
      data rgr12_3/20.966,-15.242,31.971/
      data rgr12_c/21.273,-14.661,31.731/
C    Chain Q chromophores
      real*8 bchlq1_1(3), bchlq1_2(3), bchlq1_3(3)
      real*8 bchlq2_1(3), bchlq2_2(3), bchlq2_3(3)
      real*8 bchlq3_1(3), bchlq3_2(3), bchlq3_3(3)
      real*8 bchlq4_1(3), bchlq4_2(3), bchlq4_3(3)
      real*8 bchlq5_1(3), bchlq5_2(3), bchlq5_3(3)
      real*8 bchlq6_1(3), bchlq6_2(3), bchlq6_3(3)
      data bchlq1_1/-8.224,24.660,18.174/
      data bchlq1_2/-7.968,24.764,20.211/
      data bchlq1_3/-6.237,24.672,17.893/
      data bchlq2_1/-17.053,21.083,18.132/
      data bchlq2_2/-17.219,20.630,20.134/
      data bchlq2_3/-18.944,20.532,17.825/
      data bchlq3_1/-22.150,13.605,18.174/
      data bchlq3_2/-22.021,13.849,20.211/
      data bchlq3_3/-20.636,14.891,17.893/
      data bchlq4_1/-26.615,5.190,18.132/
      data bchlq4_2/-26.450,4.737,20.134/
      data bchlq4_3/-27.707,3.552,17.825/
      data bchlq5_1/-25.712,-3.816,18.174/
      data bchlq5_2/-25.769,-3.546,20.211/
      data bchlq5_3/-25.378,-1.857,17.893/
      data bchlq6_1/-23.722,-13.132,18.132/
      data bchlq6_2/-23.305,-13.374,20.134/
      data bchlq6_3/-23.506,-15.090,17.825/

C  Loop variables, flags
      integer i, j, k, menu
      integer rxyflag, rpxyflag, ftype, unitflag
      integer r2xyflag, rp2xyflag, shiftflag, refflag
C  Cube characteristics
      integer Natoms, N1, N2, N3
      integer, allocatable, dimension(:) :: N
      real*8 cubeorig(3), step1(3), step2(3), step3(3)
      real*8, allocatable, dimension(:) :: atno
      real*8, allocatable, dimension(:,:) :: atpo
      real*8, allocatable, dimension(:,:,:) :: q
      character*80 title, dummy
C  Vectors and matrices	
      real*8 r(3), rp(3), rhat(3), rphat(3), rtmp(3), rptmp(3), rnew(3)
      real*8 r2(3), rp2(3), r2mag, rp2mag, r2new(3), rp2new(3)
      real*8 rphi, rpphi, rtheta, rptheta, rpsi, rppsi, rmag, rpmag
      real*8 phi, theta, psi, rpnew(3), rnew2(3), r2new2(3)
      real*8 charatom1(3), charatom2(3), charatom3(3)
      real*8 refatom1(3), refatom2(3), refatom3(3), transmat(3)
      real*8 rotmat(3,3), phimat(3,3), thetamat(3,3), psimat(3,3)
      real*8 shiftmat(3)
C  Initialize parameters
      menu = 0
      rxyflag = 0
      rpxyflag = 0
      r2xyflag = 0
      rp2xyflag = 0
      ftype = 0
      unitflag = 0
      shiftflag = 0
      refflag = 0
      do i = 1,3
        charatom1(i) = 0.
        charatom2(i) = 0.
        refatom1(i) = 0.
        refatom2(i) = 0.
        do j = 1,3
          rotmat(i,j) = 0.
          phimat(i,j) = 0.
          thetamat(i,j) = 0.
          psimat(i,j) = 0.
        end do
      end do
C
C  Set up menu
C
      write(6,*)'This program will transform a Gaussian transition'
      write(6,*)' density cube and its assiciated molecular coordinates'
      write(6,*)' from their coordinates used for the density calc'
      write(6,*)' to the coordinates of a reference molecule '
      write(6,*)' entered by the user.  In this way a Gaussian'
      write(6,*)' transition density cube may be calculated for one '
      write(6,*)' molecule and applied the same molecule at any other'
      write(6,*)' position.'
      write(6,*)' '
      write(6,*)'The user must enter the coordinates for three'
      write(6,*)' characteristic atoms from the Gaussian cube output'
      write(6,*)' file, and positions of the matching three atoms from'
      write(6,*)' the reference crystal structure which is the desired'
      write(6,*)' final cube position.'
      write(6,*)' '
  100 write(6,*)'The cube data is assumed to be in a file named:'
      write(6,*)'       cube.transform.start'
      write(6,*)'The transformed cube will be written to a file named:'
      write(6,*)'       cube.transform.done'
      write(6,*)' '
      write(6,*)'Note that the program always provides a default '
      write(6,*)' whenever it queries the user.  The default can '
      write(6,*)' be chosen simply by pressing enter.'
      write(6,*)' '
      write(6,*)'Input cube is of type:'
      write(6,*)'  (0) Other'
      write(6,*)'  (1) Rhodopin Glucoside'
      write(6,*)'  (2) Bacteriochlorophyll'
      call inin('      Enter 0, 1, or 2/?',ftype)
      if ((ftype < 0) .or. (ftype > 2)) goto 100
C
C  Read in the cube datafile.
C
      open (unit=16, file='cube.transform.start')
      read (16,105) title
      read (16,105) dummy
      write (6,*) title
  105 format (A80)
C  Number of atoms and origin
      read (16,*) Natoms, cubeorig(1), cubeorig(2), cubeorig(3)
C  Number of increments and size in each direction
      read (16,110) N1, step1(1), step1(2), step1(3)
      read (16,110) N2, step2(1), step2(2), step2(3)
      read (16,110) N3, step3(1), step3(2), step3(3)
  110 format(I5, 3F12.6)
      write (6,*) Natoms, N1, N2, N3
	      
C  Atomic number, charge and atomic positions (leaving a blank 
C    in which we will place a dummy center atom for the carot)
      if( ftype == 1 ) then
        Natoms = Natoms + 1
        allocate( N(Natoms), atno(Natoms), atpo(3,Natoms) )
        N(1) = 0
        atno(1) = 0.
        atpo(1,1) = 0.
        atpo(2,1) = 0.
        atpo(3,1) = 0.
        read (16,120)(N(i),atno(i),atpo(1,i),atpo(2,i),
     &   atpo(3,i), i = 2,Natoms)
      else
        allocate( N(Natoms), atno(Natoms), atpo(3,Natoms) ) 
        read (16,120)(N(i),atno(i),atpo(1,i),atpo(2,i),
     &   atpo(3,i), i=1,Natoms)
      end if
  120 format(I5, 4F12.6)
C 
C  Read in the density cube and calculate charge and dipole moment
C  (Note that units are es and Bohr.)
      allocate( q(N3, N2, N1) )
      do i=1,N1
        do j=1,N2
          read (16, '(6E13.5)') (q(k,j,i), k=1,N3)          
        end do
      end do
      close (unit=16)
C
C  Prompt for positions of cube's characteristic atoms.
C
C    First, we can set charatom positions to be the expected atomic
C      positions from the cube files
      if (ftype == 1) then  !rhodopin glucoside
        do i = 1,3
          charatom1(i) = atpo(i,14)  !C16
          charatom2(i) = atpo(i,13)  !C14
          charatom3(i) = atpo(i,31)  !C15
        end do
      else if (ftype == 2) then  !bacteriochlorophyll
        do i = 1,3
          charatom1(i) = atpo(i,1)  !Mg
          charatom2(i) = atpo(i,2)  !NA
          charatom3(i) = atpo(i,3)  !NB
        end do
      end if
      write(6,*)''
      write(6,*)'Enter cubes characteristic atom positions(in Bohr).'
      write(6,*)''
      write(6,*)'First characteristic atom defines the origin.'
      write(6,*)'  For RG this is C16, For Bchl this is Mg.'
      call indp('  Enter x/?', charatom1(1) )
      call indp('  Enter y/?', charatom1(2) )
      call indp('  Enter z/?', charatom1(3) )
      write(6,*)'Second characteristic atom defines the x axis.'
      write(6,*)'  For RG this is C14, For Bchl this is NA.'
      call indp('  Enter x/?', charatom2(1) )
      call indp('  Enter y/?', charatom2(2) )
      call indp('  Enter z/?', charatom2(3) )
      write(6,*)'Third characteristic atom defines the y axis.'
      write(6,*)'  For RG this is C15, For Bchl this is NB.'
      call indp('  Enter x/?', charatom3(1) )
      call indp('  Enter y/?', charatom3(2) )
      call indp('  Enter z/?', charatom3(3) )
C
C  Prompt user for positions of characteristic atoms of reference object.
C
      write(6,*)' '
      write(6,*)'Enter coordinates for the characteristic atoms for the'
      write(6,*)'  reference molecule. (Bohr or Angstroms)'
      write (6,*) ' '
  190 write(6,*)'Use one of the prespecified positions?'
      write(6,*)'   (0) No specified position'
      write(6,*)'   (1)  Bchl 1 position  :: Alpha B850 A'
      write(6,*)'   (2)  Bchl 2 position  :: Beta  B850 A'
      write(6,*)'   (3)  Bchl 3 position  :: Alpha B850 B'
      write(6,*)'   (4)  Bchl 4 position  :: Beta  B850 B'
      write(6,*)'   (5)  Bchl 5 position  :: Alpha B850 C'
      write(6,*)'   (6)  Bchl 6 position  :: Beta  B850 C'
      write(6,*)'   (7)  Bchl 7 position  :: B800 A'
      write(6,*)'   (8)  Bchl 8 position  :: B800 B'
      write(6,*)'   (9)  Bchl 9 position  :: B800 C'
      write(6,*)'   (10) RG1 10 position  :: RG1 A'
      write(6,*)'   (11) RG1 11 position  :: RG1 B'
      write(6,*)'   (12) RG1 12 position  :: RG1 C'
      write(6,*)'   (13) Bchl Q1 position :: Alpha B850 D'
      write(6,*)'   (14) Bchl Q2 position :: Beta  B850 D'
      write(6,*)'   (15) Bchl Q3 position :: Alpha B850 E'
      write(6,*)'   (16) Bchl Q4 position :: Beta  B850 E'
      write(6,*)'   (17) Bchl Q5 position :: Alpha B850 F'
      write(6,*)'   (18) Bchl Q6 position :: Beta  B850 F'
      write(6,*)'   (19) Bchl R1 position :: Alpha B850 G'
      write(6,*)'   (20) Bchl R2 position :: Beta  B850 G'
      write(6,*)'   (21) Bchl R3 position :: Alpha B850 H'
      write(6,*)'   (22) Bchl R4 position :: Beta  B850 H'
      write(6,*)'   (23) Bchl R5 position :: Alpha B850 I'
      write(6,*)'   (24) Bchl R6 position :: Beta  B850 I'
      write(6,*)'   (25) Bchl R9 position :: B800 I'
      write(6,*)'   (26) RG1 R12 position :: RG1 I'      
      call inin('Enter choice/?', refflag)
      unitflag = 1
      if (refflag == 0) then
        unitflag = 0
      else if (refflag == 1) then
        refatom1 = bchlp1_1
        refatom2 = bchlp1_2
        refatom3 = bchlp1_3
      else if (refflag == 2) then
        refatom1 = bchlp2_1
        refatom2 = bchlp2_2
        refatom3 = bchlp2_3
      else if (refflag == 3) then
        refatom1 = bchlp3_1
        refatom2 = bchlp3_2
        refatom3 = bchlp3_3
      else if (refflag == 4) then
        refatom1 = bchlp4_1
        refatom2 = bchlp4_2
        refatom3 = bchlp4_3
      else if (refflag == 5) then
        refatom1 = bchlp5_1
        refatom2 = bchlp5_2
        refatom3 = bchlp5_3
      else if (refflag == 6) then
        refatom1 = bchlp6_1
        refatom2 = bchlp6_2
        refatom3 = bchlp6_3
      else if (refflag == 7) then
        refatom1 = bchlp7_1
        refatom2 = bchlp7_2
        refatom3 = bchlp7_3
      else if (refflag == 8) then
        refatom1 = bchlp8_1
        refatom2 = bchlp8_2
        refatom3 = bchlp8_3
      else if (refflag == 9) then
        refatom1 = bchlp9_1
        refatom2 = bchlp9_2
        refatom3 = bchlp9_3
      else if (refflag == 10) then
        refatom1 = rgp10_1
        refatom2 = rgp10_2
        refatom3 = rgp10_3
      else if (refflag == 11) then
        refatom1 = rgp11_1
        refatom2 = rgp11_2
        refatom3 = rgp11_3
      else if (refflag == 12) then
        refatom1 = rgp12_1
        refatom2 = rgp12_2
        refatom3 = rgp12_3
      else if (refflag == 13) then
        refatom1 = bchlq1_1
        refatom2 = bchlq1_2
        refatom3 = bchlq1_3
      else if (refflag == 14) then
        refatom1 = bchlq2_1
        refatom2 = bchlq2_2
        refatom3 = bchlq2_3
      else if (refflag == 15) then
        refatom1 = bchlq3_1
        refatom2 = bchlq3_2
        refatom3 = bchlq3_3
      else if (refflag == 16) then
        refatom1 = bchlq4_1
        refatom2 = bchlq4_2
        refatom3 = bchlq4_3
      else if (refflag == 17) then
        refatom1 = bchlq5_1
        refatom2 = bchlq5_2
        refatom3 = bchlq5_3
      else if (refflag == 18) then
        refatom1 = bchlq6_1
        refatom2 = bchlq6_2
        refatom3 = bchlq6_3
      else if (refflag == 19) then
        refatom1 = bchlr1_1
        refatom2 = bchlr1_2
        refatom3 = bchlr1_3
      else if (refflag == 20) then
        refatom1 = bchlr2_1
        refatom2 = bchlr2_2
        refatom3 = bchlr2_3
      else if (refflag == 21) then
        refatom1 = bchlr3_1
        refatom2 = bchlr3_2
        refatom3 = bchlr3_3
      else if (refflag == 22) then
        refatom1 = bchlr4_1
        refatom2 = bchlr4_2
        refatom3 = bchlr4_3
      else if (refflag == 23) then
        refatom1 = bchlr5_1
        refatom2 = bchlr5_2
        refatom3 = bchlr5_3
      else if (refflag == 24) then
        refatom1 = bchlr6_1
        refatom2 = bchlr6_2
        refatom3 = bchlr6_3
      else if (refflag == 25) then
        refatom1 = bchlr9_1
        refatom2 = bchlr9_2
        refatom3 = bchlr9_3
      else if (refflag == 26) then
        refatom1 = rgr12_1
        refatom2 = rgr12_2
        refatom3 = rgr12_3
      else
        goto 190
      end if
      write (6,*) ''
      write (6,*) 'First characteristic atom defines the origin.'
      write (6,*) '  For RG this is C16, For Bchl this is Mg.'
      call indp( '  Enter x/?', refatom1(1) )
      call indp( '  Enter y/?', refatom1(2) )
      call indp( '  Enter z/?', refatom1(3) )
      write (6,*) 'Second characteristic atom defines the x axis.'
      write (6,*) '  For RG this is C14, For Bchl this is NA.'
      call indp( '  Enter x/?', refatom2(1) )
      call indp( '  Enter y/?', refatom2(2) )
      call indp( '  Enter z/?', refatom2(3) )
      write (6,*) 'Third characteristic atom defines the y axis.'
      write (6,*) '  For RG this is C15, For Bchl this is NB.'
      call indp( '  Enter x/?', refatom3(1) )
      call indp( '  Enter y/?', refatom3(2) )
      call indp( '  Enter z/?', refatom3(3) )
  200 call inin( 'Units? 1=Angstroms, 2=Bohr/?', unitflag )
      if (unitflag == 1) then
        do i = 1,3
          refatom1(i) = refatom1(i) / BohrtoAng
          refatom2(i) = refatom2(i) / BohrtoAng
          refatom3(i) = refatom3(i) / BohrtoAng
        end do
      else if (unitflag == 2) then
        continue
      else
        goto 200
      end if
      write(6,*)''
  201 write(6,*)'The program can also translate the final cube position'
      write(6,*)'  away from the reference position.'
      call inin('Translate cube? 0=no, 1=yes/?', shiftflag )
      if (shiftflag == 0) then
        do i = 1,3
            shiftmat(i) = 0.
        end do
      else if (shiftflag == 1) then
        write(6,*)'Enter translation as vector in reference atom frame.'
        call indp('  Enter x/?', shiftmat(1) )
        call indp('  Enter y/?', shiftmat(2) )
        call indp('  Enter z/?', shiftmat(3) )
  202   call inin('Units? 1=Angstroms, 2=Bohr/?', unitflag )
        if (unitflag == 1) then
          do i = 1,3
            shiftmat(i) = shiftmat(i) / BohrtoAng
          end do
        else if (unitflag == 2) then
          continue
        else
          goto 202
        end if
      else
        goto 201
      end if   
C  
C  Transform cube's characteristic atoms onto reference's
C
C    The characteristic vector1 will point from characteristic atom 1 to 2
C    Characteristic vector2 will point from characeteristic atom 1 to 3
C
      do i = 1,3
        r(i) = charatom2(i) - charatom1(i)
        rp(i) = refatom2(i) - refatom1(i)
        r2(i) = charatom3(i) - charatom1(i)
        rp2(i) = refatom3(i) - refatom1(i)
      end do
      rmag = SQRT(r(1)**2 + r(2)**2 + r(3)**2)
      rpmag = SQRT(rp(1)**2 + rp(2)**2 + rp(3)**2)
      if( (DABS(rpmag - rmag)) .gt. 1E-3 ) then
       write(6,*)'WARNING! Cube and reference vectors not same length.'
      end if
      r2mag = SQRT(r2(1)**2 + r2(2)**2 + r2(3)**2)
      rp2mag = SQRT(rp2(1)**2 + rp2(2)**2 + rp2(3)**2)
      if( (DABS(rp2mag - r2mag)) .gt. 1E-3 ) then
       write(6,*)'WARNING! Cube and reference vectors2 not same length.'
      end if
C
C  Define the rotation matrix which puts the cube characteristic vector
C    onto the reference characteristic vector.
C 
C  Set the quadrant flags to account for the inverse COS only spanning
C    180 degrees.  The flag expands this to +- 180 degrees.
C
C    For first transform (rotation within XY plane) if y is positive
C      need negative rotation and vice versa.
      if (r(2) .lt. 0) then
        rxyflag = 1
      else
        rxyflag = -1
      end if
      if (rp(2) .lt. 0) then
        rpxyflag = 1
      else
        rpxyflag = -1
      end if
C    The second transform depends on the sign of x, but the first transform
C      should have rotated both vectors to positive x so we don't need flags.
C    Third rotation depends on the first, so we need no new information.
C
C  Now the rotations...
C
C  First, we will rotate both r and rp into the XZ plane:
C    Take projection onto XY plane (z=0), unitize it, and dot into xhat.
C    This gives rotation about z axis to bring each vector into XZ plane.
C	
      rtmp(1) = r(1)
      rptmp(1) = rp(1)
      rtmp(2) = r(2)
      rptmp(2) = rp(2)
      rtmp(3) = 0.
      rptmp(3) = 0.
      call unitvec(rptmp, rphat)
      call unitvec(rtmp, rhat)
      rpphi = rpxyflag * ACOS( DOT_PRODUCT( rphat, xhat ))
      rphi = rxyflag * ACOS( DOT_PRODUCT( rhat, xhat ))
      phi = rphi - rpphi
C      write (6,*) 'First Euler rotation is phi = ',rphi
C    Create rotation matrix to rotate rp into the XZ plane.
C       (we will undo this rotation at the end for both r and rp)	
      phimat(1,1) = COS(rpphi)
      phimat(2,1) = SIN(rpphi)
      phimat(3,1) = 0.
      phimat(1,2) = -SIN(rpphi)
      phimat(2,2) = COS(rpphi)
      phimat(3,2) = 0.
      phimat(1,3) = 0.
      phimat(2,3) = 0.
      phimat(3,3) = 1.
C    Do rotation of rp into XZ plane, bringing rp2 with it
      rpnew = MATMUL(phimat, rp)
      rp2new = MATMUL(phimat, rp2)
C    Recreate rotation matrix to rotate r into the XZ plane	
      phimat(1,1) = COS(rphi)
      phimat(2,1) = SIN(rphi)
      phimat(1,2) = -SIN(rphi)
      phimat(2,2) = COS(rphi)
C    Do rotation of r into XZ plane, bringing r2 with it
      rnew = MATMUL(phimat, r)
      r2new = MATMUL(phimat, r2)
C
C  Now, both r and rp should be in XZ plane (y=0) already, so simply 
C    rotate r and rp about y so that they lie on z.  
C	
C      write(6,*) 'y component of rotated r and rp: ', rnew(2), rpnew(2)
      rtmp(1) = rnew(1)
      rptmp(1) = rpnew(1)
      rtmp(2) = rnew(2)
      rptmp(2) = rpnew(2)
      rtmp(3) = rnew(3)
      rptmp(3) = rpnew(3)
      call unitvec(rptmp, rphat)
      call unitvec(rtmp, rhat)
      rptheta = -ACOS( DOT_PRODUCT( rphat, zhat ))
      rtheta = -ACOS( DOT_PRODUCT( rhat, zhat ))
      theta = rtheta - rptheta
C      write (6,*) 'Second Euler rotation is theta = ',theta
C    Create rotation matrix	
      thetamat(1,1) = COS(rptheta)
      thetamat(2,1) = 0.
      thetamat(3,1) = -SIN(rptheta)
      thetamat(1,2) = 0.
      thetamat(2,2) = 1.
      thetamat(3,2) = 0.
      thetamat(1,3) = SIN(rptheta)
      thetamat(2,3) = 0.
      thetamat(3,3) = COS(rptheta)
C    Do second rotation of rp onto z axis and bring rp2 along
      rpnew = MATMUL(thetamat, rpnew)
      rp2new = MATMUL(thetamat, rp2new)
C    Recreate rotation matrix for r and r2
      thetamat(1,1) = COS(rtheta)
      thetamat(3,1) = -SIN(rtheta)
      thetamat(1,3) = SIN(rtheta)
      thetamat(3,3) = COS(rtheta)
C    Do second rotation of r onto z axis and bring r2 along
      rnew = MATMUL(thetamat, rnew)
      r2new = MATMUL(thetamat, r2new)
C
C  Now both r and rp should be pointing along z axis, so rotate about 
C	z to bring r2 onto rp2 without loosing r and rp overlap.
C     Again, we make use of flagging to attain +- 180 degree rotation.
C
      if (r2new(2) .lt. 0) then
        r2xyflag = 1
      else
        r2xyflag = -1
      end if
      if (rp2new(2) .lt. 0) then
        rp2xyflag = 1
      else
        rp2xyflag = -1
      end if
      rtmp(1) = r2new(1)
      rptmp(1) = rp2new(1)
      rtmp(2) = r2new(2)
      rptmp(2) = rp2new(2)
      rtmp(3) = 0.
      rptmp(3) = 0.
      call unitvec(rptmp, rphat)
      call unitvec(rtmp, rhat)
      rppsi = rp2xyflag * ACOS( DOT_PRODUCT( rphat, xhat ))
      rpsi = r2xyflag * ACOS( DOT_PRODUCT( rhat, xhat ))
      psi = rpsi - rppsi
C      write (6,*) 'Third Euler rotation is psi = ',psi
C    Create rotation matrix	
      psimat(1,1) = COS(psi)
      psimat(2,1) = SIN(psi)
      psimat(3,1) = 0.
      psimat(1,2) = -SIN(psi)
      psimat(2,2) = COS(psi)
      psimat(3,2) = 0.
      psimat(1,3) = 0.
      psimat(2,3) = 0.
      psimat(3,3) = 1.
C    Do third rotation of r2 onto rp2 and bring r along
      rnew = MATMUL(psimat, rnew)
      r2new = MATMUL(psimat, r2new)
C
C  Now let's construct rotation matrix that we have made so far.
C
      rotmat = MATMUL(psimat, thetamat)
      rotmat = MATMUL(rotmat, phimat)
C
C  Now, the pairs of vectors have been superimposed, however, we have
C    perturbed the reference vectors in order to do so.  So, rotate the
C    whole system back so that the reference vectors lie in their 
C    original orientation.
C
C    Undo theta rotation
      thetamat(1,1) = COS(-rptheta)
      thetamat(3,1) = -SIN(-rptheta)
      thetamat(1,3) = SIN(-rptheta)
      thetamat(3,3) = COS(-rptheta)
      rnew = MATMUL(thetamat, rnew)
      r2new = MATMUL(thetamat, r2new)
      rpnew = MATMUL(thetamat, rpnew)
      rp2new = MATMUL(thetamat, rp2new)
      rotmat = MATMUL(thetamat, rotmat)
C    Undo phi rotation
      phimat(1,1) = COS(-rpphi)
      phimat(2,1) = SIN(-rpphi)
      phimat(1,2) = -SIN(-rpphi)
      phimat(2,2) = COS(-rpphi)
      rnew = MATMUL(phimat, rnew)
      r2new = MATMUL(phimat, r2new)
      rpnew = MATMUL(phimat, rpnew)
      rp2new = MATMUL(phimat, rp2new)
      rotmat = MATMUL(phimat, rotmat)
C
C  Let's retransform the original r and r2 using the rotmat to double check
C
      rnew2 = MATMUL(rotmat, r)
      r2new2 = MATMUL(rotmat, r2)
C
C  We should be done, print out vectors
C
      write (6,'(" Cube char vector1:           ", 3F9.4)') r
      write (6,'(" Cube char vector2:           ", 3F9.4)') r2
      write (6,'(" Reference char vector1:      ", 3F9.4)') rp
      write (6,'(" Reference char vector2:      ", 3F9.4)') rp2
      write (6,'(" Transformed cube vector1:    ", 3F9.4)') rnew
      write (6,'(" Transformed cube vector2:    ", 3F9.4)') r2new
      write (6,'(" Transformed ref vector1:     ", 3F9.4)') rpnew
      write (6,'(" Transformed ref vector2:     ", 3F9.4)') rp2new
      write (6,'(" Matrix xform of cube vector1:", 3F9.4)') rnew2
      write (6,'(" Matrix xform of cube vector2:", 3F9.4)') r2new2
C
C  Now transform all the cube coordinates in this same manner
C
C    First, translate cube origin and all atoms such that characteristic
C	 atom 1 of cube is at the origin.
C
      transmat = charatom1
      do i = 1,3
        cubeorig(i) = cubeorig(i) - transmat(i)
        charatom1(i) = charatom1(i) - transmat(i)
        charatom2(i) = charatom2(i) - transmat(i)
        do j = 1, Natoms
          atpo(i,j) = atpo(i,j) - transmat(i)
        end do
      end do
C
C  Now, rotate all above coordinates along with cube stepping vectors
C
      cubeorig = MATMUL(rotmat, cubeorig)
      charatom1 = MATMUL(rotmat, charatom1)
      charatom2 = MATMUL(rotmat, charatom2)
      step1 = MATMUL(rotmat, step1)
      step2 = MATMUL(rotmat, step2)
      step3 = MATMUL(rotmat, step3)
      do i = 1, Natoms
        do j = 1,3
          rtmp(j) = atpo(j,i)		!need temp vector to use MATMUL
        end do
        rtmp = MATMUL(rotmat, rtmp)
        do j = 1,3
          atpo(j,i) = rtmp(j)
        end do
      end do
C
C    Now, translate cube origin and all atoms such that characteristic atom
C	 1 is at the reference atom 1 position + any shift requested by user.
C
      do i = 1,3
        cubeorig(i) = cubeorig(i) + refatom1(i) + shiftmat(i)
        charatom1(i) = charatom1(i) + refatom1(i) + shiftmat(i)
        charatom2(i) = charatom2(i) + refatom1(i) + shiftmat(i)
        do j = 1, Natoms
          atpo(i,j) = atpo(i,j) + refatom1(i) + shiftmat(i)
        end do
      end do
C
C  Now fill in dummy center atom for carotenoid positions
C
      if( refflag == 10) then
        N(1) = 2    !make the dummy center atom helium
        atno(1) = 2.
        atpo(1,1) = rgp10_c(1) / BohrtoAng + shiftmat(1)	 !convert to Bohr
        atpo(2,1) = rgp10_c(2) / BohrtoAng + shiftmat(2)  !and translate
        atpo(3,1) = rgp10_c(3) / BohrtoAng + shiftmat(3)
      else if( refflag == 11) then
        N(1) = 2  
        atno(1) = 2.
        atpo(1,1) = rgp11_c(1) / BohrtoAng + shiftmat(1)
        atpo(2,1) = rgp11_c(2) / BohrtoAng + shiftmat(2)
        atpo(3,1) = rgp11_c(3) / BohrtoAng + shiftmat(3)
      else if( refflag == 12) then
        N(1) = 2 
        atno(1) = 2.
        atpo(1,1) = rgp12_c(1) / BohrtoAng + shiftmat(1)
        atpo(2,1) = rgp12_c(2) / BohrtoAng + shiftmat(2)
        atpo(3,1) = rgp12_c(3) / BohrtoAng + shiftmat(3)
      end if
C
C  Finally, write new cube data to disk
C
      unitflag = 2
 1102 call inin('Units to write atoms to disk? 1=Ang,2=Bohr/?',unitflag)
      if (unitflag == 1) then
        do i = 1,Natoms
          do j = 1,3
            atpo(j,i) = atpo(j,i) * BohrtoAng
          end do
        end do
      else if (unitflag == 2) then
        continue
      else
        goto 1102
      end if
C
      open (unit=12, file='cube.transform.done', err=9999)
      write (12,1105) title
      write (12,1105) dummy
 1105 format (A80)
C  Number of atoms and origin
      write (12,1200) Natoms, cubeorig(1), cubeorig(2), cubeorig(3)
C  Number of increments and size in each direction
      write (12,1200) N1, step1(1), step1(2), step1(3)
      write (12,1200) N2, step2(1), step2(2), step2(3)
      write (12,1200) N3, step3(1), step3(2), step3(3)
 1200 format (I5,3F12.6)
C  Atomic number, charge and atomic positions
      write (12,'(I5,4F12.6)') (N(i),atno(i),atpo(1,i),atpo(2,i),
     & atpo(3,i), i=1,Natoms)
C  Write out the density cube (Note that units are es and Bohr.)
      do i=1,N1
        do j=1,N2
          write (12,'(6ES13.5)') (q(k,j,i), k=1,N3)          
        end do
      end do
      close (unit=12)
      deallocate( q, N, atno, atpo )
C
C  OK, we are all done
C
 2000 continue
 3000 continue
      goto 10000
 9999 write (6,*) 'ERROR!'
10000 end
C  
C  Now, define the couple of subroutines used in the main program
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C  unitvec takes a vector as input and returns the corresponding unit vector
C
      subroutine unitvec(invec, unvec)
      real*8 invec(3), unvec(3)
      real*8 mag
      mag = SQRT( invec(1)**2 + invec(2)**2 + invec(3)**2 )
      if (mag .eq. 0) then
        write (6,*) 'Warning zero magnitude vector encountered!'
        unvec(1) = 0
        unvec(2) = 0
        unvec(3) = 0
        return
      end if
      unvec(1) = invec(1)/mag
      unvec(2) = invec(2)/mag
      unvec(3) = invec(3)/mag
      mag = SQRT( unvec(1)**2 + unvec(2)**2 + unvec(3)**2 )
      if ( (DABS(mag - 1.)) .gt. 1E-8 ) then
        write (6,*) 'Warning!  Unit vector magnitude is not = 1.'
      end if 
      end
C  end of subroutine unitvec
C
C  indp and inin input numbers from the keyboard.
C    Taken from inputs.f which is part of upcvfit2
C
C************************************************************************
C	Subroutine: indp(vname,value)									  *
C	Purpose: Used to input the value of a double precision (R*8)	  *
C		variable.  Description of variable (vname) must use /? to 	  *
C	    denote its end.  (The /? will not be printed.)				  *
C************************************************************************
C
      subroutine indp(vname,value)
      integer nchars
      real*8 value,vdum
      character*15 readin
      character*80 vname
      nchars=index(vname,'/?')-1
  100 write(6,101,ADVANCE='NO') vname(1:nchars),value
  101 format(' ',A,'= ',g15.5,'  ')
      read(5,102,err=200) readin
  102 format(a15)
      if(readin .eq. '               ') return
      read(readin,103,err=200) vdum
  103 format(g15.0)
      value=vdum
      return
  200 write(6,*) '**** Error Reading Value ****'
      goto 100
      end
C  end of subroutine indp
C
C************************************************************************
C																	  *
C	Subroutine: inin(vname,ivalue)									  *
C																	  *
C	Purpose: Used to input a value (ivalue) for an integer variable	  *
C		Types variable description (vname -- must end in /?) and the  *
C         current value after which a new value can be input or the old *
C         value kept by hitting a return.								  *
C																	  *
C	Source:	mpm  8/28/85; last revision mpm 06/17/89				  *
C																	  *
C************************************************************************
C
      subroutine inin(vname,ivalue)
      integer nchars, ivalue, idum
      character*8 readin
      character*80 vname
      nchars=index(vname,'/?')-1
  100 write(6,101,ADVANCE='NO') vname(0:nchars),ivalue
  101 format(' ',A,'= ',I8,'  ')
      read(5,102,err=200) readin
  102 format(a8)
      if (readin.eq.'     ') return
      read(readin,103,err=200) idum
      ivalue=idum
  103 format(i8)
      return
  200 write(6,*) '**** Error Reading Value ****'
      goto 100
      end
C  end of subroutine inin
