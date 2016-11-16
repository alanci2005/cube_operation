C  cleanup.f   (for fortran 90)
C  Brent Krueger Nov 1997 The University of Chicago -- Fleming Group
C
C
C  This is just transform.f mashed down to take all cube elements with
C    exponents of three digits and set them to zero since these are giving
C    strange output with no E in the formatting statement.  It reads in
C    the cube 'cleanup.in' and outputs the cube 'cleanup.out'.
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
      integer i, j, k, menu, cnt
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
C
C  Read in the cube datafile.
C
      open (unit=16, file='cleanup.in')
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
C  Read in the density cube and set any very small element charge to zero
C  (Note that units are es and Bohr.)
      allocate( q(N3, N2, N1) )
      do i=1,N1
        do j=1,N2
          read (16, '(6E13.5)') (q(k,j,i), k=1,N3)          
        end do
      end do
      close (unit=16)
C
      cnt = 0
C      write(6,*) 'first element',q(1,1,1)
C      write(6,*) 'abs of first ',ABS(q(1,1,1))
C      junk = (ABS(q(1,1,1)) .LT. 1)
C      write(6,*) 'ABS < 1?',junk
C      junk = (ABS(q(1,1,1)) .LT. 1E-90)
C      write(6,*) 'ABS < 1E-90 ?',junk
C      junk = (ABS(q(1,1,1)) .LT. 1D-90)
C      write(6,*) 'ABS < 1D-90 ?',junk

      do i = 1,N1
        do j = 1,N2
          do k = 1,N3
            if ( ABS( q(k,j,i) ) < 1.001D-99 ) then 
              q(k,j,i) = 0.
              cnt = cnt+1 
            end if
          end do
        end do
      end do
      write (6,*) 'number of replacements: ',cnt
C  
C
C  Finally, write new cube data to disk
C
C
      open (unit=12, file='cleanup.out', err=9999)
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
