C  This code has been used with SGI, HP, Cray and Windows.  Some modifications
C  are necessary to use the timing function in Windows compared to UNIX.
C  These lines are commented in the text.
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

      character*80 titleA, titleB, dummy
	INTEGER MaxN
	PARAMETER (MaxN=100)
	REAL pxA(MaxN),pyA(MaxN),pzA(MaxN)
	REAL pxB(MaxN),pyB(MaxN),pzB(MaxN)
	CHARACTER*6 DUMMY6
	CHARACTER*4 DUMMY4
	CHARACTER*2 ATOMA(MAXN),SPECA(MAXN),GROUPA(MAXN)
	CHARACTER*2 ATOMb(MAXN),SPECb(MAXN),GROUPb(MAXN)
	INTEGER iDUMMY, IMOLA(MAXN), IMOLB(MAXN)
      integer IA(maxn), IB(maxn)

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
      integer i,i1,i2,i3,j1,j2,j3,i4
      integer irflag
	integer jA, jB, jamin(9), jbmin(9), idist
	real dist, distmin(9)
	logical rearrange


C  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C      expA = 13.      !RG1
C      expB = 6.13     !Bchl Qy
C      expB = 3.29     !Bchl Qx

C  The transition density cube files (from Gaussian)
C  are called cubeA and cubeB.  
C
C  Open and process cubeA
9	format(A6,I5,1x,A2,A2,A4,A2,I4,A4,3F8.3)  
3     open (unit=16, file='pdbA', err=999)      
	I=0
18    i=i+1
	read (16,9,END=19) DUMMY6,IDUMMY,ATOMA(I),SPECA(i),
     .     DUMMY4,GROUPA(I),IMOLA(i),DUMMY4,PXA(i),PYA(I),PZA(I)	
	if(atoma(i).eq.'MG') then
	  IA(I)=12
	else if(atoma(i).eq.' C') then
	  IA(I) = 6
	else if(atoma(i).eq.' N') then
	  IA(I) = 7
	else if(atoma(i).eq.' O') then
	  IA(I) = 8
	else
	  IB(I) = 0
	ENDIF
	GO TO 18
19	CLOSE(16)
	NATOMSA=I-1

      write (6,*) NatomsA
      
C  --------------------------------------------------------------
C  Now do the same for cubeB

      open (unit=14, file='pdbB', err=999)
	I=0
28    i=i+1
	read (14,9,END=29) DUMMY6,IDUMMY,ATOMb(I),SPECb(i),
     .     DUMMY4,GROUPb(I),IMOLb(i),DUMMY4,PXb(i),PYb(I),PZb(I)
	if(atomb(i).eq.'MG') then
	  IB(I)=12
	else if(atomb(i).eq.' C') then
	  IB(I) = 6
	else if(atomb(i).eq.' N') then
	  IB(I) = 7
	else if(atomb(i).eq.' O') then
	  IB(I) = 8
	else
	  IB(I) = 0
	ENDIF
	GO TO 28
29	CLOSE(14)
	NATOMSb=I-1

      write (6,*) NatomsB


      close (unit=14)


C  Cartesian coordinates (convert to PDB file using "newzmat -ixyz -opdb")
c      write (12,*) ' Cartesian coordinates for dimer'
c      write (12,200) ( IA(i),pxA(i),pyA(i),pzA(i), i=1,NatomsA )
c      write (12,200) ( IB(i),pxB(i),pyB(i),pzB(i), i=1,NatomsB )
c 200  format(2X,I5,3F12.6)
	do i=1,9
	distmin(i) = 1.e10+float(i)*1.e9
	enddo
	do jA = 1, natomsA
	  do jB = 1, natomsB
	      riA(1) = pxA(jA)
	      riA(2) = pyA(jA)
	      riA(3) = pzA(jA)
	      rjB(1) = pxB(jB)
	      rjB(2) = pyB(jB)
	      rjB(3) = pzB(jB)
	      dist=DSQRT(DOT_PRODUCT( (riA-rjB), (riA-rjB)))
	      rearrange = .false.
		  do i=9,1,-1
	        if(distmin(i).ge.dist) then
	          idist = i
	          rearrange=.true.
	        endif
	      enddo
		  if (rearrange) then
c             write(*,'(3I5,f12.5)') jA,jB,idist,dist
		    do i=9, idist+1, -1
	          distmin(i) = distmin(i-1)
	          jamin(i) = jamin(i-1)
			  jbmin(i) = jbmin(i-1)
			enddo
			distmin(idist) = dist
			jamin(idist) = jA
			jBmin(idist) = jB
	      endif		  		  
	  enddo
	enddo
	do i=1,9
        write(*,200) jAmin(i), iA(jAmin(i)), jBmin(i), iB(jBmin(i)), 
     .               distmin(i)/Bohrtoang, distmin(i)
	enddo
c	write(*,*) 'take three atoms from cubeA(1) or cubeB(2)? (2)
c	read(*,*, err=199) i
c	go to 201
199	i=2	 
201	continue
202   format(a,a,a,I5,I5, 3F12.5) 
	if (i.eq.2) then
	  j1=jbmin(1)
	  write(*,202) atomb(j1),specb(j1),groupb(j1),imolb(j1),
     .      j1, pxB(j1), pyB(j1),pzB(j1)
	      write(12,*)pxB(j1)
	      write(12,*)pyB(j1)
	      write(12,*)pzB(j1)
	  i2 = 0
	  i3 = 0
	  i4=0
	  do i1=2,9
	    j1 = jbmin(i1)
	    If(i2.le.3.and.j1.ne.jBmin(1).and.j1.ne.i3.and.j1.ne.i4) then
	      i4 = i3
		  i3=j1
	      i2=i2+1
	      write(*,202) atomb(j1),specb(j1),groupb(j1),imolb(j1),
     .      j1, pxB(j1), pyB(j1),pzB(j1)
c	      write(*,'(I3, 3F12.5)') jBmin(i1), pxB(jBmin(i1)), 
c     .           pyB(jBmin(i1)), pzB(jBmin(i1))
	      write(12,*)pxB(j1)
	      write(12,*)pyB(j1)
	      write(12,*)pzB(j1)
	    endif
	  enddo
	endif

c	do i=1,9
c        write(*,201) 
c	enddo
200	format(4I5,2f12.5)
c      close (unit=12)
C  Deallocate atomic coordinates and such

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

