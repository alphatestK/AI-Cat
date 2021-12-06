

module SSW_commsub


#ifdef MPI
  use mpi
  use mpiinfo_ssw
#endif
implicit none

#ifndef MPI
integer::cpuid=0
integer::nprocs=1
#endif
    save 

  public

  contains
!    This file contains routines adapted from 'Numerical Recipes,
!    The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky,
!    W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.


       subroutine rd_numb(x)

       double precision x
!      logical, save:: first=.true.


!      if(cpuid==0) then
!      if(first) call init_random_seed()
       call random_number(x)
!      endif

!#IFDEF MPI
!      call MPI_BARRIER(MPI_COMM_WORLD, mpier)
!       call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpier)
!#ENDIF
!      write(para%ioout,*) cpuid,x
!      first=.false.

       end subroutine



      DOUBLE PRECISION FUNCTION VOLCEL_LOC( C )

      DOUBLE PRECISION C(3,3)
      VOLCEL_LOC = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) + &
               ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) + &
               ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOLCEL_LOC = ABS( VOLCEL_LOC )
      END function

      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END function ran3


integer function elemz(element)
   implicit none
   character(len=2)::element
   if(trim(element)=='H') then
      elemz = 1
      return


   else if(trim(element)=='Li') then
      elemz = 3
      return


   else if(trim(element)=='B') then
      elemz = 5
      return

   else if(trim(element)=='N') then
      elemz = 7
      return

   else if(trim(element)=='F') then
      elemz = 9
      return

   else if(trim(element)=='Na') then
      elemz = 11
      return

   else if(trim(element)=='Al') then
      elemz = 13
      return

   else if(trim(element)=='P') then
      elemz = 15
      return

   else if(trim(element)=='Cl') then
      elemz = 17
      return

   else if(trim(element)=='K') then
      elemz = 19
      return

   else if(trim(element)=='Sc') then
      elemz = 21
      return

   else if(trim(element)=='Sc') then
      elemz = 22
      return

   else if(trim(element)=='V') then
      elemz = 23
      return

   else if(trim(element)=='Mn') then
      elemz = 25
      return

   else if(trim(element)=='Co') then
      elemz = 27
      return

   else if(trim(element)=='Cu') then
      elemz = 29
      return

   else if(trim(element)=='Ga') then
      elemz = 31
      return

   else if(trim(element)=='As') then
      elemz = 33
      return

   else if(trim(element)=='Br') then
      elemz = 35
      return

   else if(trim(element)=='Rb') then
      elemz = 37
      return

   else if(trim(element)=='Y') then
      elemz = 39
      return

   else if(trim(element)=='Nb') then
      elemz = 41
      return

   else if(trim(element)=='Tc') then
      elemz = 43
      return

   else if(trim(element)=='Rh') then
      elemz = 45
      return

   else if(trim(element)=='Ag') then
      elemz = 47
      return

   else if(trim(element)=='In') then
      elemz = 49
      return

   else if(trim(element)=='Sb') then
      elemz = 51
      return

   else if(trim(element)=='I') then
      elemz = 53
      return

   else if(trim(element)=='Cs') then
      elemz = 55
      return

   else if(trim(element)=='La') then
      elemz = 57
      return

   else if(trim(element)=='Pr') then
      elemz = 59
      return

   else if(trim(element)=='Pm') then
      elemz = 61
      return

   else if(trim(element)=='Eu') then
      elemz = 63
      return

   else if(trim(element)=='Tb') then
      elemz = 65
      return

   else if(trim(element)=='Ho') then
      elemz = 67
      return

   else if(trim(element)=='Tm') then
      elemz = 69
      return

   else if(trim(element)=='Lu') then
      elemz = 71
      return

   else if(trim(element)=='Ta') then
      elemz = 73
      return

   else if(trim(element)=='Re') then
      elemz = 75
      return

   else if(trim(element)=='Ir') then
      elemz = 77
      return

   else if(trim(element)=='Au') then
      elemz = 79
      return

   else if(trim(element)=='Tl') then
      elemz = 81
      return

   else if(trim(element)=='Bi') then
      elemz = 83
      return

   else if(trim(element)=='At') then
      elemz = 85
      return

   else if(trim(element)=='Fr') then
      elemz = 87
      return

   else if(trim(element)=='Ac') then
      elemz = 89
      return

   else if(trim(element)=='Pa') then
      elemz = 91
      return

   else if(trim(element)=='Np') then
      elemz = 93
      return

   else if(trim(element)=='Am') then
      elemz = 95
      return

   else if(trim(element)=='Bk') then
      elemz = 97
      return

   else if(trim(element)=='Es') then
      elemz = 99
      return

   else if(trim(element)=='Md') then
      elemz = 101
      return

   else if(trim(element)=='Lw') then
      elemz = 103
      return

   else if(trim(element)=='He') then
      elemz = 2
      return

   else if(trim(element)=='Be') then
      elemz = 4
      return

   else if(trim(element)=='C') then
      elemz = 6
      return

   else if(trim(element)=='O') then
      elemz = 8
      return

   else if(trim(element)=='Ne') then
      elemz = 10
      return

   else if(trim(element)=='Mg') then
      elemz = 12
      return

   else if(trim(element)=='Si') then
      elemz = 14
      return

   else if(trim(element)=='S') then
      elemz = 16
      return

   else if(trim(element)=='Ar') then
      elemz = 18
      return

   else if(trim(element)=='Ca') then
      elemz = 20
      return

   else if(trim(element)=='Ti') then
      elemz = 22
      return

   else if(trim(element)=='Cr') then
      elemz = 24
      return

   else if(trim(element)=='Fe') then
      elemz = 26
      return

   else if(trim(element)=='Ni') then
      elemz = 28
      return

   else if(trim(element)=='Zn') then
      elemz = 30
      return

   else if(trim(element)=='Ge') then
      elemz = 32
      return

   else if(trim(element)=='Se') then
      elemz = 34
      return

   else if(trim(element)=='Kr') then
      elemz = 36
      return

   else if(trim(element)=='Sr') then
      elemz = 38
      return

   else if(trim(element)=='Zr') then
      elemz = 40
      return

   else if(trim(element)=='Mo') then
      elemz = 42
      return

   else if(trim(element)=='Ru') then
      elemz = 44
      return

   else if(trim(element)=='Pd') then
      elemz = 46
      return

   else if(trim(element)=='Cd') then
      elemz = 48
      return

   else if(trim(element)=='Sn') then
      elemz = 50
      return

   else if(trim(element)=='Te') then
      elemz = 52
      return

   else if(trim(element)=='Xe') then
      elemz = 54
      return

   else if(trim(element)=='Ba') then
      elemz = 56
      return

   else if(trim(element)=='Ce') then
      elemz = 58
      return

   else if(trim(element)=='Nd') then
      elemz = 60
      return

   else if(trim(element)=='Sm') then
      elemz = 62
      return

   else if(trim(element)=='Gd') then
      elemz = 64
      return

   else if(trim(element)=='Dy') then
      elemz = 66
      return

   else if(trim(element)=='Er') then
      elemz = 68
      return

   else if(trim(element)=='Yb') then
      elemz = 70
      return

   else if(trim(element)=='Hf') then
      elemz = 72
      return

   else if(trim(element)=='W') then
      elemz = 74
      return

   else if(trim(element)=='Os') then
      elemz = 76
      return

   else if(trim(element)=='Pt') then
      elemz = 78
      return

   else if(trim(element)=='Hg') then
      elemz = 80
      return

   else if(trim(element)=='Pb') then
      elemz = 82
      return

   else if(trim(element)=='Po') then
      elemz = 84
      return

   else if(trim(element)=='Rn') then
      elemz = 86
      return

   else if(trim(element)=='Ra') then
      elemz = 88
      return

   else if(trim(element)=='Th') then
      elemz = 90
      return

   else if(trim(element)=='U') then
      elemz = 92
      return

   else if(trim(element)=='Pu') then
      elemz = 94
      return

   else if(trim(element)=='Cm') then
      elemz = 96
      return

   else if(trim(element)=='Cf') then
      elemz = 98
      return

   else if(trim(element)=='Fm') then
      elemz = 100
      return

   else if(trim(element)=='No') then
      elemz = 102
      return
   else
      elemz = 0
      return
   endif

end function


      FUNCTION SYMBOL( Z )

      character(len=2)    :: SYMBOL  ! Atomic symbol
      integer, intent(in) :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2), parameter :: NAME(NZ) =   &
              (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',  &
                'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',  &
                'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',  &
                'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',  &
                'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',  &
                'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',  &
                'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',  &
                'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',  &
                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',  &
                'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',  &
                'Md','No','Lr'/)

      IF (Z.EQ.0 .OR. Z.EQ.-100) THEN
         SYMBOL = 'BS'
      ELSE IF (ABS(Z).LE.NZ) THEN
         SYMBOL = NAME(ABS(Z))
      ELSE IF (Z.GT.200) THEN
         if(cpuid==0) write(SYMBOL,'(a1,i1)') 'S', mod(Z-200,10)
      ELSE
         if(cpuid==0) WRITE(6,*) 'SYMBOL: ERROR: No data for Z =', Z
         SYMBOL = ' '
      ENDIF

      END function symbol



      SUBROUTINE reci_latt(A,B)

!  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
!  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

      integer ::  I
      DOUBLE PRECISION A(3,3),B(3,3), c , ci
      B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
      B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
      B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
      B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      C=1.D0
      DO I=1,3
         CI=C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
         B(1,I)=B(1,I)*CI
         B(2,I)=B(2,I)*CI
         B(3,I)=B(3,I)*CI
      ENDDO
      END SUBROUTINE reci_latt



!      subroutine readarc(filename,na,iza,cell,xa)
!
!      use constants
!      use ssw_parameters,only:para
!
!      implicit          none
!
!      integer                             :: iza(na)
!      double precision                    :: xa(3,na)
!      integer                             :: unit
!      double precision                    :: cell(3,3)
!      character(*) :: filename
!
!      integer                             :: i, ia,j
!      double precision                    :: latt_abc(3)
!      double precision                    :: latt_angle(3)
!      double precision                    :: celli(3,3)
!      double precision                    :: celltran(3,3)
!      double precision                    :: angle(3)
!      double precision                    :: xfrac(3)
!      double precision                    :: xa_out(3)
!! -------------------------------------------------------------------
!      character(3)::pbc
!      character(2)::lab4
!      integer na
!      double precision::abc(3),alplp,betlp,gamlp,alp,blp,clp
!
!      unit=89128   ! random number
!      open(unit, file=filename, form = 'formatted', status='old')
!      rewind(unit)
!      read(unit,*,end=1721,err=1721)
!      read(unit,*,end=1721,err=1721)
!      read(unit,*,end=1718,err=1718)
!      !write(para%ioout,*) '11fdasf',na
!1718  continue
!      read(unit,*,end=1721,err=1721)
!      read(unit,*,end=1721,err=1721) pbc,abc(:),angle(3),angle(2),angle(1)
!      !write(para%ioout,*) '21fdasf',pbc
!      cell=0.0d0
!      alplp=angle(3)* pi/180.d0
!      betlp=angle(2)* pi/180.d0
!      gamlp=angle(1)* pi/180.d0
!      alp=abc(1)
!      blp=abc(2)
!      clp=abc(3)
!      cell(1,1) = alp
!      cell(1,2) = blp * dcos(gamlp)
!      cell(2,2) = blp * dsin(gamlp)
!      cell(1,3) = clp * dcos(betlp)
!      cell(2,3) = (clp*dcos(alplp) - clp*dcos(gamlp)*dcos(betlp))/dsin(gamlp)
!      cell(3,3) = dsqrt(clp*clp - cell(1,3)*cell(1,3) - cell(2,3)*cell(2,3))
!      do i=1,na
!          read(unit,*,end=1720,err=1720) lab4,xa_out(1:3)
!          xa(:,i) = xa_out(1:3)
!          iza(i)=elemz(lab4)
!      !   write(para%ioout,*) 'atomiza', iza(i)
!      enddo
!1720  continue
!      read(unit,*,end=1721,err=1721)
!      read(unit,*,end=1721,err=1721)
!
!      close(unit)
!      return
!1721  continue
!      if(cpuid==0) write(para%ioout,'(A20,1x,A20)') 'No structures from ',filename
!      stop
!      end subroutine readarc

      subroutine readarc(filename,na,iza,cell,xa)

      use constants
      use ssw_parameters,only:para

      implicit          none

      integer                             :: iza(na)
      double precision                    :: xa(3,na)
      integer                             :: unit
      double precision                    :: cell(3,3)
      character(*) :: filename

      integer                             :: i, ia,j
      double precision                    :: latt_abc(3)
      double precision                    :: latt_angle(3)
      double precision                    :: celli(3,3)
      double precision                    :: celltran(3,3)
      double precision                    :: angle(3)
      double precision                    :: xfrac(3)
      double precision                    :: xa_out(3)
! -------------------------------------------------------------------
      character(3)::pbc
      character(2)::lab4
      character(20)::lab20
      integer na
      double precision::abc(3),alplp,betlp,gamlp,alp,blp,clp
      character(100)::string,keyword
      logical::findPBC=.false.

      unit=89128   ! random number
      open(unit, file=filename, form = 'formatted', status='old')
      rewind(unit)

      do
        read(unit,'(A)',err=9,end=9) string
        string=adjustl(string)
        if ( string(1:1) == '#' .or. len(trim(string))==0 ) cycle
        read(string,*) keyword
        if(adjustl(trim(keyword))=="PBC") then
          findPBC=.true.
          exit
        endif
      enddo
9     continue
      if(.not.findPBC) then
        !write(ioout,*) "No PBC provided"
         goto 1720
      else
        read(string,*) pbc,abc(:),angle(3),angle(2),angle(1)
      endif
      cell=0.0d0
      alplp=angle(3)* pi/180.d0
      betlp=angle(2)* pi/180.d0
      gamlp=angle(1)* pi/180.d0
      alp=abc(1)
      blp=abc(2)
      clp=abc(3)
      cell(1,1) = alp
      cell(1,2) = blp * dcos(gamlp)
      cell(2,2) = blp * dsin(gamlp)
      cell(1,3) = clp * dcos(betlp)
      cell(2,3) = (clp*dcos(alplp) - clp*dcos(gamlp)*dcos(betlp))/dsin(gamlp)
      cell(3,3) = dsqrt(clp*clp - cell(1,3)*cell(1,3) - cell(2,3)*cell(2,3))
      do i=1,na
          read(unit,*,end=1720,err=1720) lab20,xa_out(1:3)
          lab4=lab20(1:2)
          if(lab4(2:2)=='0' .or. lab4(2:2)=='1' .or. lab4(2:2)=='2' .or. lab4(2:2)=='3' .or.  &
             lab4(2:2)=='4' .or. lab4(2:2)=='5' .or. lab4(2:2)=='6' .or. lab4(2:2)=='7' .or.  &
             lab4(2:2)=='8' .or. lab4(2:2)=='9' ) lab4=lab20(1:1)
          xa(:,i) = xa_out(1:3)
          iza(i)=elemz(lab4)
      enddo
1720  continue
      read(unit,*,end=1721,err=1721)
      read(unit,*,end=1721,err=1721)
1721  continue
      close(unit)
      return

      end subroutine readarc

      subroutine init_random_seed()
      implicit none

      integer :: i,n,clock
      integer,dimension(:),allocatable :: seed

       call random_seed(size=n)
!#IFDEF MPI
!       call MPI_BARRIER(MPI_COMM_WORLD, mpier)
!       call MPI_BCAST(n,1,MPI_integer, 0, MPI_COMM_WORLD, mpier)
!#ENDIF
       !write(para%ioout,*) 'init random',cpuid, n
      
      allocate(seed(n))
!     seed =1
!     allocate(seed(n))

!     call system_clock(count=clock)
!     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      if(cpuid==0) then
!     call random_seed(size=n)
      OPEN(89,FILE='/dev/urandom',ACCESS='stream',FORM='UNFORMATTED')
      READ(89) seed
      CLOSE(89)
      endif
!     read seed from the entropy pool
!     seed=182732
!     seed=8329

#IFDEF MPI
!      call MPI_BARRIER(MPI_COMM_WORLD, mpier)
       call MPI_BCAST(seed,n,MPI_integer, 0, MPI_COMM_WORLD, mpier)
#ENDIF
      
      call random_seed(put=seed)

!      write(para%ioout,*) 'init random',cpuid, seed

      deallocate(seed)

      end subroutine
  
    subroutine N_normal(na,N,L)

    implicit none
    integer na
    double precision N(3,na),r
    logical L

    L=.true.
    r = sum(N*N)
    if (r>1d-6) then
       N = N/dsqrt(r)
    else
       N=0d0
       L=.false.
    endif
 
    end subroutine


    subroutine SetConstraints(Nions,xa,N)

      IMPLICIT NONE
      integer Nions,m0,m1
      double precision ,DIMENSION(3,Nions) :: N,xa
      double precision ,DIMENSION(:,:),allocatable ,save :: Nx,Ny,Nz
      double precision ,DIMENSION(:,:),allocatable ,save :: RotModeX,RotModeY,RotModeZ
      double precision com(3),axisatom(3,3),V(3,3)
      integer, save :: IR
      logical,save :: first=.true.
      integer i,j
      logical L

      !write(para%ioout,*) 'set constraints',Nions
      !write(para%ioout,*) 'in SecCons--00',N

      if(first) then
         allocate(Nx(3,Nions),Ny(3,Nions),Nz(3,Nions))
         allocate(RotModeX(3,Nions),RotModeY(3,Nions),RotModeZ(3,Nions))
         Nx=0.0d0
         Ny=0.0d0
         Nz=0.0d0
         Nx(1,:)=1.0d0
         Ny(2,:)=1.0d0
         Nz(3,:)=1.0d0
!        if(Imode>=0) then
!          Nx(1,Nions-2:Nions)=0d0
!          Ny(2,Nions-2:Nions)=0d0
!          Nz(3,Nions-2:Nions)=0d0
!        endif
         Nx=Nx/sqrt(sum(Nx*Nx))
         Ny=Ny/sqrt(sum(Ny*Ny))
         Nz=Nz/sqrt(sum(Nz*Nz))
         first=.false.
      endif
      !write(para%ioout,*) 'in SecCons--10',N
      !write(para%ioout,*) 'in SecCons--1',Nx
      !write(para%ioout,*) 'in SecCons--1',Ny
      !write(para%ioout,*) 'in SecCons--1',Nz


      N=N-SUM(N*Nx)*Nx-SUM(N*Ny)*Ny-SUM(N*Nz)*Nz

!     call N_normal(Nions,N,L)
      !write(para%ioout,*) 'in SecCons',N
!     if(Imode>=0) return  ! for crystal the N direction ---only consider to remove the translational mode of atom

      !write(para%ioout,*) xa
      !write(para%ioout,*) Nions
      com(:)=0.0d0
      do i=1,Nions
      do j=1,3
        com(j)=xa(j,i)/dble(Nions)+com(j)
      enddo
      enddo
!generate three axis
      RotModeX=0d0
      RotModeY=0d0
      RotModeZ=0d0
      do j=1,3
      do i=1,3
      if(i==j)then
        axisatom(i,j)=com(i)+1.0d0
      else
        axisatom(i,j)=com(i)
      endif
      enddo
      enddo
!generate three axis end
!generate rotational mode
      do i=1,Nions
        do j=1,3
          V=0.0d0
          V(:,1)=xa(:,i)-com(:)
          V(:,2)=xa(:,i)-axisatom(:,j)
          V(1,3)=V(2,1)*V(3,2)-V(3,1)*V(2,2)
          V(2,3)=V(3,1)*V(1,2)-V(1,1)*V(3,2)
          V(3,3)=V(1,1)*V(2,2)-V(2,1)*V(1,2)
          if(j==1) RotModeX(:,i)=V(:,3)
          if(j==2) RotModeY(:,i)=V(:,3)
          if(j==3) RotModeZ(:,i)=V(:,3)
        enddo
      enddo
      if(sqrt(sum(RotModeX*RotModeX))>1d-150)then
        RotModeX=RotModeX/sqrt(sum(RotModeX*RotModeX))
      else
         RotModeX=0d0
      endif

      if(sqrt(sum(RotModeY*RotModeY))>1d-150)then
        RotModeY=RotModeY/sqrt(sum(RotModeY*RotModeY))
            else
       RotModeY=0d0
      endif

      if(sqrt(sum(RotmodeZ*RotmodeZ))>1d-150) then
        RotmodeZ=RotmodeZ/sqrt(sum(RotmodeZ*RotmodeZ))
      else
        RotModeZ=0d0
      endif

      RotModeY=RotModeY-sum(RotModeY*RotModeX )*RotModeX
      if(sqrt(sum(RotModeY*RotModeY))>1d-150) then
        RotModeY=RotModeY/sqrt(sum(RotModeY*RotModeY))
      else
        RotModeY=0d0
      endif

      RotmodeZ=RotmodeZ-sum(RotModeZ*RotModeX )*RotModeX &
                       -sum(RotModeZ*RotModeY )*RotModeY
      if(sqrt(sum(RotModeZ*RotModeZ))>1d-150) then
        RotmodeZ=RotmodeZ/sqrt(sum(RotmodeZ*RotmodeZ))
      else
        RotmodeZ=0d0
      endif

      N=N-sum(N*RotmodeX)*RotmodeX-sum(N*RotmodeY)*RotmodeY &
          -sum(N*RotmodeZ)*RotmodeZ

      !write(para%ioout,*) N
!     call N_normal(Nions,N,L)



     end subroutine SetConstraints


     subroutine neighboringlist(na,xa,atom1,Nnei,neighbor,iza,ffix)

             implicit none
             integer atom1,na
             integer iza(na),ffix(3,na)
             integer Nnei,i,j
             integer neighbor(100)
             double precision xa(3,na),rij,r,radius(2)

             neighbor(:)=0

             Nnei=0
             call species_radius(iza(atom1),radius(1))
             do j=1,na
                if(j/=atom1) then
                rij=(xa(1,atom1)-xa(1,j))**2+(xa(2,atom1)-xa(2,j))**2+(xa(3,atom1)-xa(3,j))**2
                rij=dsqrt(rij)

               call species_radius(iza(j),radius(2))
               r=radius(1)+radius(2)+0.5d0
                if(rij<r.and.ffix(1,j)/=0) then
                   Nnei=Nnei+1
                   neighbor(Nnei)=j
                endif
                endif
              enddo

      end subroutine neighboringlist

    subroutine species_radius(atom,atomicRadius)
      integer::atom
      double precision atomicRadius

      !atomicRadius = 0.0d0
      select case(atom)
        case(1)
          atomicRadius=0.310d0
        case(2)
          atomicRadius=0.280d0
        case(3)
          atomicRadius=1.280d0
        case(4)
          atomicRadius=0.960d0
        case(5)
          atomicRadius=0.850d0
        case(6)
          atomicRadius=0.760d0
        case(7)
          atomicRadius=0.710d0
        case(8)
          atomicRadius=0.660d0
        case(9)
          atomicRadius=0.570d0
        case(10)
          atomicRadius=0.580d0
        case(11)
          atomicRadius=1.660d0
        case(12)
          atomicRadius=1.410d0
        case(13)
          atomicRadius=1.210d0
        case(14)
          atomicRadius=1.110d0
        case(15)
          atomicRadius=1.070d0
        case(16)
          atomicRadius=1.050d0
        case(17)
          atomicRadius=1.020d0
        case(18)
          atomicRadius=1.060d0
        case(19)
          atomicRadius=2.030d0
        case(20)
          atomicRadius=1.760d0
        case(21)
          atomicRadius=1.700d0
        case(22)
          atomicRadius=1.600d0
        case(23)
          atomicRadius=1.530d0
        case(24)
          atomicRadius=1.390d0
        case(25)
          atomicRadius=1.390d0
        case(26)
          atomicRadius=1.320d0
        case(27)
          atomicRadius=1.260d0
        case(28)
          atomicRadius=1.240d0
        case(29)
          atomicRadius=1.320d0
        case(30)
          atomicRadius=1.220d0
        case(31)
          atomicRadius=1.220d0
        case(32)
          atomicRadius=1.200d0
        case(33)
          atomicRadius=1.190d0
        case(34)
          atomicRadius=1.200d0
        case(35)
          atomicRadius=1.200d0
        case(36)
          atomicRadius=1.160d0
        case(37)
          atomicRadius=2.200d0
        case(38)
          atomicRadius=1.950d0
        case(39)
          atomicRadius=1.900d0
        case(40)
          atomicRadius=1.750d0
        case(41)
          atomicRadius=1.640d0
        case(42)
          atomicRadius=1.540d0
        case(43)
          atomicRadius=1.470d0
        case(44)
          atomicRadius=1.460d0
        case(45)
          atomicRadius=1.420d0
        case(46)
          atomicRadius=1.390d0
        case(47)
          atomicRadius=1.450d0
        case(48)
          atomicRadius=1.440d0
        case(49)
          atomicRadius=1.420d0
        case(50)
          atomicRadius=1.390d0
        case(51)
          atomicRadius=1.390d0
        case(52)
          atomicRadius=1.380d0
        case(53)
          atomicRadius=1.390d0
        case(54)
          atomicRadius=1.400d0
        case(55)
          atomicRadius=2.440d0
        case(56)
          atomicRadius=2.150d0
        case(57)
          atomicRadius=2.070d0
        case(58)
          atomicRadius=2.040d0
        case(59)
          atomicRadius=2.030d0
        case(60)
          atomicRadius=2.010d0
        case(61)
          atomicRadius=1.990d0
        case(62)
          atomicRadius=1.980d0
        case(63)
          atomicRadius=1.980d0
        case(64)
          atomicRadius=1.960d0
        case(65)
          atomicRadius=1.940d0
        case(66)
          atomicRadius=1.920d0
        case(67)
          atomicRadius=1.920d0
        case(68)
          atomicRadius=1.890d0
        case(69)
          atomicRadius=1.900d0
        case(70)
          atomicRadius=1.870d0
        case(71)
          atomicRadius=1.870d0
        case(72)
          atomicRadius=1.750d0
        case(73)
          atomicRadius=1.700d0
        case(74)
          atomicRadius=1.620d0
        case(75)
          atomicRadius=1.510d0
        case(76)
          atomicRadius=1.440d0
        case(77)
          atomicRadius=1.410d0
        case(78)
          atomicRadius=1.360d0
        case(79)
          atomicRadius=1.360d0
        case(80)
          atomicRadius=1.320d0
        case(81)
          atomicRadius=1.450d0
        case(82)
          atomicRadius=1.460d0
        case(83)
          atomicRadius=1.480d0
        case(84)
          atomicRadius=1.400d0
        case(85)
          atomicRadius=1.500d0
        case(86)
          atomicRadius=1.500d0
        case(87)
          atomicRadius=2.600d0
        case(88)
          atomicRadius=2.210d0
        case(89)
          atomicRadius=2.150d0
        case(90)
          atomicRadius=2.060d0
        case(91)
          atomicRadius=2.000d0
        case(92)
          atomicRadius=1.960d0
        case(93)
          atomicRadius=1.900d0
        case(94)
          atomicRadius=1.870d0
        case(95)
          atomicRadius=1.800d0
        case(96)
          atomicRadius=1.690d0
        case default
          atomicRadius=2.0d0
      end select

    end subroutine species_radius



!     subroutine species_radius(atom,radius)
!
!     integer::atom
!     double precision::radius
!
!       radius=0.0d0
!       if(atom==1)   radius=0.400d0
!       if(atom==2)   radius=0.490d0
!       if(atom==3)   radius=1.400d0
!       if(atom==4)   radius=1.400d0
!       if(atom==5)   radius=0.720d0
!       if(atom==6)   radius=0.650d0
!       if(atom==7)   radius=0.620d0
!       if(atom==8)   radius=0.600d0
!       if(atom==9)   radius=0.570d0
!       if(atom==10)  radius=0.510d0
!       if(atom==11)  radius=2.230d0
!       if(atom==12)  radius=1.720d0
!       if(atom==13)  radius=1.820d0
!       if(atom==14)  radius=1.222d0
!       if(atom==15)  radius=1.230d0
!       if(atom==16)  radius=1.120d0
!       if(atom==17)  radius=1.120d0
!       if(atom==18)  radius=0.880d0
!       if(atom==19)  radius=2.770d0
!       if(atom==20)  radius=2.230d0
!       if(atom==21)  radius=2.090d0
!       if(atom==22)  radius=1.400d0
!       if(atom==23)  radius=1.920d0
!       if(atom==24)  radius=1.850d0
!       if(atom==25)  radius=1.790d0
!       if(atom==26)  radius=1.720d0
!       if(atom==27)  radius=1.670d0
!       if(atom==28)  radius=1.350d0
!       if(atom==29)  radius=1.570d0
!       if(atom==30)  radius=1.530d0
!       if(atom==31)  radius=1.810d0
!       if(atom==32)  radius=1.520d0
!       if(atom==33)  radius=1.330d0
!       if(atom==34)  radius=1.170d0
!       if(atom==35)  radius=1.120d0
!       if(atom==36)  radius=1.030d0
!       if(atom==37)  radius=2.210d0
!       if(atom==38)  radius=1.860d0
!       if(atom==39)  radius=2.270d0
!       if(atom==40)  radius=2.160d0
!       if(atom==41)  radius=2.080d0
!       if(atom==42)  radius=2.010d0
!       if(atom==43)  radius=1.950d0
!       if(atom==44)  radius=1.340d0
!       if(atom==45)  radius=1.340d0
!       if(atom==46)  radius=1.340d0
!       if(atom==47)  radius=1.340d0
!       if(atom==48)  radius=1.340d0
!       if(atom==49)  radius=2.000d0
!       if(atom==50)  radius=1.720d0
!       if(atom==51)  radius=1.530d0
!       if(atom==52)  radius=1.420d0
!       if(atom==53)  radius=1.320d0
!       if(atom==54)  radius=1.240d0
!       if(atom==55)  radius=3.340d0
!       if(atom==56)  radius=2.780d0
!       if(atom==57)  radius=2.740d0
!       if(atom==58)  radius=2.700d0
!       if(atom==59)  radius=2.670d0
!       if(atom==60)  radius=2.640d0
!       if(atom==61)  radius=2.620d0
!       if(atom==62)  radius=2.590d0
!       if(atom==63)  radius=2.560d0
!       if(atom==64)  radius=2.540d0
!       if(atom==65)  radius=2.510d0
!       if(atom==66)  radius=2.490d0
!       if(atom==67)  radius=2.470d0
!       if(atom==68)  radius=2.450d0
!       if(atom==69)  radius=2.420d0
!       if(atom==70)  radius=2.400d0
!       if(atom==71)  radius=2.250d0
!       if(atom==72)  radius=2.160d0
!       if(atom==73)  radius=2.090d0
!       if(atom==74)  radius=2.020d0
!       if(atom==75)  radius=1.370d0
!       if(atom==76)  radius=1.320d0
!       if(atom==77)  radius=1.370d0
!       if(atom==78)  radius=1.350d0
!       if(atom==79)  radius=1.350d0
!       if(atom==80)  radius=1.760d0
!       if(atom==81)  radius=2.080d0
!       if(atom==82)  radius=1.810d0
!       if(atom==83)  radius=1.630d0
!       if(atom==84)  radius=1.530d0
!       if(atom==85)  radius=1.430d0
!       if(atom==86)  radius=1.340d0
!
!     end subroutine
!

      subroutine  MetropolisMC(Energy1,Energy2,temper_loc,Laccept)


      implicit none
      logical Laccept
      double precision Energy1,Energy2,x,temper_loc

      call rd_numb(x)
      Laccept=.true.
      if((Energy2/20d0 > Energy1/20d0) .and. &
      (x>dexp(-(Energy2/20d0-Energy1/20d0)*96485d0/8.314d0/temper_loc) ) ) Laccept=.false.

      end subroutine MetropolisMC



        Subroutine Find_leastmoveAtoms(na,A,B,ffix,atompair)

!       use bpcbd, only: Imode,Lsurreac,surZ_height,surdir
!       Use ReactPattern
        use SSW_parameters
        

        implicit none
        integer na, atompair(2)
        double precision, dimension(:), allocatable :: distance, distance2
        double precision, dimension(3,na) :: A, B
        integer,dimension(3,na) :: ffix
        integer, save :: atom_save
        
        integer i, j, k, atom_min, atom_min2
        double precision r, y

!       if(first) then
        allocate (distance(na),distance2(na))
!          first=.false.
!       endif

!        A=0d0
!        B=0d0
!        do j=1, na
!         if(ffix(1,j) > 0 ) then
!            A(:,j) = IS(:,j)
!            B(:,j) = FS(:,j)
!          endif
!        enddo

        distance(:)= dsqrt( (B(1,:)-A(1,:))*(B(1,:)-A(1,:))+ &
                     (B(2,:)-A(2,:))*(B(2,:)-A(2,:))+ &
                     (B(3,:)-A(3,:))*(B(3,:)-A(3,:)) )
!       atom_min=minloc(distance(:),1,distance>0.1)    ! first least move atom
        atom_min=minloc(distance(:),1, ffix(1,:)>0)
        distance2(:)= dsqrt( (B(1,:)-B(1,atom_min))*(B(1,:)-B(1,atom_min))+ &
                     (B(2,:)-B(2,atom_min))*(B(2,:)-B(2,atom_min))+ &
                     (B(3,:)-B(3,atom_min))*(B(3,:)-B(3,atom_min)) )
        r=3d0
        if(para%Run_type<10) r=2d0   ! cluster
        do i=1,na
          if(distance2(i) < r) distance(i)=0d0
        enddo

!       atom_min2=minloc(distance(:),1,distance>0.1)     ! second least move atom
        atom_min2=minloc(distance(:),1,ffix(1,:)>0) !distance>0.1)     ! second least move atom

        distance(:)= dsqrt((B(1,:)-B(1,atom_min))*(B(1,:)-B(1,atom_min))+ &
                     (B(2,:)-B(2,atom_min))*(B(2,:)-B(2,atom_min))+ &
                     (B(3,:)-B(3,atom_min))*(B(3,:)-B(3,atom_min)))+ &
                    dsqrt( (B(1,:)-B(1,atom_min2))*(B(1,:)-B(1,atom_min2))+ &
                     (B(2,:)-B(2,atom_min2))*(B(2,:)-B(2,atom_min2))+ &
                     (B(3,:)-B(3,atom_min2))*(B(3,:)-B(3,atom_min2)) )

       atompair(1)   = maxloc(distance(:),1, ffix(1,:)>0)     ! farest atom to least-move atoms
       r =  maxval(distance(:),1, ffix(1,:)>0)
       k=0
       distance2=0d0
       do i = 1, na
          if(distance(i) > max(3d0,(r-3d0)) .and. ffix(1,i) > 0 .and. i/= atompair(1)) then
                    k=k+1
                    distance2(k)=dble(i)
          endif
       enddo
       if(k==0) then 
           atom_save=0
       else
           call rd_numb(y)
           j=int(y*dble(k))+1
           atom_save=int(distance2(j))
       endif
       !write(para%ioout,*) 'Findleast',cpuid,k
       atompair(2)   = atom_save

       deallocate (distance,distance2)

      end subroutine Find_leastmoveAtoms


    subroutine Get_atompair(na,cell,B,iza,ffix,atompair)

    use SSW_parameters

    implicit none
    integer na
    integer L
    double precision B(3,na),cell(3,3)
    integer iza(na)
    integer ffix(3,na),istep,atompair(2)

    integer atom_min, Nei,neighbor(100),i, kk
    double precision x,y,r
    double precision, dimension(:), allocatable :: distance
    integer , save :: ID,  atom_hist(50)
    data ID /0/, atom_hist /50*0/


!   write(para%ioout,*) 'iza3',iza

    ID=ID+1
    if(ID>min(50,na)) ID=1

    allocate (distance(na))

    atom_min=atompair(1)

    call rd_numb(x)

!    write(para%ioout,*) 'atom random',cpuid,x,B(1,1)
    if(x>0.5d0) then
      if(para%Run_type<10) then
       Nei=0
       call neighboringlist(na,B,atom_min,Nei,neighbor,iza,ffix)
!neighboringlist(xa,na,atom1,Nnei,neighbor,iza,ffix)
!      write(para%ioout,*) 'atom random2',Nei,x
       if(Nei>0) then
         call rd_numb(y)
         i=int(y*dble(Nei))+1
         atompair(1)=neighbor(i)
       endif
      else
       if(atompair(2)/=0)   atompair(1)= atompair(2)
      endif
    endif

    atom_min=atompair(1)
 
!   kk=0
!   do i = 1,min(50,na)
!       if (atom_hist(i)==atom_min) kk=kk+1
!   enddo
!   if(kk>min(50,na)*0.2 ) then
!     call random_number(x)
!     i=int(x*dble(na))+1
!     atompair(1) = i
!   endif
!   atom_hist(ID) = atompair(1)

    

    distance(:)= dsqrt( (B(1,:)-B(1,atom_min))*(B(1,:)-B(1,atom_min))+ &
                 (B(2,:)-B(2,atom_min))*(B(2,:)-B(2,atom_min))+ &
                 (B(3,:)-B(3,atom_min))*(B(3,:)-B(3,atom_min)) )
    r  = maxval(distance(:),1, distance>0.1) * 0.5d0
!   call random_number(x)
!   if(x>0.5d0) r=2.0d0
    if(para%Run_type<10) r=2.0d0
    if(para%Run_type>10) r=r+0.7d0

    istep=0
    do while (istep<150)
      call rd_numb(x)
      i=int(x*dble(na))+1
      if(distance(i) > r .and. ffix(1,i)/=0) then
          atompair(2) = i
          call check_forbiden(na,cell,B,atompair(1),atompair(2),1,L)
          if(L==0) cycle
          if(L==1) exit
          exit
      endif
      istep=istep+1
    enddo
!   if(istep>=150) atompair = 0

    deallocate(distance)

    end subroutine Get_atompair

    subroutine setatompair(atompair,iza1, At1, iza2, At2, Rcut,Rtype,na,cell,B,iza,ffix)


    implicit none
    integer na
    double precision B(3,na),cell(3,3),dist3,Rcut,x
    integer iza(na), iza1,At1,iza2,At2,atompair(2),Rtype
    integer ffix(3,na),k,i
    
    if(At1/=0 ) atompair(1) = At1
    if(At2/=0 ) atompair(2) = At2
    if(At1/=0 .and. At2/=0 ) return

!   write(6,*) 'atoms1', cpuid,At1, At2
    if(At1==0) then
        k=0
        do while (k< na*5)
           call rd_numb(x)
           i=int(x*dble(na))+1
           k=k+1
           if(iza(i) /=iza1 .or. ffix(1,i) ==0) cycle
           if(At2/=0) call get_dist(na,dist3,i,At2,B,cell)
           if(At2==0) then
               atompair(1)=i
               exit
           endif
           if(dist3 < 2.0d0 .and. Rtype==1 .and. iza1/=1 .and. iza2/=1) cycle
           if(dist3 < 1.6d0 .and. Rtype==1 .and. (iza1==1 .or. iza2==1) ) cycle
           if(dist3 < 1.1d0 .and. Rtype==1 .and. (iza1==1 .and. iza2==1) ) cycle
           if(dist3 < Rcut) then
                 atompair(1) = i
                 exit
           endif
!           k=k+1
         enddo

    endif

!   write(6,*) 'atoms2', cpuid,At1, At2

    if(At2==0) then
        k=0
        do while (k< na*5)
           call rd_numb(x)
           i=int(x*dble(na))+1
           k = k +1
           if(iza(i) /=iza2 .or. ffix(1,i) ==0) cycle
           if(At1/=0) then
                call get_dist(na,dist3,i,At1,B,cell)
           else
                call get_dist(na,dist3,i,atompair(1),B,cell)
            end if
           if(dist3 < 2.0d0 .and. Rtype==1 .and. iza1/=1 .and. iza2/=1) cycle
           if(dist3 < 1.6d0 .and. Rtype==1 .and. (iza1==1 .or. iza2==1) ) cycle
           if(dist3 < 1.1d0 .and. Rtype==1 .and. (iza1==1 .and. iza2==1) ) cycle
           if(dist3 < Rcut) then
                 atompair(2) = i
                 exit
           endif
!          k=k+1
         enddo
    endif

!   write(6,*) 'atoms3', At1, At2
    end subroutine setatompair

    subroutine check_forbiden(na,cell,B,atom1,atom2,RP,L)


    implicit none
    integer na,atom1,atom2,i,neig,neig2
    integer L
    double precision B(3,na),dist1,dist2,dist3,cell(3,3)
    double precision N1(3), N2(3) , coord1(3), angle
    integer RP

! make a list for the forbidden move

!   neig2=0
!    call get_dist(na,dist3,atom1,atom2,B,cell)
!    if(dist3 <3.5d0) neig2=1
!    if (dist3 <2.0d0) then
!         L= 0
!         return
!    endif
!    L = -1
!    neig=0
!    neig2=0
!    do i = 1, na
!      if (i==atom1 .or. i==atom2) cycle
!      call get_dist(na,dist1,atom1,i,B,cell)
!      call get_dist(na,dist2,atom2,i,B,cell)
!      if(  dist1 < 1.8d0 .and. dist2 < 1.80d0 .and. dist3 < dist1+dist2 ) then
!           neig=neig+1
!      endif
!      if(  dist2 < 3.5d0 ) then
!           neig2=neig2+1
!      endif
!      if ( ( dist1 < 1.25d0 .and. dist2 < 1.25d0 .and. (dist3 > dist1+dist2-0.1d0 .and. dist3 < dist1+dist2+0.1d0))  .or.  &
!           neig > 3    ) then
!          L = 0
!          exit
!      endif
!    enddo
!    if (L==0) return
!    if (neig2==0) L=1

     L = 1
     N1=  B(:,atom1) - B(:,atom2)
     do i = 1, na
       if (i==atom1 .or. i==atom2) cycle
       call get_dist(na,dist1,atom1,i,B,cell,coord1)
       if (dist1 < 2.2d0  ) then
          N2 = (B(:,atom1) - coord1(:))*dble(RP)
          angle = ACOS(SUM(N1*N2)/dsqrt(SUM(N1*N1))/dsqrt(SUM(N2*N2)))*180.0d0/3.141592653589793238d0
          if ( angle < 10d0) then
             L = 0 
             exit
          endif
          if ( angle < 60d0 .and. dist1 < 1.25d0) then
             L = 0 
             exit
          endif
          if ( angle < 30d0 .and. dist1 < 1.5d0) then
             L = 0 
             exit
          endif
       endif
     enddo
     if (L==0) return

     N1=  B(:,atom2) - B(:,atom1)
     do i = 1, na
       if (i==atom1 .or. i==atom2) cycle
       call get_dist(na,dist1,atom2,i,B,cell,coord1)
       if (dist1 < 2.2d0  ) then
          N2 = (B(:,atom2) - coord1(:)) * dble(RP)
          angle = ACOS(SUM(N1*N2)/dsqrt(SUM(N1*N1))/dsqrt(SUM(N2*N2)))*180.0d0/3.141592653589793238d0
          if ( angle < 10d0) then
             L = 0
             exit
          endif
          if ( angle < 60d0 .and. dist1 < 1.25d0) then
             L = 0 
             exit
          endif
          if ( angle < 30d0 .and. dist1 < 1.5d0) then
             L = 0 
             exit
          endif
       endif
     enddo

    endsubroutine

     SUBROUTINE SetConstraints_crystal(Nions,xa,N,fixlat)
!     use bpcbd, only : fixlat

      IMPLICIT NONE
      integer Nions
      double precision ,DIMENSION(3,Nions) :: N,xa
      double precision ,DIMENSION(:,:),allocatable ,save :: Nx,Ny,Nz
      double precision ,DIMENSION(:,:),allocatable ,save :: RotModeX,RotModeY,RotModeZ
      double precision com(3),axisatom(3,3),V(3,3),r
      integer, save :: IR
      logical,save :: first=.true.
      integer i,j,fixlat(3)




      if(first) then
         allocate(Nx(3,Nions),Ny(3,Nions),Nz(3,Nions))
         allocate(RotModeX(3,Nions),RotModeY(3,Nions),RotModeZ(3,Nions))
         Nx=0.0d0
         Ny=0.0d0
         Nz=0.0d0
         Nx(1,1:nions-3)=1.0d0
         Ny(2,1:nions-3)=1.0d0
         Nz(3,1:nions-3)=1.0d0
         Nx=Nx/sqrt(sum(Nx*Nx))
         Ny=Ny/sqrt(sum(Ny*Ny))
         Nz=Nz/sqrt(sum(Nz*Nz))
         first=.false.
      endif


      N=N-SUM(N*Nx)*Nx-SUM(N*Ny)*Ny-SUM(N*Nz)*Nz

      r=sum(N(1,nions-2:nions)**2)
      if(r<1d-6) return

      if (minval(fixlat)==0) then
          do i=Nions-2,Nions
             j=i-Nions+3
             if(fixlat(j)==0) then
                 N(j,Nions-2:Nions)=0d0
                 N(:,i)=0d0
             endif
          enddo
          return
      endif

      com(:)=0.0d0
!generate three axis
      RotModeX=0d0
      RotModeY=0d0
      RotModeZ=0d0
      do j=1,3
      do i=1,3
      if(i==j)then
        axisatom(i,j)=com(i)+1.0d0
      else
        axisatom(i,j)=com(i)
      endif
      enddo
      enddo
!generate three axis end
!generate rotational mode
      do i=nions-2,Nions
        do j=1,3
          V=0.0d0
          V(:,1)=xa(:,i)-com(:)
          V(:,2)=xa(:,i)-axisatom(:,j)
          V(1,3)=V(2,1)*V(3,2)-V(3,1)*V(2,2)
          V(2,3)=V(3,1)*V(1,2)-V(1,1)*V(3,2)
          V(3,3)=V(1,1)*V(2,2)-V(2,1)*V(1,2)
          if(j==1) RotModeX(:,i)=V(:,3)
          if(j==2) RotModeY(:,i)=V(:,3)
          if(j==3) RotModeZ(:,i)=V(:,3)
        enddo
      enddo
      if(sqrt(sum(RotModeX*RotModeX))>1d-150)then
        RotModeX=RotModeX/sqrt(sum(RotModeX*RotModeX))
      else
         RotModeX=0d0
      endif

      if(sqrt(sum(RotModeY*RotModeY))>1d-150)then
        RotModeY=RotModeY/sqrt(sum(RotModeY*RotModeY))
            else
       RotModeY=0d0
      endif

      if(sqrt(sum(RotmodeZ*RotmodeZ))>1d-150) then
        RotmodeZ=RotmodeZ/sqrt(sum(RotmodeZ*RotmodeZ))
      else
        RotModeZ=0d0
      endif

      RotModeY=RotModeY-sum(RotModeY*RotModeX )*RotModeX
      if(sqrt(sum(RotModeY*RotModeY))>1d-150) then
        RotModeY=RotModeY/sqrt(sum(RotModeY*RotModeY))
      else
        RotModeY=0d0
      endif

      RotmodeZ=RotmodeZ-sum(RotModeZ*RotModeX )*RotModeX &
                       -sum(RotModeZ*RotModeY )*RotModeY
      if(sqrt(sum(RotModeZ*RotModeZ))>1d-150) then
        RotmodeZ=RotmodeZ/sqrt(sum(RotmodeZ*RotmodeZ))
      else
        RotmodeZ=0d0
      endif

      N=N-sum(N*RotmodeX)*RotmodeX-sum(N*RotmodeY)*RotmodeY &
          -sum(N*RotmodeZ)*RotmodeZ


   END SUBROUTINE SetConstraints_crystal


   subroutine matsqrt(N,A,Asqrt)

    implicit none
    ! input number N x N matrix A , N the dimension
    ! output N x N upper trianglular matrix Asqrt

    integer :: N,info
    double precision :: A(N,N) , Asqrt(N,N)
    integer :: i,j


    Asqrt = A

    call dpotrf('U',N,Asqrt,N,info)

    if(info /= 0) print *,"matrix error! some where must be wrong!"

    do i = 1,N
        do j = 1,N
            if(i<j) Asqrt(j,i) = 0.0d0
        end do
    end do

    end subroutine





   subroutine regularcell(cell,fixlat,change, Lenlargecell)

!     use ssw_parameters

      implicit none
      double precision  :: cell(3,3),aplusb,zz
      integer    :: iv, i,j,k,fixlat(3),m
      double precision :: cellm(3),celang(3),ab(3)
      integer :: change(100,5)
      logical, optional :: Lenlargecell



      m=0
      change=0
2919     continue
      m=m+1
      if(m>100) return

        do iv = 1,3
          cellm(iv) = dot_product(cell(:,iv),cell(:,iv))
          cellm(iv) = dabs(cellm(iv))
        enddo
        i=minloc(cellm,1)
        j=maxloc(cellm,1)
        if(i==j) then
          i=1
          j=2
          k=3
        else
          k=6-i-j
        endif

        zz = dot_product(cell(:,i),cell(:,j))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=1
            change(m,2)=i
            change(m,3)=j
            change(m,5) = ceiling(dabs(zz)/cellm(i))
            goto 2919
        endif

        ab(:)=cell(:,j)-floor(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=3
            change(m,2)=i
            change(m,3)=j
            change(m,5) = floor(dabs(zz)/cellm(i))   
            goto 2919
        endif

        celang(:) = cell(:,i)+cell(:,k)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,j),celang(:))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=k
            change(m,4)=j
            change(m,5) = ceiling(dabs(zz)/aplusb)  
            goto 2919
        endif

        ab(:)=cell(:,j)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=k
            change(m,4)=j
            change(m,5) = floor(dabs(zz)/aplusb)     
            goto 2919
        endif

        celang(:) = cell(:,i)-cell(:,k)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,j),celang(:))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=-k
            change(m,4)=j
            change(m,5) = ceiling(dabs(zz)/aplusb)  
            goto 2919
        endif

        ab(:)=cell(:,j)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=-k
            change(m,4)=j
            change(m,5) = floor(dabs(zz)/aplusb)
            goto 2919
        endif

        celang(:) = cell(:,i)+cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,k),celang(:))
        ab(:)=cell(:,k)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14 .and. fixlat(k)==1)  then
            cell(:,k)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=j
            change(m,4)=k
           change(m,5) = ceiling(dabs(zz)/aplusb)
            goto 2919
        endif

        ab(:)=cell(:,k)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14 .and. fixlat(k)==1)  then
            cell(:,k)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=j
            change(m,4)=k
            change(m,5) = floor(dabs(zz)/aplusb)  
            goto 2919
        endif

        celang(:) = cell(:,i)-cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,k),celang(:))
        ab(:)=cell(:,k)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14 .and. fixlat(k)==1)  then
            cell(:,k)=ab
            change(m,1)=2
            change(m,2)=i
            change(m,3)=-j
            change(m,4)=k
             change(m,5) = ceiling(dabs(zz)/aplusb)
            goto 2919
        endif

        ab(:)=cell(:,k)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14 .and. fixlat(k)==1)  then
            cell(:,k)=ab
            change(m,1)=4
            change(m,2)=i
            change(m,3)=-j
            change(m,4)=k
            change(m,5) = floor(dabs(zz)/aplusb)
            goto 2919
        endif

        celang(:) = cell(:,k)+cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,i),celang(:))
        ab(:)=cell(:,i)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14 .and. fixlat(i)==1)  then
            cell(:,i)=ab
            change(m,1)=2
            change(m,2)=k
            change(m,3)=j
            change(m,4)=i
            change(m,5) = ceiling(dabs(zz)/aplusb)
            goto 2919
        endif

        ab(:)=cell(:,i)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14 .and. fixlat(i)==1)  then
            cell(:,i)=ab
           change(m,1)=4
            change(m,2)=k
            change(m,3)=j
            change(m,4)=i
            change(m,5) = floor(dabs(zz)/aplusb)
            goto 2919
        endif


        celang(:) = cell(:,k)-cell(:,j)
        aplusb=sum(celang*celang)
        aplusb= dabs(aplusb)
        zz = dot_product(cell(:,i),celang(:))
        ab(:)=cell(:,i)-ceiling(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14 .and. fixlat(i)==1)  then
            cell(:,i)=ab
            change(m,1)=2
            change(m,2)=k
            change(m,3)=-j
            change(m,4)=i
            change(m,5) = ceiling(dabs(zz)/aplusb)
            goto 2919
        endif

        ab(:)=cell(:,i)-floor(dabs(zz)/aplusb)*sign(1.0d0,zz)*celang(:)
        if( sqrt(sum(ab**2))<sqrt(cellm(i))-1d-14 .and. fixlat(i)==1 )  then
            cell(:,i)=ab
            change(m,1)=4
            change(m,2)=k
            change(m,3)=-j
            change(m,4)=i
            change(m,5) = floor(dabs(zz)/aplusb)
            goto 2919
        endif

        zz = dot_product(cell(:,i),cell(:,k))
        ab(:)=cell(:,k)-ceiling(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14 .and. fixlat(k)==1)  then
            cell(:,k)=ab
            change(m,1)=1
            change(m,2)=i
            change(m,3)=k
            change(m,5) = ceiling(dabs(zz)/cellm(i))
            goto 2919
        endif

        ab(:)=cell(:,k)-floor(dabs(zz)/cellm(i))*sign(1.0d0,zz)*cell(:,i)
        if( sqrt(sum(ab**2))<sqrt(cellm(k))-1d-14 .and. fixlat(k)==1)  then
            cell(:,k)=ab
            change(m,1)=3
            change(m,2)=i
            change(m,3)=k
            change(m,5) = floor(dabs(zz)/cellm(i))
            goto 2919
        endif

        zz = dot_product(cell(:,k),cell(:,j))
        ab(:)=cell(:,j)-ceiling(dabs(zz)/cellm(k))*sign(1.0d0,zz)*cell(:,k)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=1
            change(m,2)=k
            change(m,3)=j
            change(m,5) = ceiling(dabs(zz)/cellm(k))
            goto 2919
        endif

        ab(:)=cell(:,j)-floor(dabs(zz)/cellm(k))*sign(1.0d0,zz)*cell(:,k)
        if( sqrt(sum(ab**2))<sqrt(cellm(j))-1d-14 .and. fixlat(j)==1)  then
            cell(:,j)=ab
            change(m,1)=3
            change(m,2)=k
            change(m,3)=j
            change(m,5) = floor(dabs(zz)/cellm(k))
            goto 2919
        endif

!      if(present(Lenlargecell)) then
!      do i=1,3
!         if (fixlat(i)==0) then
!         cell(:,i) = cell(:,i)*1.5d0
!         endif
!      enddo
!      endif
!       if(L) then
!         call reci_latt(cell0,cell01)
!         cell_R=matmul(cell,transpose(cell01))
!         do i=1,3
!         write(para%ioout,'(A,3f8.3)') 'cell_R',cell_R(:,i)
!         enddo
!       endif


      end subroutine  regularcell



!=======================================================================
! Potential for LJ particles
      subroutine LJcorepot_old(na,iza,cell,xa0,fa,stress,E,Lperi,LJcore,LJcorePara)

      use SSW_parameters, only :  para

      implicit none
      integer na
      integer iza(na)
      double precision xa0(3,na),fa(3,na),stress(3,3),cell(3,3)
      double precision E,LJcut
      double precision celli(3,3),xa2(3)
      double precision, allocatable, save :: xacna(:,:), xa(:,:)
      double precision autoev,autoang,conf,dconf

!      parameter(autoev=   27.2113961d0)
!      parameter(autoang=0.529177249d0)
       parameter(autoev=   1d0)
       parameter(autoang=  1d0)

      double precision sig

      integer i,j,m,mp
      double precision dr,dvdr,dx,dy,dz
      double precision rij,rtmp(3),rtmp1(3)
!      double precision sig6,sig12
      double precision dr6,dr12, Aterm

      logical  , save :: first = .true.
      logical  Lperi
      integer lp,l1,l2,l3,k
      !external volcel_loc
!      data lp /10/
      integer, allocatable :: nr12(:), NN(:),Neigh_name(:,:)
      double precision, dimension(:), allocatable ::  LJcutoff , eps , sig12
      double precision, dimension(:,:), allocatable ::  neighborx,neighbory,neighborz

      integer :: LJcore 
      double precision :: LJcorePara(6,LJcore) ,vij(3),vik(3),costheta0

      save nr12,  LJcutoff , eps , sig12, NN, Neigh_name,neighborx,neighbory,neighborz

      !write(para%ioout,*) 'first',first
      if(first) then
        allocate(xacna(3,na),xa(3,na))
        allocate(nr12(LJcore), LJcutoff(LJcore) , eps(LJcore) , sig12(LJcore))
        allocate(NN(na))
        allocate(Neigh_name(na,20))
        allocate(neighborx(na,20),neighbory(na,20),neighborz(na,20))

        first=.false.
        if(cpuid==0) write(para%ioout,*) 'Total LJ Parameters',LJcore

        do i=1,LJcore
           nr12(i)=LJcorePara(3,i)
           eps(i)=LJcorePara(4,i)
           Aterm=LJcorePara(5,i)
           LJcutoff(i) = LJcorePara(6,i)
           sig12(i)=Aterm**nr12(i)
        enddo
      endif

      xa= xa0 !/autoang
  
      fa=0d0;stress=0d0;E=0d0

!     write(para%ioout,*) xa
!     write(para%ioout,*) 'cell',cell

      call reci_latt(cell, celli)
      do i=1,na
        xacna(:,i) = matmul(transpose(celli),xa0(:,i))
        xacna(:,i)= modulo(xacna(:,i)+1000.0d0,1.0d0)
      enddo
      xa = matmul(cell,xacna)


!     write(para%ioout,*) xa
!     write(para%ioout,*) 'cell',cell

      LJcut=maxval(LJcutoff(1:LJcore))

  !   do i=1,na
  !   do j=i+1,na
  !     dx=xa(1,i)-xa(1,j)
  !     dy=xa(2,i)-xa(2,j)
  !     dz=xa(3,i)-xa(3,j)
  !     rij=dsqrt(dx*dx+dy*dy+dz*dz)

  !     dr = 1.d0/rij !  (i,j)
  !     !dr6 = dr**nr6
  !     if(rij > LJcut) cycle

  !       do mp=1,LJcore
  !         if( (iza(i)==int(LJcorePara(1,mp)) .and. iza(j)==int(LJcorePara(2,mp)) ) .or. &
  !             (iza(i)==int(LJcorePara(2,mp)) .and. iza(j)==int(LJcorePara(1,mp)) ) .or. &
  !             (LJcorePara(1,mp)==0 .and. LJcorePara(2,mp)==0) ) then
! !             write(para%ioout,*) 'iza',iza(i),iza(j),int(LJcorePara(1,mp)),int(LJcorePara(2,mp))
  !             if(rij > LJcutoff(mp)) cycle
! !             write(para%ioout,*) 'inside LJ', rij,i,j

  !           call tanh3(rij,LJcutoff(mp),conf,dconf)
  !           dr12 = dr**nr12(mp) ! dr6*dr6
  !           E = E + 4.d0*eps(mp)*(sig12(mp)*dr12) *conf !-sig6*dr6)
  !           !write(para%ioout,*) 'E',E
  !           dvdr = 4.d0*eps(mp)*(-nr12(mp)*sig12(mp)*dr12*dr )*conf + 4.d0*eps(mp)*(sig12(mp)*dr12) *dconf !+ nr6*sig6*dr6*dr)
  !           fa(1,i) = fa(1,i) + dvdr*dx*dr    !dvdx(i) = dvdx(i) + dvdr*dx*dr
  !           fa(1,j) = fa(1,j) - dvdr*dx*dr
  !           fa(2,i) = fa(2,i) + dvdr*dy*dr    !dvdx(i) = dvdx(i) + dvdr*dx*dr
  !           fa(2,j) = fa(2,j) - dvdr*dy*dr
  !           fa(3,i) = fa(3,i) + dvdr*dz*dr    !dvdx(i) = dvdx(i) + dvdr*dx*dr
  !           fa(3,j) = fa(3,j) - dvdr*dz*dr
  !           stress(1,1)= stress(1,1)+ dx*dx*dr*dvdr
  !           stress(2,1)= stress(2,1)+ dy*dx*dr*dvdr
  !           stress(3,1)= stress(3,1)+ dz*dx*dr*dvdr
  !           stress(2,2)= stress(2,2)+ dy*dy*dr*dvdr
  !           stress(3,2)= stress(3,2)+ dz*dy*dr*dvdr
  !           stress(3,3)= stress(3,3)+ dz*dz*dr*dvdr
  !           exit
  !         endif
  !       enddo !mp

  !   enddo
  !   enddo

!     write(para%ioout,*) 'XXXXX'
      NN=0
      do i=1,na
!       do l1=-int(LJcut/rtmp(1))-1,int(LJcut/rtmp(1))+1
!         do l2=-int(LJcut/rtmp(2))-1,int(LJcut/rtmp(2))+1
!          do l3=-int(LJcut/rtmp(3))-1,int(LJcut/rtmp(3))+1
!       write(para%ioout,*) i,NN(i)

        do l1=-1,1
          do l2=-1,1
           do l3=-1,1
       !   l1=0; l2=0; l3=0

            do j = 1,na

!           write(para%ioout,*) 'ij',i,j
            if(l1==0 .and. l2==0 .and. l3==0 .and. i==j ) cycle

              xa2(1)=xacna(1,j)+l1
              xa2(2)=xacna(2,j)+l2
              xa2(3)=xacna(3,j)+l3
              xa2 = matmul(cell,xa2)

!             do m = 1,3
!                xa2(m) = cell(m,1) * (xacna(1,j)+l1) + &
!                         cell(m,2) * (xacna(2,j)+l2) + &
!                         cell(m,3) * (xacna(3,j)+l3)                   ! xa upper angle cartessian
!             enddo
              dx=xa(1,i)-xa2(1)
              dy=xa(2,i)-xa2(2)
              dz=xa(3,i)-xa2(3)
              rij=dsqrt(dx*dx+dy*dy+dz*dz)
!             if(l1==0 .and. l2==0 .and. l3==0) write(para%ioout,*) i,j,rij

      !       if(rij < 2.0d0) then     ! modified later
      !           NN(i) = NN(i) +1
      !           Neigh_name(i,NN(i)) = j
      !           neighborx(i,NN(i))= -dx
      !           neighbory(i,NN(i))= -dy
      !           neighborz(i,NN(i))= -dz
!     !           write(para%ioout,*) i,j,rij,NN(i)
      !       endif

              if(rij > LJcut) cycle

               dr = 1.d0/rij !  (i,j)
               do mp=1,LJcore
                 if(rij > LJcutoff(mp)) cycle
                 if( (iza(i)==int(LJcorePara(1,mp)) .and. iza(j)==int(LJcorePara(2,mp)) ) .or. &
                   (iza(i)==int(LJcorePara(2,mp)) .and. iza(j)==int(LJcorePara(1,mp)) ) .or. &
                   (int(LJcorePara(1,mp))==0 .and. int(LJcorePara(2,mp))==0) ) then
                   !write(*,*) 'inmp',rij,LJcutoff(mp)/autoang
!                   write(para%ioout,*) 'LJinside',rij,mp, LJcutoff(mp)/autoang
                   dr12 = dr**nr12(mp) ! dr6*dr6
                   call tanh3(rij,LJcutoff(mp),conf,dconf)
                   E = E + 2.d0*eps(mp)*(sig12(mp)*dr12)*conf !-sig6*dr6)
                   dvdr = 2.d0*eps(mp)*(-nr12(mp)*sig12(mp)*dr12*dr )*conf + 2.d0*eps(mp)*(sig12(mp)*dr12)*dconf   !+ nr6*sig6*dr6*dr)
                   fa(1,i) = fa(1,i) + dvdr*dx*dr    !dvdx(i) = dvdx(i) + dvdr*dx*dr
                   fa(1,j) = fa(1,j) - dvdr*dx*dr
                   fa(2,i) = fa(2,i) + dvdr*dy*dr    !dvdx(i) = dvdx(i) + dvdr*dx*dr
                   fa(2,j) = fa(2,j) - dvdr*dy*dr
                   fa(3,i) = fa(3,i) + dvdr*dz*dr    !dvdx(i) = dvdx(i) + dvdr*dx*dr
                   fa(3,j) = fa(3,j) - dvdr*dz*dr
                   stress(1,1)= stress(1,1)+ dx*dx*dr*dvdr
                   stress(2,1)= stress(2,1)+ dy*dx*dr*dvdr
                   stress(3,1)= stress(3,1)+ dz*dx*dr*dvdr
                   stress(2,2)= stress(2,2)+ dy*dy*dr*dvdr
                   stress(3,2)= stress(3,2)+ dz*dy*dr*dvdr
                   stress(3,3)= stress(3,3)+ dz*dz*dr*dvdr
                   exit
                 endif
               enddo !m

             enddo   ! j

           enddo
          enddo
        enddo

       enddo   ! i


!100   continue

!     write(para%ioout,*) 'YYYYYYY'
!     cell=cell*autoang
      fa= -fa !* autoev/autoang !dvdx=dvdx*autoev/autoang

    ! costheta0 = cos(35d0/180d0*3.14159265d0)
    ! do i=1,na
    !   if(NN(i)<2) cycle
    !   do j=1,NN(i)-1
    !      do k = j+1,NN(i)
    !          vij(1)=neighborx(i,j)
    !          vij(2)=neighbory(i,j)
    !          vij(3)=neighborz(i,j)
    !          vik(1)=neighborx(i,k)
    !          vik(2)=neighbory(i,k)
    !          vik(3)=neighborz(i,k)
    !          call repelforce (vij, vik,E, fa(:,i),fa(:,Neigh_name(i,j)), fa(:,Neigh_name(i,k)),stress, costheta0, 2.0d0)
    !      enddo
    !   enddo
    !    
    ! enddo

!     fa= -fa !* autoev/autoang !dvdx=dvdx*autoev/autoang
!     E=E !*autoev
!     stress=stress*autoev*0.5d0/volcel_loc(cell)
      stress=stress/volcel_loc(cell)
      stress(1,2)=stress(2,1)
      stress(1,3)=stress(3,1)
      stress(2,3)=stress(3,2)

      return
      end subroutine


      subroutine tanh3(d,dc,f,df)

      double precision d,dc,f,df, x,x2

      x = tanh(1.0d0-d/dc)
      x2= x*x
      f= x2*x *d**2
      !df = 1.5d0 * (1.0d0-tanh(1.0d0-d/dc)**2) *(-1.0d0/dc)
      df = -3.0d0*x2*(1.0d0-x2)/dc*d**2 + x2*x*2d0*d

      !write(para%ioout,*) d,dc,f,df



      end subroutine


!     subroutine tanh3(d,dc,f,df)
!
!     double precision d,dc,f,df
!
!
!     f = 0.5d0 * tanh(1.0d0-d/dc)**3
!     df = 1.5d0 * (1.0d0-tanh(1.0d0-d/dc)**2) *(-1.0d0/dc)
!
!     !write(para%ioout,*) d,dc,f,df
!
!
!
!     end subroutine


      subroutine readuncm(unit,na,iza,cell,xa, Lend)

      use constants
      use ssw_parameters,only:para

      implicit          none

      integer                             :: iza(na)
      double precision                    :: xa(3,na)
      integer                             :: unit
      double precision                    :: cell(3,3)

      integer                             :: i, ia,j
      double precision                    :: latt_abc(3)
      double precision                    :: latt_angle(3)
      double precision                    :: celli(3,3)
      double precision                    :: celltran(3,3)
      double precision                    :: angle(3)
      double precision                    :: xfrac(3)
      double precision                    :: xa_out(3)
! -------------------------------------------------------------------
      character(3)::pbc
      character(2)::lab4
      character(80)::xxx
      integer na
      double precision::abc(3),alplp,betlp,gamlp,alp,blp,clp
      logical Lend


      Lend = .false.

 !    write(para%ioout,*) unit   
      read(unit,*,end=1731,err=1731) xxx
 !    write(para%ioout,*) 'xxxx', xxx
      read(unit,*,end=1731,err=1731)
 !    write(para%ioout,*) 'yyyy'
      read(unit,*,end=1731,err=1731) pbc,abc(:),angle(3),angle(2),angle(1)

 !    write(para%ioout,*) 'in uncm',pbc, abc(:),angle(3),angle(2),angle(1)
  
      cell=0.0d0
      alplp=angle(3)* pi/180.d0
      betlp=angle(2)* pi/180.d0
      gamlp=angle(1)* pi/180.d0
      alp=abc(1)
      blp=abc(2)
      clp=abc(3)
      cell(1,1) = alp
      cell(1,2) = blp * dcos(gamlp)
      cell(2,2) = blp * dsin(gamlp)
      cell(1,3) = clp * dcos(betlp)
      cell(2,3) = (clp*dcos(alplp) - clp*dcos(gamlp)*dcos(betlp))/dsin(gamlp)
      cell(3,3) = dsqrt(clp*clp - cell(1,3)*cell(1,3) - cell(2,3)*cell(2,3))
     
      do i=1,na
          read(unit,*,end=1731,err=1731) lab4,xa_out(1:3)
          xa(:,i) = xa_out(1:3)
          iza(i)=elemz(lab4)
      !   write(para%ioout,*) lab4,xa_out(1:3)
      enddo
      read(unit,*,end=1731,err=1731)
      read(unit,*,end=1731,err=1731)
      Lend = .false.
      return
1731  continue
      if(cpuid==0) write(para%ioout,*) 'No more pairs from uncm.arc -- Job Done'
      Lend =.true.
      return
      end subroutine readuncm


      double precision function velo_loc(iseed,temp,mass)
      implicit none
      integer iseed
      double precision mass,temp
      double precision  arg1, arg2, gauss, med, pi,  var

      med = 0.0
      pi = acos(-1.0d0)
      var = sqrt(temp/mass)
      var = var * 0.00172309d0

      arg1 = sqrt((-2.) * log(ran3(iseed)))

      arg2 = 2.0d0 * pi * ran3(iseed)
      gauss = arg1 * cos(arg2)

      velo_loc = med + var * gauss

      end function velo_loc



      subroutine vmb2(nat,ttemp,va,ffix,x)

      implicit none
      integer nat

      double precision ttemp, va(3,nat)

! Internal variables ..................

      integer iseed, i, ix, iy, ntcon,ffix(3,nat)

      double precision &
        cell(3,3), stress(3,3), fa(3,nat), &
        massi, tempe,  vtot(3),x

      !external velo_loc

        do ix=1,3
          vtot(ix)=0.0d0
        enddo


      iseed = -int(1.0d0*10d0*(x+1.0d0))


! Loop over atoms to assing velo_loccities .................
        do i = 1,nat
          massi = 1.0d0

          if(ffix(1,i) /=0) va(1,i) = velo_loc(iseed,ttemp,massi)
          if(ffix(2,i) /=0) va(2,i) = velo_loc(iseed,ttemp,massi)
          if(ffix(3,i) /=0) va(3,i) = velo_loc(iseed,ttemp,massi)

          vtot(1) = vtot(1) + va(1,i)
          vtot(2) = vtot(2) + va(2,i)
          vtot(3) = vtot(3) + va(3,i)
        enddo
! ...............

! Impose constraint on zero center of mass velo_loccity ....................
        do i = 1,nat
          do ix=1,3
            va(ix,i) = va(ix,i) - vtot(ix)/nat
          enddo
        enddo





      return
      end subroutine

   subroutine BainStrain(h0,h1,lamda,eigv) !,lamda2,eigv2)
   use ssw_parameters,only:para


   implicit none
   double precision , intent(in) :: h0(3,3),h1(3,3)
   double precision  :: lamda(3),eigv(3,3),Mstrain(3,3) !,lamda2(3),eigv2(3,3)
   integer INFO,i
   double precision ha(3,3),hb(3,3),ALPHAR(3),ALPHAI(3),BETA(3),VL(3,3),VB(3,3),WORK(24),UNI(3,3)
   double precision hc(3,3),VR(3,3)
   integer NDIM

    NDIM=3

    hb=h1
    call reci_latt(hb,hc)
    hc=transpose(hc)
    UNI=0d0
    UNI(1,1)=1d0
    UNI(2,2)=1d0
    UNI(3,3)=1d0


    ha=transpose(h0)
    ha=matmul(ha,h0)
    !do i=1,3
    ! write(para%ioout,'(A,3F10.4)') 'FS',ha(:,i)
    !enddo


    !Mstrain=0.50d0*(matmul(matmul(transpose(hc),ha),hc)-UNI)
     Mstrain=matmul(matmul(transpose(hc),ha),hc)

       ! hb=transpose(h1)
       ! hb=matmul(hb,h1)

    !do i=1,3
    ! write(para%ioout,'(A,3F10.4)') 'IS',hb(:,i)
    !enddo

    ! call d3minv(hb,hc)
    ! ha=matmul(hc,ha)

    ha=Mstrain
    do i=1,3
     if(cpuid==0) write(para%ioout,'(A,3F10.4)') 'Bain',ha(:,i)
    enddo

!   Mstrain=0.5d0*(Mstrain-UNI)
    ! ha=Mstrain

    call DGEGV('N','V',NDIM,ha,NDIM,UNI,NDIM,ALPHAR,ALPHAI,BETA,VL,NDIM,VR,NDIM,WORK,8*NDIM,INFO)
     do i=1,3
       if(cpuid==0) write(para%ioout,'(A,F10.6)') 'Eigen',alphar(i)/beta(i)
     enddo
    do i=1,3
     VR(:,i)=VR(:,i)/dsqrt(sum(VR(:,i)**2))
      if(cpuid==0) write(para%ioout,'(A,3F10.4)') 'Eigen',VR(:,i)
    enddo
    lamda(:)=alphar(:)/beta(:)
    eigv=VR

!    call reci_latt(h1, hc)
!    ha=matmul(h0,transpose(hc))
!
!   call DGEGV('N','V',NDIM,ha,NDIM,UNI,NDIM,ALPHAR,ALPHAI,BETA,VL,NDIM,VR,NDIM,WORK,8*NDIM,INFO)
!   do i=1,3
!    VR(:,i)=VR(:,i)/dsqrt(sum(VR(:,i)**2))
!   enddo
!
!
!
!   lamda2(:)=alphar(:)/beta(:)
!   eigv2=VR



   end subroutine

   subroutine crossproduct(a,b,c)
   implicit none
   double precision, dimension(3) :: a,b,c

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

   return
   end subroutine crossproduct


      subroutine get_dist(na,dist,atom1,atom2,xa,cell,coord)

      implicit none

      integer :: na,atom1,atom2
      double precision,dimension(3,na)  :: xa
!     double precision, allocatable, save :: xa2(:,:),xacna(:,:)
      double precision  xa2(3,2),xacna(3,2)
      double precision cell(3,3),celli(3,3),dist,r,str(3)
      double precision, optional :: coord(3)

      integer                           :: i,j,i1,i2,i3,k

         call reci_latt(cell, celli)
!        do i = 1,na
!           xacna(:,i) = matmul(transpose(celli),xa(:,i))
!           xacna(:,i)= modulo(xacna(:,i)+1000.0d0,1.0d0)
!           do j = 1,3
!             xa(j,i) = cell(j,1) * (xacna(1,i)) + &
!                       cell(j,2) * (xacna(2,i)) + &
!                       cell(j,3) * (xacna(3,i))                   ! xa upper angle cartessian
!           enddo
!        enddo
         xacna(:,1)= matmul(transpose(celli),xa(:,atom1))
         xacna(:,2)= matmul(transpose(celli),xa(:,atom2))
         xacna(:,1)= modulo(xacna(:,1)+1000.0d0,1.0d0)
         xacna(:,2)= modulo(xacna(:,2)+1000.0d0,1.0d0)
         xa2(:,1)= matmul(cell,xacna(:,1))

         dist=1000000d0
         do i1=-1,1
          do i2=-1,1
           do i3=-1,1

            do j = 1,3
              xa2(j,2) = cell(j,1) * (xacna(1,2)+i1) + &
                         cell(j,2) * (xacna(2,2)+i2) + &
                         cell(j,3) * (xacna(3,2)+i3)                   ! xa upper angle cartessian
            enddo
              r= dsqrt(sum((xa2(:,1)-xa2(:,2))**2))
              if(r<dist) then
                dist=r
                if(present(coord) ) coord=xa2(:,2)
              endif
            enddo
           enddo
         enddo

      end subroutine


      subroutine check_bond(na,iza,cell,xa1,xa2,lsame)

!     use bpcbd, only : Reactatom,Reactat
      use ssw_parameters,only:para

      implicit none

      integer :: na
      integer,dimension(na)             :: iza
      double precision,dimension(3,na)  :: xa1,xa2

      double precision cell(3,3)
      integer,dimension(:,:),allocatable,save     :: conn1,conn2

      integer                           :: i,j,k,nbond0,nbond1
      logical                           :: lsame,L, onlyOH
      logical, save :: first=.true.

      if(first) then
        allocate(conn1(na,na),conn2(na,na))
        first=.false.
      endif
      conn1=0
      conn2=0

      call get_conn(na,iza,xa1,conn1,cell)
      call get_conn(na,iza,xa2,conn2,cell)

      lsame = .false.
      ! reaction only taken between O and H
      onlyOH = .true.

      do i = 1,na
!        L=.false.
!        do k=1,Reactatom
!         if(Reactat(k)==i) L=.true.
!        enddo
!        if( Reactatom/=0 .and. .not.L) cycle
         do j =i+1,na
          !  if(i==j) cycle
!           L=.false.
!           do k=1,Reactatom
!           if(Reactat(k)==i) L=.true.
!           enddo

!           if( Reactatom/=0 .and. .not.L) cycle

            if(conn1(i,j) /= conn2(i,j)) then
                 if(conn1(i,j) == 0) then
                     if(cpuid==0) write(para%ioout,*) "bond formation",i,j
                 else
                     if(cpuid==0) write(para%ioout,*) "bond break",i,j
                 end if

                 if( (iza(i) == 1 .and. iza(j) ==8) .or. &
                    (iza(j) == 1 .and. iza(i) ==8)) then
                 else
                      ! not O-H bond formation
                      onlyOH = .false.
                 end if
                 lsame = .true.   ! structure changed
            end if
         end do
      enddo
      nbond0=0
      nbond1=0
      do i = 1,na
         do j =1,na-1
            if((iza(i)==8 .and. iza(j)==1) .or.&     !OH bond
               (iza(i)==1 .and. iza(j)==8) ) then
               if(conn1(i,j)==1) nbond0=nbond0+1
               if(conn2(i,j)==1) nbond1=nbond1+1
            end if
         end do
      enddo
      if(lsame .and. nbond0==nbond1 .and. nbond0>0 .and. onlyOH) lsame=.false.


      end subroutine

      subroutine check_bond_crystal(na,iza,cell1,cell2,xa1,xa2,lsame)
      use ssw_parameters,only:para

!     use bpcbd, only : Reactatom,Reactat

      implicit none

      integer :: na
      integer,dimension(na)             :: iza
      double precision,dimension(3,na)  :: xa1,xa2

      double precision ::           cell1(3,3),cell2(3,3)
      integer,dimension(:,:),allocatable,save     :: conn1,conn2

      integer                           :: i,j,k,nbond0,nbond1
      logical                           :: lsame,L, onlyOH
      logical, save :: first=.true.

      if(first) then
        allocate(conn1(na,na),conn2(na,na))
        first=.false.
      endif
      conn1=0
      conn2=0

      call get_conn(na,iza,xa1,conn1,cell1)
      call get_conn(na,iza,xa2,conn2,cell2)

      lsame = .false.
      ! reaction only taken between O and H
      onlyOH = .true.

      do i = 1,na
         do j =i+1,na
          !  if(i==j) cycle
!           L=.false.
!           do k=1,Reactatom
!           if(Reactat(k)==i) L=.true.
!           enddo

!           if( Reactatom/=0 .and. .not.L) cycle

            if(conn1(i,j) /= conn2(i,j)) then
                 if(conn1(i,j) == 0) then
                     if(cpuid==0) write(para%ioout,*) "bond formation",i,j
                 else
                     if(cpuid==0) write(para%ioout,*) "bond break",i,j
                 end if

                 if((iza(i) /= 1 .and. iza(i) /=8) .or. &
                    (iza(j) /= 1 .and. iza(j) /=8)) then
                      ! not O-H bond formation
                      onlyOH = .false.
                 end if
                 lsame = .true.   ! structure changed
            end if
         end do
      enddo
      nbond0=0
      nbond1=0
      do i = 1,na
         do j =1,na-1
            if((iza(i)==8 .and. iza(j)==1) .or.&     !OH bond
               (iza(i)==1 .and. iza(j)==8) ) then
               if(conn1(i,j)==1) nbond0=nbond0+1
               if(conn2(i,j)==1) nbond1=nbond1+1
            end if
         end do
      enddo
      if(lsame .and. nbond0==nbond1 .and. nbond0>0 .and. onlyOH) lsame=.false.


      end subroutine


     subroutine check_slab(na,iza,cell,xa1,xa,lsame)
      use ssw_parameters,only:para

      implicit none

      integer :: na
      integer,dimension(na)             :: iza
      double precision,dimension(3,na)  :: xa1,xa
      double precision, allocatable, save :: xa2(:,:),xacna(:,:)
      double precision, allocatable, save :: xa0(:,:),xacna0(:,:)
      double precision cell(3,3),celli(3,3),distance,dtot,d

      integer,dimension(na,na)          :: conn1,conn2

      integer                           :: i,j,i1,i2,i3,k,imax
      logical                           :: lsame
      logical,save :: first=.true.

      if(first) then
        allocate(xa2(3,na),xacna(3,na))
        allocate(xa0(3,na),xacna0(3,na))
        first=.false.
      endif


         xa2=xa
         xa0=xa1
         lsame = .false.

         distance=dsqrt(sum((xa1-xa2)**2))
!         write(para%ioout,*) 'Structure diff==-1',distance
         if(distance <0.1d0) return

         call reci_latt(cell, celli)
         do i = 1,na
            xacna(:,i) = matmul(transpose(celli),xa2(:,i))
            xacna(:,i)= modulo(xacna(:,i)+100.0d0,1.0d0)
            xacna0(:,i) = matmul(transpose(celli),xa0(:,i))
            xacna0(:,i)= modulo(xacna0(:,i)+100.0d0,1.0d0)
         enddo

         do k = 1,na
           do j = 1,3
             xa0(j,k) = cell(j,1) * xacna0(1,k) + cell(j,2) * xacna0(2,k) + cell(j,3) * xacna0(3,k)
             enddo
        enddo



         dtot=0d0
         d=0d0
         do k = 1,na

           distance=1000000d0
            do i1=-1,1
             do i2=-1,1
              do i3=-1,1

               do j = 1,3
                 xa2(j,k) = cell(j,1) * (xacna(1,k)+i1) +  cell(j,2) * (xacna(2,k)+i2) + cell(j,3) * (xacna(3,k)+i3)
               enddo
              distance=min(distance, sum((xa1(:,k)-xa2(:,k))**2))
              if(dsqrt(distance)<0.2d0) distance=0d0

               enddo
              enddo
            enddo
           if(distance>0.5d0) then
               lsame=.true.
           endif
           if(d<distance) then
                imax=k
                d=distance
           endif
           dtot=dtot+distance
         enddo

          dtot=dsqrt(dtot)
          d=dsqrt(d)
          if(cpuid==0) write(para%ioout,*) 'Structure diff',dtot
          if(lsame .or. (dtot > 5d0 .and. d>0.15d0) ) then ! .and. dtot < 6d0 ) then
             if(cpuid==0) write(para%ioout,'(A,2f8.2)') 'Slab Structure Changed', dtot,d
             lsame=.true.
          endif

         ! if(lsame) then  order parameter
         ! endif

      end subroutine

      subroutine get_conn(na,iza,xa,conn,cell)
      implicit none

      integer :: na
      integer,dimension(na)             :: iza
      double precision,dimension(3,na)  :: xa
      integer,dimension(na,na)          :: conn
      integer                           :: i,j
      double precision  :: dist,co1,co2,totbond,cell(3,3)

      conn = 0
      do i = 1,na
          do j = 1,i-1
              ! dist = dsqrt(sum((xa(:,i)-xa(:,j))**2))*0.529d0
              call get_dist(na,dist,i,j,xa,cell) ! periodic system
              call fastbond(iza(i),iza(j),totbond)
              if(totbond < 0.05d0) then
                  call species_radius(iza(i),co1)
                  call species_radius(iza(j),co2)
                  totbond = co1 + co2
              end if
! define rules for counting bond
              if(dist < totbond + 0.25d0) then
                 if(iza(i)<18 .and. iza(j) <18 ) then

!             if( (dist < totbond + 0.25d0 &
!                    .and. iza(i)/=1 .and. iza(j)/=1   &
!                    .and. ((iza(i)<18 .and. iza(j) <18) &
!                    .or. (iza(i)==35 .or. iza(j)==35))) & ! count all non-H bonding, C-C,C-O
!                .or. (dist < totbond + 0.25d0  &
!                    .and. iza(i)==1 .and. iza(j)==1 ) &! count H-H bond
!                .or. (dist < totbond + 0.25d0   &
!                    .and. (iza(i)==6 .or. iza(j)==6) )  ) then ! count all C-H bond

                 ! write(para%ioout,*) "connection",i,j,dist
                  conn(i,j) = 1
                  conn(j,i) = 1
              endif
             endif
!              write(para%ioout,*) i,j,dist
          end do
      end do
      !write(*,*) "==========================================="
      end subroutine

      subroutine fastbond(zi,zj,bondlen)

      double precision,dimension(53,53),save :: bondmap
      integer ::  zi,zj,i,j
      double precision :: bondlen
      logical ,save :: init=.true.

      if(init) then
        bondmap=0d0
        bondmap(1,1) = 0.80
        bondmap(1,5) = 1.084
        bondmap(1,6) = 1.084
        bondmap(1,7) = 1.001
        bondmap(1,8) = 0.947
        bondmap(1,9) = 0.92
        bondmap(1,14) = 1.48
        bondmap(1,15) = 1.415
        bondmap(1,16) = 1.326
        bondmap(1,17) = 1.28
        bondmap(1,35) = 1.41
        bondmap(1,53) = 1.60
        bondmap(5,5) = 1.70
        bondmap(5,8) = 1.512
        bondmap(6,6) = 1.512
        bondmap(6,7) = 1.439
        bondmap(6,8) = 1.393
        bondmap(6,9) = 1.353
        bondmap(6,14) = 1.86
        bondmap(6,15) = 1.84
        bondmap(6,16) = 1.812
        bondmap(6,17) = 1.781
        bondmap(6,35) = 1.94
        bondmap(6,53) = 2.16
        bondmap(7,7) = 1.283
        bondmap(7,8) = 1.333
        bondmap(7,9) = 1.36
        bondmap(7,14) = 1.74
        bondmap(7,15) = 1.65
        bondmap(7,16) = 1.674
        bondmap(7,17) = 1.75
        bondmap(7,35) = 1.90
        bondmap(7,53) = 2.10
        bondmap(8,8) = 1.45
        bondmap(8,9) = 1.42
        bondmap(8,14) = 1.63
        bondmap(8,15) = 1.66
        bondmap(8,16) = 1.470
        bondmap(8,17) = 1.70
        bondmap(8,22) = 2.1   ! Ti-O
        bondmap(8,35) = 1.85
        bondmap(8,53) = 2.05
        bondmap(9,14) = 1.57
        bondmap(9,15) = 1.54
        bondmap(9,16) = 1.55
        bondmap(14,8) = 2.0
        bondmap(16,35) = 2.24
        bondmap(16,53) = 2.40
        bondmap(17,17) = 1.99
        bondmap(35,35) = 2.28
        bondmap(53,53) = 2.67
        init = .false.
      end if

      if(zi > 53 .or. zj > 53) then
        bondlen = 0.0d0
      else
        if(zi<zj) then
            bondlen = bondmap(zi,zj)
        else
            bondlen = bondmap(zj,zi)
        end if
      end if
      end subroutine



     subroutine regularcell_pattern(cell,change)

       implicit none
       integer m, i,j,k,iv
       double precision cell(3,3), celang(3)
       double precision zz,aplusb
       integer change(100,5)

        m=1
        do while (change(m,1)/=0)
        ! do iv = 1,3
        !   cellm(iv) = dot_product(cell(:,iv),cell(:,iv))
        !   cellm(iv) = dabs(cellm(iv))
        ! enddo

!        write(para%ioout,'(A,5I4)') 'in regularpattern', change(m,:)

         if(change(m,1)==1) then
            j=change(m,3)
            i=change(m,2)
            zz = dot_product(cell(:,i),cell(:,j))
            cell(:,j)=cell(:,j)-change(m,5)*sign(1.0d0,zz)*cell(:,i)
         elseif(change(m,1)==2) then
           i=change(m,2)
           k=change(m,3)
           j=change(m,4)
           if(k>0) then
             celang(:) = cell(:,i)+cell(:,k)
           else
             celang(:) = cell(:,i)-cell(:,-k)
           endif
           aplusb=sum(celang*celang)
           aplusb= dabs(aplusb)
           zz = dot_product(cell(:,j),celang(:))
           cell(:,j)=cell(:,j)-change(m,5)*sign(1.0d0,zz)*celang(:)
         elseif(change(m,1)==3) then
           j=change(m,3)
           i=change(m,2)
           zz = dot_product(cell(:,i),cell(:,j))
           cell(:,j)=cell(:,j)-change(m,5)*sign(1.0d0,zz)*cell(:,i)
         elseif(change(m,1)==4) then
           i=change(m,2)
           k=change(m,3)
           j=change(m,4)
           if(k>0) then
             celang(:) = cell(:,i)+cell(:,k)
           elseif(k<0) then
             celang(:) = cell(:,i)-cell(:,-k)
           endif
           aplusb=sum(celang*celang)
           aplusb= dabs(aplusb)
           zz = dot_product(cell(:,j),celang(:))
           cell(:,j)=cell(:,j)-change(m,5)*sign(1.0d0,zz)*celang(:)
         endif
         m=m+1
        enddo



     end subroutine

   subroutine  Cal_Vol(C,VOL)
!  CALCULATES THE VOLUME OF THE UNIT CELL
      DOUBLE PRECISION C(3,3), VOL

      VOL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) + &
               ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) +  &
               ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOL = ABS( VOL )
   end subroutine




  subroutine repelforce (rij, rik, f, dfdxi, dfdxj, dfdxk, dfds, costheta0, cutoff)
    implicit none
    double precision, dimension(3) :: rij, rik, rjk, drij, drik, drjk
    double precision, dimension(3) :: dfdr, dcosthetadr, dfdxi, dfdxj, dfdxk
    double precision               :: f, dfds(3,3)
    double precision               :: dij, dik, djk, dij2, dik2, djk2
    double precision               :: costheta, costheta0, cutoff
    double precision               :: fact, ftemp, dftemp, ftemp2, temp
    double precision               :: A1, dA1dr(3), A2, dA2dr(3)
    double precision               :: fcij, fcik, dfcij, dfcik
    double precision               :: S

!   f = 0.0d0; dfdx = 0.0d0; dfds = 0.0d0
    rjk = rik - rij
    dij2 = sum(rij**2); dik2 = sum(rik**2); djk2 = sum(rjk**2)
    dij  = dsqrt(dij2); dik  = dsqrt(dik2); djk  = dsqrt(djk2)

    costheta = (dij2 + dik2 - djk2) / (2.0d0 * dij * dik)
    if(costheta< costheta0) return
    dcosthetadr(1) = (dij2 + djk2 - dik2) / (2.0d0 * dij2*dik)
    dcosthetadr(2) = (dik2 + djk2 - dij2) / (2.0d0 * dik2*dij)
    dcosthetadr(3) = - djk / (dij * dik)

    drij = -rij / dij; drik = -rik /dik; drjk = -rjk /djk
    fact = 1.0d0 - costheta0
    temp = 1.0d0 - costheta

    call hsingle_fc (temp, costheta0, A1, dftemp)
    dA1dr = dftemp * (-dcosthetadr)

    call hsingle_fc (dij, cutoff, fcij, dfcij)
    call hsingle_fc (dik, cutoff, fcik, dfcik)
    
    S = 1d2

    f = A1 * fcij * fcik * S
    dfdr(1) = dA1dr(1) *  fcij * fcik + A1 * dfcij *  fcik * S
    dfdr(2) = dA1dr(2) *  fcij * fcik + A1 *  fcij * dfcik * S
    dfdr(3) = dA1dr(3) *  fcij * fcik * S

    call h2cal_ang_derivative (dfdr(1), dfdxi, dfdxj, dfds, rij, drij)
    call h2cal_ang_derivative (dfdr(2), dfdxi, dfdxk, dfds, rik, drik)
    call h2cal_ang_derivative (dfdr(3), dfdxj, dfdxk, dfds, rjk, drjk)

  end subroutine repelforce

  subroutine hsingle_fc (r, rc, f, dfdr)
    implicit none
    double precision :: r, rc, f, dfdr, temp, temp2

    temp = dtanh(1.0d0 - r/rc); temp2 = temp**2
    f = temp*temp2
    dfdr = -3.0d0*temp2 * (1.0d0 - temp2) / rc

  end subroutine hsingle_fc

  subroutine h2cal_ang_derivative (dfdr, dfdx, dfdy, dfds, x, dx)
  implicit none
  double precision, dimension(3) :: dfdx, dfdy, x, dx
  double precision, dimension(3,3) :: dfds
  double precision               :: dfdr
  integer                        :: i

  dfdx    = dfdx    + dfdr * dx
  dfdy    = dfdy    - dfdr * dx
  dfds(1,1) = dfds(1,1) - dfdr * dx(1) * x(1)
  dfds(2,1) = dfds(2,1) - dfdr * dx(2) * x(1)
  dfds(3,1) = dfds(3,1) - dfdr * dx(3) * x(1)
  dfds(2,2) = dfds(2,2) - dfdr * dx(2) * x(2)
  dfds(3,2) = dfds(3,2) - dfdr * dx(3) * x(2)
  dfds(3,3) = dfds(3,3) - dfdr * dx(3) * x(3)
  return
  end subroutine h2cal_ang_derivative

      subroutine cal_cluster_distance(na,IS0,FS0,AA,distance,distance_peratom)

      use OPTIM_MINDIST

        implicit none
        integer na,i,j,k
        double precision IS0(3,na),FS0(3,na),distance,distance_peratom,com(3),com2(3),AA(3,3)
        double precision celli(3,3)

        double precision, allocatable :: r1(:),r2(:),xac(:,:)

        allocate(r1(3*na),r2(3*na),xac(3,na))

        call reci_latt(AA, celli)

        do i = 1,na
           xac(:,i) = matmul(transpose(celli),IS0(:,i))
           xac(:,i)= modulo(xac(:,i)+100.0d0,1.0d0)
        enddo

        do k = 1,na
          do j = 1,3
            IS0(j,k) = AA(j,1) * xac(1,k) + AA(j,2) * xac(2,k) + AA(j,3) * xac(3,k)
          enddo
        enddo

        do i = 1,na
           xac(:,i) = matmul(transpose(celli),FS0(:,i))
           xac(:,i)= modulo(xac(:,i)+100.0d0,1.0d0)
        enddo

        do k = 1,na
          do j = 1,3
            FS0(j,k) = AA(j,1) * xac(1,k) + AA(j,2) * xac(2,k) + AA(j,3) * xac(3,k)
          enddo
        enddo
        com=0d0
        com2=0d0
        do i=1,na
           com(:)=com(:)+IS0(:,i)
           com2(:)=com2(:)+FS0(:,i)
        enddo
        com=com/dble(na)
        com2=com2/dble(na)
        do i=1,na
        FS0(:,i)=FS0(:,i)-com2(:)+com(:)
        enddo


        r1 = reshape(IS0,(/3*na/))
        r2 = reshape(FS0,(/3*na/))
        call MINDIST(r1,r2,na,distance,.false.,.false.,'None',.false.)

        IS0=reshape(r1,(/3,na/))
        FS0=reshape(r2,(/3,na/))

        do i=1,na
        IS0(:,i)=IS0(:,i)+com(:)
        enddo
        do i=1,na
        FS0(:,i)=FS0(:,i)+com2(:)
        enddo


        distance_peratom=0d0
        do i=1,na
        distance_peratom=max(distance_peratom,dsqrt(sum(IS0(:,i)-FS0(:,i))**2))
        enddo

        deallocate(r1,r2,xac)

      end subroutine




      subroutine check_phiT(na,cell,xa,Energy2,QphiT,Qvac)
      !note: for amorphous structures, return Qvac=-1; for crystal, return Qvac>=0
      !      the range of QphiT: [-1, 1]; for crystal, suggest QphiT: [0.8,1]; for amorphous, suggest QphiT: [:,0.8)
         implicit none
         integer                           :: i,myid,na,vac
         double precision                  :: latt_abc(3),latt_ang(3),cell(3,3), &
                                              QphiT,Qvac, inner_cut, outer_cut,Energy2
         double precision,dimension(3,na)  :: xa
         double precision,dimension(:,:) , allocatable :: coord_down
         double precision,dimension(9)     :: latt_down

          inner_cut=3.1d0
          outer_cut=3.5d0
          allocate(coord_down(3,na))

!test
!        do i = 1,na
!          write(para%ioout,*) "na: ",na
!          write(para%ioout,*) "latt_abc: "
!          write(para%ioout,*) latt_abc
!          write(para%ioout,*) "latt_ang: "
!          write(para%ioout,*) latt_ang
!          write(para%ioout,*) "pos: "
!          write(para%ioout,*) xa(:,1)
!          write(para%ioout,*) "pos: "
!          write(para%ioout,*) xa(:,na)
!        enddo
!end test

      !write(para%ioout,*)"-----init phit in ssw----"
      call upper2down_matrix(na,cell,xa,latt_down,coord_down,Energy2)
      call PHIT(na,latt_down,coord_down,inner_cut,outer_cut,Energy2,QphiT,vac)
      Qvac = float(vac)
      !write(para%ioout,*)"fiT ", QphiT,",  vac  ",Qvac
      !write(para%ioout,*)"-----end phit in ssw----"

      deallocate(coord_down)
      end subroutine



    SUBROUTINE random_pert(na,cellB,cart,xfracFS)

      integer na,j,k,i
      double precision cellB(3,3),cart(3,na),xfracFS(3,na),R(3,3),cellFS(3,3),cell(3,3)
      double precision x,y,z ,celli(3,3)
      double precision angle

       cellFS=cellB
       call reci_latt(cellFS, celli)
       do j=1,na
           xfracFS(:,j) = matmul(transpose(celli),cart(:,j))
        enddo
!     call ca2fr(StrB,xfracFS)
      call random_number(x)
      angle=0.1d0

      call random_number(x)
      call random_number(y)
      call random_number(z)

     x=(2.0d0*x-1.0d0)*angle*3.14159265d0/180d0
     y=(2.0d0*y-1.0d0)*angle*3.14159265d0/180d0
     z=(2.0d0*z-1.0d0)*angle*3.14159265d0/180d0


! a general rotation matrix
     R(1,1)=cos(y)*cos(z)
     R(2,1)=cos(y)*sin(z)
     R(3,1)=-sin(y)

     R(1,2)=-cos(x)*sin(z)+sin(x)*sin(y)*cos(z)
     R(2,2)=cos(x)*cos(z)+sin(x)*sin(y)*sin(z)
     R(3,2)=sin(x)*cos(y)

     R(1,3)=sin(x)*sin(z)+cos(x)*sin(y)*cos(z)
     R(2,3)=-sin(x)*cos(z)+cos(x)*sin(y)*sin(z)
     R(3,3)=cos(x)*cos(y)
! end rotation matrix


     cellB=matmul(R,cellFS)
     cart=matmul(cellB,xfracFS)


    END SUBROUTINE



 end module
