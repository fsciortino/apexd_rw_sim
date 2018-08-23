!C    =================================================================
      module mod_variable
!C    In this version, electric field files from Simion can be used.
!C    =================================================================
      implicit none

      double precision, parameter :: pi = 3.14159265358979323846d0;
      double precision, parameter :: xc = 2.99792458d8
      double precision, parameter :: xe = 1.6021766208d-19;
      double precision, parameter :: xm = 0.910938356d-30;
      double precision, parameter :: eom = 1.7588200236d11;
      
      double precision :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz
      integer :: nx,ny,nz,n0x,n0y,n0z,iexb,iexe,iele,ielemax,ielep,icln
      double precision,dimension(1:205) :: mt1x,mt1y,mt1z
      double precision,dimension(1:1,1:1,1:1) :: Bxm,Bym,Bzm ! to save memory for those usually not to be used 
      double precision,dimension(1:105,1:105,1:105,1:4) :: Exm,Eym,Ezm ! electric fields (unit field by electrodes) max9 
      double precision,dimension(1:12) :: Vel,Vdc,Vac ! electrode voltage
      double precision :: Fel,Phrw ! frequency and phase of RW 
      double precision :: vfp1,vfp2,vfp3,vfp4,vfp5,vfp6 
      double precision,dimension(1:52) :: ci,rc,zcl,xcl,ycl ! Coil current, loop radius, z, x, y(y is usually 0) 
      double precision :: fgrad
      end module mod_variable


!C    =================================================================
      subroutine bbr_particle(pmx,pmy,pmz,pmk1,pmk2,pma,&
           pmrwa,pmrwf,pmrwfg,pmrwd, &
           fpm1,fpm2,fpm3,fpm4,fpm5,fpm6, & 
           ipmt,dt,iend,imbk,imbk2,tE1,ci1,rc1,pid,Evers, &
           stt,stx,sty,stz,str,stgx,stgy,stgz,stgr,&
           stK,stKpr,stKpp,stps1,stps2,stmu,strt) !,vtor,r2d,ftor)
!    
!     This version was modified by FS in order to compile into a Python shared library.
!
!     Compile using a signature file:
!        f2py bbr_particle.f90 -m bbr_particle -h bbr_particle.pyf only: bbr_particle --overwrite-signature
!        --- modify .pyf file here ---
!        f2py -c bbr_particle.pyf bbr_particle.f90
!
!     Use this with particle_run.py and particle_stats.py
!        
!C    =================================================================
      use mod_variable
      implicit none
  
      double precision,dimension(1:3) :: x,v
      double precision :: ken, ken1, ken2
      double precision :: Tr
      double precision, dimension(1:3) :: u,uminus,uzero,uplus
      double precision,dimension(1:3) :: B,E,Tv,Sv
      double precision :: tini,Tsq,ddt,xu,xv, kini,xdr,ydr,zdr,kini1,kini2,xpr,ypr,zpr,rsa,rPsi,rPsipr
      double precision :: dt
!f2py intent(in) :: dt
      double precision :: xv1,xv2,xv1int,kang,bmb1,bmb2
      double precision :: xi,yi,zi,Ex,Ey,Ez,Bx,By,Bz
      integer :: imbk, imbk2
!f2py intent(in) :: imbk, imbk2
      integer :: i,iini,Nr,itf,ifile,ifile2,imbk2t,iclmap,iclbp,iclbp2,imu,ntr,ntrpr
      integer :: its,ikene,ipm,ipm1,ipm2  
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      double precision :: Hr,xgm,xbt,rnd, Btt,Psi,Psi1,Psi2,xg1,yg1,zg1,xgpr1,ygpr1,zgpr1
      double precision :: gr,xmu
      integer(8) :: iend
!f2py intent(in) :: iend
      integer(8) :: it,iend2
      integer t1, t2, t_rate, t_max, diff,iex1,iex2,iex3,iex4,icmp,ihit,ieout,ieoutc,iden,iif,irsa,irsapr
      double precision,dimension(1:100005) :: elik2,elie2
      double precision,dimension(1:1001) :: pmx,pmy,pmz,pmk1,pmk2,pma
!f2py intent(in) :: pmx,pmy,pmz,pmk1,pmk2,pma
      double precision,dimension(1:1001) :: pmk,pmrx,pmry,pmrz
      double precision,dimension(1:1001) :: pmrwa,pmrwf,pmrwfg,pmrwd
!f2py intent(in) :: pmrwa,pmrwf,pmrwfg,pmrwd      
      double precision,dimension(1:1001) :: fpm1,fpm2,fpm3,fpm4,fpm5,fpm6 
!f2py intent(in) :: fpm1,fpm2,fpm3,fpm4,fpm5,fpm6
      character filename*128,filename0*128,filenameJ*128,filenamemu*128, filetor*128
      integer dtm(8)
      character(len=10) sys_time(3)
      integer,parameter :: max_line_len1 = 100
      character(max_line_len1) cf2, Evrs
      double precision,dimension(1:105,1:10005) :: stt,stx,sty,stz,str,stgx,stgy,stgz,stgr
!f2py intent(out) :: stt,stx,sty,stz,str,stgx,stgy,stgz,stgr
      double precision, dimension(1:105,1:10005) :: stK,stKpr,stKpp,stps1,stps2,stmu,strt
!f2py intent(out) :: stK,stKpr,stKpp,stps1,stps2,stmu,strt
      !dimension(1:105,1:10005) is time step & particle number
      double precision,dimension(1:105) :: stta,stxa,stya,stza,stra,stgxa,stgya,stgza,stgra,stKa,stKpra,stKppa,stps1a,stps2a,stmua
      double precision,dimension(1:105) :: stts,stxs,stys,stzs,strs,stgxs,stgys,stgzs,stgrs,stKs,stKprs,stKpps,stps1s,stps2s,stmus
      integer, dimension(1:105) :: sttn ! averaged value and sigma (for particles)
      double precision,dimension(1:10005) :: tend 
      integer :: ist,isti
      double precision :: Vtpp1,Vtpp2,Vtpp3,Vtppn
      integer :: ipmt
!f2py intent(in) :: ipmt
      integer :: pid
!f2py intent(in) :: pid
      double precision :: tE1
!f2py intent(in) :: tE1
      double precision :: vtor,r2d,ftor
!!!!!!!!!f2py intent(out) :: vtor,r2d,ftor
      integer :: noco  !no.coils
!!!!!f2py intent(in) :: noco
      double precision :: ci1, rc1
!f2py intent(in) :: ci1,rc1      
      integer :: mucnt ! count for gyration to compute \mu
      integer :: Evers  ! version of electric field to be loaded
!f2py intent(in) :: Evers

      
      call system_clock(t1)   ! record starting time 
      
!      ipm1=1; ipm2=1; ! starting and ending line of param.txt
      !      ipmt = ipm2-ipm1+1 ! total particle number, maximum value of ipm
      !write(6,*) dt,iend, imbk, imbk2
      ipm1 = 1
      ipm2 = ipmt
      
!C    ##### time step, step number #### in BB, variable timestep doesn't current work, so its should be 0
!C    in BB, this function is omitted.
!     its=0 ! 0:fixed time step, 1:variable time step
!C    ### for 0 ###
      tini = 0.0d0
!     dt = 1.0d-12
!     iend = 200000
!     imbk = 1000 ! data is saved to a file 1 times per imbk 
!     imbk2= 10000 ! same as above, but used for statistical calculation for many particles 
      imbk2t = iend/imbk2+1
!      dt = 1.0d-11; iend = 10000000; imbk = 1000; imbk2= 10000;
!!C    ### for 1 ###
!      tsdiv = 560.0d0 ! cyclotron time is divided with this number
!      tmax = 1.0d-6 ! end of calculation
!      iend2 = 10000000000_8 ! just a large number. calculation stops when tmax<t or iend2<it

!     ##### Settings #####
      ikene = 2 ! 1: kene and direction ratio, 2: kpara and kperp
      iif=1; ! 1: save in ./orb/, 0: same folder
      ifile=1; ! 1: save orbit files for each of the particles
      icmp=1; ! 1: save additional files
      ifile2=0; ! 1: save particle statistics in a file, 0: don't compute particle statistics here
      ihit=1; ! 1: stop calculation when hitting outer electrodes

!C    ##### B coil current, radius, z, x, y (y is usually 0)
      !ci1=1.0d3; rc1=1.0d-1
      ci(1) = ci1; rc(1) = rc1; zcl(1) = 0.0d0; xcl(1) = 0.0d0; ycl(1) = 0.0d0;
      ci(2) = 1.0d3; rc(2) = 1.0d-1; zcl(2) = 0.01d0; xcl(2) = 0.0d0; ycl(2) = 0.0d0;
      ci(3) = 1.0d3; rc(3) = 1.0d-1; zcl(3) = 0.02d0; xcl(3) = 0.0d0; ycl(3) = 0.0d0;
      ci(4) = 1.0d3; rc(4) = 1.0d-1; zcl(4) = 0.03d0; xcl(4) = 0.0d0; ycl(4) = 0.0d0;
      ci(5) = 1.0d3; rc(5) = 1.0d-1; zcl(5) = 0.04d0; xcl(5) = 0.0d0; ycl(5) = 0.0d0;
      ci(6) = 1.0d3; rc(6) = 1.0d-1; zcl(6) = -0.01d0; xcl(6) = 0.0d0; ycl(6) = 0.0d0;
      ci(7) = 1.0d3; rc(7) = 1.0d-1; zcl(7) = -0.02d0; xcl(7) = 0.0d0; ycl(7) = 0.0d0;
      ci(8) = 1.0d3; rc(8) = 1.0d-1; zcl(8) = -0.03d0; xcl(8) = 0.0d0; ycl(8) = 0.0d0;
      ci(9) = 1.0d3; rc(9) = 1.0d-1; zcl(9) = -0.04d0; xcl(9) = 0.0d0; ycl(9) = 0.0d0;
      ci(10)= 1.0d3; rc(10) = 1.0d-1; zcl(10) = -0.5d0; xcl(10) = 0.0d0; ycl(10) = 0.0d0;
      ! temporary:
      noco=1;
      icln = noco; !1; ! coil number 
      iclmap = 0; ! 1:save B data 0:no
      
!C    ##### B is given by subroutine calculation or linear interpolation of external file
      iexb=0; ! 0: Biot-Savart 1: use external B field file
      
!C    ##### E is given by subroutine calculation or linear interpolation of external file
      iexe=1; ! 0: azimuthal-only field  1: use external E field file
      !if (Evrs .lt. 4) then
         ielemax = 4; ! number of electrodes
      !else
      !   ielemax = 8;
      !endif
      
!C    voltages and frequencies of electrodes. 
      ielep=1; ! when 1, param.txt is used for RW parameters in addition to the followings. 
      Vdc(1) = 0.0d0; Vdc(2) = 0.0d0; Vdc(3) = 0.0d0; Vdc(4) = 0.0d0; Vdc(5) = 0.0d0;
      Vdc(6) = 0.0d0; Vdc(7) = 0.0d0; Vdc(8) = 0.0d0; Vdc(9) = 0.0d0; Vdc(10) = 0.0d0;
      Vac(1) = 5.0d0; Vac(2) = 5.0d0; Vac(3) = 5.0d0; Vac(4) = 5.0d0; Vac(5) = 0.0d0;
      Vac(6) = 0.0d0; Vac(7) = 0.0d0; Vac(8) = 0.0d0; Vac(9) = 0.0d0; Vac(10) = 0.0d0;
      Fel = 1.0d5;
      fgrad = 5.0d9;
      
      
!C    ##### warning/stop calculation when the particle is outside of the E file region 
      ieout=2; ! 1:warning, 2: stop caculation when particle is outside E
      
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), dtm)

!      if (iif .eq. 1) then
!      write (filename0, '("./orb/log",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"bbr.txt")') dtm(1),dtm(2),dtm(3),dtm(5),dtm(6),dtm(7)
!      else
!      write (filename0, '("log",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"bbr.txt")') dtm(1),dtm(2),dtm(3),dtm(5),dtm(6),dtm(7)
!      endif
!      write(6,*) filename0
!      open (21,file=filename0)
!
!      write(21,'(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)') dtm(1),dtm(2),dtm(3),dtm(5),dtm(6),dtm(7) ! starting date and time 
!      write(21,'(i1.1)') its 
!      write(21,'(1p1d11.3,i15,1p1d11.3,i10)') dt,iend,dt*iend,imbk

!      if (its == 0) then
!      print "(a35,1p1d11.3,i15,1p1d11.3)","## fixed time step ## dt/iend/T(s):",dt,iend,dt*iend
!      endif
!      if (its == 1) then
!      print "(a39,1p2d11.3)","## variable time step ## tsdiv tmax(s):",tsdiv,tmax
!      print "(a65,1p2d11.3)","-- time step is (cyclotron period)/tsdiv at aech point --"
!      endif

      if (iexb == 1) then !C   B files to be used with interpolation. not recommended. 
         open(17, file='b3d_hd.txt')
         read(17,*) xmin,xmax;! write(6,*)  xmin,xmax
         read(17,*) ymin,ymax;! write(6,*) ymin,ymax
         read(17,*) zmin,zmax;! write(6,*) zmin,zmax
         read(17,*) nx,ny,nz;! write(6,*) nx,ny,nz
         read(17,*) dx,dy,dz;! write(6,*) dx,dy,dz
         close(17)

         do i=1,nx+1
            mt1x(i)=xmin+(i-1)*dx;! write(6,*) mt1x(i)
         end do
         do i=1,ny+1
            mt1y(i)=ymin+(i-1)*dy;! write(6,*) mt1y(i)
         end do
         do i=1,nz+1
            mt1z(i)=zmin+(i-1)*dz;! write(6,*) mt1z(i)
         end do

         open(17, file='b3d.txt')
         do iex3=1,nz+1
            do iex2=1,ny+1
               do iex1=1,nx+1
                  read (17,*) Bxm(iex1,iex2,iex3),Bym(iex1,iex2,iex3),Bzm(iex1,iex2,iex3);
               end do
            end do
         end do
         close(17)

      endif

      if (iexe == 1) then !C    E files to be used with interpolation
         write (Evrs, '("./simion_Efields/ver.", i1.1,"/e3d_hd.txt")') Evers
         open(17, file=Evrs)
         read(17,*) xmin,xmax;! write(6,*)  xmin,xmax
         read(17,*) ymin,ymax;! write(6,*) ymin,ymax
         read(17,*) zmin,zmax;! write(6,*) zmin,zmax
         read(17,*) nx,ny,nz;! write(6,*) nx,ny,nz
         read(17,*) dx,dy,dz;! write(6,*) dx,dy,dz
         close(17)

         ! write name of E-field header file
         write(6,*) "Read E-field file header: ",trim(Evrs)

         do i=1,nx+1
            mt1x(i)=xmin+(i-1)*dx;! write(6,*) mt1x(i)
         end do
         do i=1,ny+1
            mt1y(i)=ymin+(i-1)*dy;! write(6,*) mt1y(i)
         end do
         do i=1,nz+1
            mt1z(i)=zmin+(i-1)*dz;! write(6,*) mt1z(i)
         end do

         do iex4=1,ielemax
            write (cf2, '("./simion_Efields/ver.", i1.1,"/e3d", i1.1, ".txt")') Evers,iex4
            ! write E-field file that was read:
            write(6,*) "read E field file: ",trim(cf2)
            open(17, file=cf2) !C electrode
            do iex3=1,nz+1
               do iex2=1,ny+1
                  do iex1=1,nx+1
                     read (17,*) Exm(iex1,iex2,iex3,iex4),Eym(iex1,iex2,iex3,iex4),Ezm(iex1,iex2,iex3,iex4);
                  end do
               end do
            end do
            close(17)
         end do

      endif

!C    assuming x^0, v^0, v^-1/2 is caucluated by rk4

!	  ################################################################################ start
!	  ################################################################################ initial processes
!     ####### set for complete eliptic functions #######
!C    Should be 1 only for the first calculation.
      itf=0 ! 0:load and save, or 1:calculate eliptic functions.
!     ####### complete elliptic functions. usually 0. file is newly made when 1. #######

!     ##### complete eliptic functions, start ####
      if (itf.eq.1) then
!C   call efile2 only for the first time
!CC     K(k) and E(k) are calculated and saved 
      call efile2(elik2,elie2)
      OPEN (2,File="el2.txt")
      DO 101 i=1,99999
      write(2,'(1p1d12.6,1p2d25.17)') I/100000.0d0,elik2(i),elie2(i)
 101  CONTINUE
      CLOSE(2)
      endif

      if (itf.eq.0) then
!C    Load saved file (this is faster) 
!      write(6,*) "Loading complete elliptic integrals..."
      OPEN (50,File="el2.txt")
      DO 102 I=1,99999
      READ(50,*) XI,ELIK2(I),ELIE2(I)
 102  CONTINUE
      CLOSE(50)
      endif

!      The following test ensures that integrals are calculated correctly:
!      write(6,*) "Done. K(k=1) should be 1  1.5707963267948"
!      write(6,*) 1,elik2(1)

      if(ELIK2(1)-1.57079d0.gt.0.01d0) then ! Warning when E() is empty.
      write(6,*) "##### !!Check eliptic functions!! Data is null!!! #####"
      write(6,*) ""
      else
!      write(6,*) "OK. Checked that elik2(1) contains non-zero data."
!      write(6,*) ""
      endif

!    ##### complete eliptic functions, end ####

      if (iclmap == 1) then ! B calculation files are made 

!     When hd file(s) are not read, xmax, xmin, etc. are not decided. 
!     Set the calculation region of magnetic field data 
      if (xmin .gt. -1.0d-6 .and. xmin .lt. 1.0d-6) then
      if (xmax .gt. -1.0d-6 .and. xmax .lt. 1.0d-6) then
      xmin=-0.25d0; xmax=0.25d0;
      ymin=-0.25d0; ymax=0.25d0;
      zmin=-0.25d0; zmax=0.25d0;
      endif
      endif

      open (50,File="phib.txt") ! B at z=0 
            yi=0.0d0; zi=0.0d0;
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 
                  xi=(xmax/100.0d0)*iclbp;
                  call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz)
                  call bmap(elik2,elie2,xi,yi,zi,Psi)
!		write(6,*)  iclbp,xmax,xmax/10.0d0
                  write(50,'(1p3d18.10)') xi, Psi, dsqrt ( Bx*Bx+By*By+Bz*Bz )
                  end do
      close(50)

      open (50,File="ck1.txt") ! analytic solution of B at z=0, to check calculation. z=0
!	  splot [][][0:0.01] "b.txt" using 1:2:4 w l,"ck1.txt" using 1:2:3 w l   
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 
            xi=0.0d0; yi=0.0d0;
            zi=zmin+( (zmax-zmin)/100.0d0 )*iclbp; 
!		write(6,*)  iclbp,xmax,xmax/10.0d0
            write(50,'(1p3d18.10)') xi, zi, 2.0d-7*pi*ci(1)*rc(1)*rc(1) / ((rc(1)*rc(1)+zi*zi)**1.5)
            end do
      close(50)

      open (50,File="b.txt")
            yi=0.0d0
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below!
            xi=xmin+( (xmax-xmin)/100.0d0 )*iclbp;
            do iclbp2 = 0,100
            zi=zmin+( (zmax-zmin)/100.0d0 )*iclbp2;
            call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz); call bmap(elik2,elie2,xi,yi,zi,Psi);
!		write(6,*)  iclbp,xmax,xmax/10.0d0
            write(50,'(1p7d18.10)') xi,zi,Psi,dsqrt ( Bx*Bx+By*By+Bz*Bz ),Bx,By,Bz
            end do
            write(50,'(i1)')
            end do
      close(50)

      open (50,File="bi.txt")
            yi=0.0d0
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below!
            xi=xmin+( (xmax-xmin)/100.0d0 )*iclbp;
            do iclbp2 = 0,100
            zi=zmin+( (zmax-zmin)/100.0d0 )*iclbp2;
            call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz); call bmap(elik2,elie2,xi,yi,zi,Psi);
!		write(6,*)  iclbp,xmax,xmax/10.0d0
            write(50,'(1p7d18.10)') xi,zi,Psi,dsqrt ( Bx*Bx+By*By+Bz*Bz ),Bx,By,Bz
            end do
!		write(50,'(i1)')
            end do
      close(50)

      open (50,File="ck2.txt") ! BS should be close to 2 pi r Atheta near axis
!	  plot "ck2.txt" using 1:2,"ck2.txt" using 1:3 w l
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 
            yi=0.0d0; zi=0.0d0;
            xi=( (xmax-xmin)/100.0d0 )*iclbp/10.0d0; 
            call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz); call bmap(elik2,elie2,xi,yi,zi,Psi);
            write(50,'(1p3d18.10)') xi, Bz*pi*xi*xi,Psi
            end do
      close(50)
      
      endif
!	  ################################################################################ end
!	  ################################################################################ initial processes

      do ipm=ipm1,ipm2 ! loop for each particle
      ieoutc=0 ! a constant to be used for warning outside E file region 
      iden=1; ! time label for density information
      ist=1; ! time label for statistical calculation 
      
      if (ikene == 2) then !C    parallel and perpendicular K(eV) ------
      x(1)=pmx(ipm); x(2)=pmy(ipm); x(3)=pmz(ipm);
      kini1=pmk1(ipm); kini2=pmk2(ipm); kang=pma(ipm);
      endif
      
      if (ikene == 1) then !C   total K(eV) and velocity direction K(eV)
      x(1)=pmx(ipm); x(2)=pmy(ipm); x(3)=pmz(ipm);
      kini=pmk(ipm); xdr=pmrx(ipm); ydr=pmry(ipm); zdr=pmrz(ipm);
      endif
      
      if (ielep == 1) then !     RW amplitude and frequency 
      Vac(1)=pmrwa(ipm);Vac(2)=pmrwa(ipm);Vac(3)=pmrwa(ipm);Vac(4)=pmrwa(ipm);
      Fel=pmrwf(ipm);
      fgrad=pmrwfg(ipm);
      Phrw=pmrwd(ipm);
      endif

      vfp1=fpm1(ipm);vfp2=fpm2(ipm);vfp3=fpm3(ipm);vfp4=fpm4(ipm);vfp5=fpm5(ipm);vfp6=fpm6(ipm); !

      print "(a1)"," "
      print "(a32,i10,i8,i8)","---------- Particle No. out of:",ipm,ipm2-ipm1+1
!      print "(a35,1p3d12.4)","position (x,y,z) (m):   ",x(1),x(2),x(3)
!      print "(a35,1p3d12.4)","RW V, freq, dir (V) (Hz) ():   ",Vac(1),Fel,Phrw

!C    ##### initial values of E and B at x^0 ####

      if (iexb == 0) then !C    calculate B by Biot-Savart
      call bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz);
      B(1)=Bx; B(2)=By; B(3)=Bz;
!      write(6,*) "initial B with usual calculation:"
!      print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with usual calculation"
      endif

      if (iexb == 1) then !C    B is given by interpolation of external file data
      call bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz);
      print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with usual calculation"

      call bvac2(x(1),x(2),x(3),Bx,By,Bz);
      B(1)=Bx; B(2)=By; B(3)=Bz;
      print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with interpolation"
      endif
      
      if (iexe == 0) then !C    E from subroutine calculation
      call evac0(x(1),x(2),x(3),tini,Ex,Ey,Ez);
      E(1)=Ex; E(2)=Ey; E(3)=Ez;
      print "(1p3d13.5,a35)",E(1),E(2),E(3),": E (V/m) with usual calculation"
      endif

      if (iexe == 1) then !C    E is given by interpolation of external file data 
      call evac2(x(1),x(2),x(3),tini,tE1,Ex,Ey,Ez);
      E(1)=Ex; E(2)=Ey; E(3)=Ez;
!      print "(1p3d13.5,a35)",E(1),E(2),E(3),": E (V/m) with interpolation"
      endif

!C    ##### initial velocity ####

      if (ikene == 1) then
!C    total K(eV) and directions K(eV)

      ! xgm is the relativistic \gamma; xbt=\beta; xc=speed of light
         
      xgm=1.0d0+kini*1.95695119784913d-6 ! XGM=1.0D0+KENE*XE/(XME*XC*XC)
      xbt=dsqrt(1.0d0-1.0d0/(xgm*xgm)); xv=xc*xbt; xu=xgm*xv;
      u(1) = xu*xdr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))
      u(2) = xu*ydr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))
      u(3) = xu*zdr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))

      print "(a35,1p3d12.4)","four-velocity (ux,uy,uz) (m):   ",u(1),u(2),u(3)
      print "(a35,1p3d12.4)","real velocity (vx,vy,vz) (m):   ",u(1)/xgm,u(2)/xgm,u(3)/xgm
      print "(a35,1p2d18.10)","beta and gamma:   ",xbt,xgm
      print "(a35,1p2d18.10)","kinetic energy (eV):   ",kini
      
      endif

      if (ikene == 2) then
!C    parallel and perpendicular K(eV) 

      kini=kini1+kini2
      xgm=1.0d0+kini*1.95695119784913d-6 ! XGM=1.0D0+KENE*XE/(XME*XC*XC)
      xbt=dsqrt(1.0d0-1.0d0/(xgm*xgm)); xv=xc*xbt; xu=xgm*xv;
      
      xv1 = xv*dsqrt(kini1/kini); ! vpara
      xv2 = xv*dsqrt(kini2/kini); ! vperp
      print "(a35,1p3d13.5)","kinetic energy (eV), pr,pp,ttl:   ",kini1, kini2, kini
!      print "(a35,1p3d13.5)","velocity (m/s):   ",xv1, xv2, xv
      
      bmb1 = dsqrt( B(1)*B(1) + B(2)*B(2) + B(3)*B(3) ); 
      bmb2 = dsqrt( B(1)*B(1) + B(3)*B(3) )

      v(1) = xv1 * B(1) / bmb1 + xv2 * ( B(3) / bmb2 * dcos(kang) - B(1)*B(2) / (bmb1*bmb2) * dsin(kang) )
      v(2) = xv1 * B(2) / bmb1 + xv2 * ( 0.0d0 + ( B(3)*B(3)+B(1)*B(1) ) / (bmb1*bmb2) * dsin(kang) )
      v(3) = xv1 * B(3) / bmb1 + xv2 * (-B(1) / bmb2 * dcos(kang) - B(2)*B(3) / (bmb1*bmb2) * dsin(kang) )

!C    vpara component, vperp-1 (cos) component, vperp-2 (sin) component
!      write(6,*) xv1 * B(1) / bmb1, xv2 * ( B(3) / bmb2 * dcos(kang) ), - xv2 * B(1)*B(2) / (bmb1*bmb2) * dsin(kang) 
!      write(6,*) xv1 * B(2) / bmb1, 0.0d0, xv2 * (  (B(3)*B(3)+B(1)*B(1)) / (bmb1*bmb2) * dsin(kang) )
!      write(6,*) xv1 * B(3) / bmb1, xv2 * (-B(1) / bmb2 * dcos(kang)), -xv2 *  B(2)*B(3) / (bmb1*bmb2) * dsin(kang) 

      u(1) = v(1)*xgm; u(2) = v(2)*xgm; u(3) = v(3)*xgm;

      print "(a35,1p3d13.5)","four-velocity (ux,uy,uz) (m):   ",u(1),u(2),u(3)
      print "(a35,1p3d13.5)","real velocity (vx,vy,vz) (m):   ",v(1),v(2),v(3)
      print "(a35,1p2d18.10)","beta and gamma:   ",xbt,xgm

!C    Larmor radius mv/qB, v=1.75882002359937d6,B=0.01T -> 1mm
!      write(6,*) xm*v(2)/(xe*B(3))
      print "(a35,1p2d18.10)","gyroradius (m):",xv2/bmb1/eom
      print "(a35,1p2d18.10)","1/f_cyc (s):",2.0d0*pi/bmb1/eom
!      gr=xm*u(2)/(xe*B(3));
      
      endif


!C    time step in RK to calculate x^0, v^-1/2 
!     if(its == 0) then ! constant time step
      Hr = dt;
!     endif

!     if(its == 1) then ! variable time step
!     dt2 = ( 2.0d0*pi/(( dsqrt(Bx*Bx+By*By+Bz*Bz) )*eom) ) / tsdiv
!     Hr = dt2;
!     endif

!      print "(1p3d10.2,a3,1p3d10.2)",x(1),x(2),x(3)," ",u(1),u(2),u(3)

!C    ##### 1/2 back by rk4 to get v^-1/2 rk4, start ####
      iini=10000; ddt = Hr/iini;
      Nr = 6;
      Hr = ddt;
      xr(1)=x(1); xr(2)=x(2); xr(3)=x(3);
      xr(4)=u(1); xr(5)=u(2); xr(6)=u(3); ! four-vector 
!      print "(1p1d10.2,a2,1p3d10.2,a3,1p3d10.2)",tini," ",xr(1),xr(2),xr(3)," ",xr(4),xr(5),xr(6)
      open (2,file="backr.txt")
        do it = 0,iini/2-1 ! back to ^-1/2 
        Tr = tini + it*ddt      
        call rk4(elik2,elie2,Nr,Tr,Xr,K1,K2,K3,K4,YWORK,Hr,tE1)
        write(2,'(1p2d11.3,1p6d15.7)') tini-(it+1)*ddt,xr(1),xr(2),xr(3),xr(4),xr(5),xr(6)
      end do
      close(2)
      u(1)=xr(4); u(2)=xr(5); u(3)=xr(6); ! u^-1/2
!C    ##### 1/2 back by rk4 to get v^-1/2 rk4, end ####

!     x(1),x(2),x(3): x at t=0
!     u(1),u(2),u(3): u at t=-1/2

!C    'XU (four velocity), XV (real velocity)'
      xu = dsqrt( u(1)*u(1)+u(2)*u(2)+u(3)*u(3) )
      xv = xc*xu/dsqrt(xc*xc+xu*xu)
      xbt = xv/xc; xgm = 1.0d0/dsqrt(1.0d0-xbt*xbt) ! this is v^n-1/2

      v(1)=u(1)/xgm; v(2)=u(2)/xgm; v(3)=u(3)/xgm; 

!      x^n->ok for electromagnetic fields
!      v: v^n-1/2

!      if (its == 1) then ! for variable time step
!      iend = iend2
!      endif

      Hr = dt
      Tr = tini

!      write(6,*) " --- Process ID: ", pid
      if (ifile == 1) then ! orbit output files 
      
        if (iif .eq. 1) then
           write (filename, '("./orb/orb", i4.4, "_", i5.5, "bbr.txt")') ipm,pid 
           write (filenameJ, '("./orb/jp", i4.4,"_",i5.5, "bbr.txt")') ipm,pid 
           write (filenamemu, '("./orb/mu", i4.4,"_",i5.5, "bbr.txt")') ipm,pid 
        else
           write (filename, '("orb", i4.4,"_",i5.5, "bbr.txt")') ipm,pid 
           write (filenameJ, '("jp", i4.4,"_",i5.5, "bbr.txt")') ipm,pid 
           write (filenamemu, '("mu", i4.4,"_",i5.5, "bbr.txt")') ipm,pid  
        endif

        open (20,file=filename)
        open (22,file=filenameJ)
        open (24,file=filenamemu)
      
      endif
      ntr=0; ! counter for toroidal rotation 
      ntrpr=0;
      imu=0; ! saved number of mu
      xgpr1=0.0d0; ygpr1=0.0d0; zgpr1=0.0d0; xpr=0.0d0; ypr=0.0d0; zpr=0.0d0;

      ! Assign the following=0 here to avoid compiler-specific floating
      irsapr=0; rPsipr=0.0d0 ! "previous" irsa index and Psi values. 
      mucnt=0  ! count for \mu computation
      
!     ######################################################################################
!     save first, calculation of beta etc. is next. 
!C    ##### B-B method pusher, start ####
      do it = 0,iend+1
         
!C     time step for t^n,v^n-1/2 
!      if(its == 0) then ! constant time step
      Hr = dt
      Tr = dt*it;

!      gr1 = dsqrt( (x(1)-gr)*(x(1)-gr)+x(2)*x(2) ) ! use x^n x^n 

      if(ihit == 1) then ! stop calculation when hitting electrodes 
      
        if( dsqrt( x(1)*x(1)+x(2)*x(2) ) .gt. 0.3d0 ) then
        write(6,*) "positron hit the wall at r=0.3m"
        go to 10
        endif

        if( dabs( x(3) ) .gt. 0.3d0 ) then
        write(6,*) "positron hit the wall at z=\pm0.2m"
        go to 10
        endif

      endif

      if(iexe == 1) then ! only when E files are used 
      if(ieout == 1) then ! warning when particle is outside of the E file region 
      if(ieoutc == 0) then ! only once 

        if( x(1) .lt. xmin ) then
        write(6,*) "####### x < xmin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( x(1) .gt. xmax ) then
        write(6,*) "####### x > xmax of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        
        if( x(2) .lt. ymin ) then
        write(6,*) "####### y < ymin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( x(2) .gt. ymax ) then
        write(6,*) "####### y > ymax of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        
        if( x(3) .lt. ymin ) then
        write(6,*) "####### z < zmin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( x(3) .gt. ymax ) then
        write(6,*) "####### z > zmax of E file, outside of E field region, continue #######"; ieoutc=1
        endif

      endif
      endif

      if(ieout == 2) then ! stop when outside E region

        if( x(1) .lt. xmin ) then
           write(6,*) "####### outside of E field region, stop "
           go to 10;
        endif
        if( x(1) .gt. xmax ) then
        write(6,*) "####### x > xmax of E file, outside of E field region, stop #######"; go to 10;
        endif
        
        if( x(2) .lt. ymin ) then
        write(6,*) "####### y < ymin of E file, outside of E field region, stop #######"; go to 10;
        endif
        if( x(2) .gt. ymax ) then
        write(6,*) "####### y > ymax of E file, outside of E field region, stop #######"; go to 10;
        endif
        
        if( x(3) .lt. ymin ) then
        write(6,*) "####### z < zmin of E file, outside of E field region, stop #######"; go to 10;
        endif
        if( x(3) .gt. ymax ) then
        write(6,*) "####### z > zmax of E file, outside of E field region, stop #######"; go to 10;
        endif
        
      endif
      endif
   
!     Because v^n+1/2 was obtained at the end of loop, this is v^n-1/2
!     ############################################################
!     ## save data in a file, start ##
      if (ifile == 1) then
      
         if (icmp == 1) then
            ! get Kpara Kperp 
      
            ken=(xgm-1.0d0)/1.95695119784913d-6;
            xv1 = Bx*v(1)+By*v(2)+Bz*v(3) / dsqrt(Bx*Bx+By*By+Bz*Bz)
            xv2 = dsqrt( xv*xv - xv1*xv1 )
            ken1 = ken * xv1*xv1 / (xv*xv)
            ken2 = ken * xv2*xv2 / (xv*xv)
            Btt = dsqrt(Bx*Bx+By*By+Bz*Bz)
      
            call bmap(elik2,elie2,x(1),x(2),x(3),Psi); Psi1=Psi; ! \Psi 

!           guiding center calculation with real velocity 
            Vtpp1 = ( u(1)-Bx/Btt*(u(1)*Bx+u(2)*By+u(3)*Bz)/Btt )/xgm; ! get Vperp 
            Vtpp2 = ( u(2)-By/Btt*(u(1)*Bx+u(2)*By+u(3)*Bz)/Btt )/xgm;
            Vtpp3 = ( u(3)-Bz/Btt*(u(1)*Bx+u(2)*By+u(3)*Bz)/Btt )/xgm;

            Vtppn = dsqrt (Vtpp1*Vtpp1 + Vtpp2*Vtpp2 + Vtpp3*Vtpp3) ! | Vperp |
      
            gr = Vtppn*xgm/Btt/eom;
!C          XLML=XME*VPP*XGM/XE/(sqrt(BX*BX+BY*BY+BZ*BZ))

!C          guiding center position, not accurate for large rL 
            xg1 = x(1) + (Vtpp2*Bz-Vtpp3*By)/Vtppn/Btt * gr
            yg1 = x(2) + (Vtpp3*Bx-Vtpp1*Bz)/Vtppn/Btt * gr
            zg1 = x(3) + (Vtpp1*By-Vtpp2*Bx)/Vtppn/Btt * gr

!           nigr(iden,ipm) = dsqrt( nigx(iden,ipm)*nigx(iden,ipm)+nigy(iden,ipm)*nigy(iden,ipm) )

            call bmap(elik2,elie2,xg1,yg1,zg1,Psi); Psi2=Psi; ! \Psi at guiding center 
      
!C          ########## mu ########## start
            rPsi=Psi1
!           rsa = dsqrt( xr(1)*xr(1)+xr(2)*xr(2) ) - dsqrt( xpr*xpr+ypr*ypr ); ! R increase -> 1, decrease -> -1 
            rsa = rPsi - rPsipr; ! Psi increase -> 1, decrease -> -1 Psi
            if (rsa .gt. 0.0 ) then
               irsa=1;  
            endif
            if ( rsa .lt. 0.0d0 ) then
               irsa=-1; 
            endif

!           XMU=XMU+VPP*VPP/(2.0*sqrt(BX*BX+BY*BY+BZ*BZ))*XGM*XGM
            xmu = xmu + xv2*xv2/(2.0d0*Btt)*xgm*xgm  !FS: note added factor of 1/2
            ! note that this should be multiplied by *mass* (small number)
            mucnt= mucnt +1 ! count for mu since last time it was zero'ed
            
            ! initially skip, because the integration did not start from 0 
            if ( irsa*irsapr .eq. -1 .and. irsapr .eq. 1 .and. imu .eq. 0 ) then
               imu=imu+1
               xmu = 0.0d0
            endif

            ! Average value over gyration 
            if ( irsa*irsapr .eq. -1 .and. irsapr .eq. 1 .and. imu .ne. 0 ) then ! detect Psi max
               !if(mod(it,imbk).eq.0) then
               write(24,*) Tr,xmu,mucnt !'(1p3d15.7)'
               !endif
               imu=imu+1
               xmu = 0.0d0
               mucnt = 0
            endif
            
            irsapr = irsa; rPsipr = rPsi;
!C          ########## mu ########## end

!C          ########## J ##########
!C          integration along B, data output at every z=0
            xv1int = xv1int + xv1*xv1*dt  !\int \vec{v} \cdot d\vec{l}

            if (zgpr1*zg1 .lt. 0.0d0) then ! detect z=0 cross
               write(22,'(1p3d15.7)') Tr,xv1int*2.0d0,Psi2
               xv1int=0.0D0 ! write and reset the integration value
            endif
!C          ########## J ########## end

!C          # toroidal rotation number 
            if (ygpr1*yg1 .lt. 0.0d0 .and. xg1 .gt. 0.0d0) then ! detect x=0 cross with 0<y 
               ntrpr = ntr;
               ntr = ntr + 1; ! add one 
            endif

            xgpr1=xg1; ygpr1=yg1; zgpr1=zg1; ! save previous values
      
         endif ! get Kpara Kperp, end

         if(mod(it,imbk).eq.0) then
            if (icmp == 0) then
               write(20,'(1p8d15.7)') x(1),x(2),x(3),dsqrt( x(1)*x(1)+x(2)*x(2) ),v(1),v(2),v(3),Tr
            endif
            if (icmp == 1) then
               write(20,'(1p25d15.7,i6)') x(1),x(2),x(3),dsqrt( x(1)*x(1)+x(2)*x(2) ),v(1),v(2),v(3),Tr,&
                                  ken,ken1,ken2,Bx,By,Bz,Btt,Ex,Ey,Ez,Psi1,xg1,yg1,zg1,dsqrt(xg1*xg1+yg1*yg1),&
                                  Psi2,Vtppn*Vtppn/Btt,ntrpr
!     1:x 2:y 3:z 4:r 5:vx 6:vy 7:vz 8:time
!     9:Kttl 10:Kpara 11:Kperp 12:Bx 13:By 14:Bz 15:B 16:Ex 17:Ey 18:Ez
!     19 :Psi 20:xg.c. 21:yg.c. 22:zg.c. 23:rg.c. 24:Psi_g.g. 25:mu,26:t.rotation
            endif
         endif
      endif
!     ## save data, end ##
!     ############################################################


!     ############################################################
!!     ## data for statistical calculation, start ##
      !if (ifile2 == 1) then
      if(mod(it,imbk2).eq.0) then

        stt(ist,ipm) = Tr; ! (time info, particle info) 
        stx(ist,ipm) = x(1);
        sty(ist,ipm) = x(2);
        stz(ist,ipm) = x(3);
        str(ist,ipm) = dsqrt( x(1)*x(1)+x(2)*x(2) );
        stgx(ist,ipm) = xg1;
        stgy(ist,ipm) = yg1;
        stgz(ist,ipm) = zg1;
        stgr(ist,ipm) = dsqrt(xg1*xg1+yg1*yg1);
        stk(ist,ipm) = ken;
        stkpr(ist,ipm) = ken1;
        stkpp(ist,ipm) = ken2;
        stps1(ist,ipm) = Psi1;
        stps2(ist,ipm) = Psi2;
        stmu(ist,ipm) = Vtppn*Vtppn/Btt;
        strt(ist,ipm) = ntrpr;

        ! temporary debugging:
        !write(6,*) ipm,ist,stt(ist,ipm); ! (time info, particle info) 

!      write(6,*) imbk2t,ipmt
!      imbk2t ! istep/imbk2+1 total time steps 
!      ipmt = ipm2-ipm1+1 ! total particle number, or maximum value of ipm 

        ist = ist+1
      endif
      !endif
!     ## data for statistical calculation, end ##
!     ############################################################

      if (iexb == 0) then !C
      CALL bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz)
      endif

      if (iexb == 1) then !C    B interpolation
      CALL bvac2(x(1),x(2),x(3),Bx,By,Bz)
      endif
      B(1)=Bx; B(2)=By; B(3)=Bz;

      if (iexe == 0) then !C
      call evac0(x(1),x(2),x(3),it*dt,Ex,Ey,Ez)
      endif
      if (iexe == 1) then !C
      call evac2(x(1),x(2),x(3),it*dt,tE1,Ex,Ey,Ez)
      endif
      E(1)=Ex; E(2)=Ey; E(3)=Ez;

!      write(6,*) B(1),B(2),B(3),E(1),E(2),E(3)

!C    U_minus
      do i = 1,3
        uminus(i) = u(i) + eom*E(i)*Hr*0.5d0
      end do
      
!C    to calculate gamma^n, instead of gamma^n-1/2, use gamma^2 = 1+(u-/c)^2
      xgm = dsqrt ( 1 + ( uminus(1)*uminus(1)+uminus(2)*uminus(2)+uminus(3)*uminus(3) ) / (xc*xc) ) 
      
!C    T vector
      do i = 1,3
        Tv(i) = eom*B(i)*Hr*0.5d0/xgm ! this \gamma should be ^0, not ^-1/2 
        
!        write(6,*) xgmd(it)
      end do
        
      Tsq = Tv(1)*Tv(1)+Tv(2)*Tv(2)+Tv(3)*Tv(3)

!C    S vector
      do i = 1,3
        Sv(i)= 2.0d0*Tv(i)/(1.0d0+Tsq)
      end do

!C    U_zero
      uzero(1) = uminus(1) + uminus(2)*Tv(3)-uminus(3)*Tv(2)
      uzero(2) = uminus(2) + uminus(3)*Tv(1)-uminus(1)*Tv(3)
      uzero(3) = uminus(3) + uminus(1)*Tv(2)-uminus(2)*Tv(1)

!C    U_plus
      uplus(1) = uminus(1) + uzero(2)*Sv(3)-uzero(3)*Sv(2)
      uplus(2) = uminus(2) + uzero(3)*Sv(1)-uzero(1)*Sv(3)
      uplus(3) = uminus(3) + uzero(1)*Sv(2)-uzero(2)*Sv(1)

!C    U^n+1/2
      do i = 1,3
        u(i) = uplus(i) + eom*E(i)*Hr*0.5d0
      end do

!C     four-vector u(1),u(2),u(3) -> real velocity, beta, gamma. v^2=c^2 u^2 /(c^2+u^2),: ## ^n+1/2 ##
!C    'XU (four velocity), XV (real velocity)'

      xu = dsqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
      xv = xc*xu/dsqrt(xc*xc+xu*xu)
      xbt = xv/xc; xgm = 1.0d0/dsqrt(1.0d0-xbt*xbt)
      v(1)=u(1)/xgm; v(2)=u(2)/xgm; v(3)=u(3)/xgm; 
      ken=(xgm-1.0d0)/1.95695119784913d-6;
 
!C    X^n+1
        do i = 1,3
        x(i) = x(i) + u(i)*Hr/xgm
        end do

!      if(its == 1) then ! variable time step
!      Tr = Tr + Hr
!      endif

      end do
      if (ifile == 1) then
      close(20)
      close(22)
      close(24)
      close(25)
      endif

! 10   if(its == 0) then ! variable time step
10    print "(a35,i15,1p2d14.6)","step and time at stop, dt (s):",it-1,Tr-dt,dt
      
      tend(ipm)=Tr 
      !print "(a35,1p3d13.5)","(x,y,z):",x(1),x(2),x(3)
      !print "(a35,1p3d13.5)","K (eV), pr,pp,ttl:",ken1,ken2,ken
      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,x(1),x(2),x(3),ken1,ken2,ken
      !print "(a35,i6)","toroidal rotation:",ntrpr
      
!      endif

!      if(its == 1) then ! variable time step
!      print "(a35,i15,1p1d15.7)","flight step and time at stop (s):",it-1,Tr
!      endif

      end do ! come from do ipm=ipm1,ipm2      



!     ############################################################
!     minimum statistical calculation
      if (ifile2 == 1) then

      write(6,*) "# Statistical calculation done # TimeStep/ptcl:",imbk2t,ipmt

!     calculate averaged value and sigma at each of time steps 
      call avesig2(stt,stt,imbk2t,ipmt,stta,stts,sttn) 
      call avesig2(stx,stt,imbk2t,ipmt,stxa,stxs,sttn)
      call avesig2(sty,stt,imbk2t,ipmt,stya,stys,sttn)
      call avesig2(stz,stt,imbk2t,ipmt,stza,stzs,sttn)
      call avesig2(str,stt,imbk2t,ipmt,stra,strs,sttn)
      call avesig2(stgx,stt,imbk2t,ipmt,stgxa,stgxs,sttn)
      call avesig2(stgy,stt,imbk2t,ipmt,stgya,stgys,sttn)
      call avesig2(stgz,stt,imbk2t,ipmt,stgza,stgzs,sttn)
      call avesig2(stgr,stt,imbk2t,ipmt,stgra,stgrs,sttn)
      call avesig2(stK,stt,imbk2t,ipmt,stKa,stKs,sttn)
      call avesig2(stKpr,stt,imbk2t,ipmt,stKpra,stKprs,sttn)
      call avesig2(stKpp,stt,imbk2t,ipmt,stKppa,stKpps,sttn)
      call avesig2(stps1,stt,imbk2t,ipmt,stps1a,stps1s,sttn)
      call avesig2(stps2,stt,imbk2t,ipmt,stps2a,stps2s,sttn)
      call avesig2(stmu,stt,imbk2t,ipmt,stmua,stmus,sttn)

      open (51,File="orb_a0bbr.txt") ! save file 0, average
      do isti=1,imbk2t

      write(51,'(1p15d15.7,i6)')&
      stxa(isti),stya(isti),stza(isti),stra(isti),stgxa(isti),stgya(isti),stgza(isti),stgra(isti),&
      stKa(isti),stKpra(isti),stKppa(isti),stps1a(isti),stps2a(isti),stmua(isti),&
      stta(isti),sttn(isti)
!     1:x 2:y 3:z 4:r 5:gx 6:gy 7:gz 8:gr 9:K 10:Kpr 11:Kpp 12:Ps 13:gPs 14:mu 15:time

      end do
      close(51)

      open (52,File="orb_a1bbr.txt") ! save file 1, average+sigma
      do isti=1,imbk2t

      write(52,'(1p15d15.7,i6)')&
      stxa(isti)+stxs(isti),stya(isti)+stys(isti),stza(isti)+stzs(isti),stra(isti)+strs(isti),&
      stgxa(isti)+stgxs(isti),stgya(isti)+stgys(isti),stgza(isti)+stgzs(isti),stgra(isti)+stgrs(isti),&
      stKa(isti)+stKs(isti),stKpra(isti)+stKprs(isti),stKppa(isti)+stKpps(isti),&
      stps1a(isti)+stps1s(isti),stps2a(isti)+stps2s(isti),stmua(isti)+stmus(isti),&
      stta(isti),sttn(isti)

      end do
      close(52)

      open (53,File="orb_a2bbr.txt") ! save file 2, average-sigma
      do isti=1,imbk2t

      write(53,'(1p15d15.7,i6)')&
      stxa(isti)-stxs(isti),stya(isti)-stys(isti),stza(isti)-stzs(isti),stra(isti)-strs(isti),&
      stgxa(isti)-stgxs(isti),stgya(isti)-stgys(isti),stgza(isti)-stgzs(isti),stgra(isti)-stgrs(isti),&
      stKa(isti)-stKs(isti),stKpra(isti)-stKprs(isti),stKppa(isti)-stKpps(isti),&
      stps1a(isti)-stps1s(isti),stps2a(isti)-stps2s(isti),stmua(isti)-stmus(isti),&
      stta(isti),sttn(isti)

      end do
      close(53)
      
      endif
      
!     ############################################################
      
      call system_clock(t2, t_rate, t_max)   ! save ending time
       if ( t2 < t1 ) then
      diff = (t_max - t1) + t2 + 1
        else
       diff = t2 - t1
       endif

      close(21)
!      print "(A29,F10.3,F10.3)", "### calculation time (s) (h):", diff/dble(t_rate), diff/dble(t_rate)/3600.0d0

   return
   end subroutine bbr_particle



!C ###################################################
      subroutine rk4(elik2,elie2,Nr,Tr,Xr,k1,k2,k3,k4,ywork,Hr,tE1)
!     just for initial n^-1/2 value
!     ***** to get ^-1/2 from ^0 *****
!     Runge Kutta
!
!C ###################################################
      implicit none
      double precision,dimension(1:100005) :: elik2,elie2
      double precision :: Hr,Tr,F,Brt,Bzt
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      integer i,Nr
      double precision :: tE1
      
!        call bloop(elik2,elie2,0.1d0,0.2d0,0.0d0,1.0d3,Brt,Bzt)
!        write(6,*) Brt,Bzt,elik2(1000),elie2(2000)

      do 10 i=1,Nr
  10  k1(i)=Hr*F(elik2,elie2,i,Tr,Xr,Nr,tE1)
      do 11 i=1,Nr
  11  ywork(i)=Xr(i)+0.5d0*k1(i)
      do 12 i=1,Nr
  12  k2(i)=Hr*F(elik2,elie2,i,Tr+0.5d0*Hr,ywork,Nr,tE1)
      do 13 i=1,Nr
  13  YWORK(i)=Xr(i)+0.5d0*k2(i)
      do 14 i=1,Nr
  14  k3(i)=Hr*F(elik2,elie2,i,Tr+0.5d0*Hr,ywork,Nr,tE1)
      do 15 i=1,Nr
  15  ywork(i)=Xr(i)+K3(i)
      do 16 i=1,Nr
  16  k4(i)=Hr*F(elik2,elie2,i,Tr+Hr,ywork,Nr,tE1)
      do 17 i=1,Nr
  17  Xr(i)=Xr(i)+(k1(I)+2.0d0*k2(i)+2.0d0*k3(i)+k4(i))/6.0d0
      Tr=Tr+Hr
!C      write(6,*)T,H
      RETURN
      end subroutine rk4


!C ###################################################
      function F(elik2,elie2,i,t,x,n,tE1)
!     ! -t ! for backward calculation
!
!C ###################################################
      use mod_variable
      implicit none
      double precision,dimension(0:100005) :: elik2,elie2
      double precision :: t,F,XQOM,xu,xv,xbt,xgm
      double precision,dimension(1:3) :: B,E
      double precision,dimension(1:6) :: x
      integer :: n,i
      double precision :: xp,yp,zp,Ex,Ey,Ez,Bx,By,Bz
      double precision :: tE1

      xqom = 1.7588200236d11

      xp=x(1); yp=x(2); zp=x(3);

!C     four-velocity X(4),X(5),X(6) -> real velocity, v^2=c^2 u^2 /(c^2+u^2)
      xu = dsqrt(x(4)*x(4)+x(5)*x(5)+x(6)*x(6))
      xv = XC*xu/dsqrt(xc*xc+xu*xu)
!C      XV=XU/dsqrt(1.0+(XU/XC)*(XU/XC)) <-larger error 

!C      WRITE(6,*) 'XU (four velocity), XV (real velocity)'
!C      WRITE(6,602) XU, XV

      xbt=xv/xc
      xgm=1.0d0/dsqrt(1.0d0-xbt*xbt)
!C      KENE=XME*XC*XC*(XGM-1.0)


      if (iexb == 0) then !C    B by Biot-Savart 
      call bvac(elik2,elie2,XP,YP,ZP,Bx,By,Bz);
      endif
      if (iexb == 1) then !C    interpolate B in external files 
      call bvac2(XP,YP,ZP,Bx,By,Bz);
      endif
      B(1)=Bx; B(2)=By; B(3)=Bz;
      
      if (iexe == 0) then !C    E from subroutine calculation 
      call evac0(XP,YP,ZP,T,Ex,Ey,Ez);
      endif
      if (iexe == 1) then !C    interpolate E in external files 
      call evac2(XP,YP,ZP,T,tE1,Ex,Ey,Ez);
      endif
      E(1)=Ex; E(2)=Ey; E(3)=Ez;


      GO TO (11,12,13,14,15,16),I
  11  F=-X(4)/xgm
      RETURN
  12  F=-X(5)/xgm
      RETURN
  13  F=-X(6)/xgm
      RETURN
  14  F=-XQOM*( E(1) + (X(5)*B(3)-X(6)*B(2))/xgm )
      RETURN
  15  F=-XQOM*( E(2) + (X(6)*B(1)-X(4)*B(3))/xgm )
      RETURN
  16  F=-XQOM*( E(3) + (X(4)*B(2)-X(5)*B(1))/xgm )
      RETURN
 
      END

!     rk4 should include el in oder to get B
!C     =================================================================
      subroutine bvac(elik2,elie2,xp,yp,zp,Bx,By,Bz)
!C     calculate B at (x,y,z) 
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,Bx,By,Bz,xi,Bro,Bzo,ctmp,Brt,Bzt,Bx1,By1,Bz1
      double precision,dimension(1:100005) :: elik2,elie2
      integer i

!      ci(1) = 1.0d3 ! current
!      rc(1) = 1.0d-1 ! loop radius
!      zcl(1) = 0.0d0 ! loop z
!      xcl(1) = 0.0d0 ! loop offset in x direction from z axis?
!      ycl(1) = 0.0d0 ! usually zero or middle

        Bx=0.0d0; By=0.0d0; Bz=0.0d0;

        do i=1,icln ! add for each of coils 
        
        ctmp = dsqrt((xp-xcl(i))*(xp-xcl(i))+(yp-ycl(i))*(yp-ycl(i)))
        call bloop(elik2,elie2,rc(i),ctmp,zp-zcl(i),ci(i),Bro,Bzo)

!        write(6,*) ctmp

        if (ctmp .lt. 1.0d-12) then
          Bx1=0.0d0; By1=0.0d0; Bz1=Bzo; ! in order not to be divided by zero
        else
          Bx1 = Bro*(xp-xcl(i))/ctmp
          By1 = Bro*(yp-ycl(i))/ctmp
          Bz1 = Bzo
        endif

          Bx = Bx + Bx1
          By = By + By1
          Bz = Bz + Bz1

!        write(6,*) Bz1,Bz

        end do

! Homogeneous B
!        Bx = 0.0d0;
!        By = 0.0d0;
!        Bz = 0.1d0;

! Homogeneous B in finite region
!      if ( xp .gt. 0.1d0 .and. xp .lt. 0.2d0 ) then
!      if ( yp .gt. 0.1d0 .and. yp .lt. 0.2d0 ) then
!      if ( zp .gt. 0.1d0 .and. zp .lt. 0.2d0 ) then
!        Bx = 0.0d0;
!        By = 0.0d0;
!        Bz = 0.1d0;
!      endif
!      endif
!      endif

!        write(6,*) Bz

      RETURN
      END



!C     =================================================================
      subroutine bmap(elik2,elie2,xp,yp,zp,Psi)
!C     return B at (x,y,z) obtained by interpolation xyz
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,ctmp,at,Psi,Psi1
      double precision,dimension(0:100005) :: elik2,elie2
      integer i

      Psi = 0.0d0

      do i=1,icln !  add for each of coils

      ctmp = dsqrt((xp-xcl(i))*(xp-xcl(i))+(yp-ycl(i))*(yp-ycl(i)))        
      call rathe(elik2,elie2,rc(i),ctmp,zp-zcl(i),ci(i),at)
      Psi1 = 2.0d0*pi*at
      Psi = Psi + Psi1
      
      end do
      
      return
      end

!C     =================================================================
      subroutine bvac2(xp,yp,zp,Bx,By,Bz);
!C     return B at (x,y,z) obtained by interpolation 
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp
      double precision :: Bx,By,Bz
      double precision :: f111,f211,f121,f221,f112,f212,f122,f222,A1,A2,A3,A4,B1,B2
      
!      write(6,*) "initial position:"
!      write(6,*) xp,yp,zp ! serch (0,0,0) for this point 
!      write(6,*) ( (xp-xmin)/dx ),int( (xp-xmin)/dx )
!      write(6,*) ( (yp-ymin)/dy ),int( (yp-ymin)/dy )
!      write(6,*) ( (zp-ymin)/dz ),int( (zp-ymin)/dz )
      n0x=int( (xp-xmin)/dx )+1; n0y=int( (yp-ymin)/dy )+1; n0z=int( (zp-ymin)/dz )+1;

!      write(6,*) "n0x,n0y,n0z at f111:"
!      write(6,*) n0x,n0y,n0z !
!      write(6,*) "x,y,z at f111:"
!      write(6,*) mt1x(n0x),mt1y(n0y),mt1z(n0z) ! this corresponds to f111
!      write(6,*) "B at f111:"
!      write(6,*) Bxm(n0x,n0y,n0z),Bym(n0x,n0y,n0z),Bzm(n0x,n0y,n0z)
      
      f111 = Bxm(n0x,n0y,n0z); f211 = Bxm(n0x+1,n0y,n0z); f121 = Bxm(n0x,n0y+1,n0z); f221 = Bxm(n0x+1,n0y+1,n0z);
      f112 = Bxm(n0x,n0y,n0z+1); f212 = Bxm(n0x+1,n0y,n0z+1); f122 = Bxm(n0x,n0y+1,n0z+1); f222 = Bxm(n0x+1,n0y+1,n0z+1);
          A1 = f111 + ( xp-mt1x(n0x) )*( f211-f111 )/dx; A2 = f121 + ( xp-mt1x(n0x) )*( f221-f121 )/dx;
          B1 = A1 + ( yp-mt1y(n0y) )*( A2-A1 )/dy;
          A3 = f112 + ( xp-mt1x(n0x) )*( f212-f112 )/dx; A4 = f122 + ( xp-mt1x(n0x) )*( f222-f122 )/dx;
          B2 = A3 + ( yp-mt1y(n0y) )*( A4-A3 )/dy;
          Bx = B1 + ( zp-mt1z(n0z) )*( B2-B1 )/dz;

      f111 = Bym(n0x,n0y,n0z); f211 = Bym(n0x+1,n0y,n0z); f121 = Bym(n0x,n0y+1,n0z); f221 = Bym(n0x+1,n0y+1,n0z);
      f112 = Bym(n0x,n0y,n0z+1); f212 = Bym(n0x+1,n0y,n0z+1); f122 = Bym(n0x,n0y+1,n0z+1); f222 = Bym(n0x+1,n0y+1,n0z+1);
          A1 = f111 + ( xp-mt1x(n0x) )*( f211-f111 )/dx; A2 = f121 + ( xp-mt1x(n0x) )*( f221-f121 )/dx;
          B1 = A1 + ( yp-mt1y(n0y) )*( A2-A1 )/dy;
          A3 = f112 + ( xp-mt1x(n0x) )*( f212-f112 )/dx; A4 = f122 + ( xp-mt1x(n0x) )*( f222-f122 )/dx;
          B2 = A3 + ( yp-mt1y(n0y) )*( A4-A3 )/dy;
          By = B1 + ( zp-mt1z(n0z) )*( B2-B1 )/dz;

      f111 = Bzm(n0x,n0y,n0z); f211 = Bzm(n0x+1,n0y,n0z); f121 = Bzm(n0x,n0y+1,n0z); f221 = Bzm(n0x+1,n0y+1,n0z);
      f112 = Bzm(n0x,n0y,n0z+1); f212 = Bzm(n0x+1,n0y,n0z+1); f122 = Bzm(n0x,n0y+1,n0z+1); f222 = Bzm(n0x+1,n0y+1,n0z+1);
          A1 = f111 + ( xp-mt1x(n0x) )*( f211-f111 )/dx; A2 = f121 + ( xp-mt1x(n0x) )*( f221-f121 )/dx;
          B1 = A1 + ( yp-mt1y(n0y) )*( A2-A1 )/dy;
          A3 = f112 + ( xp-mt1x(n0x) )*( f212-f112 )/dx; A4 = f122 + ( xp-mt1x(n0x) )*( f222-f122 )/dx;
          B2 = A3 + ( yp-mt1y(n0y) )*( A4-A3 )/dy;
          Bz = B1 + ( zp-mt1z(n0z) )*( B2-B1 )/dz;

!      Bx=0.0d0
!      By=0.0d0
!      Bz=0.01d0

      RETURN
      END

!C     =================================================================
      subroutine evac0(xp,yp,zp,tp,Ex,Ey,Ez)
!C     return azimuthal E at (x,y,z) 
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,tp,Ex,Ey,Ez
      double precision :: rr
      
      if ( xp .gt. 0.1d0 .and. xp .lt. 0.2d0 ) then
      if ( yp .gt. 0.1d0 .and. yp .lt. 0.2d0 ) then
      if ( zp .gt. 0.1d0 .and. zp .lt. 0.2d0 ) then
        !rr = dqsrt(xp*xp + yp*yp)
        !Ex = -Vac(1) * dsin(xp/r)
        !Ey =  Vac(1) * dcos(xp/r)

        Ex = 0.0d0;
        Ey = 0.0d0;
        Ez = 0.1d0;
      endif
      endif
      endif

      RETURN
      END

!C     =================================================================
      subroutine evac(Ex,Ey,Ez)
!C     return E at (x,y,z) 
!C     =================================================================
      implicit none
      double precision :: Ex,Ey,Ez
      
! Homogeneous E
      Ex = 0.0d0
      Ey = 0.0d0
      Ez = 0.0d0

      RETURN
      END


!C     =================================================================
      subroutine evac2(xp,yp,zp,tp,tE1,Ex,Ey,Ez);
!C     return E at (x,y,z) obtained by interpolation and superposition
!C    
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,tp,tE1,Ex,Ey,Ez
      double precision :: f111,f211,f121,f221,f112,f212,f122,f222,A1,A2,A3,A4,B1,B2
      double precision :: RWamp,RWfreq,RWphase
      double precision :: fchirp
      integer :: iex1,iex2,iex3

      if (tp .gt. tE1) then
         ! Zero out E field after specified time
         Ex=0.0d0; Ey=0.0d0; Ez=0.0d0;
         
      else

      n0x=int( (xp-xmin)/dx )+1; n0y=int( (yp-ymin)/dy )+1; n0z=int( (zp-ymin)/dz )+1;
      
      ! set RW frequency to increase linearly
      fchirp = fgrad * tp + Fel
      Vel(1) = Vdc(1) + Vac(1)*dsin( 2.0d0*pi*fchirp*tp + Phrw*0.0d0*pi/2.0d0 )
      Vel(2) = Vdc(2) + Vac(2)*dsin( 2.0d0*pi*fchirp*tp + Phrw*1.0d0*pi/2.0d0 )
      Vel(3) = Vdc(3) + Vac(3)*dsin( 2.0d0*pi*fchirp*tp + Phrw*2.0d0*pi/2.0d0 )
      Vel(4) = Vdc(4) + Vac(4)*dsin( 2.0d0*pi*fchirp*tp + Phrw*3.0d0*pi/2.0d0 )

!C    Ex by interpolation 
      f111=0.0d0; f211=0.0d0; f121=0.0d0; f221=0.0d0; f112=0.0d0; f212=0.0d0; f122=0.0d0; f222=0.0d0; 
      do iele = 1,ielemax
      f111 = f111 + Vel(iele)*Exm(n0x,n0y,n0z,iele);
      f211 = f211 + Vel(iele)*Exm(n0x+1,n0y,n0z,iele);

      f121 = f121 + Vel(iele)*Exm(n0x,n0y+1,n0z,iele);
      f221 = f221 + Vel(iele)*Exm(n0x+1,n0y+1,n0z,iele);
      f112 = f112 + Vel(iele)*Exm(n0x,n0y,n0z+1,iele);
      f212 = f212 + Vel(iele)*Exm(n0x+1,n0y,n0z+1,iele);
      f122 = f122 + Vel(iele)*Exm(n0x,n0y+1,n0z+1,iele);
      f222 = f222 + Vel(iele)*Exm(n0x+1,n0y+1,n0z+1,iele);
      end do
          A1 = f111 + ( xp-mt1x(n0x) )*( f211-f111 )/dx; A2 = f121 + ( xp-mt1x(n0x) )*( f221-f121 )/dx;
          B1 = A1 + ( yp-mt1y(n0y) )*( A2-A1 )/dy;
          A3 = f112 + ( xp-mt1x(n0x) )*( f212-f112 )/dx; A4 = f122 + ( xp-mt1x(n0x) )*( f222-f122 )/dx;
          B2 = A3 + ( yp-mt1y(n0y) )*( A4-A3 )/dy;
          Ex = B1 + ( zp-mt1z(n0z) )*( B2-B1 )/dz;

!C    Ey 
      f111=0.0d0; f211=0.0d0; f121=0.0d0; f221=0.0d0; f112=0.0d0; f212=0.0d0; f122=0.0d0; f222=0.0d0; 
      do iele = 1,ielemax
      f111 = f111 + Vel(iele)*Eym(n0x,n0y,n0z,iele);
      f211 = f211 + Vel(iele)*Eym(n0x+1,n0y,n0z,iele);
      f121 = f121 + Vel(iele)*Eym(n0x,n0y+1,n0z,iele);
      f221 = f221 + Vel(iele)*Eym(n0x+1,n0y+1,n0z,iele);
      f112 = f112 + Vel(iele)*Eym(n0x,n0y,n0z+1,iele);
      f212 = f212 + Vel(iele)*Eym(n0x+1,n0y,n0z+1,iele);
      f122 = f122 + Vel(iele)*Eym(n0x,n0y+1,n0z+1,iele);
      f222 = f222 + Vel(iele)*Eym(n0x+1,n0y+1,n0z+1,iele);
      end do
          A1 = f111 + ( xp-mt1x(n0x) )*( f211-f111 )/dx; A2 = f121 + ( xp-mt1x(n0x) )*( f221-f121 )/dx;
          B1 = A1 + ( yp-mt1y(n0y) )*( A2-A1 )/dy;
          A3 = f112 + ( xp-mt1x(n0x) )*( f212-f112 )/dx; A4 = f122 + ( xp-mt1x(n0x) )*( f222-f122 )/dx;
          B2 = A3 + ( yp-mt1y(n0y) )*( A4-A3 )/dy;
          Ey = B1 + ( zp-mt1z(n0z) )*( B2-B1 )/dz;

!C    Ez 
      f111=0.0d0; f211=0.0d0; f121=0.0d0; f221=0.0d0; f112=0.0d0; f212=0.0d0; f122=0.0d0; f222=0.0d0; 
      do iele = 1,ielemax
      f111 = f111 + Vel(iele)*Ezm(n0x,n0y,n0z,iele);
      f211 = f211 + Vel(iele)*Ezm(n0x+1,n0y,n0z,iele);
      f121 = f121 + Vel(iele)*Ezm(n0x,n0y+1,n0z,iele);
      f221 = f221 + Vel(iele)*Ezm(n0x+1,n0y+1,n0z,iele);
      f112 = f112 + Vel(iele)*Ezm(n0x,n0y,n0z+1,iele);
      f212 = f212 + Vel(iele)*Ezm(n0x+1,n0y,n0z+1,iele);
      f122 = f122 + Vel(iele)*Ezm(n0x,n0y+1,n0z+1,iele);
      f222 = f222 + Vel(iele)*Ezm(n0x+1,n0y+1,n0z+1,iele);
      end do
          A1 = f111 + ( xp-mt1x(n0x) )*( f211-f111 )/dx; A2 = f121 + ( xp-mt1x(n0x) )*( f221-f121 )/dx;
          B1 = A1 + ( yp-mt1y(n0y) )*( A2-A1 )/dy;
          A3 = f112 + ( xp-mt1x(n0x) )*( f212-f112 )/dx; A4 = f122 + ( xp-mt1x(n0x) )*( f222-f122 )/dx;
          B2 = A3 + ( yp-mt1y(n0y) )*( A4-A3 )/dy;
          Ez = B1 + ( zp-mt1z(n0z) )*( B2-B1 )/dz;

       endif
      RETURN
      END



!*     =================================================================
      SUBROUTINE BLOOP(elik2,elie2,RC,R,Z,CI,BR,BZ)
!*     calculate the Br and Bz produced by a loop current
!*     RC:loop radius R,Z:position CI:coil current 
!*     =================================================================
!*                                                   Biot-Savalt formula
!*                                                   ===================
      implicit none
      double precision, dimension(0:100005) :: elik2,elie2
      DOUBLE precision :: ZK,ZZK,RC,R,Z,CI,BR,BZ,G1,FK,FE,A,G,FF,E,H,xcc,fk1,fe1,fk2,fe2
      integer :: i,i2
      
      xcc=2.0d-7
!C
      ZK=4.0D0*RC*R/((RC+R)*(RC+R)+Z*Z)
      ZZK=dsqrt(ZK)

!      write(6,*) "k",ZZK
!      write(6,*) "K",elik2(100)
!      write(6,*) "E",elie2(ZK)

      IF(ZK.GE.0.999999D0) GO TO 11

      G1=dsqrt(1.0D0-ZK)
      IF(ZK.GT.0.9998D0) FK=DLOG(4.0D0/G1)+0.25D0*(DLOG(4.0D0/G1)-1.0D0)*G1*G1
      IF(ZK.GT.0.9998D0) FE=1.0D0+0.5D0*(DLOG(4.0D0/G1)-0.5D0)*G1*G1     
      IF(ZK.GT.0.9998D0) GO TO 20

      I=IDINT(ZZK*100000.0D0)

!C    FK FE interpolation start
      i2=i+1;
      fk1=elik2(i); fe1=elie2(i);
      fk2=elik2(i2); fe2=elie2(i2);
!      write(6,'(d13.5,i7,1p3d17.9)') zzk,i,fk1,fk2,fk1+(zzk-dble(i)/1.0d5)*(fk2-fk1)
      fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
      fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);
!      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
!C    FK FE interpolation end

!C    FK FE
!      FK=elik2(I)
!      FE=elie2(I)
!      write(6,*) "I",I      
!      write(6,*) "FK,FE",FK,FE
!      write(6,*) "BR,BZ",BR,BZ

   20 A=XCC*CI/dsqrt((RC+R)*(RC+R)+Z*Z)
      G=RC*RC+R*R+Z*Z
      FF=(RC-R)*(RC-R)+Z*Z
      E=A*Z/R
      H=RC*RC-R*R-Z*Z
      BZ=A*(H*FE/FF+FK)
      BR=E*(G*FE/FF-FK)

      IF(I.EQ.0) THEN
      BR=0.0d0
      ENDIF
      
      RETURN
      
   11 BZ=0.0D0
      BR=0.0D0
      RETURN
      END


!*     =================================================================
      SUBROUTINE rathe(elik2,elie2,rc,x,z,ci,at)
!*     calculate "Atheta" produced by a loop current
!*     RC:loop radius R,Z:position CI:coil current 
!*     =================================================================
!*                                                   Biot-Savalt formula
!*                                                   ===================
      implicit none
      double precision, intent(in), dimension(0:100005) :: elik2,elie2
      double precision ZK,ZZK,A0,A1,A2,FK,FE,xcc,fk1,fe1,fk2,fe2
      double precision, intent(in) :: rc,x,z,ci
      double precision, intent(out) :: AT
      integer i,i2

      DATA XCC/2.d-7/
!C
!*     ZK is k^2 
      ZK=4.0d0*rc*x/((rc+x)*(rc+x)+z*z)
      ZZK=dsqrt(ZK)
      
      IF(ZK.GE.0.999999d0) GO TO 20
      IF(ZK.GT.0.9998d0) GO TO 12
      
      A0=2.0d0*XCC*CI*x/dsqrt(ZK)
      A1=dsqrt(RC/x)
      A2=1.0d0-0.5d0*ZK
      I=IDINT(ZZK*100000.0d0)
      
!C    FK FE interpolation start
      i2=i+1;
      fk1=elik2(i); fe1=elie2(i);
      fk2=elik2(i2); fe2=elie2(i2);
!      write(6,'(d13.5,i7,1p3d17.9)') zzk,i,fk1,fk2,fk1+(zzk-dble(i)/1.0d5)*(fk2-fk1)
      fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
      fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);
!      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
!C    FK FE interpolation end

!      FK=elik2(I)
!      FE=elie2(I)

      AT=A0*A1*(A2*FK-FE)

      IF(I.EQ.0) THEN
      AT=0.0d0
      ENDIF
      
      RETURN

 12   A1=XCC*CI*RC
      A2=dlog(8.0d0*x/dsqrt((x-RC)*(x-RC)+Z*Z))
      AT=A1*(A2-2.0d0)
      RETURN

 20   AT=0.0d0
      RETURN
      END


!     =================================================================
      SUBROUTINE EFILE2(elik2,elie2)
!     calculation speed may be adjusted by setting nmax x2
!     =================================================================
!                                                    elliptic functions
!                                                    ==================
      implicit none
      double precision,dimension(1:100005) :: elik2,elie2
      double precision ZZZAAA,FEF1,FEF2
      integer ji
      
      write(6,*) "start: elliptic fuction calculation"
      write(6,*) "... may takes several seconds ..."

      DO 5 JI=1,99999
         ZZZAAA=DBLE(JI)/100000.0d0
         ELIK2(JI)=FEF1(ZZZAAA)
         ELIE2(JI)=FEF2(ZZZAAA)

!CC       show progress
      IF ( MOD(JI,20000).EQ.0 ) THEN
      WRITE(6,*) DBLE(JI)/100000.0*100.0,'% '
      ENDIF

    5 CONTINUE
    
      WRITE(6,*) DBLE(1)/1.0*100.0,'% '
      write(6,*) "end: elliptic fuction calculation"
      
      RETURN
      END

     ! ---------------------------
      double precision function FEF1(k)  ! K(k)
        implicit none
        double precision, intent(in) :: k  ! 
        double precision :: pi, m, dt, t, tmin, tmax
        integer :: i
        integer, parameter :: nmax=2000
        double precision :: f, x

        f(m,x) = 1.0d0/dsqrt(1.0d0-(m*dsin(x))**2)

        if(k.ge.1.0d0)then
          write(*,*) "(error ! : k must 0=<k<1.)"
      return
       end if

      pi = 3.14159265358979323846d0

       tmin = 0.0d0
        tmax = pi/2.0d0
        dt = (tmax-tmin)/dble(nmax-1)

       FEF1 = 0.5d0*dt*(f(k,tmin)+f(k,tmax))
       do i=1,nmax-2
          t = tmin+dt*dble(i)
          FEF1 = FEF1+dt*f(k,t)
        end do
      
       return
      end function

      ! -------------------------------------
      double precision function FEF2(k)  ! E(k)
       implicit none
        double precision, intent(in) :: k  ! 
       double precision :: pi, m, dt, t, tmin, tmax
       integer :: i
       integer, parameter :: nmax=2000
       double precision :: f, x

        f(m,x) = dsqrt(1.0d0-(m*dsin(x))**2)

        pi = 3.14159265358979323846d0

       if(k.gt.1.0d0)then
         write(*,*) "(error) ! : k must 0=<k=<1."
          return
       end if

       tmin = 0.0d0
       tmax = pi/2.0d0
        dt = (tmax-tmin)/dble(nmax-1)

       FEF2 = 0.5d0*dt*(f(k,tmin)+f(k,tmax))
       do i=1,nmax-2
         t = tmin+dt*dble(i)
         FEF2 = FEF2+dt*f(k,t)
       end do

       return
      end function





!C    =================================================================
      subroutine avesig(data1,itt1,itt2,ave1,sig1)
!     i,j: time and particle, averaged for j.
!     ave1(i,j) is averaged for j at i ->ave2(i), sig2(i)
!     time step i:itt1, particle number j:itt2
!C    =================================================================
      implicit none
      double precision,dimension(1:105,1:1005) :: data1
      double precision,dimension(1:105) :: ave1,sig1
      integer iav1,iav2,itt1,itt2
      double precision sum1,sum2

!     ##### averaging start #####
      do iav1=1,itt1 ! for each of the time step 

         sum1 = 0.0d0
         do iav2 = 1,itt2 ! particle number 
         sum1 = sum1 + data1(iav1,iav2)
         end do

         ave1(iav1) = sum1/dble(itt2)

!      write(6,*)iav1, ave1(iav1)

      end do
!     ##### averaging end #####

!     ##### sigma start #####
      do iav1=1,itt1 ! for each of the time step 
      
         sum2 = 0.0d0
         do iav2 = 1,itt2 ! particle number 
         sum2 = sum2 + ( data1(iav1,iav2)-ave1(iav1) )*( data1(iav1,iav2)-ave1(iav1) )
         end do
         
         sig1(iav1) = dsqrt( sum2/itt2 )

!      write(6,*)iav1, sig1(iav1)

      end do
!     ##### sigma end #####

      return
      end


!C    =================================================================
      subroutine avesig2(data1,data2,itt1,itt2,ave1,sig1,flp1) 
!     When using subroutine (not this one) avesig, with ihit=1 and ieout=2, there is a problem that averaging is done even when a
!     particle is no longer in the trapping region (-> data is "0"). This results in a failure in averaging. In order to solve
!     this problem, subroutine avesig2 was added.
!     When "data2" is set to be the time data, averaging is done only when the time data is not zero.
!     
!     i,j: time and particle, averaged for j.
!     ave1(i,j) is averaged for j at i ->ave2(i), sig2(i)
!     time step i:itt1, particle number j:itt2
!C    =================================================================
      implicit none
      double precision,dimension(1:105,1:1005) :: data1,data2 ! (time info, particle info)
      double precision,dimension(1:105) :: ave1,sig1 ! (time info)
      integer,dimension(1:105) :: flp1 ! (time info)
      integer iav1,iav2,itt1,itt2,iflp
      double precision sum1,sum2

!     ##### averaging start #####
      do iav1=1,itt1 ! for each of the time step 
         sum1 = 0.0d0
         iflp = 0
         do iav2 = 1,itt2 ! for each of the particles

           if ( iav1 .eq. 1 ) then ! always do when t=0 
           sum1 = sum1 + data1(iav1,iav2)
           iflp = iflp + 1
           endif
           
!          only when the time data is not 0 ("0" means orbit calculation is already finished for that particle)         
           if ( iav1 .ne. 1 .and. data2(iav1,iav2) .ne. 0.0d0 ) then
           sum1 = sum1 + data1(iav1,iav2)
           iflp = iflp + 1
           endif

         end do

         flp1(iav1)=iflp
         ave1(iav1) = sum1/dble(iflp)
!         ave1(iav1) = sum1/dble(itt2)

!      write(6,*)iav1, ave1(iav1)

      end do
!     ##### averaging end #####

!     ##### sigma start #####
      do iav1=1,itt1 ! for each of the time step
      
         sum2 = 0.0d0
         iflp = 0
         do iav2 = 1,itt2 ! particle number 

           if ( iav1 .eq. 1 ) then ! always do when t=0  
           sum2 = sum2 + ( data1(iav1,iav2)-ave1(iav1) )*( data1(iav1,iav2)-ave1(iav1) )
           iflp = iflp + 1
           endif

!          only when the time data is not 0 ("0" means orbit calculation is finished for that particle)
           if ( iav1 .ne. 1 .and. data2(iav1,iav2) .ne. 0.0d0 ) then
           sum2 = sum2 + ( data1(iav1,iav2)-ave1(iav1) )*( data1(iav1,iav2)-ave1(iav1) )
           iflp = iflp + 1
           endif
         
         end do
         
         sig1(iav1) = dsqrt( sum2/iflp )
!         sig1(iav1) = dsqrt( sum2/itt2 )

!      write(6,*)iav1, sig1(iav1)

      end do
!     ##### sigma end #####

      return
      end
