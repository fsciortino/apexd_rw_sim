!     ======================================================================
      module mod_variable
!     Collection of variable definitions that are shared between subroutines
!    =======================================================================
      implicit none

      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, parameter :: xc = 2.99792458d8
      double precision, parameter :: xe = 1.6021766208d-19
      double precision, parameter :: xm = 0.910938356d-30
      double precision, parameter :: eom = 1.7588200236d11
      
      double precision :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz
      integer :: nx,ny,nz,n0x,n0y,n0z,iele,ielemax,icln
      double precision,dimension(1:205) :: mt1x,mt1y,mt1z
      double precision,dimension(1:1,1:1,1:1) :: Bxm,Bym,Bzm ! to save memory  
      !double precision,dimension(1:105,1:105,1:105,1:4) :: Exm,Eym,Ezm ! electric fields (unit field by electrodes)
      ! NB: electric fields are heavy on RAM. Reduce allocated memory if possible.
      double precision,dimension(1:52,1:52,1:52,1:8) :: Exm,Eym,Ezm ! electric fields (unit field by electrodes)  
      double precision,dimension(1:12) :: Vel,Vdc,Vac ! electrode voltage
      double precision :: Fel,Phrw,fgrad ! frequency,phase and gradient of RW 
      double precision,dimension(1:52) :: ci,rc,zcl,xcl,ycl ! Coil current, loop radius, z, x, y (y is usually 0)
      end module mod_variable

   
!     =================================================================
      subroutine bbr_single_sim(pmx,pmy,pmz,pmk1,pmk2,pma,pmrwa,pmrwf,pmrwfg,pmrwd, &
           dt,iend,imbk,tE1,noco,Icl_in,rcl_in,xcl_in,ycl_in,zcl_in,pid,iexe,mm,Evers,iexb)  
!    
!     This version was modified by FS in order to compile into a Python shared library.
!     It differs from other versions of the code in the fact that it is optimized to run a single particle,
!        so that it may be easier to compile with f2py and parallelize. Strictly speaking, this is a reduction in
!        capabilities and scope wrt bbr_particle.f90, but it proves convenient for many tasks, so it has been
!        made into a separate subroutine. 
!
!     Compile into a Python module with        
!          f2py -c -m apds bbr_single_particle.f90 only: bbr_single_sim        
!
!     Use this with run_single_particle.py.
!
!        Inputs:
!           (pmx,pmy,pmz) : initial particle position
!           (pmk1,pmk2) : parallel and perpendicular initial energy
!           pma : initial orientation of v_perp vector
!           (pmrwa,pmrwf,pmrfg,pmrwd) : amplitude, starting frequency, frequency gradient and direction of RW E-field
!           dt : time step [seconds]. Typical value:  1.0d-12
!           iend : maximum number of steps 
!           imbk : trajectory writing rate. Every `imbk' steps data will be saved
!           imbk2 : similar to imbk, but for the gyroaveraged \mu and bounce-averaged J
!           tE1 : end time [seconds] for RW field. Use this to interrupt the RW field at some time during the simulation
!           (Icl_in,rcl_in,xcl_in,ycl_in,zcl_in) : current, radius, (x,y,z) of the current loop center.
!                            These can be given as arrays, in which case the field from each current loop will be summer.
!           pid : Process ID (PID) for the current run
!           iexe : option to obtain E-field: (0) Use E=0 always, everywhere; 
!                                            (1) Use ideal azimuthal E-field
!                                            (2) Read E-field file from  "./simion_Efields/ver.#";
!                                            (3) Locally-rotating E-field vortex        
!           mm : azimuthal mode number to be used, only used when iexe=1
!           Evers : version of the E-field file to be read from "./simion_Efields/ver.#", only used when iexe=2
!           iexb: option to obtain B-field:  (0): Biot-Savart calculation
!                                            (1): set Bz=0.1 everywhere
!                                            (2): use external B field file (not well tested)
!        Outputs: currently written to files in ./orb/***{}_{}bbr.txtm, where *** is an identifier for the content
!                 of the file, the first {} is the particle identifier (integer) and the second {} is the process ID
!                 that produced the file. 
!           
!        Note that it is possible to modify this code to have data directly given as subroutine outputs
!        (i.e. kept in memory) but this might run into a memory limit for long/complex simulations.
!        
!    =================================================================
      use mod_variable
      implicit none
  
      double precision, INTENT(IN)        :: pmx
      double precision, INTENT(IN)        :: pmy
      double precision, INTENT(IN)        :: pmz
      double precision, INTENT(IN)        :: pmk1
      double precision, INTENT(IN)        :: pmk2
      double precision, INTENT(IN)        :: pma
      double precision, INTENT(IN)        :: pmrwa
      double precision, INTENT(IN)        :: pmrwf
      double precision, INTENT(IN)        :: pmrwfg
      double precision, INTENT(IN)        :: pmrwd
      double precision, INTENT(IN)        :: dt      
      integer, INTENT(IN)                 :: iend
      integer, INTENT(IN)                 :: imbk      
      double precision, INTENT(IN)        :: tE1
      integer, INTENT(IN)                 :: noco  
      double precision, INTENT(IN)        :: Icl_in(noco)
      double precision, INTENT(IN)        :: rcl_in(noco)
      double precision, INTENT(IN)        :: xcl_in(noco)
      double precision, INTENT(IN)        :: ycl_in(noco)
      double precision, INTENT(IN)        :: zcl_in(noco)
      integer, INTENT(IN)                 :: pid
      integer, INTENT(IN)                 :: iexe 
      integer, INTENT(IN)                 :: Evers  
      integer, INTENT(IN)                 :: mm 
      integer, INTENT(IN)                 :: iexb 
      
      ! Internal variables
      double precision,dimension(1:3) :: x,v
      double precision :: ken, ken1, ken2
      double precision :: Tr
      double precision, dimension(1:3) :: u,uminus,uzero,uplus
      double precision,dimension(1:3) :: B,E,Tv,Sv
      double precision :: tini,Tsq,ddt,xu,xv, kini,xdr,ydr,zdr,kini1,kini2,xpr,ypr,zpr,rsa,rPsi,rPsipr
      double precision :: xv1,xv2,xv1int,kang,bmb1,bmb2,tmax,tsdiv,dt2
      double precision :: xi,yi,zi,Ex,Ey,Ez,Bx,By,Bz
      integer :: i,iini,Nr,itf,ifile,iclmap,iclbp,iclbp2,imu,ntr,ntrpr
      integer :: ikene,its
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      double precision :: Hr,xgm,xbt,Btt,Psi,Psi1,Psi2,xg1,yg1,zg1,xgpr1,ygpr1,zgpr1
      double precision :: gr,xmu
      integer :: it 
      integer t1, t2, t_rate, t_max, diff,iex1,iex2,iex3,iex4,icmp,ihit,ieout,ieoutc,iif,irsa,irsapr
      double precision,dimension(1:100005) :: elik2,elie2
      double precision :: pmk,pmrx,pmry,pmrz
      character filename*128,filenameJ*128,filenamemu*128
      integer dtm(8)
      character(len=10) sys_time(3)
      integer,parameter :: max_line_len1 = 100
      character(max_line_len1) cf2, Evrs
      integer :: ist,isti
      double precision :: Vtpp1,Vtpp2,Vtpp3,Vtppn
      integer :: mucnt ! count for gyration to compute \mu
      
      call system_clock(t1)   ! record starting time 
      
!     In BB algorithm, variable timestep doesn't work, so normally keep its=0
      its=0 ! 0:fixed time step, 1:variable time step
      tini = 0.0d0
      
!     For variable time step option:
      tsdiv = 560.0d0 ! cyclotron time is divided with this number
      tmax = 1.0d-6 ! end of calculation

!     ##### Settings #####
      ikene = 2 ! 1: kene and direction ratio, 2: kpara and kperp
      iif=1; ! 1: save in ./orb/, 0: same folder
      ifile=1; ! 1: save orbit files
      icmp=1; ! 1: save extra data
      ihit=1; ! 1: stop calculation when hitting outer electrodes
      ieout=2; ! 1: warn when particle is outside of E-field region, 2: stop caculation when particle is outside E field region
      iclmap = 0; ! 1: save B data 0: don't save
      
!     ##### B coils currents, radii, (x,y,z) of coil centers
      icln = noco
      do i=1,noco
         ci(i) = Icl_in(i)
         rc(i) = rcl_in(i)
         xcl(i) = xcl_in(i)
         ycl(i) = ycl_in(i)
         zcl(i) = zcl_in(i)
      enddo
      
!     Infer number of electrodes from given file version (hard-coding to be revised)      
      if (Evers .lt. 4) then
         ielemax = 4; ! number of electrodes
      else
         ielemax = 8;
      endif
      
!     Pre-set voltages and frequencies of electrodes. (currently redundant, but useful for future upgrades) 
      Vdc(1) = 0.0d0; Vdc(2) = 0.0d0; Vdc(3) = 0.0d0; Vdc(4) = 0.0d0; Vdc(5) = 0.0d0;
      Vdc(6) = 0.0d0; Vdc(7) = 0.0d0; Vdc(8) = 0.0d0; Vdc(9) = 0.0d0; Vdc(10) = 0.0d0;
      
      Vac(1) = 0.0d0; Vac(2) = 0.0d0; Vac(3) = 0.0d0; Vac(4) = 0.0d0; Vac(5) = 0.0d0;
      Vac(6) = 0.0d0; Vac(7) = 0.0d0; Vac(8) = 0.0d0; Vac(9) = 0.0d0; Vac(10) = 0.0d0;
          
      call date_and_time(sys_time(1), sys_time(2), sys_time(3), dtm)

      if (iexb == 2) then !   B files to be used with interpolation. Not recommended. 
         open(17, file='b3d_hd.txt')
         read(17,*) xmin,xmax;
         read(17,*) ymin,ymax;
         read(17,*) zmin,zmax;
         read(17,*) nx,ny,nz;
         read(17,*) dx,dy,dz;
         close(17)

         do i=1,nx+1
            mt1x(i)=xmin+(i-1)*dx;
         end do
         do i=1,ny+1
            mt1y(i)=ymin+(i-1)*dy;
         end do
         do i=1,nz+1
            mt1z(i)=zmin+(i-1)*dz;
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

      if (iexe == 2) then !    E files to be used with interpolation (default)
         write (Evrs, '("./simion_Efields/ver.", i1.1,"/e3d_hd.txt")') Evers
         open(17, file=Evrs)
         read(17,*) xmin,xmax;
         read(17,*) ymin,ymax;
         read(17,*) zmin,zmax;
         read(17,*) nx,ny,nz;
         read(17,*) dx,dy,dz;
         close(17)

         ! write name of E-field header file
         !write(6,*) "Read E-field file header: ",trim(Evrs)

         do i=1,nx+1
            mt1x(i)=xmin+(i-1)*dx;
         end do
         do i=1,ny+1
            mt1y(i)=ymin+(i-1)*dy;
         end do
         do i=1,nz+1
            mt1z(i)=zmin+(i-1)*dz;
         end do

         do iex4=1,ielemax
            write (cf2, '("./simion_Efields/ver.", i1.1,"/e3d", i1.1, ".txt")') Evers,iex4
            ! write E-field file that was read:
            !write(6,*) "read E field file: ",trim(cf2)
            open(17, file=cf2) ! electrode-specific data
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

!     x^0, v^0, v^-1/2 is calculated by rk4

!     ################################################################################ start setup
      
!     ####### set for complete elliptic functions #######
!     Should be 1 only for the first calculation (produces .txt files in this case)
      itf=0 ! 0:load and save, or 1:calculate elliptic functions

!     ##### complete elliptic functions, start ####
      if (itf.eq.1) then
         !     Calculate K(k) and E(k) and save
      call efile2(elik2,elie2)
      OPEN (2,File="el2.txt")
      DO 101 i=1,99999
         write(2,'(1p1d12.6,1p2d25.17)') I/100000.0d0,elik2(i),elie2(i)
 101  CONTINUE
      CLOSE(2)
      endif

      if (itf.eq.0) then
!     Load saved file (faster) 
      OPEN (50,File="el2.txt")
      DO 102 I=1,99999
         READ(50,*) XI,ELIK2(I),ELIE2(I)
 102  CONTINUE
      CLOSE(50)
      endif

!      Uncomment to check that integrals are calculated correctly:
!      write(6,*) "Done. K(k=1) should be 1.5707963267948"
!      write(6,*) elik2(1)

      if(ELIK2(1)-1.57079d0.gt.0.01d0) then ! Warning when E() is empty.
         write(6,*) "##### !!Check eliptic functions!! Data is null!!! #####"
         write(6,*) ""
      endif

!    ##### complete elliptic functions, end ####

      if (iclmap == 1) then ! B calculation files are made 

!     If header ('hd') files are not read, xmax, xmin, etc. are still to be chosen at this point
!     Set the acceptable calculation region
      if (xmin .gt. -1.0d-6 .and. xmin .lt. 1.0d-6) then
         if (xmax .gt. -1.0d-6 .and. xmax .lt. 1.0d-6) then
            xmin=-0.25d0; xmax=0.25d0;
            ymin=-0.25d0; ymax=0.25d0;
            zmin=-0.25d0; zmax=0.25d0;
         endif
      endif

      write(6,*) "xmin,xmax,ymin,ymax:", xmin,xmax,ymin,ymax
      
      open (50,File="phib.txt") ! B at z=0 
            yi=0.0d0; zi=0.0d0;
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 
                  xi=(xmax/100.0d0)*iclbp;
                  call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz)
                  call bmap(elik2,elie2,xi,yi,zi,Psi)
                  write(50,'(1p3d18.10)') xi, Psi, dsqrt ( Bx*Bx+By*By+Bz*Bz )
                  end do
      close(50)

      open (50,File="ck1.txt") ! analytic solution of B at z=0, to check calculation
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 
               xi=0.0d0; yi=0.0d0;
               zi=zmin+( (zmax-zmin)/100.0d0 )*iclbp; 
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
                  write(50,'(1p7d18.10)') xi,zi,Psi,dsqrt ( Bx*Bx+By*By+Bz*Bz ),Bx,By,Bz
               end do
            end do
      close(50)

      open (50,File="ck2.txt") ! BS should be close to 2 pi r Atheta near axis
            ! plot "ck2.txt" using 1:2,"ck2.txt" using 1:3 w l
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 
               yi=0.0d0; zi=0.0d0;
               xi=( (xmax-xmin)/100.0d0 )*iclbp/10.0d0; 
               call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz); call bmap(elik2,elie2,xi,yi,zi,Psi);
               write(50,'(1p3d18.10)') xi, Bz*pi*xi*xi,Psi
            end do
      close(50)
      
      endif
!     ################################################################################ end setup

      ieoutc=0 ! constant to be used for warning outside E file region 
      
      if (ikene == 2) then !    parallel and perpendicular K(eV) ------
         x(1)=pmx; x(2)=pmy; x(3)=pmz;
         kini1=pmk1; kini2=pmk2; kang=pma;
      endif
      
      if (ikene == 1) then !   total K(eV) and velocity direction K(eV)
         x(1)=pmx; x(2)=pmy; x(3)=pmz;
         kini=pmk; xdr=pmrx; ydr=pmry; zdr=pmrz;
      endif
      
      !  RW amplitude and frequency 
      Vac(1)=pmrwa;Vac(2)=pmrwa;Vac(3)=pmrwa;Vac(4)=pmrwa;
      Fel=pmrwf;
      fgrad=pmrwfg;
      Phrw=pmrwd;

!     ##### initial values of E and B at x^0 ####

      if (iexb == 0) then ! calculate B by Biot-Savart
         call bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz);
         B(1)=Bx; B(2)=By; B(3)=Bz;
         write(6,*) "Directly computing B from Biot-Savart Law"
         ! print "(1p3d13.5,a35)",B(1),B(2),B(3),": Bx,By,Bz"
      endif
      if (iexb == 1) then ! give constant Bz=0.1
         Bx = 0.0d0;
         By = 0.0d0;
         Bz = 0.1d0;
         B(1)=Bx; B(2)=By; B(3)=Bz;
         write(6,*) "Simulation with constant B_z = 0.1"
      endif
      if (iexb == 2) then ! B is given by interpolation of external file data
         call bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz);
         print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with usual calculation"

         call bvac2(x(1),x(2),x(3),Bx,By,Bz);
         B(1)=Bx; B(2)=By; B(3)=Bz;
         print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with interpolation"
      endif

      
      if (iexe == 0) then ! set E=0 everywhere
         !call evac0(Ex,Ey,Ez);
         Ex = 0.0d0
         Ey = 0.0d0
         Ez = 0.0d0
         E(1)=Ex; E(2)=Ey; E(3)=Ez;
         write(6,*) "Simulation with fixed E=0"
         !print "(1p3d13.5,a35)",E(1),E(2),E(3),": E (V/m) set to 0"
      endif

      if (iexe == 1) then ! ideal azimuthal E-field
         call evac1(x(1),x(2),x(3),tini,tE1,mm,Ex,Ey,Ez); 
         E(1)=Ex; E(2)=Ey; E(3)=Ez;
         write(6,*) "Activated ideal azimuthal E-field"
         ! may give (correctly) all 0s depending on starting position:
         !print "(1p3d13.5,a35)",E(1),E(2),E(3),": ideally-azimuthal E (V/m) " 
      endif

      if (iexe == 2) then ! E is given by interpolation of external file data 
         call evac2(x(1),x(2),x(3),tini,tE1,Ex,Ey,Ez);
         E(1)=Ex; E(2)=Ey; E(3)=Ez;
         write(6,*) "Loading E-field data from external files"
         !print "(1p3d13.5,a35)",E(1),E(2),E(3),": E (V/m) with interpolation"
      endif

      if (iexe == 3) then ! locally-rotating E-field vortex
         call evac3(x(1),x(2),x(3),tini,tE1,mm,Ex,Ey,Ez); 
         E(1)=Ex; E(2)=Ey; E(3)=Ez;
         write(6,*) "Activated locally-rotating E-field vortex"
         !print "(1p3d13.5,a35)",E(1),E(2),E(3),": locally-rotating E vortex (V/m) "
      endif
      
!     ##### initial velocity ####

      if (ikene == 1) then
         ! total K(eV) and directions K(eV)
         ! xgm is the relativistic \gamma; xbt=\beta; xc=speed of light
         
         xgm=1.0d0+kini*1.95695119784913d-6 ! XGM=1.0D0+KENE*XE/(XME*XC*XC)
         xbt=dsqrt(1.0d0-1.0d0/(xgm*xgm)); xv=xc*xbt; xu=xgm*xv;
         u(1) = xu*xdr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))
         u(2) = xu*ydr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))
         u(3) = xu*zdr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))

         !print "(a35,1p3d12.4)","four-velocity (ux,uy,uz) (m):   ",u(1),u(2),u(3)
         !print "(a35,1p3d12.4)","real velocity (vx,vy,vz) (m):   ",u(1)/xgm,u(2)/xgm,u(3)/xgm
         !print "(a35,1p2d18.10)","beta and gamma:   ",xbt,xgm
         !print "(a35,1p2d18.10)","kinetic energy (eV):   ",kini
      
      endif

      if (ikene == 2) then
         ! parallel and perpendicular K(eV) 

         kini=kini1+kini2
         xgm=1.0d0+kini*1.95695119784913d-6 ! XGM=1.0D0+KENE*XE/(XME*XC*XC)
         xbt=dsqrt(1.0d0-1.0d0/(xgm*xgm)); xv=xc*xbt; xu=xgm*xv;
      
         xv1 = xv*dsqrt(kini1/kini); ! vpara
         xv2 = xv*dsqrt(kini2/kini); ! vperp
         !print "(a35,1p3d13.5)","kinetic energy (eV), pr,pp,ttl:   ",kini1, kini2, kini
         !print "(a35,1p3d13.5)","velocity (m/s):   ",xv1, xv2, xv
      
         bmb1 = dsqrt( B(1)*B(1) + B(2)*B(2) + B(3)*B(3) ); 
         bmb2 = dsqrt( B(1)*B(1) + B(3)*B(3) )

         v(1) = xv1 * B(1) / bmb1 + xv2 * ( B(3) / bmb2 * dcos(kang) - B(1)*B(2) / (bmb1*bmb2) * dsin(kang) )
         v(2) = xv1 * B(2) / bmb1 + xv2 * ( 0.0d0 + ( B(3)*B(3)+B(1)*B(1) ) / (bmb1*bmb2) * dsin(kang) )
         v(3) = xv1 * B(3) / bmb1 + xv2 * (-B(1) / bmb2 * dcos(kang) - B(2)*B(3) / (bmb1*bmb2) * dsin(kang) )

         ! vpara component, vperp-1 (cos) component, vperp-2 (sin) component
         ! write(6,*) xv1 * B(1) / bmb1, xv2 * ( B(3) / bmb2 * dcos(kang) ), - xv2 * B(1)*B(2) / (bmb1*bmb2) * dsin(kang) 
         ! write(6,*) xv1 * B(2) / bmb1, 0.0d0, xv2 * (  (B(3)*B(3)+B(1)*B(1)) / (bmb1*bmb2) * dsin(kang) )
         ! write(6,*) xv1 * B(3) / bmb1, xv2 * (-B(1) / bmb2 * dcos(kang)), -xv2 *  B(2)*B(3) / (bmb1*bmb2) * dsin(kang) 

         u(1) = v(1)*xgm; u(2) = v(2)*xgm; u(3) = v(3)*xgm;

         !print "(a35,1p3d13.5)","four-velocity (ux,uy,uz) (m):   ",u(1),u(2),u(3)
         !print "(a35,1p3d13.5)","real velocity (vx,vy,vz) (m):   ",v(1),v(2),v(3)
         !print "(a35,1p2d18.10)","beta and gamma:   ",xbt,xgm

         ! Larmor radius mv/qB, v=1.75882002359937d6,B=0.01T -> 1mm
         ! write(6,*) xm*v(2)/(xe*B(3))
         !print "(a35,1p2d18.10)","gyroradius (m):",xv2/bmb1/eom
         !print "(a35,1p2d18.10)","1/f_cyc (s):",2.0d0*pi/bmb1/eom
      
      endif

!    time step in RK to calculate x^0, v^-1/2
      if(its == 0) then ! constant time step
         Hr = dt;
      endif

      if(its == 1) then ! variable time step
         dt2 = ( 2.0d0*pi/(( dsqrt(Bx*Bx+By*By+Bz*Bz) )*eom) ) / tsdiv
         Hr = dt2;
      endif

!      print "(1p3d10.2,a3,1p3d10.2)",x(1),x(2),x(3)," ",u(1),u(2),u(3)

!     ##### 1/2 back by rk4 to get v^-1/2 rk4, start ####
      iini=10000; ddt = Hr/iini;
      Nr = 6;
      Hr = ddt;
      xr(1)=x(1); xr(2)=x(2); xr(3)=x(3);
      xr(4)=u(1); xr(5)=u(2); xr(6)=u(3); ! four-vector 
!      print "(1p1d10.2,a2,1p3d10.2,a3,1p3d10.2)",tini," ",xr(1),xr(2),xr(3)," ",xr(4),xr(5),xr(6)
      open (2,file="backr.txt")
        do it = 0,iini/2-1 ! back to ^-1/2 
        Tr = tini + it*ddt      
        call rk4(elik2,elie2,Nr,Tr,Xr,K1,K2,K3,K4,YWORK,Hr,iexe,tE1,mm,iexb)
        write(2,'(1p2d11.3,1p6d15.7)') tini-(it+1)*ddt,xr(1),xr(2),xr(3),xr(4),xr(5),xr(6)
      end do
      close(2)
      u(1)=xr(4); u(2)=xr(5); u(3)=xr(6); ! u^-1/2
!     ##### 1/2 back by rk4 to get v^-1/2 rk4, end ####

!     x(1),x(2),x(3): x at t=0
!     u(1),u(2),u(3): u at t=-1/2

!    'XU (four velocity), XV (real velocity)'
      xu = dsqrt( u(1)*u(1)+u(2)*u(2)+u(3)*u(3) )
      xv = xc*xu/dsqrt(xc*xc+xu*xu)
      xbt = xv/xc; xgm = 1.0d0/dsqrt(1.0d0-xbt*xbt) ! this is v^n-1/2

      v(1)=u(1)/xgm; v(2)=u(2)/xgm; v(3)=u(3)/xgm; 

!      x^n->ok for electromagnetic fields
!      v: v^n-1/2

      Hr = dt
      Tr = tini

!      write(6,*) " --- Process ID: ", pid
      if (ifile == 1) then ! orbit output files 
      
        if (iif .eq. 1) then
           write (filename, '("./orb/orb_", i5.5, "bbr.txt")') pid 
           write (filenameJ, '("./orb/jp_", i5.5, "bbr.txt")') pid 
           write (filenamemu, '("./orb/mu_", i5.5, "bbr.txt")') pid 
        else
           write (filename, '("orb_", i5.5, "bbr.txt")') pid 
           write (filenameJ, '("jp_", i5.5, "bbr.txt")') pid 
           write (filenamemu, '("mu_", i5.5, "bbr.txt")') pid  
        endif

        open (20,file=filename)
        open (22,file=filenameJ)
        open (24,file=filenamemu)
      
     endif
     
      ntr=0; ! counter for toroidal rotation 
      ntrpr=0;
      imu=0; ! saved number of mu
      xgpr1=0.0d0; ygpr1=0.0d0; zgpr1=0.0d0; xpr=0.0d0; ypr=0.0d0; zpr=0.0d0;
      irsapr=0; rPsipr=0.0d0 ! "previous" irsa index and Psi values. 
      mucnt=0  ! count for \mu computation
      
!     ######################################################################################
!     ##### B-B method pusher, start ####
      do it = 0,iend+1
         
!      time step for t^n,v^n-1/2 
      if(its == 0) then ! constant time step
         Hr = dt
         Tr = dt*it;
      endif

      if(its == 1) then ! variable time step
         if (Tr.gt.tmax) go to 10 ! stop calculation when t>tmax
         
         ! variable time step is adjusted to cyclotron frequency
         CALL bvac(elik2,elie2,xr(1),xr(2),xr(3),Bx,By,Bz)
         dt2 = ( 2.0d0*pi/(( dsqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3)) )*eom) ) / tsdiv
         Hr = dt2
      endif

      
      if(ihit == 1) then ! stop calculation when hitting electrodes 
      
        if( dsqrt( x(1)*x(1)+x(2)*x(2) ) .gt. 0.3d0 ) then
        write(6,*) " ****** positron hit the wall at r=0.3m"
        go to 10
        endif

        if( dabs( x(3) ) .gt. 0.3d0 ) then
        write(6,*) " ******* positron hit the wall at z=\pm0.2m"
        go to 10
        endif

      endif

      if(iexe == 2) then ! only when E-field files are used 
      if(ieout == 1) then ! warn when particle is outside of the E file region 
      if(ieoutc == 0) then ! stop if ieoutc is ever changed to 1 

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

        ! Use y-limits for z too (shortcut)
        if( x(3) .lt. ymin ) then
        write(6,*) "####### z < zmin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( x(3) .gt. ymax ) then
        write(6,*) "####### z > zmax of E file, outside of E field region, continue #######"; ieoutc=1
        endif

      endif
      endif

      if(ieout == 2) then ! stop when outside E-field region

        if( x(1) .lt. xmin ) then
           write(6,*) "####### outside of E field region, stop "; go to 10;
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

        ! Use y-limits for z too (shortcut)
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
!           XLML=XME*VPP*XGM/XE/(sqrt(BX*BX+BY*BY+BZ*BZ))

!           guiding center position, not accurate for large rL 
            xg1 = x(1) + (Vtpp2*Bz-Vtpp3*By)/Vtppn/Btt * gr
            yg1 = x(2) + (Vtpp3*Bx-Vtpp1*Bz)/Vtppn/Btt * gr
            zg1 = x(3) + (Vtpp1*By-Vtpp2*Bx)/Vtppn/Btt * gr

            call bmap(elik2,elie2,xg1,yg1,zg1,Psi); Psi2=Psi; ! \Psi at guiding center 
      
!           ########## mu ########## start
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
!           ########## mu ########## end

!           ########## J ##########
            ! Integration along B, data output at every z=0
            ! This will only be correct for the first value if the particle starts from z=0
            xv1int = xv1int + xv1*xv1*dt  !\int \vec{v} \cdot d\vec{l}

            if (zgpr1*zg1 .lt. 0.0d0) then ! detect z=0 cross
               write(22,'(1p3d15.7)') Tr,xv1int*2.0d0,Psi2
               xv1int=0.0D0 ! write and reset the integration value
            endif
!           ########## J ########## end

!           # toroidal rotation number 
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
!              1:x 2:y 3:z 4:r 5:vx 6:vy 7:vz 8:time
!              9:Kttl 10:Kpara 11:Kperp 12:Bx 13:By 14:Bz 15:B 16:Ex 17:Ey 18:Ez
!              19 :Psi 20:xg.c. 21:yg.c. 22:zg.c. 23:rg.c. 24:Psi_g.g. 25:mu,26:t.rotation
            endif
         endif
      endif
!     ## save data, end ##
!     ############################################################
!     Now that data is saved, push particle further

      if (iexb == 0) then
         CALL bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz)
      endif

      if (iexb == 1) then ! give constant Bz=0.1
         Bx = 0.0d0;
         By = 0.0d0;
         Bz = 0.1d0;
      endif
      
      if (iexb == 2) then !    B loading and interpolation
         CALL bvac2(x(1),x(2),x(3),Bx,By,Bz)
      endif
      
      B(1)=Bx; B(2)=By; B(3)=Bz;

      if (iexe == 0) then ! set E=0 everywhere
         Ex = 0.0d0
         Ey = 0.0d0
         Ez = 0.0d0
         !call evac0(x(1),x(2),x(3),it*dt,Ex,Ey,Ez)
      endif
      if (iexe == 1) then ! ideal azimuthal E-field
         call evac1(x(1),x(2),x(3),it*dt,tE1,mm,Ex,Ey,Ez)
      endif
      if (iexe == 2) then ! E is given by interpolation of external file data 
         call evac2(x(1),x(2),x(3),it*dt,tE1,Ex,Ey,Ez)
      endif
      if (iexe == 3) then ! locally-rotating E-field vortex
         call evac3(x(1),x(2),x(3),it*dt,tE1,mm,Ex,Ey,Ez)
      endif
      E(1)=Ex; E(2)=Ey; E(3)=Ez;
      
!      write(6,*) B(1),B(2),B(3),E(1),E(2),E(3)

!     U_minus
      do i = 1,3
        uminus(i) = u(i) + eom*E(i)*Hr*0.5d0
      end do
      
!     to calculate gamma^n, instead of gamma^n-1/2, use gamma^2 = 1+(u-/c)^2
      xgm = dsqrt ( 1 + ( uminus(1)*uminus(1)+uminus(2)*uminus(2)+uminus(3)*uminus(3) ) / (xc*xc) ) 
      
!     T vector
      do i = 1,3
        Tv(i) = eom*B(i)*Hr*0.5d0/xgm ! this \gamma should be ^0, not ^-1/2 
        
!        write(6,*) xgmd(it)
      end do
        
      Tsq = Tv(1)*Tv(1)+Tv(2)*Tv(2)+Tv(3)*Tv(3)

!     S vector
      do i = 1,3
        Sv(i)= 2.0d0*Tv(i)/(1.0d0+Tsq)
      end do

!     U_zero
      uzero(1) = uminus(1) + uminus(2)*Tv(3)-uminus(3)*Tv(2)
      uzero(2) = uminus(2) + uminus(3)*Tv(1)-uminus(1)*Tv(3)
      uzero(3) = uminus(3) + uminus(1)*Tv(2)-uminus(2)*Tv(1)

!     U_plus
      uplus(1) = uminus(1) + uzero(2)*Sv(3)-uzero(3)*Sv(2)
      uplus(2) = uminus(2) + uzero(3)*Sv(1)-uzero(1)*Sv(3)
      uplus(3) = uminus(3) + uzero(1)*Sv(2)-uzero(2)*Sv(1)

!     U^n+1/2
      do i = 1,3
        u(i) = uplus(i) + eom*E(i)*Hr*0.5d0
      end do

!      four-vector u(1),u(2),u(3) -> real velocity, beta, gamma. v^2=c^2 u^2 /(c^2+u^2),: ## ^n+1/2 ##
!     'XU (four velocity), XV (real velocity)'

      xu = dsqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
      xv = xc*xu/dsqrt(xc*xc+xu*xu)
      xbt = xv/xc; xgm = 1.0d0/dsqrt(1.0d0-xbt*xbt)
      v(1)=u(1)/xgm; v(2)=u(2)/xgm; v(3)=u(3)/xgm; 
      ken=(xgm-1.0d0)/1.95695119784913d-6;
 
!     X^n+1
        do i = 1,3
        x(i) = x(i) + u(i)*Hr/xgm
        end do

      if(its == 1) then ! variable time step
         Tr = Tr + Hr
      endif

     end do
     if (ifile == 1) then
        close(20)
        close(22)
        close(24)
        close(25)
     endif

!    Interrupted computation:
10   write(6,*) " ----  Ended computation"

     
      call system_clock(t2, t_rate, t_max)   ! save ending time
      if ( t2 < t1 ) then
         diff = (t_max - t1) + t2 + 1
      else
         diff = t2 - t1
      endif

      close(21)
      print "(A29,F10.3,F10.3)", "### calculation time (s) (h):", diff/dble(t_rate), diff/dble(t_rate)/3600.0d0

   return
   end 



!  ###################################################
   subroutine rk4(elik2,elie2,Nr,Tr,Xr,k1,k2,k3,k4,&
        ywork,Hr,iexe,tE1,mm,iexb)
!     just for initial n^-1/2 value
!     ***** to get ^-1/2 from ^0 *****
!     Runge Kutta
!
!  ###################################################
      implicit none
      double precision,dimension(1:100005) :: elik2,elie2
      double precision :: Hr,Tr,F
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      integer i,Nr
      double precision :: tE1
      integer :: mm,iexe,iexb

      do 10 i=1,Nr
  10  k1(i)=Hr*F(elik2,elie2,i,Tr,Xr,iexe,tE1,mm,iexb)
      do 11 i=1,Nr
  11  ywork(i)=Xr(i)+0.5d0*k1(i)
      do 12 i=1,Nr
  12  k2(i)=Hr*F(elik2,elie2,i,Tr+0.5d0*Hr,ywork,iexe,tE1,mm,iexb)
      do 13 i=1,Nr
  13  YWORK(i)=Xr(i)+0.5d0*k2(i)
      do 14 i=1,Nr
  14  k3(i)=Hr*F(elik2,elie2,i,Tr+0.5d0*Hr,ywork,iexe,tE1,mm,iexb)
      do 15 i=1,Nr
  15  ywork(i)=Xr(i)+K3(i)
      do 16 i=1,Nr
  16  k4(i)=Hr*F(elik2,elie2,i,Tr+Hr,ywork,iexe,tE1,mm,iexb)
      do 17 i=1,Nr
  17  Xr(i)=Xr(i)+(k1(I)+2.0d0*k2(i)+2.0d0*k3(i)+k4(i))/6.0d0
      Tr=Tr+Hr

      RETURN
      end subroutine rk4


!  ###################################################
      function F(elik2,elie2,i,t,x,iexe,tE1,mm,iexb)
!     ! -t ! for backward calculation
!
!  ###################################################
      use mod_variable
      implicit none
      double precision,dimension(0:100005) :: elik2,elie2
      double precision :: t,F,XQOM,xu,xv,xbt,xgm
      double precision,dimension(1:3) :: B,E
      double precision,dimension(1:6) :: x
      integer :: i  ! n seems redundant!!!!!!!!
      double precision :: xp,yp,zp,Ex,Ey,Ez,Bx,By,Bz
      double precision :: tE1
      integer :: mm, iexe, iexb

      xqom = 1.7588200236d11

      xp=x(1); yp=x(2); zp=x(3);

!      four-velocity X(4),X(5),X(6) -> real velocity, v^2=c^2 u^2 /(c^2+u^2)
      xu = dsqrt(x(4)*x(4)+x(5)*x(5)+x(6)*x(6))
      xv = XC*xu/dsqrt(xc*xc+xu*xu)
!      XV=XU/dsqrt(1.0+(XU/XC)*(XU/XC)) <-larger error 

      xbt=xv/xc
      xgm=1.0d0/dsqrt(1.0d0-xbt*xbt)
!      KENE=XME*XC*XC*(XGM-1.0)

      if (iexb == 0) then !    B by Biot-Savart 
         call bvac(elik2,elie2,XP,YP,ZP,Bx,By,Bz);
      endif
      if (iexb == 1) then ! give constant Bz=0.1
         Bx = 0.0d0;
         By = 0.0d0;
         Bz = 0.1d0;
      endif
      if (iexb == 2) then !    interpolate B in external files 
         call bvac2(XP,YP,ZP,Bx,By,Bz);
      endif
      

      B(1)=Bx; B(2)=By; B(3)=Bz;
      
      if (iexe == 0) then ! set E=0 everywhere
         !call evac0(XP,YP,ZP,T,Ex,Ey,Ez);
         Ex = 0.0d0
         Ey = 0.0d0
         Ez = 0.0d0
      endif
      if (iexe == 1) then ! ideal azimuthal E-field
         call evac1(XP,YP,ZP,T,tE1,mm,Ex,Ey,Ez);
      endif
      if (iexe == 2) then ! interpolate E in external files 
         call evac2(XP,YP,ZP,T,tE1,Ex,Ey,Ez);
      endif
      if (iexe == 3) then ! locally-rotating E-field vortex
         call evac3(XP,YP,ZP,T,tE1,mm,Ex,Ey,Ez); 
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

    
!     =================================================================
      subroutine bvac(elik2,elie2,xp,yp,zp,Bx,By,Bz)
!     calculate B at (x,y,z) 
!     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,Bx,By,Bz,Bro,Bzo,ctmp,Bx1,By1,Bz1
      double precision,dimension(1:100005) :: elik2,elie2
      integer i

!      ci(1) = 1.0d3 ! current
!      rc(1) = 1.0d-1 ! loop radius
!      zcl(1) = 0.0d0 ! loop z
!      xcl(1) = 0.0d0 ! loop offset in x direction from z axis
!      ycl(1) = 0.0d0 ! usually zero or middle

        Bx=0.0d0; By=0.0d0; Bz=0.0d0;

        do i=1,icln ! add for each of coils 
        
           ctmp = dsqrt((xp-xcl(i))*(xp-xcl(i))+(yp-ycl(i))*(yp-ycl(i)))
           call bloop(elik2,elie2,rc(i),ctmp,zp-zcl(i),ci(i),Bro,Bzo)

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

        end do

!      Homogeneous B in finite region
!      if ( xp .gt. 0.1d0 .and. xp .lt. 0.2d0 ) then
!      if ( yp .gt. 0.1d0 .and. yp .lt. 0.2d0 ) then
!      if ( zp .gt. 0.1d0 .and. zp .lt. 0.2d0 ) then
!        Bx = 0.0d0;
!        By = 0.0d0;
!        Bz = 0.1d0;
!      endif
!      endif
!      endif

      RETURN
      END



!     =================================================================
      subroutine bmap(elik2,elie2,xp,yp,zp,Psi)
!     return B at (x,y,z) obtained by interpolation xyz
!     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,ctmp,at,Psi,Psi1
      double precision,dimension(0:100005) :: elik2,elie2
      integer :: i

      Psi = 0.0d0

      do i=1,icln !  add fields of each coil

         ctmp = dsqrt((xp-xcl(i))*(xp-xcl(i))+(yp-ycl(i))*(yp-ycl(i)))        
         call rathe(elik2,elie2,rc(i),ctmp,zp-zcl(i),ci(i),at)
         Psi1 = 2.0d0*pi*at
         Psi = Psi + Psi1
      
      end do
      
      return
      end

!     =================================================================
      subroutine bvac2(xp,yp,zp,Bx,By,Bz);
!     return B at (x,y,z) obtained by interpolation 
!     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp
      double precision :: Bx,By,Bz
      double precision :: f111,f211,f121,f221,f112,f212,f122,f222,A1,A2,A3,A4,B1,B2
      
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

      
!     =================================================================
      subroutine evac0(Ex,Ey,Ez)
!     return homogeneous E=0 at all times and positions 
!     =================================================================
      implicit none
      double precision :: Ex,Ey,Ez
      
      ! Homogeneous E
      Ex = 0.0d0
      Ey = 0.0d0
      Ez = 0.0d0

      RETURN
      END


!     =================================================================
      subroutine evac1(xp,yp,zp,tp,tE1,mm,Ex,Ey,Ez)
!        
!     Return ideal azimuthal electric field at (x,y,z).
!     This subroutine gives the ideal electrostatic field of
!              Murakami, Sato & Hasegawa, 'Nonadiabatic behavior of the magnetic moment of a charged particle
!                               in a dipole magnetic field', Physics of Fluids B: Plasma Physics 2, 715 (1990)
!        
!     Inputs: (xp,yp,zp,tp) particle position and time;
!             tE1: time to interrupt E-field;
!             mm: azimuthal mode number
!     Outputs: (Ex,Ey,Ez): local electric field
!
!     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,tp,tE1,Ex,Ey,Ez
      double precision :: rr, ttheta
      double precision :: fchirp
      integer :: mm ! mode number for ideal azimuthal E-field

      if (tp .gt. tE1) then
         ! Zero out E field after specified time
         Ex=0.0d0; Ey=0.0d0; Ez=0.0d0;
         
      else
         ! allow RW frequency to increase linearly (set fgrad=0 to have fixed frequency)
         fchirp = fgrad * tp + Fel
         
         rr = dsqrt(xp*xp + yp*yp)
         ttheta = atan2(yp, xp) ! atan2 finds angle in correct quadrant
                 
         ! Assume that all electrodes are running at the same bias:
         Ex = - (real(mm)* Vac(1)/rr) * dsin(real(mm) * ttheta - fchirp * tp) * dsin(ttheta)
         Ey =  (real(mm) * Vac(1)/rr) * dsin(real(mm) * ttheta - fchirp * tp) * dcos(ttheta)
         Ez = 0.0d0
 
      endif

      RETURN
      END


!     =================================================================
      subroutine evac2(xp,yp,zp,tp,tE1,Ex,Ey,Ez)
!
!     Return E at (xp,yp,zp,tp) obtained by interpolation and superposition
!     of data from a finite number of electrodes   
!
!     Inputs: (xp,yp,zp,tp) particle position and time;
!             tE1: time to interrupt E-field;
!             mm: azimuthal mode number
!     Outputs: (Ex,Ey,Ez): local electric field
!    
!     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,tp,tE1,Ex,Ey,Ez
      double precision :: f111,f211,f121,f221,f112,f212,f122,f222,A1,A2,A3,A4,B1,B2
      !double precision :: RWamp,RWfreq,RWphase
      double precision :: fchirp
      !integer :: iex1,iex2,iex3

      if (tp .gt. tE1) then
         ! Zero out E field after specified time
         Ex=0.0d0; Ey=0.0d0; Ez=0.0d0;
         
      else

      n0x=int( (xp-xmin)/dx )+1; n0y=int( (yp-ymin)/dy )+1; n0z=int( (zp-ymin)/dz )+1;
      
      ! allow RW frequency to increase linearly (set fgrad=0 to have fixed frequency)
      fchirp = fgrad * tp + Fel
      !Vel(1) = Vdc(1) + Vac(1)*dsin( 2.0d0*pi*fchirp*tp + Phrw*0.0d0*pi/2.0d0 )
      !Vel(2) = Vdc(2) + Vac(2)*dsin( 2.0d0*pi*fchirp*tp + Phrw*1.0d0*pi/2.0d0 )
      !Vel(3) = Vdc(3) + Vac(3)*dsin( 2.0d0*pi*fchirp*tp + Phrw*2.0d0*pi/2.0d0 )
      !Vel(4) = Vdc(4) + Vac(4)*dsin( 2.0d0*pi*fchirp*tp + Phrw*3.0d0*pi/2.0d0 )

!     Get electric potential at chosen position from all electrodes (consider both DC and AC potentials) 
      do iele = 1,ielemax
         Vel(iele) = Vdc(iele) + Vac(iele)*dsin( 2.0d0*pi*fchirp*tp + Phrw*real(iele-1)*pi/2.0d0 )
      end do
   
!     Ex by linear interpolation 
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

!     Ey by linear interpolation
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

!     Ez by linear interpolation
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


!     =================================================================
      subroutine evac3(xp,yp,zp,tp,tE1,mm,Ex,Ey,Ez)
!        
!     return locally-rotating E-field at (x,y,z)
!     NB: The field from this subroutine is different from the one in Murakami et al. 1990.
!         It seems to produce interesting effects, but is not physically motivated.        
!        
!     Inputs: (xp,yp,zp,tp) particle position and time;
!             tE1: time to interrupt E-field;
!             mm: azimuthal mode number
!     Outputs: (Ex,Ey,Ez): local electric field
!
!     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,tp,tE1,Ex,Ey,Ez
      double precision :: rr, ttheta
      double precision :: fchirp
      integer :: mm ! mode number for ideal azimuthal E-field

      if (tp .gt. tE1) then
         ! Zero out E field after specified time
         Ex=0.0d0; Ey=0.0d0; Ez=0.0d0;
         
      else
         ! allow RW frequency to increase linearly (set fgrad=0 to have fixed frequency)
         fchirp = fgrad * tp + Fel
         
         rr = dsqrt(xp*xp + yp*yp)
         ttheta = atan2(yp, xp) ! atan2 finds angle in correct quadrant
                 
         ! Assume that all electrodes are running at the same bias:
         Ex = - (real(mm)* Vac(1)/rr) * dsin(real(mm) * ttheta - fchirp * tp)
         Ey =  (real(mm) * Vac(1)/rr) * dcos(real(mm) * ttheta - fchirp * tp)
         Ez = 0.0d0
 
      endif

      RETURN
      END


!      =================================================================
      SUBROUTINE BLOOP(elik2,elie2,RC,R,Z,CI,BR,BZ)
!      calculate the Br and Bz produced by a loop current
!      RC:loop radius R,Z:position CI:coil current 
!      =================================================================
!                                                    Biot-Savalt formula
!                                                    ===================
      implicit none
      double precision, dimension(0:100005) :: elik2,elie2
      DOUBLE precision :: ZK,ZZK,RC,R,Z,CI,BR,BZ,G1,FK,FE,A,G,FF,E,H,xcc,fk1,fe1,fk2,fe2
      integer :: i,i2
      
      xcc=2.0d-7

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

!      FK FE interpolation start
      i2=i+1;
      fk1=elik2(i); fe1=elie2(i);
      fk2=elik2(i2); fe2=elie2(i2);
!      write(6,'(d13.5,i7,1p3d17.9)') zzk,i,fk1,fk2,fk1+(zzk-dble(i)/1.0d5)*(fk2-fk1)
      fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
      fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);
!      write(21,'(i5,1p7d15.7)') Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
!      FK FE interpolation end

!      FK FE
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
!
!*     ZK is k^2 
      ZK=4.0d0*rc*x/((rc+x)*(rc+x)+z*z)
      ZZK=dsqrt(ZK)
      
      IF(ZK.GE.0.999999d0) GO TO 20
      IF(ZK.GT.0.9998d0) GO TO 12
      
      A0=2.0d0*XCC*CI*x/dsqrt(ZK)
      A1=dsqrt(RC/x)
      A2=1.0d0-0.5d0*ZK
      I=IDINT(ZZK*100000.0d0)
      
!     FK FE interpolation start
      i2=i+1;
      fk1=elik2(i); fe1=elie2(i);
      fk2=elik2(i2); fe2=elie2(i2);
!      write(6,'(d13.5,i7,1p3d17.9)') zzk,i,fk1,fk2,fk1+(zzk-dble(i)/1.0d5)*(fk2-fk1)
      fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
      fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);
!      write(21,'(i5,1p7d15.7)') Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
!     FK FE interpolation end

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

!       show progress
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
        double precision, intent(in) :: k  
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


