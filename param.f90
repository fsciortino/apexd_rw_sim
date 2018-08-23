!C    =================================================================
      program param
!C    generate parameter file
!     initial x,y,z,Epara,Eperp,angle values
!C    =================================================================
      implicit none
      double precision,dimension(1:3) :: x,u,v,uminus,uzero,uplus
      double precision,dimension(1:3) :: B,E,Tv,Sv
!      double precision,dimension(0:10000000) :: xd,yd,zd,vxd,vyd,vzd,uxd,uyd,uzd,ke,td,xgmd
      double precision :: tini,t,dt,Tsq,ddt,tt,xu,xv, kini,xdr,ydr,zdr
      double precision :: pi,xc,xe,xm,eom, xi,yi,zi,vxi,vyi,vzi, Ex,Ey,Ez,Bx,By,Bz
      integer :: it,iend,i,iini,imbk,Nr,ir,ir2,itf
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      double precision :: Hr,Tr, kene,xgm,xbt, rnd
      double precision,dimension(0:10000005) :: ene,epara,eperp,eparai,eperpi
      double precision,dimension(0:100005) :: pepara,peperp,rd1,rd2
      double precision,dimension(0:1005) :: cke,ck1,ck2 ! frequency distribution
      double precision :: Tfl,Tpara,Tperp,epst
      integer :: ip2,npt,ick,ipt,ipl,ird
      double precision,dimension(0:10000005) :: br,bwr,bwri
      double precision,dimension(0:100005) :: rd3,rd4,prds,prds1,prds2,rd5
      double precision,dimension(0:1005) :: ckr,ck3 ! frequency distribution
      double precision :: rpst,rsgm,pinix,piniy,piniz,Bro,Bzo
      double precision,dimension(0:100005) :: elik2,elie2
      integer(8) :: npst,ip
      double precision :: Vrw,frw,Drw

!     KENE,XC,XE,XGM,XBT,XV,XDRC,YDRC,ZDRC
!C    ##### physical parameters, start ####
!      pi = 3.14159265358979323846d0
      pi = 3.141592653589793238462643383279d0
      xc = 2.99792458d8
      xe = 1.6021766208d-19
      xm = 0.910938356d-30
      eom = 1.7588200236d11
!C    ##### physical parameters, end ####

      ird = 1 ! 0: simple (without random data) 1: with random data

!C    #############################################################
!C    #############################################################
      if (ird==0) then ! simple (without random data)

!     ### save data in a file

      npt=10

      open (2,file="param.txt")
        do ip = 1,npt
        Vrw = 5.0d0
        frw = 50.0d3 + ip*10.0d3
        Drw = 1.0d0 ! RW direction, 1 or -1
        write(2,'(1p6d16.8,1p3d13.5)') 0.2d0,0.0d0,0.0d0,5.0d0,1.0d0,pi/4.0d0,Vrw,frw,Drw
        ipl=ipl+1
        end do
      close(2)
      
      endif
!C    #############################################################
!C    #############################################################







!C    #############################################################
!C    #############################################################
      if (ird==1) then ! with random data
      
      Vrw = 5.0d0 ! rw amplitude
      frw = 100.0d3 ! rw frequency
      Drw = 1 ! rw direction
      
      ir=0 ! 1: save many files 0: no
      ir2=0 ! save large filesÅids1.txt, ds11.txtÅj

!C      #############################################################
!C    ##### make parameter file #####

      write(6,*)"########### parameter calculation mode ###########"  

      ipl = 1; ! particle (line) number
      do ipt=1,1
      
!     ### kinetic energies, particle number
      Tfl=5.0d0; Tpara=1.0d0; Tperp=1.0d0; npt=100;

!     ### spatial spread and center position
      rsgm=1.0d-2; pinix=0.2d0; piniy=0.0d0; piniz=0.0d0;

!C    ##### velocity space, start #####

!     ### Energy step and span used for mapping
!     npst*epst gives maximum energy. This value should be large enough.
      npst=1000000_8; epst=0.00001d0
      write(6,'(a30,1p2d20.12)') "energy range is (eV)",npst*epst

!     ### calculate the values of distribution function at E
      do ip = 0,npst
        ene(ip) = ip*epst;
        epara(ip) = dexp( -(ene(ip)-Tfl)**2/(2.0d0*Tpara**2) )
        eperp(ip) = dexp( -ene(ip)/Tperp )
      end do
      
!     ### error function by integration, normalized for [0:1]
      eparai(0)=0.0d0; eperpi(0)=0.0d0;
      do ip = 0,npst-1
        eparai(ip+1) = eparai(ip) + epara(ip)
        eperpi(ip+1) = eperpi(ip) + eperp(ip)
      end do
      do ip = 0,npst
        eparai(ip) = eparai(ip) / eparai(npst)
        eperpi(ip) = eperpi(ip) / eperpi(npst)
      end do

!     ### write file; energy distribution and error function
      if(ir2.eq.1) then
      open (2,file="ds1.txt")
        do ip = 0,npst
        write(2,'(1p5d15.7)') ene(ip),epara(ip),eperp(ip),eparai(ip),eperpi(ip)
        end do
      close(2)
      endif

!     ### random function to be used for mapping. ok without seed.
      do ip = 1,npt
      call random_number(rnd)
      rd1(ip) = rnd
      call random_number(rnd)
      rd2(ip) = rnd
      end do

!     ### write file; random function to be used for mapping
      if(ir2.eq.1) then
      open (2,file="ds2.txt")
        do ip = 1,npt
        write(2,'(1p2d15.7)') rd1(ip),rd2(ip)
        end do
      close(2)
      endif
      
!     ### mapping according to the error function
      do ip = 1,npt
        do ip2 = 0,npst-1
          if ( eparai(ip2) .le. rd1(ip) .and. rd1(ip) .lt. eparai(ip2+1) ) then
            pepara(ip) = ene(ip2)
!                write(6,*) rd1(ip),eparai(ip2),eparai(ip2+1)
          endif
          if ( eperpi(ip2) .le. rd2(ip) .and. rd2(ip) .lt. eperpi(ip2+1) ) then
            peperp(ip) = ene(ip2)
          endif
        end do
      end do

!     ### write file; kinetic energy of particles
      if(ir.eq.1) then
      open (2,file="ds3.txt")
        do ip = 1,npt
        write(2,'(1p2d15.7)') pepara(ip),peperp(ip)
        end do
      close(2)
      endif
      
!     ### check calculation by using frequency distribution
      do ip = 0,20
      cke(ip)=ip*0.5d0;
      ck1(ip)=0.0d0; ck2(ip)=0.0d0;
      end do

      ick=0
      do ip = 1,npt
!      write(6,*) pepara(ip)
        do ip2 = 0,19
          if ( cke(ip2) .le. pepara(ip) .and. pepara(ip) .lt. cke(ip2+1) ) then
            ck1(ip2) = ck1(ip2) + 1.0d0
            ick = ick + 1
!                write(6,*) pepara(ip),ip2,ck1(ip)
          endif
          if ( cke(ip2) .le. peperp(ip) .and. peperp(ip) .lt. cke(ip2+1) ) then
            ck2(ip2) = ck2(ip2) + 1.0d0
          endif
        end do
      end do

      if (ick.eq.npt) then
      write(6,*) "checked that particle number is",ick,"."
      endif
      
!     ### save file; frequency distribution
      if(ir.eq.1) then
      open (2,file="ds4.txt")
        do ip = 0,19
        write(2,'(1p3d15.7)') cke(ip),ck1(ip),ck2(ip)
        end do
      close(2)
      endif
      
!C    ##### velocity space, end #####



!C    ##### real space, start #####

!     ### r step and span used for mapping
!     npst*epst gives maximum r. This value should be large enough.
      npst=1000000_8; rpst=0.00001d0

!     ### distribution function at each step
      do ip = 0,npst
        br(ip) = ip*rpst;
        bwr(ip) = dexp ( -2.0d0 * br(ip)*br(ip) / rsgm**2 )
      end do
      
!     ### error function by integration, normalized for [0:1]
      bwri(0)=0.0d0;
      do ip = 0,npst-1
        bwri(ip+1) = bwri(ip) + bwr(ip)
      end do
      do ip = 0,npst
        bwri(ip) = bwri(ip) / bwri(npst)
      end do

!     ### write file; spatial distribution and error function
      if(ir2.eq.1) then
      open (2,file="ds11.txt")
        do ip = 0,npst
        write(2,'(1p3d15.7)') br(ip),bwr(ip),bwri(ip)
        end do
      close(2)
      endif
      
!     ### random function to be used for mapping. ok without seed.
      do ip = 1,npt
      call random_number(rnd)
      rd3(ip) = rnd
      call random_number(rnd)
      rd4(ip) = rnd
      end do

!     ### write file; random function to be used for mapping
      if(ir.eq.1) then
      open (2,file="ds12.txt")
        do ip = 1,npt
        write(2,'(1p1d15.7)') rd3(ip)
        end do
      close(2)
      endif
      
!     ### mapping according to the error function
      do ip = 1,npt
        do ip2 = 0,npst-1
          if ( bwri(ip2) .le. rd3(ip) .and. rd3(ip) .lt. bwri(ip2+1) ) then
            prds(ip) = br(ip2)
!                write(6,*) rd3(ip),bwri(ip2)
          endif
        end do
      end do

      do ip = 1,npt
      prds1(ip) = prds(ip)*dcos(rd4(ip)*2.0d0*pi)
      prds2(ip) = prds(ip)*dsin(rd4(ip)*2.0d0*pi)
      end do

!     ### save file; positions of particles, reflecting the spatial profile
      if(ir.eq.1) then
      open (2,file="ds13.txt")
        do ip = 1,npt
        write(2,'(1p3d15.7)') prds(ip),prds1(ip),prds2(ip)
        end do
      close(2)
      endif
      
!     ### check calculation by using frequency distribution
      do ip = 0,20
      ckr(ip)=ip*0.001d0;
      ck3(ip)=0.0d0;
      end do

      do ip = 1,npt
!      write(6,*) pepara(ip)
        do ip2 = 0,19
          if ( ckr(ip2) .le. prds(ip) .and. prds(ip) .lt. ckr(ip2+1) ) then
            ck3(ip2) = ck3(ip2) + 1.0d0
!                write(6,*) pepara(ip),ip2,ck1(ip)
          endif
        end do
      end do

!     ### write file; frequency distribution
      if(ir.eq.1) then
      open (2,file="ds14.txt")
        do ip = 0,19
        write(2,'(1p2d15.7)') ckr(ip),ck3(ip)
        end do
      close(2)
      endif

      print "(a20,1p3d12.4)","## parameters ##"
      print "(a35,i10)","particle number:   ",npt
      print "(a35,1p3d12.4)","Tfl,Tpara,Tperp (eV):   ",Tfl,Tpara,Tperp
!      print "(a20,1p3d12.4)","## parameters ##"
!      print "(a35,i10)","particle number:   ",npt
      print "(a35,1p3d12.4)","r_sigma (m):   ",rsgm
      print "(a40)","Output file of Tz and T para: ds3.txt"
      print "(a40)","Output file of r,x,y: ds13.txt"
      print "(a70)","f(E) and f(r) can be checked with ds4.txt and ds14.txt:"
      print "(a70)","plot 'ds4.txt' using 1:2 w l,'ds4.txt' using 1:3 w l : Tpara and Tperp"
      print "(a70)","plot 'ds14.txt' using 1:2 w l"
      print "(a70)","plot 'ds13.txt' using 2:3 w p"
!      print "(a50)","Distribution functions can be checked with ds14.txt"

!      write(6,*) dcos(60.0d0*pi/180.0d0)
      end do

!C    ##### real space, end #####

!     ### randum function for mapping
      do ip = 1,npt
      call random_number(rnd)
      rd5(ip) = rnd * 2.0d0 * pi
      end do
      

!     ### save file; parameter file
      open (2,file="param.txt")
        do ip = 1,npt
        write(2,'(1p9d16.8)') pinix+prds1(ip),piniy,piniz+prds2(ip),pepara(ip),peperp(ip),rd5(ip),Vrw,frw,Drw
        ipl=ipl+1
        end do
      close(2)

!     initial (x,y,z)ÅCKpara, Kperp, Kperp direction
!     RW amplitude, frequency, direction

      endif
!C    #############################################################
!C    #############################################################

      end


      
      
      
      
      
