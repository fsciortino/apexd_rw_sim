!C    =================================================================
      module mod_variable
!C    In this version, electric field files from simion can be used.
!C    Calculation flow is,
!C     1. make electric field files by simion, and run efile.f90 (see comments in efile.f90).
!C     2. copy e3d_hd.txt(grid information) and e3d*.txt into the same folder as rkr.f90.
!C     3. compile and run rkr.90 with param.txt, which can be made by param.f90.
!C     4. plot data or statistical treatment should be done.
!     
!     Output files:
!     orb****rk.txt:
!     1:x 2:y 3:z 4:r 5:vx 6:vy 7:vz 8:time
!     9:Kttl 10:Kpara 11:Kperp 12:Bx 13:By 14:Bz 15:B 16:Ex 17:Ey 18:Ez 19:Psi
!     20:xg.c. 21:yg.c. 22:zg.c. 23:rg.c. 24:Psi_g.g. 25:mu 26:t.rotation
!     splot "orb0001rk.txt" w l
!
!     mu:24mu,J:22jp,Psi:(23)
!
!     orb_a0rk.txt,orb_a1rk.txt,orb_a2rk.txt: averatge, average+sigma, average-sigma
!     1:iden 2:ipm 3:time 4:Psi 5:r 6:x 7:y 8:z 9:K 10:K1 11:K2   g.c.->12:Psi 13:r 14:x 15:y 16:z
!     1:x 2:y 3:z 4:r 5:gx 6:gy 7:gz 8:gr 9:K 10:Kpr 11:Kpp 12:Ps 13:gPs 14:mu 15:time
!
!C    =================================================================
      implicit none
!C    ##### physical parameters, start ####
      double precision, parameter :: pi = 3.14159265358979323846d0;
      double precision, parameter :: xc = 2.99792458d8
      double precision, parameter :: xe = 1.6021766208d-19;
      double precision, parameter :: xm = 0.910938356d-30;
      double precision, parameter :: eom = 1.7588200236d11;
      
      double precision :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz
      integer :: nx,ny,nz,n0x,n0y,n0z,iexb,iexe,iele,ielemax,ielep,icln
      double precision,dimension(1:205) :: mt1x,mt1y,mt1z
      double precision,dimension(1:1,1:1,1:1) :: Bxm,Bym,Bzm ! to save memory for those usually not to be used 磁場配列は通常使わないので，メモリ節約
      double precision,dimension(1:205,1:205,1:205,1:8) :: Exm,Eym,Ezm ! electric fields (unit field by electrodes) max9 電場ファイル (電極による単位電場) max9
      double precision,dimension(1:12) :: Vel,Vdc,Vac ! electrode voltage 電極電位
      double precision :: Fel,Phrw ! frequency and phase of RW 回転電場（等）の周波数,位相
      double precision :: vfp1,vfp2,vfp3,vfp4,vfp5,vfp6 ! for param.txt の最後の6個は，グローバル変数として与えられる．
      double precision,dimension(1:52) :: ci,rc,zcl,xcl,ycl ! Coil current, loop radius, z, x, y(y is usulally 0) コイル電流，ループ半径，z位置，x位置，y位置（通常ゼロ）

      end module mod_variable

!C    =================================================================
      program rkr
!C    Integrator with relativistic Runge-Kutta
!C    *ring coil current*  *variable time step*
!C    ##### use param.txt in this version ####
!C    =================================================================
      use mod_variable
      implicit none
      double precision,dimension(1:3) :: x,v,u,vminus,vzero,vplus
      double precision,dimension(1:3) :: B,E,Tv,Sv
      double precision :: tini,t,dt,Tsq,ddt,tt,xv,xvs,xgm,xbt,xu,ken,tmax,tsdiv,dt2,xpr,ypr,zpr,rsa,rPsi,rPsipr
      double precision :: kini,xdr,ydr,zdr,kini1,kini2,xv1,xv2,xv1int,kang,bmb1,bmb2,bmb3
      double precision :: xi,yi,zi,vxi,vyi,vzi,ken1,ken2,xmu
      double precision :: Ex,Ey,Ez,Bx,By,Bz,Btt,Psi,Psi1,Psi2,xg1,yg1,zg1,xgpr1,ygpr1,zgpr1
      integer :: i,iini,imbk,Nr,itf,its,ikene,iparam,ipm,ipm1,ipm2,ipmt,ifile,ifile2,imbk2,imbk2t,iclmap,iclbp,iclbp2,iif,ntr,ntrpr
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      double precision :: Hr,Tr,gr,gr1
      integer(8) :: it,iend,iend2
      integer t1, t2, t_rate, t_max, diff,iex1,iex2,iex3,iex4,icmp,ihit,izero,ieout,ieoutc,iden,imu,irsa,irsapr
      double precision,dimension(0:100005) :: elik2,elie2
      double precision,dimension(0:1001) :: pmx,pmy,pmz,pmk1,pmk2,pma,pmk,pmrx,pmry,pmrz
      double precision,dimension(0:1001) :: pmrwa,pmrwf,pmrwd
      double precision,dimension(0:1001) :: fpm1,fpm2,fpm3,fpm4,fpm5,fpm6 ! for param.txt can be used as vfp1,vfp2,vfp3,vfp4,vfp5,vfp6 as global parameters
      character filename*128,filename0*128,filenameJ*128,filenameP*128,filenamemu*128
      integer dtm(8)
      character(len=10) sys_time(3)
      integer,parameter :: max_line_len1 = 100
      character(max_line_len1) cf1,cf2
!      double precision,dimension(1:105,1:1005) :: nit,niPsi,nir,nix,niy,niz,niK,niK1,niK2
!      double precision,dimension(1:105,1:1005) :: nigPsi,nigr,nigx,nigy,nigz ! for density calculation
      double precision,dimension(1:105,1:10005) :: stt,stx,sty,stz,str,stgx,stgy,stgz,stgr,stK,stKpr,stKpp,stps1,stps2,stmu
      integer, dimension(1:105,1:10005) :: strt ! time step & particle number（to be averaged with the latter）これらは時間番号，粒子番号（後者で平均を取る）
      double precision,dimension(1:105) :: stta,stxa,stya,stza,stra,stgxa,stgya,stgza,stgra,stKa,stKpra,stKppa,stps1a,stps2a,stmua
      double precision,dimension(1:105) :: stts,stxs,stys,stzs,strs,stgxs,stgys,stgzs,stgrs,stKs,stKprs,stKpps,stps1s,stps2s,stmus
      integer, dimension(1:105) :: strta,strts ! averaged value and sigma (for particles) 粒子番号（後者）での平均値と標準偏差
      integer :: ist,isti
      double precision :: Vtpp1,Vtpp2,Vtpp3,Vtppn
      call system_clock(t1)   ! record starting time 開始時を記録

!C    ##### line numbers in param.txt to be used in the calculation #### パラメータファイルの行数指定
      ipm1=1; ipm2=2; ! starting and ending line of param.txt
      ipmt = ipm2-ipm1+1 ! total particle number, maximum value of ipm トータル粒子数，ipmの最大値
      
!C    ##### time step, step number ####
      its=0 ! 0:fixed time step, 1:variable time step
!C    ### for its=0 ###
      tini = 0.0d0
      dt = 1.0d-11
      iend = 1000000
      imbk = 100 ! data is saved to a file 1 times per imbk ファイル記録用の間引き
      imbk2= 10000 ! same as above, but used for statistical calculation for many particles 複数粒子の統計処理用の間引き
      imbk2t = iend/imbk2+1
!      dt = 1.0d-11; iend = 10000000; imbk = 1000; imbk2= 10000;
!C    ### for its=1 ###
      tsdiv = 56.0d0 ! cyclotron time is divided with this number
      tmax = 1.0d-6 ! end of calculation
      iend2 = 10000000000_8 ! just a large number. calculation stops when tmax<t or iend2<it
!C    ##### just for a test, to apply v=0 as initial value 初期速度を0で与えるテスト用 ####
      izero=0; ! 1: start from x=0 with v=0
!C    ##### how to apply energy. エネルギーの与え方 2 is convenient to give random energy ####
      ikene = 2 ! 1: kene and direction ratio, 2: kpara and kperp
!C    ##### save files in an independent folder 独立したフォルダにファイルを保存するか ####
      iif=0; ! 1: save in ./orb/, 0: same folder
!C    ##### save obt files for each of the particles obtファイルを保存するか ####
      ifile=1; ! 1: save orbit files
!C    ##### save additional files 追加的な計算も保存 ####
      icmp=1; ! 1:
!C    ##### save files for statitics 密度計算用のファイルを保存するか ####
      ifile2=1; ! 1: save orbit files
!C    ##### stop calculation by hitting electrodes 障害物に衝突で計算終了 ####
      ihit=1; ! 1: stop caculation when hitting something
      ci(1) = 1.0d3; rc(1) = 1.0d-1; zcl(1) = 0.0d0; xcl(1) = 0.0d0; ycl(1) = 0.0d0;
      ci(2) = 1.0d3; rc(2) = 1.0d-1; zcl(2) = 0.1d0; xcl(2) = 0.0d0; ycl(2) = 0.0d0;
      ci(3) = 1.0d3; rc(3) = 1.0d-1; zcl(3) = 0.2d0; xcl(3) = 0.0d0; ycl(3) = 0.0d0;
      ci(4) = 1.0d3; rc(4) = 1.0d-1; zcl(4) = 0.3d0; xcl(4) = 0.0d0; ycl(4) = 0.0d0;
      ci(5) = 1.0d3; rc(5) = 1.0d-1; zcl(5) = 0.4d0; xcl(5) = 0.0d0; ycl(5) = 0.0d0;
      ci(6) = 1.0d3; rc(6) = 1.0d-1; zcl(6) = -0.1d0; xcl(6) = 0.0d0; ycl(6) = 0.0d0;
      ci(7) = 1.0d3; rc(7) = 1.0d-1; zcl(7) = -0.2d0; xcl(7) = 0.0d0; ycl(7) = 0.0d0;
      ci(8) = 1.0d3; rc(8) = 1.0d-1; zcl(8) = -0.3d0; xcl(8) = 0.0d0; ycl(8) = 0.0d0;
      ci(9) = 1.0d3; rc(9) = 1.0d-1; zcl(9) = -0.4d0; xcl(9) = 0.0d0; ycl(9) = 0.0d0;
      ci(10) = 1.0d3; rc(10) = 1.0d-1; zcl(10) = -0.5d0; xcl(10) = 0.0d0; ycl(10) = 0.0d0;
      icln = 1; ! coil number コイルの数 1からiclnまでが使われる
      iclmap = 1; ! 1:save B data 0:no 磁場の計算結果を保存 0:しない
!C    ##### B is geven by linear interpolation (not recommended) #### need to modify arrays in modules
      iexb=0; ! 0:Biot-Savart 1: use external B field file 磁場データを補完して使用 あまり使わない，使う場合はmoduleで磁場の配列を変更
!C    ##### E is geven by linear interpolation 電場データを補完して使用 ####
      iexe=0; ! 0:given as constant values  1: use external E field file
      ielemax = 4; ! number of electrodes 電極の数
!C    voltages and frequencies of electrodes. 各電極の電位と周波数．パラメータファイルを使う場合，ielep=1として調整．
      ielep=1; ! when 1, param.txt is used for RW parameters in addition to the followings. 1の場合，param.txtから読むので，下記だけでは設定しない．
      Vdc(1) = 0.0d0; Vdc(2) = 0.0d0; Vdc(3) = 0.0d0; Vdc(4) = 0.0d0; Vdc(5) = 0.0d0;
      Vdc(6) = 0.0d0; Vdc(7) = 0.0d0; Vdc(8) = 0.0d0; Vdc(9) = 0.0d0; Vdc(10) = 0.0d0;
      Vac(1) = 5.0d0; Vac(2) = 5.0d0; Vac(3) = 5.0d0; Vac(4) = 5.0d0; Vac(5) = 0.0d0;
      Vac(6) = 0.0d0; Vac(7) = 0.0d0; Vac(8) = 0.0d0; Vac(9) = 0.0d0; Vac(10) = 0.0d0;
      Fel = 1.0d5;
!C    ##### warning/stop calculation when the particle is outside of the E file region 電場ファイルの領域を出た場合に警告/計算停止 ####
      ieout=2; ! 1:warning, 2: stop caculation when particle is outside E

      call date_and_time(sys_time(1), sys_time(2), sys_time(3), dtm)

      if (iif .eq. 1) then
      write (filename0, '("./orb/log",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"rkr.txt")') dtm(1),dtm(2),dtm(3),dtm(5),dtm(6),dtm(7)
      else
      write (filename0, '("log",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"rkr.txt")') dtm(1),dtm(2),dtm(3),dtm(5),dtm(6),dtm(7)
      endif
      write(6,*) filename0
      open (21,file=filename0)

      write(21,'(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)') dtm(1),dtm(2),dtm(3),dtm(5),dtm(6),dtm(7) ! starting date and time 開始日時
      write(21,'(i1.1)') its ! constant or variable time step 時間ステップの取り方
      if (its==0) then
      write(21,'(1p1d11.3,i15,1p1d11.3,i10)') dt,iend,dt*iend,imbk
      endif
      if (its==1) then
      write(21,'(1p2d15.7)') tsdiv,tmax
      endif


      if (its == 0) then
      print "(a35,1p1d11.3,i15,1p1d11.3)","## fixed time step ## dt/iend/T(s):",dt,iend,dt*iend
      endif
      if (its == 1) then
      print "(a39,1p2d11.3)","## variable time step ## tsdiv tmax(s):",tsdiv,tmax
      print "(a65,1p2d11.3)","-- time step is (cyclotron period)/tsdiv at aech point --"
      endif

      if (ikene == 2) then !C parallel and perpendicular K(eV) 縦と横のK(eV)で与える
!C    param.txt(x,y,z,Kpara,Kperp,angle) 外部ファイルで与えるパラメータを読む．通常こちら
!C     use external file to read parameters
        open (50,File="param.txt")
          do i=1,ipm2
          read(50,*) pmx(i),pmy(i),pmz(i),pmk1(i),pmk2(i),pma(i),& ! initial value x,y,z,K1,K2,K2angle 初期位置
          pmrwa(i),pmrwf(i),pmrwd(i),& ! Vrw,Vfreq,direction 回転電場
          fpm1(i),fpm2(i),fpm3(i),fpm4(i),fpm5(i),fpm6(i) ! to be used
!          print "(1p10d11.3)", pmx(i),pmy(i),pmz(i),pmk1(i),pmk2(i),pma(i),pmrwa(i),pmrwf(i),pmrwd(i),fpm1(i) 
          end do
        close(50)
      endif

      if (ikene == 1) then !C    total K(eV) and ratio of velocity directions K(eV)と向きの比率から与える
!C    param.txt(x,y,z,K,Vx:Vy:Vz)
!C     use external file to read parameters
        open (50,File="param.txt")
          do i=1,ipm2
          read(50,*) pmx(i),pmy(i),pmz(i),pmk(i),pmrx(i),pmry(i),pmrz(i)
          write(6,*) pmx(i),pmy(i),pmz(i),pmk(i),pmrx(i),pmry(i),pmrz(i)
          end do
        close(50)
      endif

      if (iexb == 1) then !C   B files to be used with interpolation. not recommended. 磁場を外部ファイルで読み込み．あまり使わない．

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
!          read (17,*) i,i,i,mtx(iex1,iex2,iex3),mty(iex1,iex2,iex3),mtz(iex1,iex2,iex3),&
!            Bxm(iex1,iex2,iex3),Bym(iex1,iex2,iex3),Bzm(iex1,iex2,iex3);
!          write(6,*) mtx(iex1,iex2,iex3),mty(iex1,iex2,iex3),mtz(iex1,iex2,iex3),&
!            Bxm(iex1,iex2,iex3),Bym(iex1,iex2,iex3),Bzm(iex1,iex2,iex3);
          end do
        end do
      end do
      close(17)

      endif

      if (iexe == 1) then !C    E files to be used with interpolation 電場を外部ファイルで読み込み．
      open(17, file='e3d_hd.txt')
        read(17,*) xmin,xmax;! write(6,*)  xmin,xmax
        read(17,*) ymin,ymax;! write(6,*) ymin,ymax
        read(17,*) zmin,zmax;! write(6,*) zmin,zmax
        read(17,*) nx,ny,nz;! write(6,*) nx,ny,nz
        read(17,*) dx,dy,dz;! write(6,*) dx,dy,dz
      close(17)

      write(6,*) "E field header file ### e3d_hd.txt ### was used."

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
      write (cf2, '("e3d", i1.1, ".txt")') iex4
      write(6,*) "read E field file: ",trim(cf2)
      open(17, file=cf2) !C electrode1 電極1
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

!	  ################################################################################
!	  ################################################################################ start
!	  ################################################################################ initial processes
!     ####### set for complete eliptic functions #######
!C    Should be 1 only for the first calculation.
      itf=1 ! 0:load and save, or 1:calculate eliptic functions.
!     ####### complete elliptic functions. usually 0. file is newly made when 1. ここは完全楕円関数，通常は0  #######

!     ##### complete eliptic functions, start ####
      if (itf.eq.1) then
!C   call efile2 only for the first time
!CC     K(k) and E(k) are calculated and saved 楕円関数の計算結果を計算して書き込み
      call efile2(elik2,elie2)
      OPEN (2,File="el2.txt")
      DO 101 i=0,99999
      write(2,'(1p1d12.6,1p2d25.17)') I/100000.0d0,elik2(i),elie2(i)
 101  CONTINUE
      CLOSE(2)
      endif

      if (itf.eq.0) then
!C    Load saved file (this is faster) 楕円関数の計算結果を読み込み（こちらが早い）
!      write(6,*) "Loading complete elliptic integrals..."
      OPEN (50,File="el2.txt")
      DO 102 I=0,99999
      READ(50,*) XI,ELIK2(I),ELIE2(I)
 102  CONTINUE
      CLOSE(50)
      endif

!      write(6,*) "Done. K(k=0) should be 0  1.5707963267948"
!      write(6,*) 0,elik2(0)

      if(ELIK2(1)-1.57079d0.gt.0.01d0) then ! Warning when E() is empty.
      write(6,*) "##### !!Check eliptic functions!! Data is null!!! #####"
      write(6,*) ""
      else
!      write(6,*) "OK. Checked that elik2(1) contains non-zero data."
!      write(6,*) ""
      endif

!    ##### complete eliptic functions, end ####

      if (iclmap == 1) then ! B calculation files are made 磁場の計算結果を出力

!     When hd file(s) are not read, xmax, xmin, etc. are not decided. グリッドデータを使わない時，xmax等が未定．
!     Set the calculation region of magnetic field data ここで範囲を指定する．
      if (xmin .gt. -1.0d-6 .and. xmin .lt. 1.0d-6) then
      if (xmax .gt. -1.0d-6 .and. xmax .lt. 1.0d-6) then
      xmin=-0.25d0; xmax=0.25d0;
      ymin=-0.25d0; ymax=0.25d0;
      zmin=-0.25d0; zmax=0.25d0;
      endif
      endif

      open (50,File="phib.txt") ! B at z=0 z=0での磁場
            yi=0.0d0; zi=0.0d0;
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 細かさを変える時はxmax/10.0d0 等の 10も変えること！
                  xi=(xmax/100.0d0)*iclbp;
                  call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz)
                  call bmap(elik2,elie2,xi,yi,zi,Psi)
!		write(6,*)  iclbp,xmax,xmax/10.0d0
                  write(50,'(1p3d18.10)') xi, Psi, dsqrt ( Bx*Bx+By*By+Bz*Bz )
                  end do
      close(50)

      open (50,File="ck1.txt") ! analytic solution of B at z=0, to check calculation. z=0での磁場の解析解．計算チェック用．
!	  splot [][][0:0.01] "b.txt" using 1:2:4 w l,"ck1.txt" using 1:2:3 w l   
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 細かさを変える時はxmax/10.0d0 等の 10も変えること！
            xi=0.0d0; yi=0.0d0;
            zi=zmin+( (zmax-zmin)/100.0d0 )*iclbp; 
!		write(6,*)  iclbp,xmax,xmax/10.0d0
            write(50,'(1p3d18.10)') xi, zi, 2.0d-7*pi*ci(1)*rc(1)*rc(1) / ((rc(1)*rc(1)+zi*zi)**1.5)
            end do
      close(50)

      open (50,File="b.txt")
            yi=0.0d0
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 細かさを変える時はxmax/10.0d0 等の 10も変えること！
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
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 細かさを変える時はxmax/10.0d0 等の 10も変えること！
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

      open (50,File="ck2.txt") ! BS should be close to 2 pi r Atheta near axis 計算チェック用．
!	  plot "ck2.txt" using 1:2,"ck2.txt" using 1:3 w l
            do iclbp = 0,100 ! when changing here, for example 100, xmax/10.0d0 etc. should be changed below! 細かさを変える時はxmax/10.0d0 等の 10も変えること！
            yi=0.0d0; zi=0.0d0;
            xi=( (xmax-xmin)/100.0d0 )*iclbp/10.0d0; 
            call bvac(elik2,elie2,xi,yi,zi,Bx,By,Bz); call bmap(elik2,elie2,xi,yi,zi,Psi);
            write(50,'(1p3d18.10)') xi, Bz*pi*xi*xi,Psi
            end do
      close(50)
      
      endif
!	  ################################################################################
!	  ################################################################################ end
!	  ################################################################################ initial processes




      do ipm=ipm1,ipm2 ! loop for each particle 各粒子についてのループ
      ieoutc=0 ! a constant to be used for warning outside E file region 電場領域の外に出る警告を与えるための定数
      iden=1; ! time label for density information 密度情報の時間ラベル
      ist=1; ! time label for statistical calculation 統計処理用の時間ラベル

      if (ikene == 2) then !C    parallel and perpendicular K(eV) 縦と横のK(eV)で与える
      x(1)=pmx(ipm); x(2)=pmy(ipm); x(3)=pmz(ipm);
      kini1=pmk1(ipm); kini2=pmk2(ipm); kang=pma(ipm);
      endif

      if (ikene == 1) then !C    total K(eV) and velocity direction K(eV)と向きの比率から与える
      x(1)=pmx(ipm); x(2)=pmy(ipm); x(3)=pmz(ipm);
      kini=pmk(ipm); xdr=pmrx(ipm); ydr=pmry(ipm); zdr=pmrz(ipm);
!!      kini=1.7588200236042485d1; xdr=0.0d0; ydr=1.0d0; zdr=1.0d0; ! -> rL=0.1mm
      endif

      if (ielep == 1) then !     RW amplitude and frequency 回転電場の電位と周波数 サブルーチン内で使用するためにグローバル変数に与える
      Vac(1)=pmrwa(ipm);Vac(2)=pmrwa(ipm);Vac(3)=pmrwa(ipm);Vac(4)=pmrwa(ipm);
      Fel=pmrwf(ipm);
      Phrw=pmrwd(ipm);
      endif

      vfp1=fpm1(ipm);vfp2=fpm2(ipm);vfp3=fpm3(ipm);vfp4=fpm4(ipm);vfp5=fpm5(ipm);vfp6=fpm6(ipm); ! サブルーチン内で使用するためにグローバル変数に与える

      print "(a1)"," "
      print "(a25,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2,a27)","calculation started at ",&
      dtm(1),".",dtm(2),".",dtm(3)," ",dtm(5),":",dtm(6),".",dtm(7), " # RK4 relativistic #"
      print "(a38,i10,i8,i8)","## initial values ## ptcl No. out of:",ipm,ipm2-ipm1+1
      print "(a35,1p3d12.4)","position (x,y,z) (m):   ",x(1),x(2),x(3)
      print "(a35,1p3d12.4)","RW V, freq, dir (V) (Hz) ():   ",Vac(1),Fel,Phrw
      
!C    ##### initial values of E and B at x^0 ####

      if (iexb == 0) then !C    calculate B by Biot-Savart 磁場を計算．
      call bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz);
      B(1)=Bx; B(2)=By; B(3)=Bz;
      print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with usual calculation"
      endif

      if (iexb == 1) then !C    B is given by interpolation of external file data 磁場を外部ファイルで読み込み補間する．
      call bvac(elik2,elie2,x(1),x(2),x(3),Bx,By,Bz);
      print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with usual calculation"

      call bvac2(x(1),x(2),x(3),Bx,By,Bz);
      B(1)=Bx; B(2)=By; B(3)=Bz;
      print "(1p3d13.5,a35)",B(1),B(2),B(3),": B (T) with interpolation"
      endif
      
      if (iexe == 0) then !C    E from usual subroutine 電場は通常のサブルーチンで．
      call evac(x(1),x(2),x(3),tini,Ex,Ey,Ez);
      E(1)=Ex; E(2)=Ey; E(3)=Ez;
      print "(1p3d13.5,a35)",E(1),E(2),E(3),": E (V/m) with usual calculation"
      endif

      if (iexe == 1) then !C    E is given by interpolation of external file data 電場を外部ファイルで読み込み補間する．
      call evac2(x(1),x(2),x(3),tini,Ex,Ey,Ez);
      E(1)=Ex; E(2)=Ey; E(3)=Ez;
      print "(1p3d13.5,a35)",E(1),E(2),E(3),": E (V/m) with interpolation"
      endif


!C    ##### initial velocity ####

      if (ikene == 1) then
!C    total K(eV) and directions K(eV)と向きの比率から与える

      xgm=1.0d0+kini*1.95695119784913d-6 ! XGM=1.0D0+KENE*XE/(XME*XC*XC)
      xbt=dsqrt(1.0d0-1.0d0/(xgm*xgm)); xv=xc*xbt; xu=xgm*xv;
      u(1) = xu*xdr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))
      u(2) = xu*ydr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))
      u(3) = xu*zdr/(dsqrt(xdr*xdr+ydr*ydr+zdr*zdr))

!      write(6,*) u(1),u(2),u(3)
      print "(a35,1p3d12.4)","four-velocity (ux,uy,uz) (m):   ",u(1),u(2),u(3)
      print "(a35,1p3d12.4)","real velocity (vx,vy,vz) (m):   ",u(1)/xgm,u(2)/xgm,u(3)/xgm
      print "(a35,1p2d18.10)","beta and gamma:   ",xbt,xgm
      print "(a35,1p2d18.10)","kinetic energy (eV):   ",kini
      
      endif


      if (ikene == 2) then
!C    parallen and perpendicular K(eV) 縦と横のK(eV)で与える
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
!C    vpara成分の寄与，vperpのAからの寄与(cos成分)，vperpのCからの寄与(sin成分)をプロット
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

      if (izero == 1) then !C    !C    ##### initial v=0 for test 初期速度を0で与える ####
      x(1)=0.0d0; x(2)=0.0d0; x(3)=0.0d0; u(1)=0.0d0; u(2)=0.0d0; u(3)=0.0d0; 
      kini=0.0d0; kini1=0.0d0; kini2=0.0d0;
      write(6,*) "########## This is a test run with initial x=v=0   ##########"
      endif

      Nr = 6;
      Hr = dt;
      xr(1)=x(1); xr(2)=x(2); xr(3)=x(3);
      xr(4)=u(1); xr(5)=u(2); xr(6)=u(3);

      if (its == 1) then ! for variable time step
      iend = iend2
      endif

      Tr = tini

      if (ifile == 1) then ! orbit output files 軌道計算結果のファイルを作成
      
        if (iif .eq. 1) then
        write (filename, '("./orb/orb", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
        write (filenameJ, '("./orb/jp", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
        write (filenamemu, '("./orb/mu", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
!        write (filenameP, '("./orb/p", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
        else
        write (filename, '("orb", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
        write (filenameJ, '("jp", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
        write (filenamemu, '("mu", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
!        write (filenameP, '("p", i4.4, "rkr.txt")') ipm; !      write(6,*) filename
        endif

      open (20,file=filename)
      open (22,file=filenameJ)
!      open (23,file=filenameP)
      open (24,file=filenamemu)
      endif
      ntr=0; ! counter for toroidal rotation トロイダル方向の周回数のカウンタ
      ntrpr=0
      imu=0; ! saved number of mu muの保存回数
      xgpr1=0.0d0; ygpr1=0.0d0; zgpr1=0.0d0; xpr=0.0d0; ypr=0.0d0; zpr=0.0d0; irsapr=0;

!     ######################################################################################
!C    ##### R-K method pusher, start ####
      do it = 0,iend

      if(its == 0) then ! constant time step
      Tr = dt*it;
      endif

      if(its == 1) then ! variable time step
      if (Tr.gt.tmax) go to 10 ! stop calculation when t>tmax tmaxを超えている時は計算終了
      CALL bvac(elik2,elie2,xr(1),xr(2),xr(3),Bx,By,Bz)
      dt2 = ( 2.0d0*pi/(( dsqrt(Bx*Bx+By*By+Bz*Bz) )*eom) ) / tsdiv
      Hr = dt2
      endif ! rk4 calculates new Tr rk4は新しいTrを計算する

      if(iexe == 1) then ! only when E files are used 電場ファイルを使う時のみ
      if(ieout == 1) then ! warning when particle is outside of the E file region 電場領域の外に出た時警告
      if(ieoutc == 0) then ! only once 1回のみ

        if( xr(1) .lt. xmin ) then
        write(6,*) "####### x < xmin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( xr(1) .gt. xmax ) then
        write(6,*) "####### x > xmax of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        
        if( xr(2) .lt. ymin ) then
        write(6,*) "####### y < ymin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( xr(2) .gt. ymax ) then
        write(6,*) "####### y > ymax of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        
        if( xr(3) .lt. ymin ) then
        write(6,*) "####### z < zmin of E file, outside of E field region, continue #######"; ieoutc=1
        endif
        if( xr(3) .gt. ymax ) then
        write(6,*) "####### z > zmax of E file, outside of E field region, continue #######"; ieoutc=1
        endif

      endif
      endif

      if(ieout == 2) then ! stop when outside E region 電場領域の外に出た時終了

        if( xr(1) .lt. xmin ) then
        write(6,*) "####### x < xmin of E file, outside of E field region, stop #######"; go to 10;
        endif
        if( xr(1) .gt. xmax ) then
        write(6,*) "####### x > xmax of E file, outside of E field region, stop #######"; go to 10;
        endif
        
        if( xr(2) .lt. ymin ) then
        write(6,*) "####### y < ymin of E file, outside of E field region, stop #######"; go to 10;
        endif
        if( xr(2) .gt. ymax ) then
        write(6,*) "####### y > ymax of E file, outside of E field region, stop #######"; go to 10;
        endif
        
        if( xr(3) .lt. ymin ) then
        write(6,*) "####### z < zmin of E file, outside of E field region, stop #######"; go to 10;
        endif
        if( xr(3) .gt. ymax ) then
        write(6,*) "####### z > zmax of E file, outside of E field region, stop #######"; go to 10;
        endif
        
      endif
      endif

      if(ihit == 1) then ! stop when hitting electrodes 障害物に衝突した場合に計算を終了
      
        if( dsqrt( xr(1)*xr(1)+xr(2)*xr(2) ) .gt. 0.3d0 ) then
        write(6,*) "positron hit the wall at r=0.3m"
        go to 10
        endif

        if( dabs( xr(3) ) .gt. 0.3d0 ) then
        write(6,*) "positron hit the wall at z=\pm0.2m"
        go to 10
        endif

      endif

!C     四元速度 u(1),u(2),u(3)から実速度やbeta gammaを出す．v^2=c^2 u^2 /(c^2+u^2),: push前
!C    'XU (four velocity), XV (real velocity)'
      xu = dsqrt( xr(4)*xr(4)+xr(5)*xr(5)+xr(6)*xr(6) );
      xv = xc*xu/dsqrt( xc*xc+xu*xu );
      xbt = xv/xc; xgm = 1.0d0/dsqrt(1.0d0-xbt*xbt);
!!        gr1 = dsqrt( (xr(1)-1.0d-4)*(xr(1)-1.0d-4)+xr(2)*xr(2) )
!!        gr1 = dsqrt( vxd(it)*vxd(it)+vyd(it)*vyd(it) )

!     ############################################################
!     ## save data in a file ファイルに値を記録, start ##
      if (ifile == 1) then

      if (icmp == 1) then ! get Kpara Kperp 位置情報に加えてKpara Kperpも計算する
      
      call bvac(elik2,elie2,xr(1),xr(2),xr(3),Bx,By,Bz) ! otherwise no B values これが無いと磁場を計算していない
      ken = (xgm-1.0d0)/1.95695119784913d-6
      xv1 = Bx*xr(4)/xgm+By*xr(5)/xgm+Bz*xr(6)/xgm / dsqrt(Bx*Bx+By*By+Bz*Bz)
      xv2 = dsqrt( xv*xv - xv1*xv1 )
      ken1 = ken * xv1*xv1 / (xv*xv)
      ken2 = ken * xv2*xv2 / (xv*xv)
      Btt = dsqrt(Bx*Bx+By*By+Bz*Bz)
      
      call bmap(elik2,elie2,xr(1),xr(2),xr(3),Psi); Psi1=Psi; ! \Psi 磁気面関数の計算

!C     ########## mu ########## start

      rPsi=Psi1
!      rsa = dsqrt( xr(1)*xr(1)+xr(2)*xr(2) ) - dsqrt( xpr*xpr+ypr*ypr ); ! R increase -> 1, decrease -> -1 Rが増加した時は1を，減少した時は-1を返す．
      rsa = rPsi - rPsipr; ! Psi increase -> 1, decrease -> -1 Psiが増加した時は1を，減少した時は-1を返す．
      if (rsa .gt. 0.0 ) then
      irsa=1; ! write(6,*) irsa;
      endif
      if ( rsa .lt. 0.0d0 ) then
      irsa=-1; ! write(6,*) irsa;
      endif

!           XMU2=XMU2+VPP*VPP/(2.0*sqrt(BX*BX+BY*BY+BZ*BZ))*XGM*XGM
      xmu = xmu + xv2*xv2/(1.0d0*Btt)*xgm*xgm

!      write(6,'(1p3d15.7)') rPsi,rPsipr,rsa
      if ( irsa*irsapr .eq. -1 .and. irsapr .eq. 1 .and. imu .ne. 0 ) then ! detect Psi max Psiの極大点を検出
!      write(6,'(i5,1p3d15.7,i5)') -1, Tr, Psi1, xmu, imu;
      write(24,'(1p3d15.7)') Tr,xmu
      xmu = 0.0d0
      imu=imu+1
      endif
      if ( irsa*irsapr .eq. -1 .and. irsapr .eq. 1 .and. imu .eq. 0 ) then ! initially skip, because the integration did not start from 0 初回は積分の開始が中途半端なのでskip
      imu=imu+1
      xmu = 0.0d0
      endif
      
      irsapr = irsa; rPsipr = rPsi;
!      xpr=xr(1); ypr=xr(2); zpr=xr(3); 

!C     ########## mu ########## end

!     guiding center position - real velocity 案内中心の座標を計算 real velocityで
      Vtpp1 = ( xr(4)-Bx/Btt*(xr(4)*Bx+xr(5)*By+xr(6)*Bz)/Btt )/xgm; ! get Vperp 磁場成分のpara成分を引く->Vperp
      Vtpp2 = ( xr(5)-By/Btt*(xr(4)*Bx+xr(5)*By+xr(6)*Bz)/Btt )/xgm;
      Vtpp3 = ( xr(6)-Bz/Btt*(xr(4)*Bx+xr(5)*By+xr(6)*Bz)/Btt )/xgm;

      Vtppn = dsqrt (Vtpp1*Vtpp1 + Vtpp2*Vtpp2 + Vtpp3*Vtpp3) ! | Vperp | Vperpの大きさ
      gr = Vtppn*xgm/Btt/eom;
!C      XLML=XME*VPP*XGM/XE/(sqrt(BX*BX+BY*BY+BZ*BZ))

!C    guiding center position, not accurate for large rL 案内中心の座標
      xg1 = xr(1) + (Vtpp2*Bz-Vtpp3*By)/Vtppn/Btt * gr
      yg1 = xr(2) + (Vtpp3*Bx-Vtpp1*Bz)/Vtppn/Btt * gr
      zg1 = xr(3) + (Vtpp1*By-Vtpp2*Bx)/Vtppn/Btt * gr

      call bmap(elik2,elie2,xg1,yg1,zg1,Psi); Psi2=Psi; ! \Psi at guiding center 磁気面関数の計算，案内中心で
      
!C     ########## J ##########
!C     integration along B, data out put at every z=0 Jを計算するために，磁力線に沿った積分を行う．赤道面に帰るたびに出力．
      xv1int = xv1int + xv1*xv1*dt
!      xv1int = xv1int + dabs(xv1*dt) !test, to calculate field line length これは磁力線長のテスト
      if (zgpr1*zg1 .lt. 0.0d0) then ! detect z=0 cross 赤道面を横切る場合を検出．
      write(22,'(1p3d15.7)') Tr,xv1int*2.0d0,Psi2
!      write(6,'(1p3d15.7)') Tr,xv1int*2.0d0,Psi2
      xv1int=0.0D0 ! reset the integration value 積分計算を一旦リセットする．
      endif
!C     ########## J ##########

!C     ########## \Psi ########## toroidal rotation number トロイダル方向の周回数を記録
      if (ygpr1*yg1 .lt. 0.0d0 .and. xg1 .gt. 0.0d0) then ! detect x=0 cross with 0<y 0<yでx=0の面を横切る場合を検出
!      write(23,'(1p4d15.7,i4)') Tr,xg1,yg1,zg1,ntr
!      write(6,'(1p4d15.7,i4)') Tr,xg1,yg1,zg1,ntr;
!      write(6,*) "test",ygpr1*yg1
      ntrpr = ntr;
      ntr = ntr + 1; ! add one カウンタを1追加
      endif
!C     ########## \Psi ##########

      xgpr1=xg1; ygpr1=yg1; zgpr1=zg1; ! save previous values 前回の値を記録しておく
      
      endif ! get Kpara Kperp, 位置情報に加えてKpara Kperpも計算する end

        if(mod(it,imbk).eq.0) then
          if (icmp == 0) then
          write(20,'(1p8d15.7)') xr(1),xr(2),xr(3),dsqrt( xr(1)*xr(1)+xr(2)*xr(2) ),xr(4)/xgm,xr(5)/xgm,xr(6)/xgm,Tr
          endif
          if (icmp == 1) then
          write(20,'(1p25d15.7,i6)') xr(1),xr(2),xr(3),dsqrt( xr(1)*xr(1)+xr(2)*xr(2) ),xr(4)/xgm,xr(5)/xgm,xr(6)/xgm,Tr,&
                                  ken,ken1,ken2,Bx,By,Bz,Btt,Ex,Ey,Ez,Psi1,xg1,yg1,zg1,dsqrt(xg1*xg1+yg1*yg1),&
                                  Psi2,Vtppn*Vtppn/Btt, ntrpr
!     1:x 2:y 3:z 4:r 5:vx 6:vy 7:vz 8:time
!     9:Kttl 10:Kpara 11:Kperp 12:Bx 13:By 14:Bz 15:B 16:Ex 17:Ey 18:Ez 19:Psi 20:xg.c. 21:yg.c. 22:zg.c. 23:rg.c. 24:Psi_g.g. 25:mu,26:t.rotation


          endif
        endif
      endif
!     ## save data ファイルに値を記録, end ##
!     ############################################################


!     ############################################################
!!     ## data for statistical calculation 統計処理用の情報を記録, start ##
      if (ifile2 == 1) then
        if(mod(it,imbk2).eq.0) then

        stt(ist,ipm) = Tr; ! (time info, particle info) これらは時間番号，粒子番号
        stx(ist,ipm) = xr(1);
        sty(ist,ipm) = xr(2);
        stz(ist,ipm) = xr(3);
        str(ist,ipm) = dsqrt( xr(1)*xr(1)+xr(2)*xr(2) );
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
        
!        write(6,*) ipm,ist,stt(ist,ipm); ! (time info, particle info) 時間番号，粒子番号

!      write(6,*) imbk2t,ipmt
!      imbk2t ! istep/imbk2+1 total time steps トータルの時間ステップ数
!      ipmt = ipm2-ipm1+1 ! total particle number, or maximum value of ipm トータル粒子数，ipmの最大値

        ist = ist+1
        endif
      endif
!     ## data for statistical calculation, 統計処理用の情報を記録 end ##
!     ############################################################


      call rk4(elik2,elie2,Nr,Tr,Xr,K1,K2,K3,K4,YWORK,Hr)
!     Hr is added to Tr in rk4.

      end do
      if (ifile == 1) then
      close(20)
      close(22)
      close(24)
      endif

 10   if(its == 0) then ! variable time step
      print "(a35,i15,1p2d14.6)","step and time at stop, dt (s):",it-1,Tr-dt,dt
      print "(a35,1p3d13.5)","(x,y,z):",xr(1),xr(2),xr(3)
      print "(a35,1p3d13.5)","K (eV), pr,pp,ttl:",ken1,ken2,ken
      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
      print "(a35,i12)","toroidal rotation:",ntrpr
      endif

      if(its == 1) then ! variable time step
      print "(a35,i12,1p1d15.7)","step and time at stop (s):",it-1,Tr-dt
      print "(a35,1p3d13.5)","(x,y,z):",xr(1),xr(2),xr(3)
      print "(a35,1p3d13.5)","K (eV), pr.pp,ttl:",ken1,ken2,ken
      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
      print "(a35,i12)","toroidal rotation:",ntrpr
      endif

      end do ! come from do ipm=ipm1,ipm2


!     ############################################################
!     minimum statistical calculation 統計処理をここで
      if (ifile2 == 1) then
      
!      double precision,dimension(1:105) :: stta,stxa,stya,stza,stra,stgxa,stgya,stgza,stgra,stKa,stKpra,stKppa,stps1a,stps2a,stmua
!      double precision,dimension(1:105) :: stts,stxs,stys,stzs,strs,stgxs,stgys,stgzs,stgrs,stKs,stKprs,stKpps,stps1s,stps2s,stmus
!      integer, dimension(1:105) :: strta,strts ! averaged value and sigma for many particles 粒子番号（後者）での平均値と標準偏差
!
!      call avesig(stt,imbk2t,ipmt,stt2,stt3)
!      imbk2t ! istep/imbk2+1 total time steps トータルの時間ステップ数
!      ipmt = ipm2-ipm1+1 ! total particle number, max of ipm トータル粒子数，ipmの最大値
      write(6,*) "# Statistical calculation done # TimeStep/ptcl:",imbk2t,ipmt

!     calculate averaged value and sigma at each of time steps 各時間ステップで平均値と標準偏差を計算する．
      call avesig(stt,imbk2t,ipmt,stta,stts)
      call avesig(stx,imbk2t,ipmt,stxa,stxs)
      call avesig(sty,imbk2t,ipmt,stya,stys)
      call avesig(stz,imbk2t,ipmt,stza,stzs)
      call avesig(str,imbk2t,ipmt,stra,strs)
      call avesig(stgx,imbk2t,ipmt,stgxa,stgxs)
      call avesig(stgy,imbk2t,ipmt,stgya,stgys)
      call avesig(stgz,imbk2t,ipmt,stgza,stgzs)
      call avesig(stgr,imbk2t,ipmt,stgra,stgrs)
      call avesig(stK,imbk2t,ipmt,stKa,stKs)
      call avesig(stKpr,imbk2t,ipmt,stKpra,stKprs)
      call avesig(stKpp,imbk2t,ipmt,stKppa,stKpps)
      call avesig(stps1,imbk2t,ipmt,stps1a,stps1s)
      call avesig(stps2,imbk2t,ipmt,stps2a,stps2s)
      call avesig(stmu,imbk2t,ipmt,stmua,stmus)
      
!      print "(a5,1p3d12.4,a10,1p2d12.4)","1,t",stt(1,1),stt(1,2),stt(1,3),"ave/sig:",stta(1),stts(1)
!      print "(a5,1p3d12.4,a10,1p2d12.4)","2,t",stt(2,1),stt(2,2),stt(2,3),"ave/sig:",stta(2),stts(2)
!      print "(a5,1p3d12.4,a10,1p2d12.4)","3,t",stt(3,1),stt(3,2),stt(3,3),"ave/sig:",stta(3),stts(3)

!      print "(a5,1p3d12.4,a10,1p2d12.4)","1,x",stx(1,1),stx(1,2),stx(1,3),"ave/sig:",stxa(1),stxs(1)
!      print "(a5,1p3d12.4,a10,1p2d12.4)","2,x",stx(2,1),stx(2,2),stx(2,3),"ave/sig:",stxa(2),stxs(2)
!      print "(a5,1p3d12.4,a10,1p2d12.4)","3,x",stx(3,1),stx(3,2),stx(3,3),"ave/sig:",stxa(3),stxs(3)

      open (51,File="orb_a0rkr.txt") ! save file 0, average 記録ファイル
      do isti=1,imbk2t

      write(51,'(1p15d15.7)')&
      stxa(isti),stya(isti),stza(isti),stra(isti),stgxa(isti),stgya(isti),stgza(isti),stgra(isti),&
      stKa(isti),stKpra(isti),stKppa(isti),stps1a(isti),stps2a(isti),stmua(isti),&
      stta(isti)
!     1:x 2:y 3:z 4:r 5:gx 6:gy 7:gz 8:gr 9:K 10:Kpr 11:Kpp 12:Ps 13:gPs 14:mu 15:time

      end do
      close(51)

      open (52,File="orb_a1rkr.txt") ! save file 1, average+sigma 記録ファイル
      do isti=1,imbk2t

      write(52,'(1p15d15.7)')&
      stxa(isti)+stxs(isti),stya(isti)+stys(isti),stza(isti)+stzs(isti),stra(isti)+strs(isti),&
      stgxa(isti)+stgxs(isti),stgya(isti)+stgys(isti),stgza(isti)+stgzs(isti),stgra(isti)+stgrs(isti),&
      stKa(isti)+stKs(isti),stKpra(isti)+stKprs(isti),stKppa(isti)+stKpps(isti),&
      stps1a(isti)+stps1s(isti),stps2a(isti)+stps2s(isti),stmua(isti)+stmus(isti),&
      stta(isti)

      end do
      close(52)

      
      open (53,File="orb_a2rkr.txt") ! save file 2, average-sigma 記録ファイル
      do isti=1,imbk2t

      write(53,'(1p15d15.7)')&
      stxa(isti)-stxs(isti),stya(isti)-stys(isti),stza(isti)-stzs(isti),stra(isti)-strs(isti),&
      stgxa(isti)-stgxs(isti),stgya(isti)-stgys(isti),stgza(isti)-stgzs(isti),stgra(isti)-stgrs(isti),&
      stKa(isti)-stKs(isti),stKpra(isti)-stKprs(isti),stKppa(isti)-stKpps(isti),&
      stps1a(isti)-stps1s(isti),stps2a(isti)-stps2s(isti),stmua(isti)-stmus(isti),&
      stta(isti)

      end do
      close(53)
      
      endif
      
!     ############################################################

     call system_clock(t2, t_rate, t_max)   ! save ending time 終了時を記録
       if ( t2 < t1 ) then
      diff = (t_max - t1) + t2 + 1
        else
       diff = t2 - t1
       endif

      close(21)

      print "(A23,F10.3,F10.3)", "### End of calculation."
      print "(A29,F10.3,F10.3)", "### calculation time (s) (h):", diff/dble(t_rate), diff/dble(t_rate)/3600.0d0
      
     
      end




!C ###################################################
      subroutine rk4(elik2,elie2,Nr,Tr,Xr,k1,k2,k3,k4,ywork,Hr)
!C ###################################################
      implicit none
      double precision,dimension(0:100005) :: elik2,elie2
      double precision :: Hr,Tr,F,Brt,Bzt
      double precision,dimension(1:6) :: Xr,k1,k2,k3,k4,ywork
      integer i,Nr

!        call bloop(elik2,elie2,0.1d0,0.2d0,0.0d0,1.0d3,Brt,Bzt)
!        write(6,*) Brt,Bzt,elik2(1000),elie2(2000)

      do 10 i=1,Nr
  10  k1(i)=Hr*F(elik2,elie2,i,Tr,Xr,Nr)
      do 11 i=1,Nr
  11  ywork(i)=Xr(i)+0.5d0*k1(i)
      do 12 i=1,Nr
  12  k2(i)=Hr*F(elik2,elie2,i,Tr+0.5d0*Hr,ywork,Nr)
      do 13 i=1,Nr
  13  YWORK(i)=Xr(i)+0.5d0*k2(i)
      do 14 i=1,Nr
  14  k3(i)=Hr*F(elik2,elie2,i,Tr+0.5d0*Hr,ywork,Nr)
      do 15 i=1,Nr
  15  ywork(i)=Xr(i)+K3(i)
      do 16 i=1,Nr
  16  k4(i)=Hr*F(elik2,elie2,i,Tr+Hr,ywork,Nr)
      do 17 i=1,Nr
  17  Xr(i)=Xr(i)+(k1(I)+2.0d0*k2(i)+2.0d0*k3(i)+k4(i))/6.0d0
      Tr=Tr+Hr
!C      write(6,*)T,H
      RETURN
      END



!C ###################################################
      function F(elik2,elie2,i,t,x,n)
!C ###################################################
      use mod_variable
      implicit none
      double precision,dimension(0:100005) :: elik2,elie2
      double precision :: t,F,XQOM,xu,xv,xbt,xgm,Brt,Bzt
      double precision,dimension(1:3) :: B,E
      double precision,dimension(1:6) :: x
      integer :: n,i
      double precision :: xp,yp,zp,Ex,Ey,Ez,Bx,By,Bz

      XQOM = 1.7588200236d11

      XP=X(1)
      YP=X(2)
      ZP=X(3)

!C     four-velocity X(4),X(5),X(6) -> real velocity, v^2=c^2 u^2 /(c^2+u^2)
!C     四元速度 X(4),X(5),X(6)から，実速度を出す．v^2=c^2 u^2 /(c^2+u^2)
      xu = dsqrt( x(4)*x(4)+x(5)*x(5)+x(6)*x(6) )
      xv = xc * xu / dsqrt( xc*xc+xu*xu )
!C      XV=XU/sqrt(1.0+(XU/XC)*(XU/XC)) <-larger error こちらの方が誤差が大きい

!C      WRITE(6,*) 'XU (four velocity), XV (real velocity)'
!C      WRITE(6,602) XU, XV

      xbt = xv/xc
      xgm = 1.0d0/dsqrt(1.0d0-xbt*xbt)

      if (iexb == 0) then !C B by Biot-Savart
      CALL bvac(elik2,elie2,xp,yp,zp,Bx,By,Bz)
      endif

      if (iexb == 1) then !C    interpolate B in external files 磁場を外部ファイルで読み込み補間する．
      CALL bvac2(xp,yp,zp,Bx,By,Bz)
      endif

      B(1)=Bx
      B(2)=By
      B(3)=Bz

      if (iexe == 0) then !C    interpolate E in external files 電場は通常のサブルーチンで．
      CALL evac(XP,YP,ZP,T,Ex,Ey,Ez)
      endif

      if (iexe == 1) then !C    interpolate E in external files 電場を外部ファイルで読み込み補間する．
      call evac2(XP,YP,ZP,T,Ex,Ey,Ez);
      endif
      
      E(1)=Ex; E(2)=Ey; E(3)=Ez;
      

      GO TO (11,12,13,14,15,16),I
  11  F=X(4)/xgm
      RETURN
  12  F=X(5)/xgm
      RETURN
  13  F=X(6)/xgm
      RETURN
  14  F=XQOM*( E(1) + (X(5)*B(3)-X(6)*B(2))/xgm )
      RETURN
  15  F=XQOM*( E(2) + (X(6)*B(1)-X(4)*B(3))/xgm )
      RETURN
  16  F=XQOM*( E(3) + (X(4)*B(2)-X(5)*B(1))/xgm )
      RETURN
 
      END


!     rk4 should include el in oder to get B
!     bvacでbloopを呼ぶと磁場が計算されない elが渡されていなかった -> rk4に変数追加
!C     =================================================================
      subroutine bvac(elik2,elie2,xp,yp,zp,Bx,By,Bz)
!C     calculate B at (x,y,z) xyz座標の点での磁場を計算．各粒子位置で毎回計算する．
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,Bx,By,Bz,xi,Bro,Bzo,ctmp,Brt,Bzt,Bx1,By1,Bz1
      double precision,dimension(0:100005) :: elik2,elie2
      integer i,itf

!      ci(1) = 1.0d3 ! current
!      rc(1) = 1.0d-1 ! loop radius
!      zcl(1) = 0.0d0 ! loop z
!      xcl(1) = 0.0d0 ! loop offset in x direction from z axis?
!      ycl(1) = 0.0d0 ! usually zero or middle

        Bx=0.0d0; By=0.0d0; Bz=0.0d0;

        do i=1,icln ! add for each of coils 各コイルについて足し合わせる 

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
!C     return \psi=2*pi*r*Atheta at (x,y,z) xyz座標の点で \psi=2*pi*r*Athetaを計算して返す
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,ctmp,at,Psi,Psi1
      double precision,dimension(0:100005) :: elik2,elie2
      integer i,itf

      Psi = 0.0d0

      do i=1,icln !  add for each of coils 各コイルについて足し合わせる

      ctmp = dsqrt((xp-xcl(i))*(xp-xcl(i))+(yp-ycl(i))*(yp-ycl(i)))        
      call rathe(elik2,elie2,rc(i),ctmp,zp-zcl(i),ci(i),at)
      Psi1 = 2.0d0*pi*at
      Psi = Psi + Psi1
      
      end do
      
      return
      end

!C     =================================================================
      subroutine bvac2(xp,yp,zp,Bx,By,Bz);
!C     return B at (x,y,z) obtained by interpolation xyz座標の点での磁場を計算．読み込んだ磁場データを使用
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp
      double precision :: Bx,By,Bz
      double precision :: f111,f211,f121,f221,f112,f212,f122,f222,A1,A2,A3,A4,B1,B2

!      write(6,*) "initial position:"
!      write(6,*) xp,yp,zp ! serch (0,0,0) for this point このポイントに対しての(0,0,0)を探す
!      write(6,*) ( (xp-xmin)/dx ),int( (xp-xmin)/dx )
!      write(6,*) ( (yp-ymin)/dy ),int( (yp-ymin)/dy )
!      write(6,*) ( (zp-ymin)/dz ),int( (zp-ymin)/dz )
      n0x=int( (xp-xmin)/dx )+1; n0y=int( (yp-ymin)/dy )+1; n0z=int( (zp-ymin)/dz )+1;

!      write(6,*) "n0x,n0y,n0z at f111:"
!      write(6,*) n0x,n0y,n0z !
!      write(6,*) "x,y,z at f111:"
!      write(6,*) mt1x(n0x),mt1y(n0y),mt1z(n0z) ! this corresponds to f111 これがf111に相当する場所
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
      subroutine evac(xp,yp,zp,tp,Ex,Ey,Ez)
!C     retuen E at (x,y,z) 電場を各点で指定する．
!C     =================================================================
      implicit none
      double precision :: xp,yp,zp,tp,Ex,Ey,Ez
      
! Homogeneous E
      Ex = 0.0d0
      Ey = 0.0d0
      Ez = 0.0d0

! Homogeneous B in finite region
!      if ( xp .gt. 0.1d0 .and. xp .lt. 0.2d0 ) then
!      if ( yp .gt. 0.1d0 .and. yp .lt. 0.2d0 ) then
!      if ( zp .gt. 0.1d0 .and. zp .lt. 0.2d0 ) then
!        Ex = 0.0d0;
!        Ey = 0.0d0;
!        Ez = 0.1d0;
!      endif
!      endif
!      endif

      RETURN
      END



!C     =================================================================
      subroutine evac2(xp,yp,zp,tp,Ex,Ey,Ez);
!C     return E at (x,y,z) obtained by interpolation and superposition
!C     電場を各点で指定する： 場所，時間，各単位電場に対して返す．
!C     点を囲む格子点上の8点を出し，その8点での電場情報から線形補間する．
!C     =================================================================
      use mod_variable
      implicit none
      double precision :: xp,yp,zp,tp,Ex,Ey,Ez
      double precision :: f111,f211,f121,f221,f112,f212,f122,f222,A1,A2,A3,A4,B1,B2
      double precision :: RWamp,RWfreq,RWphase
      integer :: iex1,iex2,iex3

!      write(6,*) tp
!      write(6,*) "initial position:"
!      write(6,*) xp,yp,zp ! find (0,0,0) for this point このポイントに対しての(0,0,0)を探す
!      write(6,*) ( (xp-xmin)/dx ),int( (xp-xmin)/dx )
!      write(6,*) ( (yp-ymin)/dy ),int( (yp-ymin)/dy )
!      write(6,*) ( (zp-ymin)/dz ),int( (zp-ymin)/dz )
      n0x=int( (xp-xmin)/dx )+1; n0y=int( (yp-ymin)/dy )+1; n0z=int( (zp-ymin)/dz )+1;

!      write(6,*) "n0x,n0y,n0z at f111:"
!      write(6,*) n0x,n0y,n0z !
!      write(6,*) "x,y,z at f111:"
!      write(6,*) mt1x(n0x),mt1y(n0y),mt1z(n0z) ! this corresponds to f111 これがf111に相当する場所
!      write(6,*) "B at f111:"
!      write(6,*) Bxm(n0x,n0y,n0z),Bym(n0x,n0y,n0z),Bzm(n0x,n0y,n0z)

!     At 8 lattice points surrounding (x,y,z), E by each of the lectrodes are cauculated and suporposed.

!      Amp(1) = 1.0d0; Amp(2) = 0.0d0; Amp(3) = 0.0d0; Amp(4) = 0.0d0; Amp(5) = 0.0d0; 
!      Amp(6) = 0.0d0; Amp(7) = 0.0d0; Amp(8) = 0.0d0; Amp(9) = 0.0d0; Amp(10) = 0.0d0
!     RWamp = 100.0d0; RWfreq = 100.0d3; RWphase = 0.0d0*pi
!     Amp3 = RWamp * dsin ( 2.0d0 * pi * RWfreq * tp - RWphase + 0.0d0*pi/2.0d0 )
!     Amp4 = RWamp * dsin ( 2.0d0 * pi * RWfreq * tp - RWphase + 1.0d0*pi/2.0d0 )
!     Amp5 = RWamp * dsin ( 2.0d0 * pi * RWfreq * tp - RWphase + 2.0d0*pi/2.0d0)
!     Amp6 = RWamp * dsin ( 2.0d0 * pi * RWfreq * tp - RWphase + 3.0d0*pi/2.0d0)

      Vel(1) = Vdc(1) + Vac(1)*dsin( 2.0d0*pi*Fel*tp + Phrw*0.0d0*pi/2.0d0 )
      Vel(2) = Vdc(2) + Vac(2)*dsin( 2.0d0*pi*Fel*tp + Phrw*1.0d0*pi/2.0d0 )
      Vel(3) = Vdc(3) + Vac(3)*dsin( 2.0d0*pi*Fel*tp + Phrw*2.0d0*pi/2.0d0 )
      Vel(4) = Vdc(4) + Vac(4)*dsin( 2.0d0*pi*Fel*tp + Phrw*3.0d0*pi/2.0d0 )

!      write(6,*) Fel,Vac(1)

!      do iele = 1,ielemax
!      Vel(iele) = Vdc(iele) + Vac(iele)*dsin( 2.0d0*pi*Fel*tp )
!!      Vel(iele) = Vdc(iele) ! fast for static E 変動電場を使わない時はこちらが高速
!      end do

!      write(6,*) Vel(1),tp

!C    Ex by interpolation Exの計算，各電極が8点で作るExを重み付きで足し合わせ，線形補間する．
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

!C    Ey Eyの計算，各電極が8点で作るEyを重み付きで足し合わせ，線形補間する．
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

!C    Ez Ezの計算，各電極が8点で作るEzを重み付きで足し合わせ，線形補間する．
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

      RETURN
      END


!*     =================================================================
      SUBROUTINE BLOOP(elik2,elie2,RC,R,Z,CI,BR,BZ)
!*     calculate the Br and Bz produced by a loop current
!*     RC:loop radius R,Z:position CI:coil current 円形電流の作る磁場を計算。
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

      IF(ZK.GE.0.999999D0) GO TO 10

      G1=dsqrt(1.0D0-ZK)
      IF(ZK.GT.0.9998D0) FK=DLOG(4.0D0/G1)+0.25D0*(DLOG(4.0D0/G1)-1.0D0)*G1*G1
      IF(ZK.GT.0.9998D0) FE=1.0D0+0.5D0*(DLOG(4.0D0/G1)-0.5D0)*G1*G1     
      IF(ZK.GT.0.9998D0) GO TO 20

      I=IDINT(ZZK*100000.0D0)
      
!C    FK FE interpolation FKとFEを線形補間する start
      i2=i+1;
      fk1=elik2(i); fe1=elie2(i);
      fk2=elik2(i2); fe2=elie2(i2);
!      write(6,'(d13.5,i7,1p3d17.9)') zzk,i,fk1,fk2,fk1+(zzk-dble(i)/1.0d5)*(fk2-fk1)
      fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
      fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);
!      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
!C    FK FE interpolation end

!C    従来通りのFK FE
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
      
   10 BZ=0.0D0
      BR=0.0D0
      RETURN
      END


!*     =================================================================
      SUBROUTINE rathe(elik2,elie2,rc,x,z,ci,at)
!*     calculate "Atheta" produced by a loop current
!*     RC:loop radius R,Z:position CI:coil current 円形電流の作る磁気面を計算。rA_thetaを計算する。
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
!*     ZK is k^2 ZK とは k^2 のこと.
      ZK=4.0d0*rc*x/((rc+x)*(rc+x)+z*z)
      ZZK=dsqrt(ZK)
      
      IF(ZK.GE.0.999999d0) GO TO 20
      IF(ZK.GT.0.9998d0) GO TO 10
      
      A0=2.0d0*XCC*CI*x/dsqrt(ZK)
      A1=dsqrt(RC/x)
      A2=1.0d0-0.5d0*ZK
      I=IDINT(ZZK*100000.0d0)
      
!C    FK FE interpolation FKとFEを線形補間する start
      i2=i+1;
      fk1=elik2(i); fe1=elie2(i);
      fk2=elik2(i2); fe2=elie2(i2);
!      write(6,'(d13.5,i7,1p3d17.9)') zzk,i,fk1,fk2,fk1+(zzk-dble(i)/1.0d5)*(fk2-fk1)
      fk=fk1+(zzk*1.0d5-dble(i))*(fk2-fk1);
      fe=fe1+(zzk*1.0d5-dble(i))*(fe2-fe1);
!      write(21,'(i5,1p7d15.7)') ipm,Tr-dt,xr(1),xr(2),xr(3),ken1,ken2,ken
!C    FK FE interpolation end

!     従来通りのFK FE
!      FK=elik2(I)
!      FE=elie2(I)

      AT=A0*A1*(A2*FK-FE)

      IF(I.EQ.0) THEN
      AT=0.0d0
      ENDIF
      
      RETURN

 10   A1=XCC*CI*RC
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
      double precision,dimension(0:100005) :: elik2,elie2
      double precision ZZZAAA,FEF1,FEF2
      integer ji
      
      write(6,*) "start: elliptic fuction calculation"
      write(6,*) "... may takes several seconds ..."

      DO 5 JI=0,99999
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
!     複数粒子の時系列データの平均値と標準偏差を計算するサブルーチン．
!     i,jは時間番号，粒子番号（後者に対して平均化）．
!     ave1(i,j)をjについて平均化したave2(i)，標準偏差sig2(i)を返す．
!     数値の個数は 時間のステップ数i:itt1, 粒子数j:itt2
!C    =================================================================
      implicit none
      double precision,dimension(1:105,1:1005) :: data1
      double precision,dimension(1:105) :: ave1,sig1
      integer iav1,iav2,itt1,itt2
      double precision sum1,sum2

!     ##### averaging 平均値 start #####
      do iav1=1,itt1 ! for each of the time step 各時間ステップ数に対して計算．

         sum1 = 0.0d0
         do iav2 = 1,itt2 ! particle number 粒子数
         sum1 = sum1 + data1(iav1,iav2)
         end do

         ave1(iav1) = sum1/dble(itt2)

!      write(6,*)iav1, ave1(iav1)

      end do
!     ##### averaging 平均値 end #####

!     ##### sigma 標準偏差 start #####
      do iav1=1,itt1 ! for each of the time step 各時間ステップ数に対して計算．
      
         sum2 = 0.0d0
         do iav2 = 1,itt2 ! particle number 粒子数
         sum2 = sum2 + ( data1(iav1,iav2)-ave1(iav1) )*( data1(iav1,iav2)-ave1(iav1) )
         end do
         
         sig1(iav1) = dsqrt( sum2/itt2 )

!      write(6,*)iav1, sig1(iav1)

      end do
!     ##### sigma 標準偏差 end #####

      return
      end

