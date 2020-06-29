program insolation4_tst
!----------------------------------------------------------------------------------------
! sample program showing the use of the "twostars_m" module to compute the top of the
! atmosphere insolation on a circumbinary planet 
! 
! written by Siegfried Eggl 20140610
!
! modified by Siegfried Eggl 20180409 
!   > modified output to csv
!
! modified by Anna Zuckerman 20190814
!   > modified input to include stellar spectra
!   > modified output to include spectra, weight factors, and distances
!
! last modified by Siegfried Eggl 20200629
!   > added support for non-zero initial eccentricity of the orbit of the planet
!
! dependencies: module twostars_m
!
! Input: file "insolation4.inn"
!        file "input_spectra1.csv"
!        file "input_spectra2.csv"
!
! Output: file "twostars_out_general.csv" 
!         file "twostars_out_bndflux1.csv"
!         file "twostars_out_bndflux2.csv"
!         file "twostars_out_bndfluxnet.csv"
!         file "twostars_out_distances.csv"
!         file "twostars_out_weights.csv"
!----------------------------------------------------------------------------------------
use twostars_m 

implicit none

!define
integer::i
integer, parameter :: dp = selected_real_kind(15, 307)         !precision parameter
integer::nbands                                                !number of spectral bands
real(kind=dp),parameter::deg2rad=0.0174532925199432957692369076849_dp ![1/deg] PI / 180
real(kind=dp),parameter::rad2deg=1._dp/deg2rad                 !180/PI
real(kind=dp)::e1                                              !eccentricity of the stellar binary
real(kind=dp),dimension(1:2)::teff,rs,semia,man                !stellar effective temperatures, stellar radii, semimajor axes of binary and planetary orbits, initial mean longitudes 
real(kind=dp),dimension(1:3)::mass                             !system masses: m(1:2) stars, m(3) planet
real(kind=dp),dimension(1:4)::spin                             !spin quantities: (1) obliquity [deg] (2) angle of precession [deg] (3) rotation period [D] (4) initial hour angle [deg] (=initial position of continents wrt x-axis) 
real(kind=dp),dimension(:,:,:),allocatable::bands              !spectral band information (upper and lower limits per band per star for each band)
real(kind=dp)::t,tend,dt                                       !time, end time, output timestep
real(kind=dp),dimension(:,:),allocatable::spect1               !Col 1 is spectrum of star1, col 2 is lower band limit, col 3 is upper band limit. ADZ ADD 8/1/19
real(kind=dp),dimension(:,:),allocatable::spect2               !Col 1 is spectrum of star2, col 2 is lower band limit, col 3 is upper band limit. ADZ ADD 8/1/19
real(kind=dp),dimension(:,:),allocatable::spectra              !Spectra of stars 1, 2 (row -> star, cols-> spectra)(flux at each wavelength band) ADZ ADD 8/1/19

real(kind=dp) :: e2t0, e2t                                      !initial & evolved orbital eccentricity of the circumbinary planet
real(kind=dp):: w1, w2                                          ! argument of pericenter of binary and planetary orbit [deg]
!results
integer::ok                                                     !configuration check parameter: 0-> everything is ok, else...
real(kind=dp),dimension(1:2)::decl,phi                          !declination and hour angle wrt each star
real(kind=dp),dimension(:,:),allocatable::bndflux               !band fluxes for each star and each band
real(kind=dp),dimension(1:2)::d2,dr1,dr2                        !distances to stars  !ADZ CHANGE 7/17/19
real(kind=dp),dimension(1:2)::weights                           !scaling factor by which to weight each stars spectrum at each time step (due to distance to star and area blocked during eclipses) !ADZ ADD 8/4/19

character(len=200)::fmt_out1 !changed 20180409 
character(len=200)::fmt_out2 !ADZ CHANGE 7/25/19
character(len=200)::fmt_out3 !ADZ CHANGE 7/25/19
character(len=200)::fmt_out4 !ADZ CHANGE 8/8/19
character(len=25)::str_nbands !ADZ CHANGE 8/4/19 string representation of nbands for output formatting

real(kind=dp):: aq(1:2)



 !-----------------------------------------------------------------
 !read input file
 !----------------------------------------------------------------
 open(21,file="insolation4.inn",status='old')
 open(22,file="input_spectra1.csv", status = 'old') !ADZ ADD 7/29/19  
 open(23,file="input_spectra2.csv", status = 'old') !ADZ ADD 8/1/19
 read(21,*)tend,dt        !end time, output timestep
 read(21,*)mass           !masses of star1,star2,planet [Msun]
 read(21,*)rs             !stellar radii [Rsun]
 read(21,*)semia          !binary and planetary semimajor axes [au]
 read(21,*)e1, e2t0       !binary eccentricity / initial planetary orbital eccentricity
 read(21,*)man            !initial mean longitudes [rad] (starting positions of the planet and the stars on their orbits)
 read(21,*)spin           !planetary spin quantities: (1) obliquity [deg] (2) angle of precession [deg] (3) rotation period [D] (4) initial hour angle [deg] (=initial position of continents wrt x-axis) 
 read(21,*)nbands         !number of spectral bands        !ADZ CHANGE 7/29/29 
 
allocate(spect1(nbands,3),spect2(nbands,3)) ! ADZ ADD 8/1/19
 read(22,*)
 read(23,*)
 do i=1,nbands
    read(22,*) spect1(i,:)
    read(23,*) spect2(i,:)
 end do

allocate(bands(2,nbands,2),bndflux(2,nbands),spectra(2,nbands)) !ADZ ADD 8/1/19
 bands(1,:,1) = spect1(:,2) ! lower band limits star 1
 bands(1,:,2) = spect1(:,3) ! upper band limits star 1
 bands(2,:,1) = spect2(:,2) ! lower band limits star 2
 bands(2,:,2) = spect2(:,3) ! upper band limits star 2
 spectra(1,:) = spect1(:,1) ! star 1 spectra
 spectra(2,:) = spect2(:,1) ! star 2 spectra


 !-----------------------------------------------------------------
 !check whether the dynamical configuration is possible (stellar radii and effective temperatures are not checked)
 !----------------------------------------------------------------
 aq(1)=semia(1)
 aq(2)=semia(2)*(1._dp-e2t0) 
 call checkconf(mass,aq,e1,ok)
 write(*,*)'configuration ok?',ok
 if(ok.gt.0.and.ok.le.4) STOP

 !----------------------------------------------------------------
 ! Prepare output files
 !------------------------------------------------------------------
 fmt_out1="(5 (F21.15, ', '),2 (ES13.6, ', '))"  !changed 20180409 ADZ CHANGE 7/23/19 and  8/5/19 
 write(str_nbands,*) nbands !get nbands as a string for format string
 fmt_out2="("//trim(str_nbands)//" (E13.6, ', '))"   !ADZ CHANGE 7/25/19
 fmt_out3="(4 (F21.15, ', '))"   !ADZ CHANGE 7/25/19
 fmt_out4="(2 (F21.15, ', '))"   !ADZ CHANGE 8/8/19

 !produce the output files with the insolation3 subroutine
 open(24,file="twostars3_out_general.csv",status='replace')      !changed 20180409
 !ADZ CHANGE 7/23/19: output files for distances, spectra, and weights 
 open(25, file = 'twostars_out_distances.csv')                  !time dep distances to stars from planet (x,y coords of stars)  
 open(26, file="twostars_out_bndflux1.csv", status='replace')   !time dep spectra for star 1 reaching planet    
 open(27, file="twostars_out_bndflux2.csv", status='replace')   !time dep spectra for star 2 
 open(28, file="twostars_out_bndfluxnet.csv", status='replace') !time dep net spectra reaching planet
 open(29, file = 'twostars_out_weights.csv', status='replace')  !scaling factor for flux from each star over time (due to distance, eclipses)
 open(30, file = 'eccentricity.dat', status='replace')           !eccentricity and argument of pericenter over time for planet and binary  
 
 !write header into ouput files
 write(24,*)'1:  time [D], 2: declination wrt star1 [rad], 3: declination wrt star2 [rad],',&
 '4: hour angle wrt star 1 [rad], 5: hour angle wrt star2 [rad],', &
 '6: sum flux star1 [W/m^2], 7: sum flux star2 [W/m^2]'
 write(25,*) 'Distance to star 1 [AU],', '' , 'Distance to star 2 [AU] '             
 write(25,*) 'x coord, y coord, x coord, y coord'
 write(26,*) 'Star 1 Spectra at Planet [W/m^2])  NOTE: cols-> wavelengths  rows-> time steps  row 1-> spectral bin lower limit'
 write(26,fmt_out2) bands(1,:,1)	 
 write(27,*) 'Star 2 Spectra at Planet  [W/m^2]  NOTE: cols-> wavelengths  rows-> time steps  row 1-> spectral bin lower limit'
 write(27,fmt_out2) bands(2,:,1) 

 write(28,*) 'Net Flux [W/m^2]  NOTE: cols -> wavelengths  rows-> time steps  row 1-> spectral bin lower limit'
 if(all(bands(1,:,1).EQ.bands(2,:,1))) then
    write(28,fmt_out2) bands(1,:,1)
 else
    write(28,*) 'SPECTRAL BANDS MUST MATCH FOR BOTH STELLAR SPECTRA'
 end if
 write(29,*) 'Scaling factor for each star''s spectra due to distance to stars and area blocked during eclipses' 

 !--------------------------------------------------------------------
 ! Calculations start here (time evolution of the system)
 !--------------------------------------------------------------------
 t=0.d0
 do while(t.lt.tend)
  !call module function
  call insolation4(spectra,mass,rs,semia,e1,e2t0,man,spin,t,nbands,bands,bndflux,decl,phi,d2,dr1,dr2,weights,e2t,w1,w2)!ADZ ADD 7/29/19 added spectra and weights, remove teff 8/15/19
  !write output at each timestep
  write(24,fmt=fmt_out1)t,decl(:),phi(:),sum(bndflux(1,:)),sum(bndflux(2,:))     !changed 20180409   !ADZ CHANGE 7/17/19
  write(25,fmt=fmt_out3) dr1, dr2                                                !ADZ CHANGE 7/23/19
  write(26,fmt=fmt_out2) bndflux(1,1:nbands)                                     !ADZ CHANGE 7/23/19  
  write(27,fmt=fmt_out2) bndflux(2,1:nbands)                                     !ADZ CHANGE 7/24/19 
  write(28,fmt=fmt_out2) bndflux(1,1:nbands)+bndflux(2,1:nbands)                 !ADZ CHANGE 7/24/19 
  write(29,fmt=fmt_out4) weights                                                 !ADZ CHANGE 8/8/19
  write(30,*)t,e2t,w2*rad2deg,e1,w1*rad2deg  
  t=t+dt 
 end do

! wrapup
close(22)
close(23)
close(24)
close(25)
close(26)
close(27)
close(28)
close(29)
close(30)
end program 
