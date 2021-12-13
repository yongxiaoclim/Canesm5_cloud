;********************************************************
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/shea_util.ncl"

; This script is to calculated low cloud fraction below 700hpa.
; the low cloud fraction is is under a random-overlap assumption [Manabe and Strickler, 1964; 
;Ramanathan et al., 1983; Stephens, 1984].

begin

;prepare lat lon a for interpolation in the following calculation
lat=fspan(-90, 90, 90)
lon=fspan(0, 359, 180)

; loop across model members  

do z = 0, dimsizes(tune)-1

    ncdf     = addfile(tune(z)+"_cl_layer.nc" ,"c")  ; create a netCDF file
    print(tune(z))  ;print the loop step
    
    ;open file seperately by 'tune' , 'tune' is an array consisted by all available tunning member
    fls_z = systemfunc("ls "+path+"CanESM5-tune-esm001-"+tune(z)+"/pdControl/r"+inv(z)+"i1p2f1/Amon/cl/gn/v20190429/"+"cl_Amon_CanESM5-tune-esm001-"+tune(z)+"*")
    fz=addfiles(fls_z,"r")
    ListSetType(fz,"join")
    
    ; read lat lon 
    lat_old=fz[0] ->lat
    lon_old=fz[0] ->lon
    
    ; the time change change, the follwoing time selcted to be same as 
    ; oberved quantity described in Zhai 2014 paper
    yrStrt = 1980
    yrLast = 2004
    year  = ispan(yrStrt, yrLast, 1)
    TIME   = fz[0]->time
    ymStrt = yrStrt*100 +  1
    ymLast = yrLast*100 + 12
    YYYY   = cd_calendar(TIME,-1)/100                 
    iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    ; read cl and ps
    cl = fz[:]->cl(:,iYYYY,:,:,:)
    PS = fz[:]->ps(:,iYYYY,:,:) 
    
    ;delete TIME and YYYY after use up
    delete(TIME)
    delete(YYYY)
    
    ; define cloud fraction in monthly mean calculation
    cl_mon=new((/12,dimsizes(cl(:,0,0,0,0)),dimsizes(cl(0,0,:,0,0)),dimsizes(lat_old),dimsizes(lon_old)/), float, cl@_FillValue)

    ; define cl_new for new array of cloud fraction; define air_p for the following air pressure calculation  
    cl_new=new((/12,dimsizes(cl(:,0,0,0,0)),dimsizes(cl(0,0,:,0,0)),dimsizes(lat),dimsizes(lon)/), typeof(cl), cl@_FillValue)
    ; air_p=new((/12,dimsizes(cl(:,0,0,0,0)),dimsizes(cl(0,0,:,0,0)),dimsizes(lat),dimsizes(lon)/), float, cl@_FillValue)
    
    do m = 0, dimsizes(cl(:,0,0,0,0))-1 ; m represent a realization in multiple initial condition for a specific model member, this is mainly designed for CMIP6 models with multiple ensembles. 
      
      do l = 0, dimsizes(cl(0,0,:,0,0))-1 ; l represent a atmosperic layer. due to the three dimensional limit for clmMonTLL function  
        cl_mon(:,m,l,:,:)=clmMonTLL(cl(time|:,ncl_join|m,lev|l,lat|:,lon|:))
      end do
       
       cl_new(:,m,:,:,:)=linint2_Wrap (lon_old,lat_old,cl_mon(:,m,:,:,:), True, lon,lat, 0) ; remap cloud fraction 

    end do

  ;define low cloud fraction and will will asigned quantity in the following loop
  lcc=new((/12,dimsizes(cl(:,0,0,0,0)),dimsizes(lat),dimsizes(lon)/), float, cl@_FillValue)

  ;loop for all atmospheric layers and detect low cloud as cloud layer below 700hPa
  ;and low cloud fraction is under a random-overlap assumption [Manabe and Strickler, 1964; 
  ;Ramanathan et al., 1983; Stephens, 1984].

  do m = 0, dimsizes(cl(:,0,0,0,0))-1 
      do t=0,12-1 ;loop 12 months
        do i = 0, dimsizes(lat)-1
           do j=0,dimsizes(lon)-1
              aa1=1
              aa=tofloat(aa1)
              do n = 0, dimsizes(hybm)-1
                 if (.not.ismissing(air_p(t,m,n,i,j)).and.air_p(t,m,n,i,j).gt.70000) then
                    aa=(1-cl_new(t,m,n,i,j)/100)*aa
                 end if
              end do
              lcc(t,m,i,j)=1-aa 
          end do
      end do
    end do
  end do


  ; assgin Coordinate to low cloud fraction calculated above 
  lat@units= "degrees_north"                                                                        
  lon@units= "degrees_east"  
  lcc!0="month" 
  lcc!1="runs" 
  lcc!2="lat"
  lcc!3="lon"
  lcc&lat=lat
  lcc&lon=lon

  ; write out the variable    
  ncdf->lcc =lcc  

end do

end