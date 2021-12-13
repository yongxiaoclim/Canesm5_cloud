;********************************************************
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/tsanta/liang/code/ncl/ncl/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
a    = addfile("/tsanta/liang/code/ncl/cloud_obs/landsea.nc","r")
lsdata = a->LSMASK


;prepare lat lon a for interpolation in the following calculation
lat=fspan(-90, 90, 90)
lon=fspan(0, 359, 180)

; loop across model members 
do i = 0, (dimsizes(model_name)-1)

    ncdf     = addfile(model_name(i)+"_cloud.nc" ,"c")  ;  create a netCDF file
    print(model_name(i)) ;print the loop step
    
    ; read cloud fraction at 850 and 950 hPa
    fz850=addfile(path+model_name(i)+"_cl_layer.nc","r")
    fz950=addfile(path+model_name(i)+"_cl_layer.nc","r")

    
    ; read wap whcih has been interpolated and the timeslice has been specified as cloud fraction
    fw=addfile("/net/venus/kenes/data/cmip6/CanESM5_tune/omega/wget/"+model_name(i)+"_cl_layer.nc","r")
    
    ; the time change change, the follwoing time selcted to be same as 
    ; description in Brient 2016 paper
    yrStrt = 2006
    yrLast = 2012
    TIME   = fz1[0]->time
    ymStrt = yrStrt*100 +  1
    ymLast = yrLast*100 + 12
    YYYY   = cd_calendar(TIME,-1)/100                 
    iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    
    clt850 =dim_avg_n_Wrap(fz850->lcc,0)
    clt950 =dim_avg_n_Wrap(fz950->lcc,0)

    ; wap should be in units of mb/day
    wap = fw->wap_remap
    wap_enmean_clim=dim_avg_n_Wrap(wap, (/0,1/))*36*24

    ; detect weak subsidence region using wap
    do j = 0, dimsizes(clt850(:,0,0))-1
       clt850(j,:,:)=where(wap_enmean_clim.gt.0.and.wap_enmean_clim.lt.10 ,clt850(j,:,:),clt850@_FillValue)
       clt950(j,:,:)=where(wap_enmean_clim.gt.0.and.wap_enmean_clim.lt.10 ,clt950(j,:,:),clt850@_FillValue)
 
     end do
    lat@units= "degrees_north"
    lat!0="lat" 
    lat&lat=lat
    
    ;calculate BCS metric by area mean in tropical region 
    rad    = 4.0*atan(1.0)/180.0
    lats=lat({-30:30})
    lcc_850_ave = wgt_areaave_Wrap(clt850(:,{-30:30},:), cos(lats*rad), 1.0, 1)
    lcc_950_ave = wgt_areaave_Wrap(clt950(:,{-30:30},:), cos(lats*rad), 1.0, 1)
    
    ;calculate BCS metric which which quantifies the ratio of 
    ;cloud fraction below the 900-hPa level (∼1 km), CF950,
    ;to the total low-cloud fraction below the 
    ;800-hPa level (∼2 km), CF850 + CF950. 
    index=lcc_950_ave/(lcc_850_ave+lcc_950_ave)

    ncdf->index=index
    
    delete(fls_z1)
    delete(fz1)     
    delete(lat_old)
    delete(lon_old)
    delete(iYYYY)
    delete(YYYY)
    delete(clt850)
    delete(clt950)
    delete(wap)
    delete(TIME)
    delete(lcc_850_ave)
    delete(lcc_950_ave)
    delete(index)

 end do
end