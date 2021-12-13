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

do i = 0, (dimsizes(model_name)-1)

    ncdf     = addfile(model_name(i)+"_cloud.nc" ,"c")  ; open output netCDF file
    ; print(model_name(i)) ;print the loop step

    ; read low cloud fraction have been calculated
    fz=addfile("/net/venus/kenes/data/cmip6/CanESM5_tune/lcc/wget/open/"+model_name(i)+"_cl_layer.nc","r")
    
    ; load ts file
    fls_z1 = systemfunc("ls "+"/net/venus/kenes/data/cmip6/CanESM5_tune/data/pd_control/ftp.cccma.ec.gc.ca/CCCMA/products/TUNING/CMIP6/CMIP/CCCma/"+"CanESM5-tune-esm001-"+model_name(i)+"/pdControl/r"+inv(i)+"i1p2f1/Amon/ts/gn/v20190429/"+"ts_Amon_CanESM5-tune-esm001-"+model_name(i)+"*")
    fz1=addfiles(fls_z1,"r")
    ListSetType(fz1,"join")
    
    ; load omega, units is pa/s
    fw=addfile("/net/venus/kenes/data/cmip6/CanESM5_tune/omega/wget/"+model_name(i)+"_cl_layer.nc","r")
    
    ; load inital lat lon to be interpolated
    lat_old=fz1[0] ->lat
    lon_old=fz1[0] ->lon

    ; select time 
    yrStrt = 2007
    yrLast = 2014
    TIME   = fz1[0]->time
    ymStrt = yrStrt*100 +  1
    ymLast = yrLast*100 + 12
    YYYY   = cd_calendar(TIME,-1)/100                 
    iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    
    ; variable should be well prepared, cloud fraction 'clt', surface temperature 'ts', 'wap'
    clt_new = fz->lcc
    ts = fz1[:]->ts(:,iYYYY,:,:)
    wap = fw->wap_remap
    wap_enmean=dim_avg_n_Wrap(wap, 1)
    
    ; calculate climatology monthly mean surface temperature
    ts_mon=new((/12,dimsizes(ts(:,0,0,0)),dimsizes(ts(0,0,:,0)),dimsizes(ts(0,0,0,:))/), float, ts@_FillValue)
    do m = 0, dimsizes(ts(:,0,0,0))-1
      ts_mon(:,m,:,:)=clmMonTLL(ts(time|:,ncl_join|m,lat|:,lon|:))
    end do
    
    ts_new=linint2_Wrap (lon_old,lat_old,ts_mon, True, lon,lat, 0)
    
    ; mask land 
    lsm = landsea_mask(lsdata,clt_new&lat,clt_new&lon)
    clt_new = mask(clt_new,lsm.eq.0,True)*100 ; multiple 100 to get %
    ts_new = mask(ts_new,lsm.eq.0,True)
   
   ; get subsidence region detected by vertical volocity
   do j = 0, dimsizes(clt_new(0,:,0,0))-1
      clt_new(:,j,:,:)=mask(clt_new(:,j,:,:), wap_enmean.gt.0, True)
      ts_new(:,j,:,:)=mask(ts_new(:,j,:,:), wap_enmean.gt.0, True)
   end do

    clt_new!1="runs"
    ts_new!1="runs"
    copy_VarMeta(ts_new, clt_new)
    
    ; area mean of low cloud fraction and surface temperature
    rad    = 4.0*atan(1.0)/180.0
    lat@units= "degrees_north"
    lat!0="lat" 
    lat&lat=lat
    lats=lat({-40:-20})
    latn=lat({20:40})
    lcc_n=wgt_areaave_Wrap(clt_new(:,:,{20:40},:), cos(latn*rad), 1.0, 1)
    lcc_s=wgt_areaave_Wrap(clt_new(:,:,{-40:-20},:), cos(lats*rad), 1.0, 1)

    ts_n=wgt_areaave_Wrap(ts_new(:,:,{20:40},:), cos(latn*rad), 1.0, 1)
    ts_s=wgt_areaave_Wrap(ts_new(:,:,{-40:-20},:), cos(lats*rad), 1.0, 1)
    
    ; regress low cloud fraction change to surafce 
    ;temperature change, and get the  linear regression coefficient
    lcc_sst_s1=regCoef(ts_s(runs|:,month|:),lcc_s(runs|:,month|:))
    lcc_sst_n1=regCoef(ts_n(runs|:,month|:),lcc_n(runs|:,month|:))
   
    lcc_sst_g1=0.5*(lcc_sst_n1 + lcc_sst_s1)
    
    ;MBLC metric is lcc_sst_g
    lcc_sst_g= rm_single_dims(lcc_sst_g1)
    lcc_sst_s= rm_single_dims(lcc_sst_s1)
    lcc_sst_n= rm_single_dims(lcc_sst_n1)
    ; print(lcc_sst_g)
    
    ncdf->lcc_sst_g=lcc_sst_g
    ; ncdf->lcc_sst_s=lcc_sst_s
    ; ncdf->lcc_sst_n=lcc_sst_n

    delete(TIME)
    delete(fls_z1)
    delete(fz1)     
    delete(lat_old)
    delete(lon_old)
    delete(iYYYY)
    delete(YYYY)
    delete(ts)
    delete(clt)
    delete(ts_mon)
    delete(clt_new)
    delete(ts_new)
    ; delete(lcc_sst)
    delete(wap)
    delete(lcc_n)
    delete(lcc_s)
    delete(ts_n)
    delete(ts_s)
    delete(lcc_sst_s)
    delete(lcc_sst_n)
    delete(lcc_sst_g)
    delete(lcc_sst_g1)
    delete(lcc_sst_s1)
    delete(lcc_sst_n1)

   end do
end