#!/usr/bin/env python
# coding: utf-8

# In[6]:


def search_cldsatfile(initial_time,day_of_year):
    cldlon=[]
    cldlat=[]
    cldtime=[]
    FU=[]
    FU_NA=[]
    cldsat_fnames=sorted(glob.glob(cldsat_dir+year+str(day_of_year).zfill(3)+'*'))
    Ncldsat=len(cldsat_fnames)
    if (Ncldsat > 0):
        cldsat_initial_time=np.zeros((Ncldsat),'f')
        for fi in np.arange(Ncldsat):
            cldsat_initial_time[fi]=float(cldsat_fnames[fi][-57:-55])+            float(cldsat_fnames[fi][-55:-53])/60.0+float(cldsat_fnames[fi][-53:-51])/3600.0
        
            timediff=abs(initial_time-cldsat_initial_time)
            fileloc=np.argmin(timediff)
        
            cldsat_fname=cldsat_fnames[fileloc]
            cldsat_initial_time1=cldsat_initial_time[fileloc]  
        
            #== read cloudsat data ===
            f=HDF(cldsat_fname,SDC.READ)
            vs=f.vstart()
            atcldlat=vs.attach('Latitude')
            cldlat=np.array(atcldlat[:])
            atcldlon=vs.attach('Longitude')
            cldlon=np.array(atcldlon[:])
            atcldtime=vs.attach('Profile_time')
            cldtime=np.array(atcldtime[:])/3600.0+cldsat_initial_time1
            ind=np.where(cldtime >=24)[0]
            cldtime[ind]=cldtime[ind]-24
            atFD = vs.attach('FD_TOA_IncomingSolar')
            FD = np.array(atFD[:])/10.0
            scene_st=vs.attach('Scene_status')
            scene_status=np.array(scene_st[:])
            scene_st.detach()
            atFD.detach()
            atcldlat.detach()
            atcldlon.detach()
            atcldtime.detach()
            hdf=SD(cldsat_fname,SDC.READ)
            FU = hdf.select('FU_TOA')[:]/10.0
            FU_NA = hdf.select('FU_NA_TOA')[:]/10.0
            TOACRE = hdf.select('TOACRE')[:]
            
    #print('rematch',cldsat_initial_time,cldsat_fname)        
    return cldlon,cldlat,cldtime,FU,FU_NA,TOACRE,scene_status
    
def vfm_feature_flag(val):
    feature_type=0
    feature_type_qa = 0
    ice_water_phase = 0
    ice_water_phase_qa = 0
    feature_subtype = 0
    cloud_aerosol_psc_type_qa = 0
    horizontal_averaging = 0
    aerotyp_flag=0
    
    for bi in np.arange(16):
        if (val%2 != 0):
            #print('bit ',bi+1)
        
            if (bi+1 == 1):
                feature_type=feature_type+1
            elif (bi+1 == 2):
                feature_type=feature_type+2
            elif (bi+1 == 3):
                feature_type=feature_type+4
            elif (bi+1 == 4): 
                feature_type_qa=feature_type_qa+1
            elif (bi+1 == 5):
                feature_type_qa=feature_type_qa+2  
            elif (bi+1 == 6):
                ice_water_phase=ice_water_phase+1
            elif (bi+1 == 7): 
                ice_water_phase=ice_water_phase+2
            elif (bi+1 == 8):
                ice_water_phase_qa=ice_water_phase_qa+1
            elif (bi+1 == 9):
                ice_water_phase_qa=ice_water_phase_qa+2
            elif (bi+1 == 10):
                feature_subtype = feature_subtype+1
            elif (bi+1 == 11):
                feature_subtype = feature_subtype+2
            elif (bi+1 == 12):
                feature_subtype = feature_subtype+4
            elif (bi+1 == 13):
                cloud_aerosol_psc_type_qa = cloud_aerosol_psc_type_qa + 1
            elif (bi+1 == 14):
                horizontal_averaging = horizontal_averaging + 1
            elif (bi+1 == 15):
                horizontal_averaging = horizontal_averaging + 2
            else:
                horizontal_averaging = horizontal_averaging + 4
        val=int(val/2)
    return feature_type, ice_water_phase   
                         

# this script is to colocate CALIPSO and 2B-FLXHR-LIDAR data to obtain cirrus impact on aerosol DRE
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime, date, timedelta
import numpy as np
import math
from pyhdf.SD import SD, SDC
from pyhdf.HDF import *
from pyhdf.VS import *
import h5py
import glob

print('start', datetime.now())

year='2007'

wfname='aerosol_radiative_effect_cldsat_calipso_'+year+'.h5'
cldsat_dir='/data/keeling/a/yulanh/c/BW_backup/CloudSat/2B-FLXHR-LIDAR/'+year+'/'
fnames=sorted(glob.glob('/data/gdi/f/yulanh/CALIPSO/MLay/'+year+'/'+'*ZD.hdf'))
lonres=6
latres=5
Ntau=501 # binres 0.01
Nci =471 #binres 0.015
dimx=int(360/lonres)
dimy=int(180/latres)
write_lon=np.arange(dimx)*lonres-180+lonres/2.0
write_lat=np.arange(dimy)*latres-180+latres/2.0

obsnum=np.zeros((dimx,dimy),'i')
aeroonly_taupdf=np.zeros((dimx,dimy,Ntau),'i')
aeroci_taupdf=np.zeros((dimx,dimy,Ntau,Nci),'i')
aeroonly_are = np.zeros((dimx,dimy,Ntau),'f')
aeroci_are = np.zeros((dimx,dimy,Ntau,Nci),'f')

#fnames=['/data/gdi/f/yulanh/CALIPSO/MLay/2007/CAL_LID_L2_05kmMLay-Standard-V4-20.2007-01-01T11-02-31ZD.hdf']

for fname in fnames:
    #print(fname)
    # === get the julday ===
    month=fname[-20:-18]
    day = fname[-17:-15]
    initial_hh= int(fname[-14:-12])
    initial_mm= int(fname[-11:-9])
    initial_time=initial_hh+initial_mm/60.0
    day_of_year=date(int(year),int(month),int(day)).timetuple().tm_yday
        
    cldsat_fnames=sorted(glob.glob(cldsat_dir+year+str(day_of_year).zfill(3)+'*'))
    Ncldsat=len(cldsat_fnames)
    #===get inital time for each cloudsat fname ===
    if (Ncldsat > 0):
        cldsat_initial_time=np.zeros((Ncldsat),'f')
        for fi in np.arange(Ncldsat):
            cldsat_initial_time[fi]=float(cldsat_fnames[fi][-57:-55])+                float(cldsat_fnames[fi][-55:-53])/60.0+float(cldsat_fnames[fi][-53:-51])/3600.0
        
        timediff=abs(initial_time-cldsat_initial_time)
        fileloc=np.argmin(timediff)
        
        cldsat_fname=cldsat_fnames[fileloc]
        cldsat_initial_time1=cldsat_initial_time[fileloc]
        
        if (initial_time < cldsat_initial_time[fileloc]):
            cldsat_fname=cldsat_fnames[fileloc-1]
            cldsat_initial_time1=cldsat_initial_time[fileloc-1]
        
        #print('frist match',cldsat_fname)
        
        #== read cloudsat data ===
        f=HDF(cldsat_fname,SDC.READ)
        vs=f.vstart()
        atcldlat=vs.attach('Latitude')
        cldlat=np.array(atcldlat[:])
        atcldlon=vs.attach('Longitude')
        cldlon=np.array(atcldlon[:])
        atcldtime=vs.attach('Profile_time')
        cldtime=np.array(atcldtime[:])/3600.0+cldsat_initial_time1
        ind=np.where(cldtime >=24)[0]
        cldtime[ind]=cldtime[ind]-24
        atFD = vs.attach('FD_TOA_IncomingSolar')
        FD = np.array(atFD[:])/10.0
        scene_st=vs.attach('Scene_status')
        scene_status=np.array(scene_st[:])
        scene_st.detach()
        atFD.detach()
        atcldlat.detach()
        atcldlon.detach()
        atcldtime.detach()

        hdf=SD(cldsat_fname,SDC.READ)
        FU = hdf.select('FU_TOA')[:]/10.0
        FU_NA = hdf.select('FU_NA_TOA')[:]/10.0
        TOACRE = hdf.select('TOACRE')[:]
  
        #==== read aerosol data ===
        hdf=SD(fname)
        lidlat=hdf.select('Latitude')
        lidlat=lidlat[:,:]
        lidlon=hdf.select('Longitude')
        lidlon=lidlon[:,:]
        lidtime=hdf.select('Profile_UTC_Time')[:,:]
        cad_score=hdf.select('CAD_Score')
        cad_score=cad_score[:,:]
        extinc_qc=hdf.select('ExtinctionQC_532')[:,:]
        feature_flag=hdf.select('Feature_Classification_Flags')[:,:]
        cod=hdf.select('Column_Optical_Depth_Cloud_532')[:]
        cod=cod.reshape(len(cod))
        aod=hdf.select('Column_Optical_Depth_Tropospheric_Aerosols_532')[:]
        aod=aod.reshape(len(aod))
        layer_top=hdf.select('Layer_Top_Altitude')[:,:]
        layer_base=hdf.select('Layer_Base_Altitude')[:,:]

        matchdata_flag=0 # used for mark if rematch the data
        #aodind=np.where(aod > 0)[0]
        
        for si in range(len(aod)):
            #print(si,aodind[si])
            tplidlat=lidlat[si,1]
            tplidlon=lidlon[si,1]
            tplidtime=(lidtime[si,1]-int(lidtime[si,1]))*24
            tpcad_score=cad_score[si,:]
            tpext = extinc_qc[si,:]
            tpcod = cod[si]
            tpaod = aod[si]
            tplaytop=layer_top[si,:]
            tplaybase=layer_base[si,:]
            tpfeature=feature_flag[si,:]

            # record obsnum
            lon_res= int(np.round((tplidlon+180-lonres/2.0)/lonres))
            lat_res= int(np.round((tplidlat+90-latres/2.0)/latres))
            if (lon_res < 0):
                lon_res=0
            if (lon_res > dimx-1):
                lon_res = dimx-1
                
            if (lat_res < 0):
                lat_res=0
            if (lat_res > dimy-1):
                lat_res=dimy-1
                
            obsnum[lon_res,lat_res]=obsnum[lon_res,lat_res]+1
           
            tphour=int(tplidtime/3600)
            residual=((tplidtime/3600-tphour)*60)
            tpminute=int(residual)
            tpsecond=int((residual-tpminute)*60)

            # === get clear sky and below-cirrus aerosol flag === 
          
            altind=np.where((tplaytop > 0))[0]
            Nh=len(altind)
            lowcad_flag=0
            extqc_flag=0
            cirrus_flag=0
            liquid_flag=0
            aero_flag=0
            aero_topmost=0
            cirrus_baselow =0
            #print(tpcad_score,tplaybase,tplaytop,tpext,tpfeature)
            for hi in range(Nh):
                tplaybase_1=tplaybase[altind[hi]]
                tplaytop_1 = tplaytop[altind[hi]]
                #print('cadscore',tpcad_score[altind[hi]])
                if ((abs(tpcad_score[altind[hi]]) < 20) | (abs(tpcad_score[altind[hi]]) > 100)):
                    lowcad_flag =1
                   
                tpext_1=tpext[altind[hi]]
               
                if ((tpext_1 != 0) & (tpext_1 != 1) & (tpext_1 != 16) & (tpext_1 != 18)):
                    if (extqc_flag ==0):
                        extqc_flag=1
             
                if ((abs(tpcad_score[altind[hi]]) >= 20) & (abs(tpcad_score[altind[hi]]) <= 100) & (extqc_flag ==0)):
                    feature_type,cloud_phase=vfm_feature_flag(tpfeature[altind[hi]])
                    # dealing with ice clouds
                    if ((cloud_phase == 1) | (cloud_phase == 3)):
                        cirrus_flag=1
                        if (cirrus_baselow ==0):
                            cirrus_baselow=tplaybase_1
                        if ((cirrus_baselow !=0) & (cirrus_baselow > tplaybase_1)):
                            cirrus_baselow=tplaybase_1
                            
                    if (cloud_phase ==2):
                        liquid_flag=1
                        
                    # dealing with aerosols
                    if ((feature_type == 3) | (feature_type ==4)):
                        aero_flag=1
                        if (aero_topmost ==0):
                            aero_topmost=tplaytop_1
                        if ((aero_topmost != 0) & (aero_topmost < tplaytop_1)):
                            aero_topmost=tplaytop_1
       
                    
            #=== end deal with high confidence aero and cloud ======
            
            # start dealing with aerosol only ====
            if ((aero_flag == 1) & (tpcod == 0) & (extqc_flag ==0) & (lowcad_flag ==0) & (tpaod >0)):
                #=== match to cloudsat =====
                swath_locdiff=np.sqrt((tplidlon-cldlon)**2+(tplidlat-cldlat)**2)
                swath_loc=np.argmin(swath_locdiff)
                
                #print(si,tplidtime,cldtime[swath_loc],swath_locdiff[swath_loc])
                if ((swath_locdiff[swath_loc] > 0.035) & (matchdata_flag ==0)):
                    #print('not match',si,aodind[si],swath_locdiff[swath_loc],fname)
                    # read new file
                    lidyymmdd=str(int(lidtime[si,1]))
                    month=lidyymmdd[1:3]
                    day=lidyymmdd[3:5]
                    day_of_year=date(int(year),int(month),int(day)).timetuple().tm_yday
                    initial_time=(lidtime[si,1]-int(lidtime[si,1]))*24
                    cldlon,cldlat,cldtime,FU,FU_NA,TOACRE,scene_status=search_cldsatfile(initial_time,day_of_year)
                    matchdata_flag=1
                    swath_locdiff=np.sqrt((tplidlon-cldlon)**2+(tplidlat-cldlat)**2)
                    swath_loc=np.argmin(swath_locdiff)
                
                scene_status_flag=bin(scene_status[swath_loc][0]).replace("0b","")
                scene_status_flag=scene_status_flag.zfill(16)
                subpixel_cloud_flag=1
                                
                if ((scene_status_flag[-2] == '0') & (scene_status_flag[-3] == '0') & (scene_status_flag[-4] == '0') &
                     (scene_status_flag[-5] == '0') & (scene_status_flag[-6] == '0')):
                    subpixel_cloud_flag=0
                    
                if ((swath_locdiff[swath_loc] <= 0.035) & (subpixel_cloud_flag==0) & (TOACRE[0, swath_loc] == 0)):
                    aerotau_scp=int((np.log10(tpaod)+4)/0.01)
                    if (aerotau_scp < 0):
                        aerotau_scp=0
                    if (aerotau_scp > Ntau-1):
                        aerotau_scp=Ntau-1
                    aeroonly_taupdf[lon_res,lat_res,aerotau_scp]=                    aeroonly_taupdf[lon_res,lat_res,aerotau_scp]+1
                    
                    aeroonly_are[lon_res,lat_res,aerotau_scp]=                    aeroonly_are[lon_res,lat_res,aerotau_scp]+(FU_NA[0,swath_loc]-FU[0,swath_loc])
                    #print(si,'aeroonly',tpaod,FU_NA[0,swath_loc]-FU[0,swath_loc], scene_status_flag, TOACRE[0,swath_loc])
                   
                    
            #print('flag',si,extqc_flag,lowcad_flag,aero_flag,cirrus_flag,tpaod,tpcod,swath_locdiff[swath_loc])    
            if ((extqc_flag ==0) & (lowcad_flag ==0) & (aero_flag ==1) & (cirrus_flag ==1) &  
               (aero_topmost < cirrus_baselow) & (tpcod > 0) & (tpaod >0) & (liquid_flag ==0)): 
                #print(si,tpcod,tpaod,aero_topmost, cirrus_baselow)
                #=== match to cloudsat =====
                swath_locdiff=np.sqrt((tplidlon-cldlon)**2+(tplidlat-cldlat)**2)
                swath_loc=np.argmin(swath_locdiff)
                #print(si,tplidtime,cldtime[swath_loc],swath_locdiff[swath_loc])
                if ((swath_locdiff[swath_loc] > 0.035) & (matchdata_flag==0)):
                    lidyymmdd=str(int(lidtime[si,1]))
                    month=lidyymmdd[1:3]
                    day=lidyymmdd[3:5]
                    day_of_year=date(int(year),int(month),int(day)).timetuple().tm_yday
                    initial_time=(lidtime[si,1]-int(lidtime[si,1]))*24
                    cldlon,cldlat,cldtime,FU,FU_NA,TOACRE,scene_status=search_cldsatfile(initial_time,day_of_year)
                    matchdata_flag=1
                    swath_locdiff=np.sqrt((tplidlon-cldlon)**2+(tplidlat-cldlat)**2)
                    swath_loc=np.argmin(swath_locdiff)
               
                scene_status_flag=bin(scene_status[swath_loc][0]).replace("0b","")
                scene_status_flag=scene_status_flag.zfill(16)
                if ((scene_status_flag[-3] == '0') & (scene_status_flag[-4] == '1') & (scene_status_flag[-5] == '0') & 
                    (scene_status_flag[-6] == '0')):
                    subpixel_cloud_flag=0
                    
                if ((swath_locdiff[swath_loc] <= 0.035) & (subpixel_cloud_flag == 0) & (abs(TOACRE[0,swath_loc]) > 0)):
                    aerotau_scp=int((np.log10(tpaod)+4)/0.01)
                    if (aerotau_scp < 0):
                        aerotau_scp=0
                    if (aerotau_scp > Ntau-1):
                        aerotau_scp=Ntau-1
                    if (np.isnan(tpcod)):
                        print(fname,tpcad_score)
                        
                    citau_scp=int((np.log10(tpcod)+4)/0.015)
                    if (citau_scp < 0):
                        citau_scp=0
                    if (citau_scp > Nci-1):
                        citau_scp=Nci-1
                    
                    aeroci_taupdf[lon_res,lat_res,aerotau_scp,citau_scp]=                    aeroci_taupdf[lon_res,lat_res,aerotau_scp,citau_scp]+1
                    
                    aeroci_are[lon_res,lat_res,aerotau_scp,citau_scp]=                    aeroci_are[lon_res,lat_res,aerotau_scp,citau_scp]+(FU_NA[0,swath_loc]-FU[0,swath_loc])
                    #print(si,'cirrus_aero',tpcod,tpaod,FU_NA[0,swath_loc]-FU[0,swath_loc],scene_status_flag,TOACRE[0,swath_loc])
                    
                       
                    
hf=h5py.File(wfname,'w')
hf.create_dataset('lon',data=write_lon)
hf.create_dataset('lat',data=write_lat)
hf.create_dataset('aeroci_taupdf',data=aeroci_taupdf)
hf.create_dataset('aeroci_are',data=aeroci_are)
hf.create_dataset('aeroonly_taupdf',data=aeroonly_taupdf)
hf.create_dataset('aeroonly_are',data=aeroonly_are)
hf.create_dataset('obsnum',data=obsnum)
hf.close()

print('finish', datetime.now())


