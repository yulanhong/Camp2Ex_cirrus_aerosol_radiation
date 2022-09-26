#!/usr/bin/env python
# coding: utf-8

# In[1]:


# this script test simulations only for cirrus-aerosol sky for four simulations, cirrus-aero, cirrus-only, aeroonly, clean

def write_swinput(SZA, cirrus_od, aero_tau,swoutfile, day_of_year):
    swfname = open('swinput1', 'w')
    swfname.write(
        'atmosphere_file /data/keeling/a/yulanh/c/software_install/libradtran/2.0.2/share/libRadtran/data/atmmod/afglt.dat \n')
    swfname.write(
        'source solar /data/keeling/a/yulanh/c/software_install/libradtran/2.0.2/share/libRadtran/data/solar_flux/kurudz_1.0nm.dat \n')
    swfname.write('wavelength 250 4000\n')
    swfname.write('mol_abs_param fu\n')
    swfname.write('radiosonde atm_file_ecmwf.dat H2O MMR O3 MMR\n')
    swfname.write('rte_solver twostr\n')
    # swfname.write('sur_temperature '+str(surf_temp)+'\n')
    swfname.write('sza '+str(SZA[0])+'\n')

    if (cirrus_od > 0):
        swfname.write('ic_properties yang2013 interpolate \n')
        #swfname.write('ic_properties fu \n')
        swfname.write('ic_habit_yang2013 column_8elements smooth \n')
        swfname.write('ic_file 1D ICE_1H.dat \n')
        swfname.write('ic_modify tau set '+str(float(cirrus_od))+'\n')
        
    if (aero_tau > 0):
        swfname.write('aerosol_default\n')
        #swfname.write('aerosol_species_library OPAC\n')
        swfname.write('aerosol_file tau aero_ext.dat\n')
        swfname.write('aerosol_file ssa aero_ssa.dat\n')
        swfname.write('aerosol_file gg aero_gg.dat\n')

    #swfname.write('albedo '+str(surf_albedo[0])+' \n')
    swfname.write('albedo_file albedo_file.dat \n')
    swfname.write('day_of_year '+str(day_of_year) +
                  '\n')  # 265 for 9.21, 260 for 9.16
    swfname.write('zout TOA\n')
    #swfname.write('zout '+str(float(flight_alt))+'\n')
    swfname.write('output_user eglo eup \n')
    swfname.write('output_process sum\n')
    swfname.write('output_file  ' + swoutfile+'\n')
    swfname.write('quiet\n')
    swfname.close()
    
    
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
    return feature_type, feature_subtype, ice_water_phase   
                         
def read_text(fname):
    fo=open(fname,'r') #create file object
    #headstr=fo.readline()
    data=fo.readlines()
    #print (data)
    fo.close()  # close object
    return data
    
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
import os

print('start', datetime.now())

year='2008'

cldsat_dir='/data/keeling/a/yulanh/c/BW_backup/CloudSat/2B-FLXHR-LIDAR/'+year+'/'
fnames=sorted(glob.glob('/data/gdi/f/yulanh/CALIPSO/MLay/'+year+'/'+'CAL_LID_L2_05kmMLay-Standard-V4-20.2008-*'))

#fnames=['/data/gdi/f/yulanh/CALIPSO/MLay/2007/CAL_LID_L2_05kmMLay-Standard-V4-20.2007-01-01T11-02-31ZD.hdf']
#==== getting albedo =====
albedo_wv=['0.47','0.555','0.659','0.858','1.24','1.64','2.13']
albedo_day=np.array([   1,   17,   33,   49, 65,   81,   97,  113,  129,  145,              161,  177,  193,  209,  225,  241,  257,  273,  289,  305,            321,337,353])
albedo_day_str=['001','017','033','049','065','081','097','113','129','145',            '161','177','193','209','225','241','257','273','289','305',            '321','337','353']


albedo_dir='/data/keeling/a/yulanh/b/albedo/albedo/'


# ===============================================

for fname in fnames[0:95]:   
    
    # === get the julday ===
    month=fname[-20:-18]
    day = fname[-17:-15]

    day_of_year=date(int(year),int(month),int(day)).timetuple().tm_yday
    
    dnflag=fname[-6:-4]
    

    #===get inital time for each cloudsat fname ===
    if ((dnflag == 'ZD')):
        
        albedo_spec=np.zeros((360,720,len(albedo_wv)),'f')
        diff_day=abs(day_of_year-albedo_day)
        day_scp=diff_day.argmin()
        if (albedo_day[day_scp] > day_of_year):
            day_scp=day_scp-1

        for wvi in np.arange(len(albedo_wv)):
            albedo_fname=albedo_dir+'landsea_albedo_'+albedo_wv[wvi]+'/albedo_landsea_'+albedo_wv[wvi]+'_'+albedo_day_str[day_scp]+'.hdf'
            hdf=SD(albedo_fname)
            albedo_lat=hdf.select('lat')
            albedo_lat=albedo_lat[:]
            albedo_lon=hdf.select('lon')
            albedo_lon=albedo_lon[:]
            albedo = hdf.select('albedo_'+albedo_wv[wvi])
            albedo=albedo[:]/1000.0
            albedo_spec[:,:,wvi]=albedo
        
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
        feature_od=hdf.select('Feature_Optical_Depth_532')[:,:]
        SZA = hdf.select('Solar_Zenith_Angle')[:]

        matchdata_flag=0 # used for mark if rematch the data
        #aodind=np.where(aod > 0)[0]
        #=== initialize writing ==================
        wfname='simulate_calipsoonly_aerosol_radiation_'+fname[-25:-4]+'.h5'
        Ns=len(aod)
        calipso_cod=np.zeros((Ns),'f')
        calipso_aod=np.zeros((Ns),'f')
        waerosol_type=[''  for waerosol_type in range(len(aod))]
        wflux_wciup = np.zeros((Ns),'f')
        wflux_wcina_up=np.zeros((Ns),'f')
        wflux_nciup = np.zeros((Ns),'f')
        wflux_ncina_up=np.zeros((Ns),'f')
        
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
            tpfeature_od=feature_od[si,:]
            tpsza=SZA[si]
           
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
            cloud_flag=0
            aero_flag=0
            aero_topmost=0
            aero_baselow=0
            cirrus_baselow =0
            cirrus_topmost=0
            tpaerotype=0

            #print(tpcad_score,tplaybase,tplaytop,tpext,tpfeature)
            #=== write aerosol file ====
            ssa_fname=open('aero_ssa.dat','w')
            gg_fname =open('aero_gg.dat','w')
            ext_fname=open('aero_ext.dat','w')
            cirrus_fname=open('ICE_1H.dat','w')
            
            mean_ssa =0.0
            mean_gg =0.0
            Naero_lay=0
                       
            for hi in range(Nh-1,-1,-1):
                tplaybase_1=tplaybase[altind[hi]]
                tplaytop_1 = tplaytop[altind[hi]]
                #print('cadscore',tpcad_score[altind[hi]],lowcad_flag)
                if (tplaybase_1 < 0):
                    tplaybase_1=0
                
                if ((abs(tpcad_score[altind[hi]]) < 20) | (abs(tpcad_score[altind[hi]]) > 100)):
                    if (lowcad_flag == 0):
                        lowcad_flag=1
                   
                tpext_1=tpext[altind[hi]]
               
                if ((tpext_1 != 0) & (tpext_1 != 1) & (tpext_1 != 16) & (tpext_1 != 18)):
                    if (extqc_flag ==0):
                        extqc_flag=1
                #print('extqc',si,extqc_flag,tpext_1,tpext)
                if ((abs(tpcad_score[altind[hi]]) >= 20) & (abs(tpcad_score[altind[hi]]) <= 100) & (extqc_flag ==0)):
                    feature_type, feature_subtype,cloud_phase=vfm_feature_flag(tpfeature[altind[hi]])
                    # dealing with ice clouds
                    if ((cloud_phase == 1) | (cloud_phase == 3)):
                        cirrus_flag=1
                        if (cirrus_baselow ==0):
                            cirrus_baselow=tplaybase_1
                        if ((cirrus_baselow !=0) & (cirrus_baselow > tplaybase_1)):
                            cirrus_baselow=tplaybase_1
                        if (cirrus_topmost ==0):
                            cirrus_topmost=tplaytop_1
                        if ((cirrus_topmost != 0) & (cirrus_topmost < tplaytop_1)):
                            cirrus_topmost=tplaytop_1
                            
                    if (cloud_phase ==2):
                        liquid_flag=1
                        
                    # dealing with aerosols
                    aero_ssa=0.0
                    aero_gg =0.0
                
                    if ((feature_type == 3) | (feature_type ==4)):
                        aero_flag=1
                        if (aero_topmost ==0):
                            aero_topmost=tplaytop_1
                        if ((aero_topmost != 0) & (aero_topmost < tplaytop_1)):
                            aero_topmost=tplaytop_1
                        if ((aero_baselow == 0) & (tplaybase_1 >0)):
                            aero_baselow=tplaybase_1
                        if ((aero_baselow !=0) & (aero_baselow > tplaybase_1) & (tplaybase_1 >0)):
                            aero_baselow=tplaybase_1
                            
                            #add ssa and g same as Tabl2 2 in 2B-FLXHR-LIDAR document
                            
                        if (feature_subtype == 1):
                            aero_ssa=0.96
                            aero_gg =0.77
                            tpaerotype=tpaerotype+1
                        elif (feature_subtype ==2):
                            aero_ssa=0.93
                            aero_gg=0.74
                            tpaerotype=tpaerotype+10
                        elif (feature_subtype==3):
                            aero_ssa=0.99
                            aero_gg=0.76
                            tpaerotype=tpaerotype+100
                        elif (feature_subtype==4):
                            aero_ssa=0.99
                            aero_gg=0.74
                            tpaerotype=tpaerotype+1000
                        elif (feature_subtype==5):
                            aero_ssa=0.94
                            aero_gg = 0.70
                            tpaerotype=tpaerotype+10000
                        elif (feature_subtype ==6):
                            aero_ssa=0.83
                            aero_gg=0.66
                            tpaerotype=tpaerotype+100000
                        elif (feature_subtype ==7):
                            aero_ssa=0.94
                            aero_gg=0.71
                            tpaerotype=tpaerotype+1000000
                        else:
                            aero_ssa=0.95
                            aero_gg=0.70
                            tpaerotype=tpaerotype+10000000


                        mean_ssa=mean_ssa+aero_ssa
                        mean_gg =mean_gg +aero_gg
                        Naero_lay=Naero_lay+1
                        
                    if (feature_type == 2):
                        cloud_flag=1
            
                #==write aerosol file ===
                
            if (aero_flag ==1):
                ext_fname.write('{} {} \n'.format(aero_topmost,0.0))
                ssa_fname.write('{} {} \n'.format(aero_topmost,0.0))
                gg_fname.write('{} {} \n'.format(aero_topmost,0.0))


                ext_fname.write('{} {} \n'.format(aero_baselow,tpaod))#/(1000*(tplaytop_1-tplaybase_1))))
                ssa_fname.write('{} {} \n'.format(aero_baselow,mean_ssa/Naero_lay))
                gg_fname.write('{} {} \n'.format(aero_baselow,mean_gg/Naero_lay))    

                ext_fname.close()
                ssa_fname.close()
                gg_fname.close()
            if (cirrus_flag ==1):
                cirrus_fname.write('{} {} {} \n'.format(cirrus_topmost, 0.0, 0.0))
                cirrus_fname.write('{} {} {} \n'.format(cirrus_baselow,0.01816,20 ))
                cirrus_fname.close()
            #=== end deal with high confidence aero and cloud ======
          
            if ((extqc_flag ==0) & (lowcad_flag ==0) & (aero_flag ==1) & (cirrus_flag ==1) &                 (aero_topmost < cirrus_baselow) & (tpcod > 0) & (tpaod >0) & (liquid_flag ==0) &                (tpsza < 90)):
                
            #==== to write albedo file ===
                albedo_latscp=int((90-tplidlat-0.25)/0.5)
                albedo_lonscp=int((tplidlon+180-0.25)/0.5)
       
                tpspec_albedo=albedo_spec[albedo_latscp,albedo_lonscp,:]
                tpspec_albedo=np.insert(tpspec_albedo,0,0)
                tpspec_albedo=np.insert(tpspec_albedo,8,0)
                dat=np.array([[250,470,555,659,858,1240,1640,2130,5000],tpspec_albedo])
                dat=dat.T
                np.savetxt('albedo_file.dat',dat,delimiter=' ')
               
       
                #print(si,tpcod)
                #== write the script to run ===
                write_swinput(tpsza,tpcod, tpaod,'aero_swout.txt', day_of_year)
                os.system('uvspec <swinput1> test')
                #os.system('cat swinput1')
                #=== read output
                data=read_text('aero_swout.txt')
                splitcol=data[0].split(' ')
                simu_swout=np.zeros((2),'f')
                k=0
                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] !='\n':
                        simu_swout[k]=float(splitcol[j])
                        k=k+1
                    
                calipso_aod[si] = tpaod
                calipso_cod[si] = tpcod
                waerosol_type[si]=str(tpaerotype).zfill(8)
                wflux_wciup[si]=simu_swout[1]
                
                write_swinput(tpsza,tpcod, 0.0,'cinoaero_swout.txt', day_of_year)
                os.system('uvspec <swinput1> test')
                #=== read output
                data=read_text('cinoaero_swout.txt')
                splitcol=data[0].split(' ')
                simu_swout=np.zeros((2),'f')
                k=0
                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] !='\n':
                        simu_swout[k]=float(splitcol[j])
                        k=k+1
                    
                wflux_wcina_up[si]=simu_swout[1]
                
                write_swinput(tpsza,0.0, tpaod,'nociaero_swout.txt', day_of_year)
                os.system('uvspec <swinput1> test')
                #=== read output
                data=read_text('nociaero_swout.txt')
                splitcol=data[0].split(' ')
                simu_swout=np.zeros((2),'f')
                k=0
                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] !='\n':
                        simu_swout[k]=float(splitcol[j])
                        k=k+1
                    
                wflux_nciup[si]=simu_swout[1]
                
    
                write_swinput(tpsza,0.0,0.0,'nocinoaero_swout.txt', day_of_year)
                os.system('uvspec <swinput1> test')
                    #=== read output
                data=read_text('nocinoaero_swout.txt')
                splitcol=data[0].split(' ')
                simu_swout=np.zeros((2),'f')
                k=0
                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] !='\n':
                        simu_swout[k]=float(splitcol[j])
                        k=k+1
                    
                wflux_ncina_up[si]=simu_swout[1]
                    
                                  
                    
        hf=h5py.File(wfname,'w')
        hf.create_dataset('calipso_cod',data=calipso_cod)
        hf.create_dataset('calipso_aod',data=calipso_aod)
        hf.create_dataset('calipso_aeortype',data=waerosol_type)
        hf.create_dataset('simulation_ciaero_uptoa',data=wflux_wciup)
        hf.create_dataset('simulation_cinoaero_uptoa',data=wflux_wcina_up)
        hf.create_dataset('simulation_nociaero_uptoa',data=wflux_nciup)
        hf.create_dataset('simulation_nocinoaero_uptoa',data=wflux_ncina_up)
        hf.close()

print('finish', datetime.now())
