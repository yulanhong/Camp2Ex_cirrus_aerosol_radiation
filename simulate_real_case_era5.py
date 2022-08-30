import os
import matplotlib as mpl
import glob
import netCDF4 as nc
import h5py
from pyhdf.SD import SD, SDC
import pandas as pd
import math
import numpy as np
from datetime import datetime, date, timedelta
import matplotlib.pyplot as plt


def read_text(fname):
    fo = open(fname, 'r')  # create file object
    # headstr=fo.readline()
    data = fo.readlines()
    #print (data)
    fo.close()  # close object
    return data

# combine cod from beer mat
surf_albedo = 0.06
timerange = [25.5, 25.6]
era5_hour = 1
yyyymmdd = '20190827'
day_of_year = date(2019, 8, 27).timetuple().tm_yday
print('day of year', day_of_year)
#swlw_outflight_level = open(
#    'swlw_outflight_level_'+yyyymmdd+'_28_beermatcod.txt', 'w')
hy5file='simulation_'+yyyymmdd+'_25_beermatcod_albedo.h5'

era5_fn = 'CAMP2Ex_201908270006_ERA5.nc'

hsrl_fname = 'CAMP2EX-HSRL2_P3B_'+yyyymmdd+'_R0.h5'
fid = h5py.File(hsrl_fname, 'r')
hsrl_time = fid['Nav_Data/gps_time']
hsrl_lon = fid['Nav_Data/gps_lon']
hsrl_lat = fid['Nav_Data/gps_lat']
hsrl_alt = fid['Nav_Data/gps_alt']
hsrl_sza = fid['Nav_Data/solar_zenith']
hsrl_ext = fid['DataProducts/532_ext']
hsrl_aeroID = fid['DataProducts/Aerosol_ID']
hsrl_tauhi = fid['DataProducts/532_AOT_hi']
hsrl_tauac = fid['DataProducts/532_AOT_above_cloud']
hsrl_vert = fid['DataProducts/Altitude']
hsrl_press = fid['State/Pressure']
hsrl_ozone = fid['State/O3']
hsrl_temp = fid['State/Temperature']
hsrl_rh = fid['State/Relative_Humidity']
hsrl_sza = fid['Nav_Data/solar_zenith']
hsrl_Ans = fid['DataProducts/Angstrom_532_355']
hsrl_tauuv = fid['DataProducts/355_AOT_hi']
hsrl_aot_alts = fid['DataProducts/532_AOT_alts']

hsrl_time = hsrl_time[:]
hsrl_ext = hsrl_ext[:]
hsrl_aeroID = hsrl_aeroID[:]
hsrl_alt = hsrl_alt[:]/1000.0
hsrl_vert = hsrl_vert[:]/1000.0
hsrl_press = hsrl_press[:]
hsrl_ozone = hsrl_ozone[:]
hsrl_temp = hsrl_temp[:]
hsrl_rh = hsrl_rh[:]
hsrl_sza = hsrl_sza[:]
hsrl_Ans = hsrl_Ans[:]
hsrl_tauuv = hsrl_tauuv[:]
hsrl_aot_alts = hsrl_aot_alts[:]

ind = np.where((hsrl_time < 10.0))
ind = ind[0]
hsrl_time[ind] = hsrl_time[ind]+24
ind = np.where((hsrl_time > timerange[0]) & (hsrl_time < timerange[1]))
ind = ind[0]
hsrl_time1 = hsrl_time[ind]
hsrl_lon1 = hsrl_lon[ind]
hsrl_lat1 = hsrl_lat[ind]
hsrl_ext1 = hsrl_ext[ind, :]
hsrl_aeroID1 = hsrl_aeroID[ind, :]
hsrl_press1 = hsrl_press[ind, :]
hsrl_ozone1 = hsrl_ozone[ind, :]
hsrl_temp1 = hsrl_temp[ind, :]
hsrl_rh1 = hsrl_rh[ind, :]
hsrl_tauhi1 = hsrl_tauhi[ind]
hsrl_sza1 = hsrl_sza[ind]
hsrl_alt1 = hsrl_alt[ind]
hsrl_Ans1 = hsrl_Ans[ind, :]
hsrl_tauuv1 = hsrl_tauuv[ind]
hsrl_aot_alts1 = hsrl_aot_alts[ind, :]/1000.0
# print(hsrl_aot_alts1.shape)
dimension = hsrl_ext1.shape
Nx = dimension[0]
Ny = dimension[1]


# hsrl 2hz data
hsrl_2hz_fname = '/data/keeling/a/yulanh/c/camp2ex/HSRL_CTH_2Hz/' + \
    yyyymmdd+'_F1_special.h5'
fid = h5py.File(hsrl_2hz_fname, 'r')
CTH_2hz = fid['DataProducts/cloud_top_height']
CTH_2hz = CTH_2hz[:]/1000.0
time_2hz = fid['Nav_Data/gps_time']
time_2hz = time_2hz[:]
ind = np.where(time_2hz < 10)
time_2hz[ind] = time_2hz[ind]+24
ind = np.where((time_2hz > timerange[0]) & (time_2hz < timerange[1]))
ind = ind[0]
CTH_2hz_1 = CTH_2hz[ind]
time_2hz_1 = time_2hz[ind]

# === era5 ===
ds = nc.Dataset(era5_fn)
era_lat = ds['latitude'][:]
era_lon = ds['longitude'][:]
era_pres = ds['level'][:]
era_o3 = ds['o3'][:]
era_t = ds['t'][:]
era_r = ds['r'][:]

# === plot transmittance ===
spn_tran = 'spn_transmittance_'+yyyymmdd+'_hiac.txt'
data = read_text(spn_tran)
splitcol = data[0].split(' ')
Ncol = len(splitcol)-splitcol.count('')
Nrow = len(data)
dataT = np.zeros((Nrow, Ncol), 'f')
for i in range(Nrow):
    splitcol = data[i].split(' ')
    k = 0
    for j in range(len(splitcol)):
        if splitcol[j] != '' and splitcol[j] != '\n':
            dataT[i, k] = float(splitcol[j])
            k = k+1
spn_time = dataT[:, 0]
spn_tran = dataT[:, 1]
ind=np.where(spn_time < 10)
spn_time[ind]=spn_time[ind]+24

# === plot optical depth ===
spn_fname = 'spn_optical_depth_'+yyyymmdd+'_hiac.txt'
data = read_text(spn_fname)
splitcol = data[0].split(' ')
Ncol = len(splitcol)-splitcol.count('')
Nrow = len(data)
dataT = np.zeros((Nrow, Ncol), 'f')
for i in range(Nrow):
    splitcol = data[i].split(' ')
    k = 0
    for j in range(len(splitcol)):
        if splitcol[j] != '' and splitcol[j] != '\n':
            dataT[i, k] = float(splitcol[j])
            k = k+1
# spn_time=dataT[:,0]
spn_tau = dataT[:, 1]

# read matt cod
matt_fname = glob.glob('/data/keeling/a/yulanh/c/BW_backup/mydata/HSRL_out/CAMP2Ex/SPN_cirrus/SPN_COD_Matt/' + yyyymmdd+'*')
matt_fname = matt_fname[0]
fid = h5py.File(matt_fname, 'r')
SPNmatt_time = fid['UTC']
matt_tau = fid['COD_870nm']
matt_dr = fid['DR_870nm']

# === read original SPN data ====
#spn_irrad_fname = '/data/gdi/f/yulanh/Simulate_camp2ex/SPN/CAMP2EX-SPNS_P3B_'+yyyymmdd+'_R0.h5'
#fid = h5py.File(spn_irrad_fname, 'r')
#spn_irad_alt = fid['Data/Altitude']
#spn_irad_time = fid['Data/HR_UTC']
#spn_irad_diff = fid['Data/Irradiance_diffuse']
#spn_irad_dire = fid['Data/Irradiance_direct']
#spn_irad_tot = fid['Data/Irradiance_total']
#spn_wave = fid['Data/Wavelengths']


# == plot optical depth ====
rsp_fname = 'CAMP2EX-RSP1-AER_P3B_'+yyyymmdd+'_R2.h5'
fid = h5py.File(rsp_fname, 'r')
rsp_tau = fid['aerosol_tau_total_555']
rsp_ssa = fid['aerosol_ssa_total_555']
rsp_time = fid['rsp_time']
rsp_time = rsp_time[:]
rsp_tau = rsp_tau[:]
rsp_ssa = rsp_ssa[:]
ind = np.where(rsp_time < 10)
rsp_time[ind] = rsp_time[ind]+24
ind = np.where((rsp_time > timerange[0]) & (rsp_time < timerange[1]))
tprsp = rsp_ssa[ind]

# == read plane attitude ====
metdir = '/data/keeling/a/yulanh/c/BW_backup/mydata/HSRL_out/CAMP2Ex/Metnav/'
met_file = metdir+'CAMP2EX-MetNav_P3B_'+yyyymmdd+'_R0.ict'
df = pd.read_csv(met_file, delimiter=',', header=73)
met_time = df['Time_Start']/3600.0
met_hgt = df['GPS_Altitude']/1000.0
met_sza = df['Solar_Zenith_Angle']
met_roll = df['Roll_Angle']
met_pitch = df['Pitch_Angle']
met_sam = df['Sun_Azimuth']
met_lat = df['Latitude']
met_lon = df['Longitude']
met_head = df['True_Heading']
met_drift = df['Drift_Angle']

tind = np.where(met_time < 10)[0]
met_time[tind] = met_time[tind]+24.0

# ===== plot BBR data
#BBR_fname = 'CAMP2EX-BBR_P3B_'+yyyymmdd+'_R2.ict'
#fo = open(BBR_fname, 'r')  # create file object
#data = fo.readlines()[54:]
#fo.close()
#splitcol = data[0].split(',')
#Ncol = len(splitcol)-splitcol.count('')
#Nrow = len(data)
#dataT = np.zeros((Nrow, Ncol), 'f')
#for i in range(Nrow):
#    splitcol = data[i].split(',')
#    k = 0
#    for j in range(len(splitcol)):
#        if splitcol[j] != '' and splitcol[j] != '\n':
#            dataT[i, k] = float(splitcol[j])
#            k = k+1


#BBR_time = dataT[:, 0]/3600.0
# BBR_dnlw=dataT[:,1]
# BBR_uplw=dataT[:,2]
#BBR_dnsw = dataT[:, 3]
#BBR_upsw = dataT[:, 4]
#BBRspn1_dnsw = dataT[:, 5]
#BBRspn1_upsw = dataT[:, 6]
#ind = np.where(BBR_time < 10)
#BBR_time[ind] = BBR_time[ind]+24
#ind = np.where(BBR_dnsw < 0.0)[0]
#BBR_dnsw[ind] = 'nan'
#ind = np.where(BBR_upsw < 0.0)[0]
#BBR_upsw[ind] = 'nan'

# read SSFR and get albedo
wv_len = [415, 440, 500, 550, 675, 870, 990, 1020, 1064, 1250, 1650, 2100]
SSFR_dir = '/data/gdi/f/yulanh/Simulate_camp2ex/SSFR/'
SSFR_fname = SSFR_dir+'CAMP2EX-SSFR-Partial_P3B_'+yyyymmdd+'_R0.ict'

df1 = pd.read_csv(SSFR_fname, delimiter=',', header=70)
data = df1.to_numpy()

SSFR_time = data[:, 0]/3600.0

ind = np.where(SSFR_time < 10.0)[0]
SSFR_time[ind] = SSFR_time[ind]+24

dn_wv = np.zeros((len(SSFR_time), len(wv_len)), 'f')
up_wv = np.zeros((len(SSFR_time), len(wv_len)), 'f')

# dn_wv[:,0]=df['DN440']
dn_wv[:, 0] = data[:, 1]  # 415
dn_wv[:, 1] = data[:, 2]  # 440
dn_wv[:, 2] = data[:, 3]  # 500
dn_wv[:, 3] = data[:, 4]  # 550
dn_wv[:, 4] = data[:, 5]  # 675
dn_wv[:, 5] = data[:, 6]  # 870
dn_wv[:, 6] = data[:, 7]  # 990
dn_wv[:, 7] = data[:, 8]  # 1020
dn_wv[:, 8] = data[:, 9]  # 1064
dn_wv[:, 9] = data[:, 10]  # 1250
dn_wv[:, 10] = data[:, 11]  # 1650
dn_wv[:, 11] = data[:, 12]  # 2100

up_wv[:, 0] = data[:, 13]  # 415
up_wv[:, 1] = data[:, 14]  # 440
up_wv[:, 2] = data[:, 15]  # 500
up_wv[:, 3] = data[:, 16]  # 550
up_wv[:, 4] = data[:, 17]  # 675
up_wv[:, 5] = data[:, 18]  # 870
up_wv[:, 6] = data[:, 19]  # 990
up_wv[:, 7] = data[:, 20]  # 1020
up_wv[:, 8] = data[:, 21]  # 1064
up_wv[:, 9] = data[:, 22]  # 1250
up_wv[:, 10] = data[:, 23]  # 1650
up_wv[:, 11] = data[:, 24]  # 2100

ssfr_lon = data[:, 25]
ssfr_lat = data[:, 26]

dn_wv[dn_wv == -999.0] = 'nan'

up_wv[dn_wv == -999.0] = 'nan'

# SSFR whole spectral
SSFR_fname1 = SSFR_dir+'CAMP2EX-SSFR_P3B_'+yyyymmdd+'_R0.h5'
#fid = h5py.File(SSFR_fname1, 'r')
#ssfr_nadflux = fid['nad_flux'][...]
#ssfr_nadwvl = fid['nad_wvl'][...]
#ssfr_zenflux = fid['zen_flux'][...]
#ssfr_zenwvl = fid['zen_wvl'][...]
#ssfr_time1 = fid['tmhr'][...]

# =============================


def write_swinput(SZA, ci_tau, aero_tau, ssa, gg, flight_alt, swoutfile, surf_albedo, day_of_year):
    swfname = open('swinput1', 'w')
    swfname.write(
        'atmosphere_file /data/keeling/a/yulanh/c/software_install/libradtran/2.0.2/share/libRadtran/data/atmmod/afglt.dat \n')
    swfname.write(
        'source solar /data/keeling/a/yulanh/c/software_install/libradtran/2.0.2/share/libRadtran/data/solar_flux/kurudz_1.0nm.dat \n')
    swfname.write('wavelength 250 3600\n')
    swfname.write('mol_abs_param reptran coarse\n')
    swfname.write('radiosonde atm_file_era.dat H2O RH O3 MMR\n')
    swfname.write('rte_solver disort\n')
    # swfname.write('sur_temperature '+str(surf_temp)+'\n')
    swfname.write('sza '+str(SZA[0])+'\n')

    if (ci_tau > 0):
        swfname.write('ic_properties yang2013 interpolate \n')
        #swfname.write('ic_properties fu \n')
        swfname.write('ic_habit_yang2013 column_8elements smooth \n')
        swfname.write('ic_file 1D ICE_1H.dat\n')
        swfname.write('ic_modify tau set '+str(float(ci_tau))+'\n')

    if (aero_tau > 0):
        swfname.write('aerosol_default\n')
        swfname.write('aerosol_species_library OPAC\n')
        swfname.write('aerosol_file tau aero_file_full.dat\n')
        swfname.write('aerosol_modify ssa set '+str(ssa)+'\n')
        swfname.write('aerosol_modify gg set '+str(gg)+'\n')

   # swfname.write('albedo '+str(surf_albedo)+' \n')
    swfname.write('albedo_file albedo_file.dat \n')
    swfname.write('day_of_year '+str(day_of_year) +
                  '\n')  # 265 for 9.21, 260 for 9.16
    swfname.write('zout 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20\n')
    #swfname.write('zout '+str(float(flight_alt))+'\n')
    swfname.write('output_user lambda zout sza albedo eglo eup edn heat\n')
    swfname.write('output_process integrate\n')
    swfname.write('output_file  ' + swoutfile+'\n')
    swfname.write('quiet\n')
    swfname.close()
    os.system('uvspec <swinput1> test')

# apply for simulation
# data needed HSRL:hsrl_ext1, hsrl_aeroID1, hsrl_vert,hsrl_time1,CTH_2hz_1,time_2hz_1, hsrl_press1,hsrl_ozone1,hsrl_rh1,hsrl_temp1
# BBR_uplw,BBR_dnlw,BBR_upsw,BBR_dnsw,BBR_time
# RSP_time,RSP_ssa
# SPN tau

print('start', datetime.now())
Ns = len(hsrl_time1)
hsrl_vert = np.squeeze(hsrl_vert)
hsrl_time1 = np.squeeze(hsrl_time1)
Nh = hsrl_vert.shape
Nh = Nh[0]
Nh_p = len(era_pres)
# one standard atmosphere

hsrl_ext1[np.isnan(hsrl_ext1)] = 0.0
ind = np.where(abs(hsrl_vert) == abs(hsrl_vert).min())
ind = ind[0]
surf_bin_scp = ind
# ==== obtain the mean ssa
ind = np.where((rsp_time > timerange[0]) & (rsp_time < timerange[1]))[0]
ssa = 0.0
if (len(ind) > 0):
    tprsp = rsp_ssa[ind]
    ssa = np.nanmean(tprsp)

if (ssa == 0.0):
    ssa = 0.95

gg = 0.55  # from

#= = obtain albedo database ===
#=== obtain surface albedo===

albedo_wv=['0.47','0.555','0.659','0.858','1.24','1.64','2.13']
albedo_day=np.array([   1,   17,   33,   49, 65,   81,   97,  113,  129,  145,\
              161,  177,  193,  209,  225,  241,  257,  273,  289,  305,\
            321,337,353])
albedo_day_str=['001','017','033','049','065','081','097','113','129','145',\
            '161','177','193','209','225','241','257','273','289','305',\
            '321','337','353']
albedo_spec=np.zeros((360,720,len(albedo_wv)),'f')

albedo_dir='/data/keeling/a/yulanh/b/albedo/albedo/'

diff_day=abs(day_of_year-albedo_day)
day_scp=diff_day.argmin()
if (albedo_day[day_scp] > day_of_year):
    day_scp=day_scp-1

for wvi in np.arange(len(albedo_wv)):
    albedo_fname=albedo_dir+'landsea_albedo_'+albedo_wv[wvi]+'/albedo_landsea_'+albedo_wv[wvi]+'_'+albedo_day_str[day_scp]+'.hdf'
    print(albedo_fname)
    hdf=SD(albedo_fname)
    albedo_lat=hdf.select('lat')
    albedo_lat=albedo_lat[:]
    albedo_lon=hdf.select('lon')
    albedo_lon=albedo_lon[:]
    albedo = hdf.select('albedo_'+albedo_wv[wvi])
    albedo=albedo[:]/1000.0
    albedo_spec[:,:,wvi]=albedo
    

swup_output = np.zeros((Ns,41,4), 'f')
swdn_output = np.zeros((Ns,41,4), 'f') #0-ice_aero,1-ice_noaero,2-noice_aero,3-noice_noaero
swheat_output = np.zeros((Ns,41,4), 'f')
cod_tau_output = np.zeros((Ns),'f')
aero_tau_output=np.zeros((Ns),'f')
cod_tau_matt_output=np.zeros((Ns),'f')

for si in np.arange(Ns):

    # write input files
    # aerosol tau file, atmosphere file(pressure, temperature, water vapor, ozone)
    # print(hsrl_tauhi1[si],np.nansum(hsrl_ext1[si,:]*0.015))
    # exclude the situation when high resolution hsrl CTH has valid values
    tphsrl_time = hsrl_time1[si]
    tphsrl_lat = hsrl_lat1[si]
    tphsrl_lon = hsrl_lon1[si]
    # ===== obtain reanalysis grid ====
    time_scp = era5_hour
    lonind = (tphsrl_lon-era_lon[0])/0.25
    lonscp = round(lonind[0])
    latind = (era_lat[0]-tphsrl_lat)/0.25
    latscp = round(latind[0])
    tpera_temp = era_t[time_scp, :, latscp, lonscp]
    tpera_ozone = era_o3[time_scp, :, latscp, lonscp]
    tpera_rh = era_r[time_scp, :, latscp, lonscp]

    tpera_rh[tpera_rh < 0] = 0.0
    # print(tphsrl_lon,lonind,era_lon[lonscp],tphsrl_lat,latind,era_lat[latscp])

    ind = np.where((time_2hz_1 > (tphsrl_time - 0.00138)) & (time_2hz_1 < (tphsrl_time + 0.00139)))
    ind = ind[0]
    # print(tphsrl_time,time_2hz_1[ind])
    # print(np.nanmean(CTH_2hz_1[ind]),'mean',CTH_2hz_1[ind])
    ind1 = np.where(CTH_2hz_1[ind] > 0)
    ind1 = ind1[0]
    CTH_num = len(ind1)
    # surf_temp=float(hsrl_temp1[si,0])
# == check roll and pitch angle
    ind = np.where((met_time > (tphsrl_time - 0.00138)) & (met_time < (tphsrl_time + 0.00139)))[0]
    tproll_angle = np.nanmean(met_roll[ind])
    tppitch_angle = np.nanmean(met_pitch[ind])

    # CTH_num is used to control whether there are low clouds
    if ((hsrl_tauhi1[si] > 0) & (CTH_num <= 5) & (np.abs(tproll_angle) < 3.0) & (np.abs(tppitch_angle) < 3.0)):
        
        aero_fname = open('aero_file_full.dat', 'w')
        # down scale the resolution of aerosol, from surface to 6 km
        if (math.isnan(hsrl_aot_alts1[si,0])):
            hsrl_aot_alts1[si,0]=0.0
        if (math.isnan(hsrl_aot_alts1[si,1])):
            hsrl_aot_alts1[si,1]=6.5

        for hi in np.arange(14):
            hgt_top = 6.5-hi*0.5
            hgt_base = hgt_top-0.5
            ind = np.where((hsrl_vert > hgt_base) & (hsrl_vert <= hgt_top))
            ind = ind[0]
            tpmean_hgt = np.nanmean(hsrl_vert[ind])
            tpmean_aero_tau = np.nansum(hsrl_ext1[si, ind]*0.0149)
            if (tpmean_aero_tau < 0):
                tpmean_aero_tau=0.0
            
            if ((hgt_base > hsrl_aot_alts1[si, 0]) & (hgt_top < hsrl_aot_alts1[si, 1])):
                #aero_fname.write('{} {}\n'.format(tpmean_hgt, tpmean_aero_tau))
                aero_fname.write(str(tpmean_hgt)+'  '+str(tpmean_aero_tau)+'\n')

        aero_fname.close()
        


        dat=np.array([era_pres,tpera_temp,tpera_rh,tpera_ozone])
        dat=dat.T
        np.savetxt('atm_file_era.dat',dat,delimiter=' ')

        # add albedo file
        albedo_latscp=int((90-tphsrl_lat-0.25)/0.5)
        albedo_lonscp=int((tphsrl_lon+180-0.25)/0.5)
        tpspec_albedo=albedo_spec[albedo_latscp,albedo_lonscp,:]
        tpspec_albedo=np.insert(tpspec_albedo,0,0)
        tpspec_albedo=np.insert(tpspec_albedo,8,0)
        dat=np.array([[250,470,555,659,858,1240,1640,2130,4000],tpspec_albedo])
        dat=dat.T
        np.savetxt('albedo_file.dat',dat,delimiter=' ')
 
    
        # albedo_fname.close
        # ==== finish writing atmosphere and aerosol profile
        # == surface temperature
        # surf_temp=float(hsrl_temp1[si,surf_bin_scp])
        # === solar zenith angle ====

        # === SPN tau ====
        # print(spn_time.shape,hsrl_time1.shape)
        ind=np.where((SPNmatt_time > (hsrl_time1[si]-0.0014)) & (SPNmatt_time < (hsrl_time1[si]+0.0014)))
        ind=ind[0]
        tpspn_mean1=np.nanmean(matt_tau[ind])
        #tpspn_dr = np.nanmean(matt_dr[ind])

        ind = np.where((spn_time > (hsrl_time1[si]-0.0014)) & (spn_time < (hsrl_time1[si]+0.0014)))
        tpspn_mean2 = np.nanmean(spn_tau[ind[0]])

        tpspn_mean=0.0 # combination of two cod
        if (tpspn_mean2 > 0) :
            tpspn_mean = tpspn_mean2
        if ((tpspn_mean2 <= 0) & (tpspn_mean1 > 0)):
            tpspn_mean = tpspn_mean1
            cod_tau_matt_output[si]=tpspn_mean

        #tpspn_dr = 0
        tpspn_time = np.nanmean(spn_time[ind[0]])
        cod_tau_output[si] =tpspn_mean
        aero_tau_output[si]=hsrl_tauhi1[si]
        print(si,Ns,tpspn_mean,hsrl_tauhi1[si])

        if ((tpspn_mean > 0.0) & (hsrl_tauhi1[si] > 0)):
            swoutfile1 = 'swirradiance_out_icaero'
            write_swinput(hsrl_sza1[si], tpspn_mean, hsrl_tauhi1[si], ssa, gg, hsrl_alt1[si], swoutfile1, surf_albedo, day_of_year)

            # === record outputs =====
            data = read_text(swoutfile1)
            splitcol = data[0].split(' ')
            Ncol = len(splitcol)-splitcol.count('')
            Nrow = len(data)
            swic_dataT = np.zeros((Nrow, Ncol), 'f')
            for i in range(Nrow):
                splitcol = data[i].split(' ')
                k = 0

                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] != '\n':
                        swic_dataT[i, k] = float(splitcol[j])
                        k = k+1

            swic_dataT = np.squeeze(swic_dataT)
            input_albedo = swic_dataT[2]
            simuwv = swic_dataT[:, 0]
            tpswdn_output = swic_dataT[:, 4]/1000.0
            tpswup_output = swic_dataT[:, 5]/1000.0
            tpswht_output = swic_dataT[:, 7]/1000.0

#            print('a',tpswdn_output.shape)
            swup_output[si,:,0] = tpswup_output
            swdn_output[si,:,0] = tpswdn_output
            swheat_output[si,:,0]=tpswht_output
        
            swoutfile2 = 'swirradiance_out_noicaero'
            write_swinput(hsrl_sza1[si], 0.0, hsrl_tauhi1[si], ssa, gg, hsrl_alt1[si], swoutfile2, surf_albedo, day_of_year)
            # === record outputs =====
            data = read_text(swoutfile2)
            splitcol = data[0].split(' ')
            Ncol = len(splitcol)-splitcol.count('')
            Nrow = len(data)
            swic_dataT = np.zeros((Nrow, Ncol), 'f')
            for i in range(Nrow):
                splitcol = data[i].split(' ')
                k = 0

                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] != '\n':
                        swic_dataT[i, k] = float(splitcol[j])
                        k = k+1

            swic_dataT = np.squeeze(swic_dataT)
            input_albedo = swic_dataT[2]
            simuwv = swic_dataT[:, 0]
            tpswdn_output = swic_dataT[:, 4]/1000.0
            tpswup_output = swic_dataT[:, 5]/1000.0
            tpswht_output = swic_dataT[:, 7]/1000.0

            swup_output[si,:,1] = tpswup_output
            swdn_output[si,:,1] = tpswdn_output
            swheat_output[si,:,1]=tpswht_output

            swoutfile3 = 'swirradiance_out_icnoaero'
            write_swinput(hsrl_sza1[si], tpspn_mean, 0.0, ssa, gg, hsrl_alt1[si], swoutfile3, surf_albedo, day_of_year)
            # === record outputs =====
            data = read_text(swoutfile3)
            splitcol = data[0].split(' ')
            Ncol = len(splitcol)-splitcol.count('')
            Nrow = len(data)
            swic_dataT = np.zeros((Nrow, Ncol), 'f')
            for i in range(Nrow):
                splitcol = data[i].split(' ')
                k = 0

                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] != '\n':
                        swic_dataT[i, k] = float(splitcol[j])
                        k = k+1

            swic_dataT = np.squeeze(swic_dataT)
            input_albedo = swic_dataT[2]
            simuwv = swic_dataT[:, 0]
            tpswdn_output = swic_dataT[:, 4]/1000.0
            tpswup_output = swic_dataT[:, 5]/1000.0
            tpswht_output = swic_dataT[:, 7]/1000.0

            swup_output[si,:,2] = tpswup_output
            swdn_output[si,:,2] = tpswdn_output
            swheat_output[si,:,2]=tpswht_output

            swoutfile4 = 'swirradiance_out_noicnoaero'
            write_swinput(hsrl_sza1[si],0.0, 0.0, ssa, gg, hsrl_alt1[si], swoutfile4, surf_albedo, day_of_year)
            # === record outputs =====
            data = read_text(swoutfile4)
            splitcol = data[0].split(' ')
            Ncol = len(splitcol)-splitcol.count('')
            Nrow = len(data)
            swic_dataT = np.zeros((Nrow, Ncol), 'f')
            for i in range(Nrow):
                splitcol = data[i].split(' ')
                k = 0
                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] != '\n':
                        swic_dataT[i, k] = float(splitcol[j])
                        k = k+1

            swic_dataT = np.squeeze(swic_dataT)
            input_albedo = swic_dataT[2]
            simuwv = swic_dataT[:, 0]
            tpswdn_output = swic_dataT[:, 4]/1000.0
            tpswup_output = swic_dataT[:, 5]/1000.0
            tpswht_output = swic_dataT[:, 7]/1000.0

            swup_output[si,:,3] = tpswup_output
            swdn_output[si,:,3] = tpswdn_output
            swheat_output[si,:,3]=tpswht_output
           
            os.system('rm '+swoutfile1)
            os.system('rm '+swoutfile2)
            os.system('rm '+swoutfile3)
            os.system('rm '+swoutfile4)

        if ((tpspn_mean == 0.0) & (hsrl_tauhi1[si] > 0)):
            swoutfile1 = 'swirradiance_out_noicaero'
            write_swinput(hsrl_sza1[si], 0.0, hsrl_tauhi1[si], ssa, gg, hsrl_alt1[si], swoutfile1, surf_albedo, day_of_year)

            # === record outputs =====
            data = read_text(swoutfile1)
            splitcol = data[0].split(' ')
            Ncol = len(splitcol)-splitcol.count('')
            Nrow = len(data)
            swic_dataT = np.zeros((Nrow, Ncol), 'f')
            for i in range(Nrow):
                splitcol = data[i].split(' ')
                k = 0

                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] != '\n':
                        swic_dataT[i, k] = float(splitcol[j])
                        k = k+1

            swic_dataT = np.squeeze(swic_dataT)
            input_albedo = swic_dataT[2]
            simuwv = swic_dataT[:, 0]
            tpswdn_output = swic_dataT[:, 4]/1000.0
            tpswup_output = swic_dataT[:, 5]/1000.0
            tpswht_output = swic_dataT[:, 7]/1000.0

            swup_output[si,:,1] = tpswup_output
            swdn_output[si,:,1] = tpswdn_output
            swheat_output[si,:,1]=tpswht_output


            swoutfile2 = 'swirradiance_out_noicnoaero'
            write_swinput(hsrl_sza1[si], 0.0, 0.0, ssa, gg, hsrl_alt1[si], swoutfile2, surf_albedo, day_of_year)
            # === record outputs =====
            data = read_text(swoutfile2)
            splitcol = data[0].split(' ')
            Ncol = len(splitcol)-splitcol.count('')
            Nrow = len(data)
            swic_dataT = np.zeros((Nrow, Ncol), 'f')
            for i in range(Nrow):
                splitcol = data[i].split(' ')
                k = 0

                for j in range(len(splitcol)):
                    if len(splitcol[j]) != 0 and splitcol[j] != '\n':
                        swic_dataT[i, k] = float(splitcol[j])
                        k = k+1

            swic_dataT = np.squeeze(swic_dataT)
            input_albedo = swic_dataT[2]
            simuwv = swic_dataT[:, 0]
            tpswdn_output = swic_dataT[:, 4]/1000.0
            tpswup_output = swic_dataT[:, 5]/1000.0
            tpswht_output = swic_dataT[:, 7]/1000.0
 
            swup_output[si,:,3] = tpswup_output
            swdn_output[si,:,3] = tpswdn_output
            swheat_output[si,:,3]=tpswht_output

            os.system('rm '+swoutfile1)
            os.system('rm '+swoutfile2)

hf=h5py.File(hy5file,'w')
hf.create_dataset('time',data=hsrl_time1)
hf.create_dataset('swup_output',data=swup_output)
hf.create_dataset('swdn_output',data=swdn_output)
hf.create_dataset('swht_output',data=swheat_output)
hf.create_dataset('cod_tau_output',data=cod_tau_output)
hf.create_dataset('cod_tau_matt_output',data=cod_tau_matt_output)
hf.create_dataset('aero_tau_output',data=aero_tau_output)
hf.close()
print('finish', datetime.now())
