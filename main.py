from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from numpy import zeros,where,array,mean,std,ones,nan,linspace,meshgrid,isnan,log
from numpy.ma import is_masked
from osgeo import gdal,gdalconst,osr,ogr
import calendar
import csv
import gzip
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy.ma as ma
import os
import rasterio

def point_to_cell(point_x, point_y, cellx, celly, xmin, ymax):
    col = int((point_x - xmin) / cellx)
    row = int((point_y - ymax) / -celly)
    return row,col


def cell_to_coord(col, row, cellx, celly, xmin, ymax):
	lon = cellx*col + xmin + cellx/2.0
	lat = ymax-celly*row - celly/2.0
	return lat,lon




###famine data, directly from https://s3.amazonaws.com/shapefiles.fews.net/ALL_HFIC.zip https://fews.net/data/acute-food-insecurity
os.mkdir('./OUTPUT/rasters')

input_files = [i for i in os.listdir('./INPUT/ALL_HFIC/all_lay') if i[-6:]=='CS.shp']
for input_file in input_files:
    InputVector = './INPUT/ALL_HFIC/all_lay/'+input_file
    OutputImage = './OUTPUT/rasters/'+input_file[:-3]+'tif'
    RefImage = './INPUT/ref_raster.tif'
    gdalformat = 'GTiff'
    datatype = gdal.GDT_Byte
    ##########################################################
    # Get projection info from reference image
    Image = gdal.Open(RefImage, gdal.GA_ReadOnly)
    # Open Shapefile
    Shapefile = ogr.Open(InputVector)
    Shapefile_layer = Shapefile.GetLayer()
    # Rasterise
    Output = gdal.GetDriverByName(gdalformat).Create(OutputImage, Image.RasterXSize, Image.RasterYSize, 2, datatype, options=['COMPRESS=DEFLATE'])
    Output.SetProjection(Image.GetProjectionRef())
    Output.SetGeoTransform(Image.GetGeoTransform())
    # Write data to band 1
    Band = Output.GetRasterBand(1)
    Band.SetNoDataValue(0)
    gdal.RasterizeLayer(Output, [1], Shapefile_layer, options=['ATTRIBUTE=CS'])#burn_values=[burnVal])
    Band = Output.GetRasterBand(1)
    Band.SetNoDataValue(0)
    gdal.RasterizeLayer(Output, [2], Shapefile_layer, options=['ATTRIBUTE=HA0'])#burn_values=[burnVal])
    # Close datasets
    Band = None
    Output = None
    Image = None
    Shapefile = None
    print (input_file)




####prepare GDP data from https://zenodo.org/records/7898409
###both at 0.5 and 1 degree
fff = os.listdir('./INPUT/GDP/025d')
target_dirs = ['./OUTPUT/GDP_05_deg/','./OUTPUT/GDP_1_deg/']
referencefiles = ['./INPUT/ref_raster.tif','./INPUT/ref_raster_1deg.tif']
for i in range(2):
    os.mkdir(target_dirs[i])
    referencefile = referencefiles[i]
    reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
    referenceProj = reference.GetProjection()
    referenceTrans = reference.GetGeoTransform()
    x = reference.RasterXSize
    y = reference.RasterYSize
    for ff in fff:
        inputfile = './INPUT/GDP/025d/' + ff
        input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
        inputProj = input.GetProjection()
        inputTrans = input.GetGeoTransform()
        outputfile = target_dirs[i] + ff  # Path to output file
        driver = gdal.GetDriverByName('GTiff')
        bandreference = input.GetRasterBand(1)
        output = driver.Create(outputfile, x, y, 1, bandreference.DataType, options=['COMPRESS=DEFLATE'])
        output.SetGeoTransform(referenceTrans)
        output.SetProjection(referenceProj)
        gdal.ReprojectImage(input, output, inputProj, referenceProj, gdalconst.GRA_Sum)
        del output
        print(ff)


###climate data historical
###temperature
dataset = Dataset("./INPUT/NOOA_climate/temp/air.mon.mean.nc")
#1948/01 to 2023/07
os.mkdir('./OUTPUT/NOOA_climate_temp_rasters')
meta = rasterio.open("./INPUT/land_05_wgs84.tif").meta
meta['nodata'] = -99999
meta['dtype'] = 'float32'
dataset = Dataset("./INPUT/NOOA_climate/temp/air.mon.mean.nc")
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
var_layer = dataset.variables['air'][:]
n_mon = len(var_layer)
for mon in range(n_mon):
    year = 1948+int(mon/12)
    month = mon%12 + 1
    m = zeros([len(lat),len(lon)])
    for i in range(len(lat)):
        for j in range(len(lon)):
            if j >= 360:
                x = j-360
            else:
                x = j+360
            if not(is_masked(var_layer[mon][i][j])):
                m[i][x] = var_layer[mon][i][j]
            else:
                m[i][x] = -99999
    out = rasterio.open('./OUTPUT/NOOA_climate_temp_rasters/'+str(year) +'_'+ str(month)+ '.tif', 'w', **meta)
    out.write(m.astype(rasterio.float64), 1)
    out.close()
    print (str(year) +'_'+ str(month))


####PREC
#1891/01 to 2019/12
os.mkdir('./OUTPUT/NOOA_climate_prec_rasters')
meta = rasterio.open("./INPUT/land_05_wgs84.tif").meta
meta['nodata'] = -99999
meta['dtype'] = 'float32'
dataset = Dataset("./INPUT/NOOA_climate/prec/precip.mon.total.0.5x0.5.v2020.nc")
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
var_layer = dataset.variables['precip'][:]
n_mon = len(var_layer)

start_year = 1980
start_month = (start_year-1891)*12
for mon in range(start_month,n_mon):
    year = 1891+int(mon/12)
    month = mon%12 + 1
    m = zeros([len(lat),len(lon)])
    for i in range(len(lat)):
        for j in range(len(lon)):
            if j >= 360:
                x = j-360
            else:
                x = j+360
            if not(is_masked(var_layer[mon][i][j])):
                m[i][x] = var_layer[mon][i][j]
            else:
                m[i][x] = -99999
    out = rasterio.open('./OUTPUT/NOOA_climate_prec_rasters/'+str(year) +'_'+ str(month)+ '.tif', 'w', **meta)
    out.write(m.astype(rasterio.float64), 1)
    out.close()
    print (str(year) +'_'+ str(month))


###Get precipitation data from 2020 to now (need to run outside the JRC network, of course is blocked here).
###from: https://opendata.dwd.de/climate_environment/GPCC/monitoring_v2022/

os.mkdir('./OUTPUT/NOOA_climate_prec_2020_2023_rasters')

meta = rasterio.open("./INPUT/land_1_wgs84.tif").meta
meta['nodata'] = -99999
meta['dtype'] = 'float32'

for year in range(2020,2024):
    for month in range(1,13):
        if year<2023 or month<6:
            local_file = "./INPUT/NOOA_climate/prec/2020_2023/monitoring_v2022_10_"+str(year)+"_"+"{:02d}".format(month)+".nc.gz"
            dataset = Dataset('dummy', mode='r', memory=gzip.open(local_file).read())
            lat = dataset.variables['lat'][:]
            lon = dataset.variables['lon'][:]
            var_layer = dataset.variables['p'][0][:]
            m = zeros([len(lat),len(lon)])
            for i in range(len(lat)):
                for j in range(len(lon)):
                    # if j >= 180:
                    #     x = j-180
                    # else:
                    #     x = j+180
                    x = j
                    if not(is_masked(var_layer[i][j])):
                        m[i][x] = var_layer[i][j]
                    else:
                        m[i][x] = -99999
            out = rasterio.open('./OUTPUT/NOOA_climate_prec_2020_2023_rasters/'+str(year) +'_'+ str(month)+ '.tif', 'w', **meta)
            out.write(m.astype(rasterio.float64), 1)
            out.close()
            print (str(year) +'_'+ str(month))


###downscale recent prec data (2020-2023) to 0.5 x 0.5

referencefile = './INPUT/ref_raster_1deg.tif'#Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize

#prec
fff = [i for i in os.listdir('./OUTPUT/NOOA_climate_prec_rasters') if i[-4:]=='.tif']
os.mkdir('./OUTPUT/NOOA_climate_prec_rasters_1/')

for ff in fff:
    inputfile = './OUTPUT/NOOA_climate_prec_rasters/'+ff
    input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
    inputProj = input.GetProjection()
    inputTrans = input.GetGeoTransform()
    outputfile = './OUTPUT/NOOA_climate_prec_rasters_1/'+ff #Path to output file
    driver= gdal.GetDriverByName('GTiff')
    bandreference = input.GetRasterBand(1)
    output = driver.Create(outputfile,x,y,1,bandreference.DataType,options=['COMPRESS=DEFLATE'])
    output.SetGeoTransform(referenceTrans)
    output.SetProjection(referenceProj)
    gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Bilinear)
    del output
    print (ff)
    

###same with population data from 2000 to 2022 for the historical scenario
for year in range(2000,2023):
    meta = rasterio.open("./INPUT/pop_density/2000-2022/landscan-global-"+str(year)+"-assets/landscan-global-"+str(year)+".tif").meta
    meta['nodata'] = 0
    inputfile = "./INPUT/pop_density/2000-2022/landscan-global-"+str(year)+"-assets/landscan-global-"+str(year)+".tif"
    outputfile = "./INPUT/pop_density/2000-2022/landscan-global-"+str(year)+"-assets/landscan-global-"+str(year)+"_ok.tif"
    m = rasterio.open(inputfile).read(1)
    m[where(m<=0)] = 0
    out = rasterio.open(outputfile, 'w', **meta)
    out.write(m, 1)
    out.close()
    print (year)




referencefile = './INPUT/pop_density/2000-2022/ref_pop_1x1.tif'#Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize


for year in range(2000,2023):
    inputfile = "./INPUT/pop_density/2000-2022/landscan-global-"+str(year)+"-assets/landscan-global-"+str(year)+"_ok.tif"
    input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
    inputProj = input.GetProjection()
    inputTrans = input.GetGeoTransform()
    outputfile = "./INPUT/pop_density/2000-2022/"+str(year)+".tif"
    driver= gdal.GetDriverByName('GTiff')
    bandreference = input.GetRasterBand(1)
    output = driver.Create(outputfile,x,y,1,bandreference.DataType)#,options=['COMPRESS=DEFLATE'])
    output.SetGeoTransform(referenceTrans)
    output.SetProjection(referenceProj)
    gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Sum)
    del input
    del output
    print (year)
    
    
        

os.mkdir("./OUTPUT/pop_density_2000_2022")
meta = rasterio.open('./INPUT/pop_density/2000-2022/ref_pop_1x1.tif').meta
meta['nodata'] = 0
meta['dtype'] = 'float64'
for year in range(2000,2023):
    inputfile = "./INPUT/pop_density/2000-2022/"+str(year)+".tif"
    outputfile = "./OUTPUT/pop_density_2000_2022/"+str(year)+"_float.tif"
    m = rasterio.open(inputfile).read(1)
    m[where(m<=0)] = 0
    m[:, 359] = 0
    m_float = m.astype('float64')
    out = rasterio.open(outputfile, 'w', **meta)
    out.write(m_float, 1)
    out.close()
    print (year,m_float.sum()/10**9)





######Future ensemble


meta = rasterio.open("./INPUT/land_1_wgs84.tif").meta
meta['nodata'] = -99999
meta['dtype'] = 'float32'


rcps = ['historical','ssp126','ssp245','ssp370','ssp460','ssp585']
vars = ['pr','tas']
fold_names = ['historical','ssp1_rcp26','ssp2_rcp45','ssp3_rcp7','ssp4_rcp6','ssp5_rcp85']

os.mkdir("./OUTPUT/future_clim_layers")
for rcp in range(len(rcps)):
    try:
        os.mkdir('./OUTPUT/future_clim_layers/'+fold_names[rcp])
    except:
        pass
    for var in vars:
        dataset = Dataset('./INPUT/future_clim_layers/ensemble_nc_files/'+var+'_Global_ensemble_'+rcps[rcp]+'_r1i1p1f1_mean.nc')
        lat = dataset.variables['lat'][:]
        lon = dataset.variables['lon'][:]
        if rcps[rcp]=='historical':
            var_layer = dataset.variables[var][:] #use [-180:] for only from 2000 to 2015 (i.e. 01/2000 to 12/2014)
            y = 1850
        else:
            var_layer = dataset.variables[var][:]
            y = 2015
        n_mon = len(var_layer)
        for mon in range(n_mon):
            m = var_layer[mon]
            m[:,range(180,360)],m[:,range(180)] = m[:,range(180)],m[:,range(180,360)] 
            m = m[::-1,:]
            #convert from kg/m2/sec to monthly total
            if var=='pr':
                n_days = calendar.monthrange(y,mon%12+1)[1]
                m*=24*60*60*n_days
            if mon>0 and mon%12 == 0:
                y+=1
            if rcps[rcp]=='historical':
                out = rasterio.open('./OUTPUT/future_clim_layers/'+fold_names[rcp]+'/'+var+'_'+str(y)+'_'+str(1+mon%12)+'.tif', 'w', **meta)
            else:
                out = rasterio.open('./OUTPUT/future_clim_layers/'+fold_names[rcp]+'/'+var+str(mon)+'.tif', 'w', **meta)
            out.write(m.astype(rasterio.float64), 1)
            out.close()
            print (var+str(y)+' '+str(mon))                
            

 
####################
meta = rasterio.open("./INPUT/land_1_wgs84.tif").meta
meta['nodata'] = -99999
meta['dtype'] = 'float32'

###
referencefile = './OUTPUT/NOOA_climate_prec_rasters/1980_1.tif' #Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize

#prec
fff = [i for i in os.listdir('./OUTPUT/NOOA_climate_prec_2020_2023_rasters') if i[-4:]=='.tif']
for ff in fff:
    inputfile = './OUTPUT/NOOA_climate_prec_2020_2023_rasters/'+ff
    input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
    inputProj = input.GetProjection()
    inputTrans = input.GetGeoTransform()
    outputfile = './OUTPUT/NOOA_climate_prec_rasters/'+ff #Path to output file
    driver= gdal.GetDriverByName('GTiff')
    bandreference = input.GetRasterBand(1)
    output = driver.Create(outputfile,x,y,1,bandreference.DataType,options=['COMPRESS=DEFLATE'])
    output.SetGeoTransform(referenceTrans)
    output.SetProjection(referenceProj)
    gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Bilinear)
    del output
    print (ff)
    

###temp
fff = [i for i in os.listdir('./OUTPUT/NOOA_climate_temp_rasters') if i[-4:]=='.tif']
try:
    os.mkdir('./OUTPUT/NOOA_climate_temp_rasters_1/')
except:
    pass

for ff in fff:
    inputfile = './OUTPUT/NOOA_climate_temp_rasters/'+ff
    input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
    inputProj = input.GetProjection()
    inputTrans = input.GetGeoTransform()
    outputfile = './OUTPUT/NOOA_climate_temp_rasters_1/'+ff #Path to output file
    driver= gdal.GetDriverByName('GTiff')
    bandreference = input.GetRasterBand(1)
    output = driver.Create(outputfile,x,y,1,bandreference.DataType,options=['COMPRESS=DEFLATE'])
    output.SetGeoTransform(referenceTrans)
    output.SetProjection(referenceProj)
    gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Bilinear)
    del output
    print (ff)
    
    
###population density
referencefile = './INPUT/ref_raster_1deg.tif'#Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize
import shutil
for rcp in ['ssp3_rcp7','ssp1_rcp26','ssp2_rcp45','ssp4_rcp6']:
    try:
        os.mkdir('./OUTPUT/pop_density_world_pop_'+rcp+'_1deg/')
    except:
        pass
    try:
        os.mkdir('./OUTPUT/pop_density_world_pop_' + rcp + '_05deg/')
    except:
        pass
    fff = [i for i in os.listdir('./INPUT/pop_density/world_pop_'+rcp+'_05deg/') if i[-4:]=='.tif']
    for ff in fff:
        inputfile = './INPUT/pop_density/world_pop_'+rcp+'_05deg/'+ff
        input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
        inputProj = input.GetProjection()
        inputTrans = input.GetGeoTransform()
        outputfile = './OUTPUT/pop_density_world_pop_'+rcp+'_1deg/'+ff #Path to output file
        driver= gdal.GetDriverByName('GTiff')
        bandreference = input.GetRasterBand(1)
        output = driver.Create(outputfile,x,y,1,bandreference.DataType,options=['COMPRESS=DEFLATE'])
        output.SetGeoTransform(referenceTrans)
        output.SetProjection(referenceProj)
        gdal.ReprojectImage(input,output,inputProj,referenceProj,gdalconst.GRA_Sum)
        del output
        shutil.copyfile(inputfile,'./OUTPUT/pop_density_world_pop_'+rcp+'_05deg/'+ff)
        print (ff)
 



###rasters
def cell_to_coord(col, row, cellx, celly, xmin, ymax):
	lon = cellx*col + xmin + cellx/2.0
	lat = ymax-celly*row - celly/2.0
	return lat,lon


#########anomaly and year of event instead of month
#ref climate
temp_ref,prec_ref = [],[]
for month in range(1,13):
    mon_t = zeros([360, 720])
    mon_p = zeros([360, 720])
    for year in range(1981,2011): #mean of
        mon_t+=rasterio.open('./OUTPUT/NOOA_climate_temp_rasters/' + str(year) + '_' + str(month) + '.tif').read(1)
        mon_p += rasterio.open('./OUTPUT/NOOA_climate_prec_rasters/' + str(year) + '_' + str(month) + '.tif').read(1)
    temp_ref.append(mon_t/30)
    prec_ref.append(mon_p/30)


famine_rasters = os.listdir('./OUTPUT/rasters')
fam_locs = sorted(list(set([i.split('_')[0] for i in famine_rasters])))
fam_ym = sorted(list(set([(int(i.split('_')[1][:4]),int(i.split('_')[1][4:])) for i in famine_rasters])))

all_cells = zeros([360, 720])
famine_rasters_merged_and_sorted = []
for year,month in fam_ym:
    rast_complete = zeros([360, 720])
    rast_complete_ha = zeros([360, 720])
    for fam_loc in fam_locs:
        try:
            famine_raster = rasterio.open('./OUTPUT/rasters/' + fam_loc +'_'+str(year)+str(f"{month:02}")+'_CS.tif').read(1)
            famine_raster[where(famine_raster > 6)] = 0
            famine_raster_ha = rasterio.open('./OUTPUT/rasters/' + fam_loc +'_'+str(year)+str(f"{month:02}")+'_CS.tif').read(2)
            famine_raster_ha[where(famine_raster_ha > 6)] = 0
            rast_complete+=famine_raster + famine_raster_ha
            all_cells+=rast_complete
        except:
            pass
            #print (year,month,fam_loc)
    rast_complete[where(rast_complete > 5)] = 5
    famine_rasters_merged_and_sorted.append(rast_complete)
    print (year,month,rast_complete.sum())


to_remove = where(all_cells==0)
acc = rasterio.open('./INPUT/accessibility_05.tif').read(1)
###note that accessibility was not used in the analyses. I am keeping it here anyway
cris_freq = zeros([360,720])
matched_climate_crisis = []
for ym_n in (range(3,len(fam_ym)-1)):###from 2010 to 2020
    year,month = fam_ym[ym_n]
    pop_dens = rasterio.open('./OUTPUT/pop_density_world_pop_ssp1_rcp26_05deg/CMPI6wp_grid_pop_count' + str(year) + '_SSP1_RCP2_6.tif').read(1)
    if year<2021:
        gdp = rasterio.open('./OUTPUT/GDP_05_deg/GDP' + str(year) + '.tif').read(1) #https://zenodo.org/records/7898409
    else:
        gdp = rasterio.open('./OUTPUT/GDP_05_deg/GDP2020.tif').read(1)
    famine_raster = famine_rasters_merged_and_sorted[ym_n]
    pre = famine_rasters_merged_and_sorted[ym_n-1]
    post = famine_rasters_merged_and_sorted[ym_n+1]
    crisis = zeros(famine_raster.shape)
    crisis[where(famine_raster>2)] = 1
    crisis[where(post<3)] = 0 #if the alert is only for one period, no crisis
    crisis[where(pre>2)] = 0 #if previous period was already a crisis, then need to be removed.
    crisis[where((1*(pre>2)*(famine_raster>2))==1)] = -99999 ##but the cells where the crisis continued are removed from the model as not informative!
    crisis[to_remove] = -99999
    crisis_continued = zeros(famine_raster.shape)
    crisis_continued[where(famine_raster>2)] = 1
    crisis_continued[where(pre<3)] = 0
    crisis_continued[to_remove] = -99999
    crisis_famine = zeros(famine_raster.shape)
    crisis_famine[where(famine_raster>2)] = 1
    crisis_famine[to_remove] = -99999
    cris_freq[where(crisis==1)]+=1
    obs_wind = 3
    y_start = int(year - obs_wind)
    m_start = 1#month
    t_range = []
    for i in range(int(obs_wind*12)):
        t_range.append([y_start, m_start])
        m_start += 1
        if m_start == 13:
            m_start = 1
            y_start += 1
    temp_lays = [rasterio.open('./OUTPUT/NOOA_climate_temp_rasters/' + str(i[0]) + '_' + str(i[1]) + '.tif').read(1)-temp_ref[i[1]-1] for i in t_range]
    prec_lays = [rasterio.open('./OUTPUT/NOOA_climate_prec_rasters/' + str(i[0]) + '_' + str(i[1]) + '.tif').read(1)-prec_ref[i[1]-1] for i in t_range]
    cris_n = crisis.sum()
    R,C = famine_raster.shape
    for row in range(R):
        for col in range(C):
            crisis_val,crisis_continued_val,crisis_famine_val = crisis[row][col],crisis_continued[row][col],crisis_famine[row][col]
            if any([val != -99999 for val in [crisis_val,crisis_continued_val,crisis_famine_val]]):
                lat, lon = cell_to_coord(col, row, 0.5, 0.5, -180, 90)
                loc_temp = [temp_lay[row][col] for temp_lay in temp_lays]
                loc_prec = [prec_lay[row][col] for prec_lay in prec_lays]
                matched_climate_crisis.append([lat, lon, crisis_val, year, month] + loc_temp + loc_prec + [acc[row][col], pop_dens[row][col],gdp[row][col],crisis_continued_val,crisis_famine_val])
                print(matched_climate_crisis[-1])




meta = rasterio.open("./INPUT/land_05_wgs84.tif").meta
meta['nodata'] = -99999
meta['dtype'] = 'float32'
out = rasterio.open('./OUTPUT/cris_freq_prol.tif', 'w', **meta)
out.write(cris_freq.astype(rasterio.float64), 1)
out.close()


out = open('./OUTPUT/matched_climate_crisis.csv','w')
out.write('row,col,fam_lev,year,month,'+','.join(['t'+str(i) for i in range(int(obs_wind*12))]+['p'+str(i) for i in range(int(obs_wind*12))]+['acc','pop_density','gdp','crisis_continued_val','crisis_famine_val'])+'\n')
for i in matched_climate_crisis:
    out.write(','.join(map(str,i))+'\n')


out.close()


##########rescale the cris frequency to 1x1
referencefile = './INPUT/ref_raster_1deg.tif'  # Path to reference file
reference = gdal.Open(referencefile, gdalconst.GA_ReadOnly)
referenceProj = reference.GetProjection()
referenceTrans = reference.GetGeoTransform()
x = reference.RasterXSize
y = reference.RasterYSize
inputfile = './OUTPUT/cris_freq_prol.tif'
input = gdal.Open(inputfile, gdalconst.GA_ReadOnly)
inputProj = input.GetProjection()
inputTrans = input.GetGeoTransform()
outputfile = './OUTPUT/cris_freq_prol_1x1.tif' # Path to output file
driver = gdal.GetDriverByName('GTiff')
bandreference = input.GetRasterBand(1)
output = driver.Create(outputfile, x, y, 1, bandreference.DataType, options=['COMPRESS=DEFLATE'])
output.SetGeoTransform(referenceTrans)
output.SetProjection(referenceProj)
gdal.ReprojectImage(input, output, inputProj, referenceProj, gdalconst.GRA_Sum)
del output

###############
###make future projections
###years from 2015 to 31/12/2100
##projections start from 2020; ref is target year - 10year
#target month
#8/2023


##used approach: compute anomalies from historical ssp model
temp_ref_,prec_ref_ = [],[]
for month in range(1,13):
    mon_t = zeros([180, 360])
    mon_p = zeros([180, 360])
    for year in range(1981,2011): #mean of
        mon_t+=rasterio.open('./OUTPUT/future_clim_layers/historical/tas_' + str(year) + '_' + str(month) + '.tif').read(1)
        mon_p += rasterio.open('./OUTPUT/future_clim_layers/historical/pr_' + str(year) + '_' + str(month) + '.tif').read(1)
    temp_ref_.append(mon_t/30.0)
    prec_ref_.append(mon_p/30.0)



land = rasterio.open("./INPUT/land_1_wgs84.tif").read(1)
ok_cells = where(land==1)
try:
    os.mkdir('./OUTPUT/future_data_for_projections')
except:
    pass


for rcp in ['ssp1_rcp26','ssp2_rcp45','ssp3_rcp7','ssp4_rcp6']:
    out_temp = open('./OUTPUT/future_data_for_projections/future_temp_' + rcp + '.csv', 'w')
    out_prec = open('./OUTPUT/future_data_for_projections/future_prec_' + rcp + '.csv', 'w')
    out_pop = open('./OUTPUT/future_data_for_projections/future_pop_' + rcp + '.csv', 'w')
    out_gdp = open('./OUTPUT/future_data_for_projections/future_gdp_' + rcp + '.csv', 'w')
    for target_year in range(2015,2100):
        pop_dens = rasterio.open('./OUTPUT/pop_density_world_pop_'+rcp+'_1deg/CMPI6wp_grid_pop_count' + str(target_year) + '_'+rcp.upper().split('RCP')[0]+'RCP'+'_'.join(list(map(str,rcp.split('rcp')[1])))+'.tif').read(1)
        if target_year<2021:
            gdp = rasterio.open('./OUTPUT/GDP_1_deg/GDP'+str(target_year)+'.tif').read(1)
        else:
            gdp = rasterio.open('./OUTPUT/GDP_1_deg/GDP'+str(target_year-target_year%5+5)+'_'+rcp.split('_')[0]+'.tif').read(1)
        for month in range(12):
            target_month = (target_year-2015)*12+month
            temp = rasterio.open('./OUTPUT/future_clim_layers/'+rcp+'/tas' + str(target_month) + '.tif').read(1)-temp_ref[month]
            prec = rasterio.open('./OUTPUT/future_clim_layers/'+rcp+'/pr' + str(target_month) + '.tif').read(1)-prec_ref[month] ##doublecheck the x4
            out_temp.write(','.join(map(str, temp[ok_cells])) + '\n')
            out_prec.write(','.join(map(str, prec[ok_cells])) + '\n')
            out_pop.write(','.join(map(str, pop_dens[ok_cells])) + '\n')
            out_gdp.write(','.join(map(str, gdp[ok_cells])) + '\n')
        print(target_year,rcp)
    out_temp.close()
    out_prec.close()
    out_pop.close()
    out_gdp.close()



####make projections from 2010 to 2023 based on current climate
###anomalies are computed using NOOA climate

temp_ref,prec_ref = [],[]
for month in range(1,13):
    mon_t = zeros([180, 360])
    mon_p = zeros([180, 360])
    for year in range(1981,2011): #mean of
        mon_t+=rasterio.open('./OUTPUT/NOOA_climate_temp_rasters_1/' + str(year) + '_' + str(month) + '.tif').read(1)
        mon_p += rasterio.open('./OUTPUT/NOOA_climate_prec_rasters_1/' + str(year) + '_' + str(month) + '.tif').read(1)
    temp_ref.append(mon_t/30.0)
    prec_ref.append(mon_p/30.0)


land = rasterio.open("./INPUT/land_1_wgs84.tif").read(1)
ok_cells = where(land==1)

out_temp = open('./OUTPUT/future_data_for_projections/future_temp_cur.csv', 'w')
out_prec = open('./OUTPUT/future_data_for_projections/future_prec_cur.csv', 'w')
out_pop = open('./OUTPUT/future_data_for_projections/future_pop_cur.csv', 'w')
out_gdp = open('./OUTPUT/future_data_for_projections/future_gdp_cur.csv', 'w')

for target_year in range(2000,2023):
    pop_dens = rasterio.open('./OUTPUT/pop_density_2000_2022/'+str(target_year) + '_float.tif').read(1)
    if target_year < 2021:
        gdp = rasterio.open('./OUTPUT/GDP_1_deg/GDP' + str(target_year) + '.tif').read(1)
    else:
        gdp = zeros([180, 360])
        for rcp in ['ssp1_rcp26', 'ssp2_rcp45', 'ssp3_rcp7','ssp4_rcp6']:
            gdp += rasterio.open('./OUTPUT_GDP_1_deg/GDP' + str(target_year - target_year % 5 + 5) + '_' + rcp.split('_')[0] + '.tif').read(1)
        gdp/=4.0
    for month in range(12):
        temp = rasterio.open('./OUTPUT/NOOA_climate_temp_rasters_1/' + str(target_year) + '_' + str(month+1) + '.tif').read(1)-temp_ref[month]
        prec = rasterio.open('./OUTPUT/NOOA_climate_prec_rasters_1/' + str(target_year) + '_' + str(month+1) + '.tif').read(1)-prec_ref[month]
        out_temp.write(','.join(map(str, temp[ok_cells])) + '\n')
        out_prec.write(','.join(map(str, prec[ok_cells])) + '\n')
        out_pop.write(','.join(map(str, pop_dens[ok_cells])) + '\n')
        out_gdp.write(','.join(map(str, gdp[ok_cells])) + '\n')
    print(target_year)


out_temp.close()
out_prec.close()
out_pop.close()
out_gdp.close()


land = rasterio.open("./INPUT/land_1_wgs84.tif").read(1)
ok_cells = where(land==1)

def cell_to_coord(col, row, cellx, celly, xmin, ymax):
	lon = cellx*col + xmin + cellx/2.0
	lat = ymax-celly*row - celly/2.0
	return lat,lon


def deg_to_km(lat,res,R=6371.0072):
	return (math.sin(math.radians(lat+(res/2.0)))-math.sin(math.radians(lat-(res/2.0))))*math.radians(res)*(R**2)


out = open('./OUTPUT/lat_lon.csv','w')
out.write('latitude,longitude,area\n')
for i in range(len(ok_cells[0])):
    row,col = ok_cells[0][i],ok_cells[1][i]
    lat,lon = cell_to_coord(col,row,1,1,-180,90)
    area = deg_to_km(lat-1, 1)
    out.write(','.join(map(str,[lat,lon,area]))+'\n')


out.close()