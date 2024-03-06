import os
import re
import copy
import glob
import cdsapi
import cfgrib
import datetime
import subprocess
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from tqdm import tqdm
from .fstrings import *
from .regex_strings import *
from multiprocessing.pool import ThreadPool
from cartopy.feature import NaturalEarthFeature






def is_valid_date(year:int, month:int, day:int)->bool:
    """Function by 'Anon' that checks a date to see if it exists
    https://stackoverflow.com/a/51981596 
    """    
    day_count_for_month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    if year%4==0 and (year%100 != 0 or year%400==0):
        day_count_for_month[2] = 29
    
    return (1 <= month <= 12 and 1 <= day <= day_count_for_month[month])


class YesNo:
    """Class to associate various spellings of "yes" and "no" with True and False
    """
    
    def __init__(self,val:str|float):
        
        if val in ['YES','Yes','yes']:
            self.bool = True
            self.val = val
            
        elif val in ['NO','No','no']:
            self.bool = False
            self.val = val
            
        elif val == True:
            self.bool = True
            self.val = 'Yes'
            
        elif val == False:
            self.bool = False
            self.val = 'No'
            
        # the class needs to be able to be initialised from another
        # instance for easy program flow
        elif type(val) == type(self):
            self.bool = val.bool
            self.val = val.val
            
        else:
            raise TypeError(val," not one of Yes, No, True, False")
    
    def __repr__(self):
        # the value of the YesNo object that is inserted into the fstring
        return(self.val)
    
    def __str__(self):
        # the value of the YesNo object that is inserted into the fstring
        return(self.val)
    
    def __bool__(self):
        # the value of the YesNo object that is used for logical operations
        return(self.bool)
    


class OnOff:
    """Class to associate various spellings of "On" and "Off" with True and False
    """
    
    def __init__(self,val:str|float):
        """
        """
        
        if val in ['ON','On','on']:
            self.bool = True
            self.val = val
            
        elif val in ['OFF','Off','off']:
            self.bool = False
            self.val = val
            
        elif val == True:
            self.bool = True
            self.val = 'On'
            
        elif val == False:
            self.bool = False
            self.val = 'Off'
            
        # the class needs to be able to be initialised from another
        # instance for easy program flow
        elif type(val) == type(self):
            self.bool = val.bool
            self.val = val.val
            
        else:
            raise TypeError(val,"not one of On, Off, True, False")
    
    def __repr__(self):
        """
        """
        # the value of the OnOff object that is inserted into the fstring
        return(self.val)
    
    def __str__(self):
        """
        """
        # the value of the YesNo object that is inserted into the fstring
        return(self.val)
    
    def __bool__(self):
        """
        """
        # the value of the OnOff object that is used for logical operations
        return(self.bool)

class Longitude:
    """Class to remove annoying, hard to debug longitude definition errors
    """

    def __init__(self, value:float, type:int):

        
        if type == 360:
            if value<0:
                raise TypeError("Longitude cannot be negative when defined 0-360")
            else:
                self._360 = value
                self._180 = (value + 180) % 360 - 180

        if type == 180:
            if value>180:
                raise TypeError("Longitude cannot be > 180 when defined -180 to 180")
            else:
                self._360 = value %360
                self._180 = value


def Lon360(value:float)->Longitude:
    """Quicker to type than Longitude(value,360)
    """
    return(Longitude(value, 360))
    

def Lon180(value:float)->Longitude:
    """Quicker to type than Longitude(value, 180)
    """
    return(Longitude(value, 180))





class Section:
    """Generic class for a Fall3D input file section. 
    All sections need to be able to do 1 of 3 things:
    
    (1) create new from string from valid Fall3D input file
    (2) write to a string that is a valid part of a Fall3D input file
    (3) update values while performing all initialization type and value checks all over again
    
    plus domain specific visualisations, etc.., which are implemented in the child classes
    
    """
    @classmethod    
    def from_string(cls,string:str):
        """
        """
        
        # get variables from string
        variables = re.search(cls.regex, string,re.DOTALL).groupdict()
        
        # initialise a dict to hold the variables after we turn them into
        # the appropriate types
        variables_typed = {}
        
        # 
        for key, item in variables.items():
        
            variables_typed[key] = cls.types[key](item)
        
        cls = cls(**variables_typed)
        
        return(cls)
    
    def to_string(self):
        """
        """
        
        variables = {}
        
        for key in self.types.keys():
            
            variables[key] = getattr(self,key)
        
        string = self.fstring.format(**variables)
        
        return(string)


    def __repr__(self):
        """
        """
        return(self.to_string())

    def __str__(self):
        """
        """
        return(self.to_string())


    def update(self,new_values:dict):
        """Update witgh values in dict
        """
        # createb a dict that wil be used to reinitialise theb object
        params = {}

        for key, item in self.types.items():
            
            params[key] = getattr(self,key)

        for key, item in new_values.items():

            params[key] = item

        self.__init__(**params)


        




class TimeUTC(Section):
    
    fstring = fstring_time_utc
    
    regex = regex_time_utc
    
    types = {
        'year':int,
        'month':int,
        'day':int,
        'run_start':int,
        'run_end':int,
        'initial_condition':str,
        'restart_file':str,
        'restart_ensemble_basepath':str
    }
    
    
    
    def __init__(self, 
                year:int, 
                month:int, 
                day:int, 
                run_start:int, 
                run_end:int, 
                initial_condition:str, 
                restart_file:str, 
                restart_ensemble_basepath:str
                ):
        
        #---------------------------------------------------------------------------------
        # TEST INPUTS
        #---------------------------------------------------------------------------------
        

        # tests for year
        if (year<1970) or (year>2024):
            raise ValueError('Year must be >= 1970 and <= 2024')

        if type(year) != int:
            raise TypeError("Year must be type int")


        # tests for month
        if (type(month) !=int):
            raise TypeError("Month must be type int")

        if (month<1) or (month>12):
            raise TypeError("Month must be >= 1 and <= 12")
        
        # tests for day
        if (type(day) !=int):
            raise TypeError("Day must be type int")

        if (day<1) or (day>31):
            raise TypeError('Day must be between 1 and 31')

        # test for whole date
        assert (is_valid_date(year,month,day)), ["/"].join(str(year),str(month),str(day))+" is not a valid date."
        
        # tests for initial_condition
        assert (type(initial_condition)==str), "initial_conditions must be type str"
        assert (initial_condition in ['NONE', 'RESTART', 'INSERTION']), "initial_condition must be one of NONE, RESTART, INSERTION"
        
        
        # tests for restart_file
        assert (type(restart_file)==str), "restart_file must be type string"
        
        # tests for restart_ensemble_basepath
        assert (type(restart_ensemble_basepath)==str), "restart_ensemble_basepath must be type string"
        
        #---------------------------------------------------------------------------------       
        # OPTIONS THAT HAVE BECOME ACTIVATED
        #---------------------------------------------------------------------------------

        if initial_condition in ['RESTART']:
            print("RESTART_FILE in use as INITIAL_CONDITION = RESTART  )")
            
        self.year = year
        self.month = month
        self.day = day
        self.run_start = run_start
        self.run_end = run_end
        self.initial_condition = initial_condition
        self.restart_file = restart_file
        self.restart_ensemble_basepath = restart_ensemble_basepath  
        




class InsertionData(Section):
    
    fstring = fstring_insertion_data
    
    regex = regex_insertion_data
    
    types = {
        'insertion_file':str,
        'insertion_dictionary_file':str,
        'insertion_time_slab':int,
        'diameter_cut_off':int
    }
    
    def __init__(self,
                insertion_file:str, 
                insertion_dictionary_file:str, 
                insertion_time_slab:int,
                diameter_cut_off:int
                ):
        
        # TEST INPUTS
        
        # tests for insertion_file
        assert (type(insertion_file)==str), "insertion_file must be string"
        
        # tests for insertion_dictionary_file
        assert (type(insertion_dictionary_file)==str), "insertion_dictionary_file must be string"
        
        
        # tests for insertion_time_slab
        assert (type(insertion_time_slab)==int), "insertion_time_slab must be int"
        
        # tests for diameter_cut_off_mic
        assert (type(diameter_cut_off)==int)
        
        self.insertion_file = insertion_file
        self.insertion_dictionary_file = insertion_dictionary_file
        self.insertion_time_slab = insertion_time_slab
        self.diameter_cut_off = diameter_cut_off
        


class MeteoData(Section):
    
    
    fstring = fstring_meteo
    
    regex = regex_meteo
    
    types={
        
        'meteo_data_format':str,
        'meteo_data_dictionary_file':str,
        'meteo_data_file':str,
        'meteo_ensemble_basepath':str,
        'meteo_levels_file':str,
        'dbs_begin_meteo_data':int,
        'dbs_end_meteo_data':int,
        'meteo_coupling_interval':int,
        'memory_chunk_size':int
    }
    
    def __init__(self, 
                meteo_data_format:str, 
                meteo_data_dictionary_file:str,
                meteo_data_file:str,
                meteo_ensemble_basepath:str,
                meteo_levels_file:str,
                dbs_begin_meteo_data:int,
                dbs_end_meteo_data:int,
                meteo_coupling_interval:int,
                memory_chunk_size:int
                ):

        # TEST INPUTS
        
        # test meteo_data_format
        assert(type(meteo_data_format)==str)
        assert( meteo_data_format in ['WRF', 'GFS' , 'ERA5' , 'ERA5ML' , 'IFS' , 'CARRA'])
        
        # test meteo_data_dictionary
        assert(type(meteo_data_dictionary_file)==str)
        
        # test meteo_data_file
        assert(type(meteo_data_file)==str)
        
        # test meteo_ensemble_basepath
        assert(type(meteo_ensemble_basepath)==str)
        
        # test meteo_levels_file
        assert(type(meteo_levels_file)==str)
        
        # test dbs_begin_meteo_data
        assert(type(dbs_begin_meteo_data)==int)
        assert(dbs_begin_meteo_data>=0)
        
        # test dbs_end_meteo_data
        assert(type(dbs_end_meteo_data)==int)
        assert(dbs_end_meteo_data>=0)
        assert(dbs_end_meteo_data>dbs_begin_meteo_data)
        
        # test meteo_coupling_interval
        assert(type(meteo_coupling_interval)==int)
        assert(meteo_coupling_interval>0)
        
        # test memory_chunk_size
        assert(type(memory_chunk_size)==int)
        assert(memory_chunk_size>0)
        
        # print options that have become activated
        if meteo_data_format in ['ERA5ML', 'IFS']:
            print("METEO_DATA_FORMAT has been set to ERA5ML or IFS so METEO_LEVELS_FILE is in use.")
        
        self.meteo_data_format = meteo_data_format 
        self.meteo_data_dictionary_file = meteo_data_dictionary_file  
        self.meteo_data_file = meteo_data_file 
        self.meteo_ensemble_basepath = meteo_ensemble_basepath 
        self.meteo_levels_file = meteo_levels_file
        self.dbs_begin_meteo_data = dbs_begin_meteo_data
        self.dbs_end_meteo_data = dbs_end_meteo_data
        self.meteo_coupling_interval = meteo_coupling_interval
        self.memory_chunk_size = memory_chunk_size

    def plot_on_map(self,ax=None):

        if not os.path.exists(self.meteo_data_file):
            raise OSError("Meteo data file does not exist.")

        ds = xr.open_dataset(self.meteo_data_file)

        # get outline of meteo data

        west_lat = ds.latitude.sel(x=0).values
        west_lon = ds.longitude.sel(x=0).values

        east_lat = ds.latitude.sel(x=-1).values
        east_lon = ds.longitude.sel(x=-1).values

        south_lat = ds.latitude.sel(y=0).values
        south_lon = ds.longitude.sel(y=0).values

        north_lat = ds.latitude.sel(y=-1).values
        north_lon = ds.longitude.sel(y=-1).values
        
        extent = [
                    ds.longitude.min(), 
                    ds.longitude.max(),
                    ds.latitude.min(), 
                    ds.latitude.max()
                ]

        

        if ax is None:
            plt.figure("Test Map")
            crs = ccrs.PlateCarree()
            ax = plt.subplot(111, projection=crs)
            
            ax.set_extent(extent, crs=crs)
            
            #ax.add_feature(NaturalEarthFeature('physical', 'ocean', '50m'))
            ax.coastlines(resolution='10m',color='blue')
    
            ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)

        ax.plot(west_lon-360, west_lat, color='blue')
        ax.plot(east_lon-360, east_lat, color='blue')
        ax.plot(north_lon-360, north_lat, color='blue')
        ax.plot(south_lon-360, south_lat, color='blue')

        return(extent)
        
        #plt.show()
        





class Grid(Section):
    
    
    fstring = fstring_grid

    
    regex = regex_grid
    
    types={
        'horizontal_mapping':str,
        'vertical_mapping':str,
        'lonmin':float,
        'lonmax':float,
        'latmin':float,
        'latmax':float,
        'nx':int,
        'ny':int,
        'nz':int,
        'zmax':int,
        #'sigma_values':str,
    }
    
    def __init__(self,
                horizontal_mapping:str,
                vertical_mapping:str, 
                lonmin:float, 
                lonmax:float, 
                latmin:float,
                latmax:float, 
                nx:int, 
                ny:int, 
                nz:int, 
                zmax:int#, 
                #sigma_values:str
                ):

        print("WARNING sigma_values_currently_disabled")
        # TEST INPUTS
        
        # test horizontal_mapping
        assert(type(horizontal_mapping)==str)
        assert( horizontal_mapping in ['CARTESIAN', 'SPHERICAL' ])
        
        # test vertical_mapping
        assert(type(vertical_mapping)==str)
        assert( vertical_mapping in ['SIGMA_NO_DECAY', 'SIGMA_LINEAR_DECAY', 'SIGMA_EXPONENTIAL_DECAY'])
        
        # test lonmin
        assert(type(lonmin)==float)
        assert(lonmin>-180)
        assert(lonmin<180)
        
        # test lonmax
        assert(type(lonmax)==float)
        assert(lonmin>-180)
        assert(lonmin<180)       
        
        # check lonmax > lonmin
        assert(lonmax>lonmin)
        
        # test latmin
        assert(type(latmin)==float)
        assert(latmin>-90)
        assert(latmin<90)
                
        
        # test latmax
        assert(type(latmax)==float)
        assert(latmax>-90)
        assert(latmax<90)
        
        # check latmin < latmax
        assert(latmin < latmax)
        
        # test nx
        assert(type(nx)==int)
        assert(nx>1)
        
        #test ny
        assert(type(ny)==int)
        assert(ny>1)
        
        #test nz
        assert(type(nz)==int)
        assert(nz>1)
        
        #test zmax
        assert(type(zmax)==int)
        assert(zmax>0)
        
        #test sigma_values
        #assert(type(sigma_values)==str)
        
        
        
        
        
        
        self.horizontal_mapping = horizontal_mapping
        self.vertical_mapping = vertical_mapping 
        self.lonmin = lonmin 
        self.lonmax = lonmax 
        self.latmin = latmin 
        self.latmax = latmax 
        self.nx = nx 
        self.ny = ny 
        self.nz = nz 
        self.zmax = zmax
        #self.sigma_values =sigma_values

    def plot_on_map(self, ax=None):

        extent = [self.lonmin-1, self.lonmax+1,self.latmin-1, self.latmax+1]

        bbox = np.array([
                    [self.lonmin, self.latmin], # lower left
                    [self.lonmin, self.latmax], # upper left
                    [self.lonmax, self.latmax], # upper right
                    [self.lonmax, self.latmin], # lower right
                    [self.lonmin, self.latmin] # back to lower left again
        ])


        if ax == None:
            plt.figure("Test Map")
            ax = plt.subplot(111, projection=ccrs.PlateCarree())
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            
            #ax.add_feature(NaturalEarthFeature('physical', 'ocean', '50m'))
            ax.coastlines(resolution='10m',color='blue')
    
            ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)

        ax.plot(*bbox.T, color='red')
        
        #plt.show()





class Species(Section):
    
    fstring = fstring_species
    
    regex = regex_species
    
    types={
        
        'tephra':OnOff,
        'dust':OnOff,
        'h2o':OnOff,
        'so2':OnOff,
        'cs134':OnOff,
        'cs137':OnOff,
        'i131':OnOff,
        'sr90':OnOff,
        'y90':OnOff,
        'h2o_mass_fraction':float,
        'so2_mass_fraction':float,
        'cs134_mass_fraction':float,
        'cs137_mass_fraction':float,
        'i131_mass_fraction':float,
        'sr90_mass_fraction':float,
        'y90_mass_fraction':float
    }
    
    def __init__(self, 
                tephra:OnOff, 
                dust:OnOff, 
                h2o:OnOff, 
                so2:OnOff, 
                cs134:OnOff, 
                cs137:OnOff, 
                i131:OnOff, 
                sr90:OnOff, 
                y90:OnOff, 
                h2o_mass_fraction:float,
                so2_mass_fraction:float, 
                cs134_mass_fraction:float, 
                cs137_mass_fraction:float,
                i131_mass_fraction:float, 
                sr90_mass_fraction:float, 
                y90_mass_fraction:float
                ):
        
        

        assert((type(tephra)==bool)or(type(tephra)==OnOff))
        assert((type(dust)==bool)or(type(dust)==OnOff))
        assert((type(h2o)==bool)or(type(h2o)==OnOff))
        assert(type(h2o_mass_fraction)==float)
        assert((type(so2)==bool)or(type(so2)==OnOff))
        assert(type(so2_mass_fraction)==float)
        assert((type(cs137)==bool)or(type(cs137)==OnOff))
        assert(type(cs137_mass_fraction)==float)
        assert((type(sr90)==bool)or(type(sr90)==OnOff))
        assert(type(sr90_mass_fraction)==float)
        assert((type(y90)==bool)or(type(y90)==OnOff))
        assert(type(y90_mass_fraction)==float)
        
        # check constraints
        # (1) only one of tephra or dust is allowed
        assert(not (tephra and dust))
        
        # (2) if tephra = Off, then mass fraction of aeraosols must sum to 1 (100% apparently)
        if not tephra:
            assert(h2o_mass_fraction+so2_mass_fraction==100)
        
        # (3) cannot run PARTICLES (tephra, dust) and RADIONUCLEIDES at the same time
        assert(
            not (
                    any([tephra, dust]) and any([cs134, cs137, i131, sr90,  y90])
            )
        )
        
        # (4) mass fraction of radionucleides must sum to 1 (100%, apparently)
        if (cs134 or cs137 or i131 or sr90 or y90):
            assert(
                (cs134_mass_fraction + cs137_mass_fraction + 
                 i131_mass_fraction + sr90_mass_fraction + 
                 y90_mass_fraction) == 100
            )
        
        self.tephra = OnOff(tephra)
        self.dust = OnOff(dust) 
        self.h2o = OnOff(h2o)
        self.h2o_mass_fraction = h2o_mass_fraction
        self.so2 = OnOff(so2)
        self.so2_mass_fraction = so2_mass_fraction
        self.cs134 = OnOff(cs134)
        self.cs134_mass_fraction = cs134_mass_fraction
        self.cs137 = OnOff(cs137)
        self.cs137_mass_fraction = cs137_mass_fraction
        self.i131 = OnOff(i131)
        self.i131_mass_fraction = i131_mass_fraction
        self.sr90 = OnOff(sr90)
        self.sr90_mass_fraction = sr90_mass_fraction
        self.y90 = OnOff(y90)
        self.y90_mass_fraction = y90_mass_fraction
               


class TephraTgsd(Section):
    
    fstring = fstring_tephra_tgsd
    
    regex = regex_tephra_tgsd
    
    types = {
        'number_of_bins':int,
        'fi_range1':int, 
        'fi_range2':int,
        'density_range1':int, 
        'density_range2':int,
        'sphericity_range1':float, 
        'sphericity_range2':float,
        'distribution':str,
        'gaussian_fi_mean':float, 
        'gaussian_fi_disp':float,
        'bigaussian_fi_mean1':float, 
        'bigaussian_fi_mean2':float,
        'bigaussian_fi_disp1':float, 
        'bigaussian_fi_disp2':float,
        'bigaussian_mixing_factor':float,
        'weibull_fi_scale':float, 
        'weibull_w_shape':float,
        'biweibull_fi_scale1':float, 
        'biweibull_fi_scale2':float,
        'biweibull_w_shape1':float, 
        'biweibull_w_shape2':float,
        'biweibull_mixing_factor':float,
        'custom_file':str,
        'estimate_viscosity':str,
        'estimate_height_above_vent':str
    }
    
    def __init__(self,
                number_of_bins:int,
                fi_range1:int, 
                fi_range2:int,
                density_range1:int, 
                density_range2:int,
                sphericity_range1:float, 
                sphericity_range2:float,
                distribution:str,
                gaussian_fi_mean:float, 
                gaussian_fi_disp:float,
                bigaussian_fi_mean1:float, 
                bigaussian_fi_mean2:float,
                bigaussian_fi_disp1:float, 
                bigaussian_fi_disp2:float,
                bigaussian_mixing_factor:float,
                weibull_fi_scale:float, 
                weibull_w_shape:float,
                biweibull_fi_scale1:float, 
                biweibull_fi_scale2:float,
                biweibull_w_shape1:float, 
                biweibull_w_shape2:float,
                biweibull_mixing_factor:float,
                custom_file:str,
                estimate_viscosity:str,
                estimate_height_above_vent:str
                ):
        
        assert(type(number_of_bins)==int)
        assert(type(fi_range1)==int)
        assert(type(fi_range2)==int)
        assert(type(density_range1)==int) 
        assert(type(density_range2)==int)
        assert(type(sphericity_range1)==float) 
        assert(type(sphericity_range2)==float)
        assert(type(distribution)==str)
        assert(type(gaussian_fi_mean)==float) 
        assert(type(gaussian_fi_disp)==float)
        assert(type(bigaussian_fi_mean1)==float) 
        assert(type(bigaussian_fi_mean2)==float)
        assert(type(bigaussian_fi_disp1)==float) 
        assert(type(bigaussian_fi_disp2)==float)
        assert(type(bigaussian_mixing_factor)==float)
        assert(type(weibull_fi_scale)==float)
        assert(type(weibull_w_shape)==float)
        assert(type(biweibull_fi_scale1)==float) 
        assert(type(biweibull_fi_scale2)==float)
        assert(type(biweibull_w_shape1)==float)
        assert(type(biweibull_w_shape2)==float)
        assert(type(biweibull_mixing_factor)==float)
        assert(type(custom_file)==str)
        assert(type(estimate_viscosity)==str)
        assert(type(estimate_height_above_vent)==str)
        
        self.number_of_bins = number_of_bins
        self.fi_range1 = fi_range1 
        self.fi_range2 = fi_range2
        self.density_range1 = density_range1
        self.density_range2 = density_range2
        self.sphericity_range1 = sphericity_range1
        self.sphericity_range2 = sphericity_range2
        self.distribution = distribution
        self.gaussian_fi_mean = gaussian_fi_mean  
        self.gaussian_fi_disp = gaussian_fi_disp
        self.bigaussian_fi_mean1 = bigaussian_fi_mean1 
        self.bigaussian_fi_mean2 = bigaussian_fi_mean2
        self.bigaussian_fi_disp1 = bigaussian_fi_disp1 
        self.bigaussian_fi_disp2 = bigaussian_fi_disp2
        self.bigaussian_mixing_factor = bigaussian_mixing_factor
        self.weibull_fi_scale = weibull_fi_scale 
        self.weibull_w_shape = weibull_w_shape
        self.biweibull_fi_scale1 = biweibull_fi_scale1 
        self.biweibull_fi_scale2 = biweibull_fi_scale2
        self.biweibull_w_shape1 = biweibull_w_shape1 
        self.biweibull_w_shape2 = biweibull_w_shape2 
        self.biweibull_mixing_factor = biweibull_mixing_factor
        self.custom_file = custom_file
        self.estimate_viscosity = estimate_viscosity 
        self.estimate_height_above_vent = estimate_height_above_vent
        
        
        
        
        
        
class RadionucleidesTgsd(Section):
    
    fstring = fstring_radionucleides_tgsd

    regex = regex_radionucleides_tgsd
        
    types = {'number_of_bins':int,
            'fi_range1':int,
            'fi_range2':int,
            'density_range1':int,
            'density_range2':int,
            'sphericity_range1':float,
            'sphericity_range2':float,
            'distribution':str,
            'gaussian_fi_mean':float,
            'gaussian_fi_disp':float,
            'bigaussian_fi_mean1':float,
            'bigaussian_fi_mean2':float,
            'bigaussian_fi_disp1':float,
            'bigaussian_fi_disp2':float,
            'bigaussian_mixing_factor':float,
            'weibull_fi_scale':float,
            'weibull_w_shape':float,
            'biweibull_fi_scale1':float,
            'biweibull_fi_scale2':float,
            'biweibull_w_shape1':float,
            'biweibull_w_shape2':float,
            'biweibull_mixing_factor':float,
            'custom_file':str,
            'estimate_viscosity':str,
            'estimate_height_above_vent':str}
    
    def __init__(self,
                number_of_bins:int,
                fi_range1:int,
                fi_range2:int,
                density_range1:int,
                density_range2:int,
                sphericity_range1:float,
                sphericity_range2:float,
                distribution:str,
                gaussian_fi_mean:float,
                gaussian_fi_disp:float,
                bigaussian_fi_mean1:float,
                bigaussian_fi_mean2:float,
                bigaussian_fi_disp1:float,
                bigaussian_fi_disp2:float,
                bigaussian_mixing_factor:float,
                weibull_fi_scale:float,
                weibull_w_shape:float,
                biweibull_fi_scale1:float,
                biweibull_fi_scale2:float,
                biweibull_w_shape1:float,
                biweibull_w_shape2:float,
                biweibull_mixing_factor:float,
                custom_file:str,
                estimate_viscosity:str,
                estimate_height_above_vent:str,
                ):
    
        assert(type(number_of_bins)==int)
        assert(type(fi_range1)==int)
        assert(type(fi_range2)==int)
        assert(type(density_range1)==int)
        assert(type(density_range2)==int)
        assert(type(sphericity_range1)==float)
        assert(type(sphericity_range2)==float)
        assert(type(distribution)==str)
        assert(type(gaussian_fi_mean)==float)
        assert(type(gaussian_fi_disp)==float)
        assert(type(bigaussian_fi_mean1)==float)
        assert(type(bigaussian_fi_mean2)==float)
        assert(type(bigaussian_fi_disp1)==float)
        assert(type(bigaussian_fi_disp2)==float)
        assert(type(bigaussian_mixing_factor)==float)
        assert(type(weibull_fi_scale)==float)
        assert(type(weibull_w_shape)==float)
        assert(type(biweibull_fi_scale1)==float)
        assert(type(biweibull_fi_scale2)==float)
        assert(type(biweibull_w_shape1)==float)
        assert(type(biweibull_w_shape2)==float)
        assert(type(biweibull_mixing_factor)==float)
        assert(type(custom_file)==str)
        assert(type(estimate_viscosity)==str)
        assert(type(estimate_height_above_vent)==str)
    
        self.number_of_bins = number_of_bins
        self.fi_range1 = fi_range1
        self.fi_range2 = fi_range2 
        self.density_range1 = density_range1
        self.density_range2 = density_range2 
        self.sphericity_range1 = sphericity_range1 
        self.sphericity_range2 = sphericity_range2
        self.distribution = distribution 
        self.gaussian_fi_mean = gaussian_fi_mean
        self.gaussian_fi_disp = gaussian_fi_disp
        self.bigaussian_fi_mean1 = bigaussian_fi_mean1 
        self.bigaussian_fi_mean2 = bigaussian_fi_mean2 
        self.bigaussian_fi_disp1 = bigaussian_fi_disp1
        self.bigaussian_fi_disp2 = bigaussian_fi_disp2 
        self.bigaussian_mixing_factor = bigaussian_mixing_factor 
        self.weibull_fi_scale =  weibull_fi_scale
        self.weibull_w_shape = weibull_w_shape 
        self.biweibull_fi_scale1 = biweibull_fi_scale1 
        self.biweibull_fi_scale2 = biweibull_fi_scale2 
        self.biweibull_w_shape1 = biweibull_w_shape1 
        self.biweibull_w_shape2 = biweibull_w_shape2
        self.biweibull_mixing_factor = biweibull_mixing_factor 
        self.custom_file = custom_file 
        self.estimate_viscosity = estimate_viscosity 
        self.estimate_height_above_vent = estimate_height_above_vent
    
    
    
    
class ParticleAggregation(Section):
    
    fstring = fstring_particle_aggregation
    
    regex = regex_particle_aggregation
    
    types = {
        'particle_cut_off':str,
        'aggregation_model':str,
        'number_of_aggregate_bins':int,
        'diameter_aggregates1':float,
        'diameter_aggregates2':float,
        'density_aggregates1':float,
        'density_aggregates2':float,
        'percentage1':float,
        'percentage2':float,
        'vset_factor':float,
        'fractal_exponent':float
    }
    
    def __init__(self,
                particle_cut_off:str,
                aggregation_model:str,
                number_of_aggregate_bins:int,
                diameter_aggregates1:float,
                diameter_aggregates2:float,
                density_aggregates1:float,
                density_aggregates2:float,
                percentage1:float,
                percentage2:float,
                vset_factor:float,
                fractal_exponent:float
                ):
        
        
        assert(type(particle_cut_off)==str)
        assert(type(aggregation_model)==str)
        assert(type(number_of_aggregate_bins)==int)
        assert(type(diameter_aggregates1)==float)
        assert(type(diameter_aggregates2)==float)
        assert(type(density_aggregates1)==float)
        assert(type(density_aggregates2)==float)
        assert(type(percentage1)==float)
        assert(type(percentage2)==float)
        assert(type(vset_factor)==float)
        assert(type(fractal_exponent)==float)
        


        self.particle_cut_off = particle_cut_off
        self.aggregation_model = aggregation_model
        self.number_of_aggregate_bins = number_of_aggregate_bins
        self.diameter_aggregates1 = diameter_aggregates1
        self.diameter_aggregates2 = diameter_aggregates2
        self.density_aggregates1 = density_aggregates1
        self.density_aggregates2 = density_aggregates2 
        self.percentage1 = percentage1 
        self.percentage2 = percentage2 
        self.vset_factor = vset_factor 
        self.fractal_exponent = fractal_exponent
        
        
        
        
        
        
        
     
    
class Source(Section):
    
    fstring = fstring_source
    
    regex = regex_source
    
    types = {
        'source_type':str,
        'source_start':int,
        'source_end':int,
        'lon_vent':float,
        'lat_vent':float,
        'vent_height':float,
        'height_above_vent':float,
        'mass_flow_rate':float,
        'alfa_plume':float,
        'beta_plume':float,
        'exit_temperature':float,
        'exit_water_fraction':float,
        'a1':float,
        'a2':float,
        'l':float,
        'thickness':float,
        'solve_plume_for':str,
        'mfr_search_range1':int,
        'mfr_search_range2':int,
        'exit_velocity':float,
        'exit_gas_water_temperature':float,
        'exit_liquid_water_temperature':float,
        'exit_solid_water_temperature':float,
        'exit_gas_water_fraction':float,
        'exit_liquid_water_fraction':float,
        'exit_solid_water_fraction':float,
        'wind_coupling':YesNo,
        'air_moisture':YesNo,
        'latent_heat':YesNo,
        'reentrainment':YesNo,
        'bursik_factor':float,
        'z_min_wind':float,
        'c_umbrella':float,
        'a_s':str,
        'a_v':str      
    }
    
    def __init__(self,
                source_type:str,
                source_start:int,
                source_end:int,
                lon_vent:float,
                lat_vent:float,
                vent_height:float,
                height_above_vent:float,
                mass_flow_rate:float,
                alfa_plume:float,
                beta_plume:float,
                exit_temperature:float,
                exit_water_fraction:float,
                a1:float,
                a2:float,
                l:float,
                thickness:float,
                solve_plume_for:str,
                mfr_search_range1:int,
                mfr_search_range2:int,
                exit_velocity:float,
                exit_gas_water_temperature:float,
                exit_liquid_water_temperature:float,
                exit_solid_water_temperature:float,
                exit_gas_water_fraction:float,
                exit_liquid_water_fraction:float,
                exit_solid_water_fraction:float,
                wind_coupling:YesNo,
                air_moisture:YesNo,
                latent_heat:YesNo,
                reentrainment:YesNo,
                bursik_factor:float,
                z_min_wind:float,
                c_umbrella:float,
                a_s:str,
                a_v:str      
                ):
    
        assert(type(source_type)==str)
        assert(type(source_start)==int)
        assert(type(source_end)==int)
        assert(type(lon_vent)==float)
        assert(type(lat_vent)==float)
        assert(type(vent_height)==float)
        assert(type(height_above_vent)==float)
        assert(type(mass_flow_rate)==float)
        assert(type(alfa_plume)==float)
        assert(type(beta_plume)==float)
        assert(type(exit_temperature)==float)
        assert(type(exit_water_fraction)==float)
        assert(type(a1)==float)
        assert(type(a2)==float)
        assert(type(l)==float)
        assert(type(thickness)==float)
        assert(type(solve_plume_for)==str)
        assert(type(mfr_search_range1)==int)
        assert(type(mfr_search_range2)==int)
        assert(type(exit_velocity)==float)
        assert(type(exit_gas_water_temperature)==float)
        assert(type(exit_liquid_water_temperature)==float)
        assert(type(exit_solid_water_temperature)==float)
        assert(type(exit_gas_water_fraction)==float)
        assert(type(exit_liquid_water_fraction)==float)
        assert(type(exit_solid_water_fraction)==float)
        assert(type(wind_coupling)==YesNo)
        assert(type(air_moisture)==YesNo)
        assert(type(latent_heat)==YesNo)
        assert(type(reentrainment)==YesNo)
        assert(type(bursik_factor)==float)
        assert(type(z_min_wind)==float)
        assert(type(c_umbrella)==float)
        assert(type(a_s)==str)
        assert(type(a_v)==str)     
        
  
        self.source_type = source_type
        self.source_start = source_start
        self.source_end = source_end
        self.lon_vent = lon_vent
        self.lat_vent = lat_vent
        self.vent_height = vent_height
        self.height_above_vent = height_above_vent
        self.mass_flow_rate = mass_flow_rate
        self.alfa_plume = alfa_plume
        self.beta_plume = beta_plume
        self.exit_temperature = exit_temperature
        self.exit_water_fraction = exit_water_fraction
        self.a1 = a1
        self.a2 = a2
        self.l = l
        self.thickness = thickness
        self.solve_plume_for = solve_plume_for
        self.mfr_search_range1 = mfr_search_range1
        self.mfr_search_range2 = mfr_search_range2
        self.exit_velocity = exit_velocity
        self.exit_gas_water_temperature = exit_gas_water_temperature
        self.exit_liquid_water_temperature = exit_liquid_water_temperature
        self.exit_solid_water_temperature = exit_solid_water_temperature
        self.exit_gas_water_fraction = exit_gas_water_fraction
        self.exit_liquid_water_fraction = exit_liquid_water_fraction
        self.exit_solid_water_fraction =  exit_solid_water_fraction
        self.wind_coupling = YesNo(wind_coupling)
        self.air_moisture = YesNo(air_moisture)
        self.latent_heat = YesNo(latent_heat) 
        self.reentrainment = YesNo(reentrainment) 
        self.bursik_factor = bursik_factor 
        self.z_min_wind = z_min_wind
        self.c_umbrella = c_umbrella
        self.a_s = a_s
        self.a_v = a_v  

    def plot_on_map(self,ax=None):


        if  ax==None:
            plt.figure("Test Map")
            ccrs.PlateCarree()
            crs = ccrs.PlateCarree()
            extent = [self.lon_vent-1, self.lon_vent+1,self.lat_vent-1, self.lat_vent+1]

            ax = plt.subplot(111, projection=crs)
            ax.set_extent(extent, crs=crs)
            
            #ax.add_feature(NaturalEarthFeature('physical', 'ocean', '50m'))
            ax.coastlines(resolution='10m',color='blue')
    
            ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)

        ax.plot([self.lon_vent], [self.lat_vent], marker='o', color='red')
            
        
        
        
        
        
class Ensemble(Section):
    
    fstring = fstring_ensemble
    
    regex = regex_ensemble
    
    types = {
        'random_numbers_from_file':YesNo,
        'perturbate_column_height':str,
        'column_height_perturbation_range':int,
        'column_height_pdf':str,
        'perturbate_mass_flow_rate':str,
        'mass_flow_rate_perturbation_range':int,
        'mass_flow_rate_pdf':str,
        'perturbate_source_start':str,
        'perturbate_source_start_range':int,
        'perturbate_source_start_pdf':str,
        'perturbate_source_duration':str,
        'perturbate_source_duration_range':int,
        'perturbate_source_duration_pdf':str,
        'perturbate_top_hat_thickness':str,
        'perturbate_top_hat_thickness_range':int,
        'perturbate_top_hat_pdf':str,
        'perturbate_suzuki_a':str,
        'perturbate_suzuki_a_range':int,
        'perturbate_suzuki_a_pdf':str,
        'perturbate_suzuki_l':str,
        'perturbate_suzuki_l_range':int,
        'perturbate_suzuki_l_pdf':str,
        'perturbate_wind':str,
        'perturbate_wind_range':int,
        'perturbate_wind_pdf':str,
        'perturbate_data_insertion_cloud_height':str,
        'perturbate_data_insertion_cloud_height_range':int,
        'perturbate_data_insertion_cloud_height_pdf':str,
        'perturbate_data_insertion_cloud_thickness':str,
        'perturbate_data_insertion_cloud_thickness_range':int,
        'perturbate_data_insertion_cloud_thickness_pdf':str,
        'perturbate_fi_mean':str,
        'perturbate_fi_range':int,
        'perturbate_fi_pdf':str,
        'perturbate_diamater_aggregates':YesNo,
        'perturbate_diamater_aggregates_range':int,
        'perturbate_diamater_aggregates_pdf':str,
        'perturbate_density_aggregates':str,
        'perturbate_density_aggregates_range':int,
        'perturbate_density_aggregates_pdf':str,
    }
    
    def __init__(self,
                random_numbers_from_file:YesNo,
                perturbate_column_height:str,
                column_height_perturbation_range:int,
                column_height_pdf:str,
                perturbate_mass_flow_rate:str,
                mass_flow_rate_perturbation_range:int,
                mass_flow_rate_pdf:str,
                perturbate_source_start:str,
                perturbate_source_start_range:int,
                perturbate_source_start_pdf:str,
                perturbate_source_duration:str,
                perturbate_source_duration_range:int,
                perturbate_source_duration_pdf:str,
                perturbate_top_hat_thickness:str,
                perturbate_top_hat_thickness_range:int,
                perturbate_top_hat_pdf:str,
                perturbate_suzuki_a:str,
                perturbate_suzuki_a_range:int,
                perturbate_suzuki_a_pdf:str,
                perturbate_suzuki_l:str,
                perturbate_suzuki_l_range:int,
                perturbate_suzuki_l_pdf:str,
                perturbate_wind:str,
                perturbate_wind_range:int,
                perturbate_wind_pdf:str,
                perturbate_data_insertion_cloud_height:str,
                perturbate_data_insertion_cloud_height_range:int,
                perturbate_data_insertion_cloud_height_pdf:str,
                perturbate_data_insertion_cloud_thickness:str,
                perturbate_data_insertion_cloud_thickness_range:int,
                perturbate_data_insertion_cloud_thickness_pdf:str,
                perturbate_fi_mean:str,
                perturbate_fi_range:int,
                perturbate_fi_pdf:str,
                perturbate_diamater_aggregates:YesNo,
                perturbate_diamater_aggregates_range:int,
                perturbate_diamater_aggregates_pdf:str,
                perturbate_density_aggregates:str,
                perturbate_density_aggregates_range:int,
                perturbate_density_aggregates_pdf:str,
                ):
        
        
        
        assert(type(random_numbers_from_file)==YesNo)
        assert(type(perturbate_column_height)==str)
        assert(type(column_height_perturbation_range)==int)
        assert(type(column_height_pdf)==str)
        assert(type(perturbate_mass_flow_rate)==str)
        assert(type(mass_flow_rate_perturbation_range)==int)
        assert(type(mass_flow_rate_pdf)==str)
        assert(type(perturbate_source_start)==str)
        assert(type(perturbate_source_start_range)==int)
        assert(type(perturbate_source_start_pdf)==str)
        assert(type(perturbate_source_duration)==str)
        assert(type(perturbate_source_duration_range)==int)
        assert(type(perturbate_source_duration_pdf)==str)
        assert(type(perturbate_top_hat_thickness)==str)
        assert(type(perturbate_top_hat_thickness_range)==int)
        assert(type(perturbate_top_hat_pdf)==str)
        assert(type(perturbate_suzuki_a)==str)
        assert(type(perturbate_suzuki_a_range)==int)
        assert(type(perturbate_suzuki_a_pdf)==str)
        assert(type(perturbate_suzuki_l)==str)
        assert(type(perturbate_suzuki_l_range)==int)
        assert(type(perturbate_suzuki_l_pdf)==str)
        assert(type(perturbate_wind)==str)
        assert(type(perturbate_wind_range)==int)
        assert(type(perturbate_wind_pdf)==str)
        assert(type(perturbate_data_insertion_cloud_height)==str)
        assert(type(perturbate_data_insertion_cloud_height_range)==int)
        assert(type(perturbate_data_insertion_cloud_height_pdf)==str)
        assert(type(perturbate_data_insertion_cloud_thickness)==str)
        assert(type(perturbate_data_insertion_cloud_thickness_range)==int)
        assert(type(perturbate_data_insertion_cloud_thickness_pdf)==str)
        assert(type(perturbate_fi_mean)==str)
        assert(type(perturbate_fi_range)==int)
        assert(type(perturbate_fi_pdf)==str)
        assert(type(perturbate_diamater_aggregates)==YesNo)
        assert(type(perturbate_diamater_aggregates_range)==int)
        assert(type(perturbate_diamater_aggregates_pdf)==str)
        assert(type(perturbate_density_aggregates)==str)
        assert(type(perturbate_density_aggregates_range)==int)
        assert(type(perturbate_density_aggregates_pdf)==str)
        
        self.random_numbers_from_file = YesNo(random_numbers_from_file)
        self.perturbate_column_height = perturbate_column_height
        self.column_height_perturbation_range  = column_height_perturbation_range
        self.column_height_pdf  = column_height_pdf
        self.perturbate_mass_flow_rate  = perturbate_mass_flow_rate
        self.mass_flow_rate_perturbation_range  = mass_flow_rate_perturbation_range
        self.mass_flow_rate_pdf  = mass_flow_rate_pdf
        self.perturbate_source_start  = perturbate_source_start
        self.perturbate_source_start_range  = perturbate_source_start_range
        self.perturbate_source_start_pdf  = perturbate_source_start_pdf
        self.perturbate_source_duration  = perturbate_source_duration
        self.perturbate_source_duration_range  = perturbate_source_duration_range
        self.perturbate_source_duration_pdf  = perturbate_source_duration_pdf
        self.perturbate_top_hat_thickness  = perturbate_top_hat_thickness
        self.perturbate_top_hat_thickness_range  = perturbate_top_hat_thickness_range
        self.perturbate_top_hat_pdf  = perturbate_top_hat_pdf
        self.perturbate_suzuki_a  = perturbate_suzuki_a
        self.perturbate_suzuki_a_range  = perturbate_suzuki_a_range
        self.perturbate_suzuki_a_pdf = perturbate_suzuki_a_pdf
        self.perturbate_suzuki_l = perturbate_suzuki_l
        self.perturbate_suzuki_l_range  = perturbate_suzuki_l_range
        self.perturbate_suzuki_l_pdf  = perturbate_suzuki_l_pdf
        self.perturbate_wind  = perturbate_wind
        self.perturbate_wind_range  = perturbate_wind_range
        self.perturbate_wind_pdf  = perturbate_wind_pdf
        self.perturbate_data_insertion_cloud_height  = perturbate_data_insertion_cloud_height
        self.perturbate_data_insertion_cloud_height_range  = perturbate_data_insertion_cloud_height_range
        self.perturbate_data_insertion_cloud_height_pdf  = perturbate_data_insertion_cloud_height_pdf
        self.perturbate_data_insertion_cloud_thickness =  perturbate_data_insertion_cloud_thickness
        self.perturbate_data_insertion_cloud_thickness_range  = perturbate_data_insertion_cloud_thickness_range
        self.perturbate_data_insertion_cloud_thickness_pdf  = perturbate_data_insertion_cloud_thickness_pdf
        self.perturbate_fi_mean  = perturbate_fi_mean
        self.perturbate_fi_range  = perturbate_fi_range
        self.perturbate_fi_pdf  = perturbate_fi_pdf
        self.perturbate_diamater_aggregates  = YesNo(perturbate_diamater_aggregates)
        self.perturbate_diamater_aggregates_range  = perturbate_diamater_aggregates_range
        self.perturbate_diamater_aggregates_pdf  = perturbate_diamater_aggregates_pdf
        self.perturbate_density_aggregates  = perturbate_density_aggregates
        self.perturbate_density_aggregates_range  = perturbate_density_aggregates_range
        self.perturbate_density_aggregates_pdf  = perturbate_density_aggregates_pdf

class TimeUTC(Section):
    
    fstring = fstring_time_utc

    regex = regex_time_utc
    
    types = {
        'year':int,
        'month':int,
        'day':int,
        'run_start':int,
        'run_end':int,
        'initial_condition':str,
        'restart_file':str,
        'restart_ensemble_basepath':str
    }
    
    
    
    def __init__(self, 
                year:int, 
                month:int, 
                day:int, 
                run_start:int, 
                run_end:int, 
                initial_condition:str, 
                restart_file:str, 
                restart_ensemble_basepath:str
                ):
        
        #---------------------------------------------------------------------------------
        # TEST INPUTS
        #---------------------------------------------------------------------------------
        
        # tests for year
        assert(type(year)==int)
        assert(year>=1970)
        assert(year<=2024)
        
        # tests for month
        assert(type(month)==int)
        assert(month>=1)
        assert(month<=12)
        
        # tests for day
        assert(type(day)==int)
        assert(day>=1)
        assert(day<=31)
        
        # tests for initial_condition
        assert(type(initial_condition)==str)
        assert(initial_condition in ['NONE', 'RESTART', 'INSERTION'])
        
        
        # tests for restart_file
        assert(type(restart_file)==str)
        
        # tests for restart_ensemble_basepath
        assert(type(restart_ensemble_basepath)==str)
        
        #---------------------------------------------------------------------------------       
        # OPTIONS THAT HAVE BECOME ACTIVATED
        #---------------------------------------------------------------------------------

        if initial_condition in ['RESTART']:
            print("RESTART_FILE in use as INITIAL_CONDITION = RESTART  )")
            
        self.year = year
        self.month = month
        self.day = day
        self.run_start = run_start
        self.run_end = run_end
        self.initial_condition = initial_condition
        self.restart_file = restart_file
        self.restart_ensemble_basepath = restart_ensemble_basepath  
        


class EnsemblePostprocess(Section):

    fstring = fstring_ensemble_postprocess
    
    regex = regex_ensemble_postprocess

    types = {
        'postprocess_members':YesNo,
        'postprocess_mean':YesNo,
        'postprocess_logmean':YesNo,
        'postprocess_median':YesNo,
        'postprocess_standard_dev':YesNo,
        'postprocess_probability':YesNo,
        'postprocess_percentiles':YesNo,
        'postprocess_probability_concentration_thresholds':int,
        'postprocess_probability_column_mass_thresholds_gm2':int,
        'postprocess_probability_column_mass_thresholds_du':int,
        'postprocess_probability_ground_load_thresholds':int ,
        'postprocess_percentiles_percentile_values':int  
    }


    def __init__(self,
                postprocess_members:YesNo,
                postprocess_mean:YesNo,
                postprocess_logmean:YesNo,
                postprocess_median:YesNo,
                postprocess_standard_dev:YesNo,
                postprocess_probability:YesNo,
                postprocess_percentiles:YesNo,
                postprocess_probability_concentration_thresholds:int,
                postprocess_probability_column_mass_thresholds_gm2:int,
                postprocess_probability_column_mass_thresholds_du:int,
                postprocess_probability_ground_load_thresholds:int ,
                postprocess_percentiles_percentile_values:int  
                ):

        assert(type(postprocess_members)==YesNo)
        assert(type(postprocess_mean)==YesNo)
        assert(type(postprocess_logmean)==YesNo)
        assert(type(postprocess_median)==YesNo)
        assert(type(postprocess_standard_dev)==YesNo)
        assert(type(postprocess_probability)==YesNo)
        assert(type(postprocess_percentiles)==YesNo)
        assert(type(postprocess_probability_concentration_thresholds)==int)
        assert(type(postprocess_probability_column_mass_thresholds_gm2)==int)
        assert(type(postprocess_probability_column_mass_thresholds_du)==int)
        assert(type(postprocess_probability_ground_load_thresholds)==int )
        assert(type(postprocess_percentiles_percentile_values)==int)  

        self.postprocess_members = YesNo(postprocess_members)
        self.postprocess_mean = YesNo(postprocess_mean)
        self.postprocess_logmean = YesNo(postprocess_logmean)
        self.postprocess_median = YesNo(postprocess_median)
        self.postprocess_standard_dev = YesNo(postprocess_standard_dev)
        self.postprocess_probability = YesNo(postprocess_probability)
        self.postprocess_percentiles = YesNo(postprocess_percentiles)
        self.postprocess_probability_concentration_thresholds = postprocess_probability_concentration_thresholds
        self.postprocess_probability_column_mass_thresholds_gm2 = postprocess_probability_column_mass_thresholds_gm2 
        self.postprocess_probability_column_mass_thresholds_du = postprocess_probability_column_mass_thresholds_du
        self.postprocess_probability_ground_load_thresholds = postprocess_probability_ground_load_thresholds 
        self.postprocess_percentiles_percentile_values = postprocess_percentiles_percentile_values

class ModelPhysics(Section):

    fstring = fstring_model_physics
    
    regex = regex_model_physics

    types = {
        'limiter':str,
        'time_marching':str,
        'cfl_criterion':str,
        'cfl_safety_factor':float,
        'terminal_velocity_model':str,
        'horizontal_turbulence_model':str,
        'vertical_turbulence_model':float,
        'rams_cs':float,
        'wet_deposition':YesNo,
        'dry_deposition':YesNo,
        'gravity_current':YesNo,
        'c_flow_rate':str,
        'lambda_grav':float,
        'k_entrain':float,
        'brunt_vaisala':float,
        'gc_start':int,
        'gc_end':int
    }

    def __init__(self,
                limiter:str,
                time_marching:str,
                cfl_criterion:str,
                cfl_safety_factor:float,
                terminal_velocity_model:str,
                horizontal_turbulence_model:str,
                vertical_turbulence_model:float,
                rams_cs:float,
                wet_deposition:YesNo,
                dry_deposition:YesNo,
                gravity_current:YesNo,
                c_flow_rate:str,
                lambda_grav:float,
                k_entrain:float,
                brunt_vaisala:float,
                gc_start:int,
                gc_end:int
                ):
        
        assert(type(limiter)==str)
        assert(type(time_marching)==str)
        assert(type(cfl_criterion)==str)
        assert(type(cfl_safety_factor)==float)
        assert(type(terminal_velocity_model)==str)
        assert(type(horizontal_turbulence_model)==str)
        assert(type(vertical_turbulence_model)==float)
        assert(type(rams_cs)==float)
        assert(type(wet_deposition)==YesNo)
        assert(type(dry_deposition)==YesNo)
        assert(type(gravity_current)==YesNo)
        assert(type(c_flow_rate)==str)
        assert(type(lambda_grav)==float)
        assert(type(k_entrain)==float)
        assert(type(brunt_vaisala)==float)
        assert(type(gc_start)==int)
        assert(type(gc_end)==int)


        self.limiter = limiter
        self.time_marching = time_marching
        self.cfl_criterion = cfl_criterion
        self.cfl_safety_factor = cfl_safety_factor
        self.terminal_velocity_model = terminal_velocity_model
        self.horizontal_turbulence_model = horizontal_turbulence_model
        self.vertical_turbulence_model = vertical_turbulence_model
        self.rams_cs = rams_cs
        self.wet_deposition = YesNo(wet_deposition)
        self.dry_deposition = YesNo(dry_deposition)
        self.gravity_current = YesNo(gravity_current)
        self.c_flow_rate = c_flow_rate
        self.lambda_grav = lambda_grav
        self.k_entrain = k_entrain
        self.brunt_vaisala = brunt_vaisala
        self.gc_start = gc_start 
        self.gc_end = gc_end




class ModelOutput(Section):

    fstring = fstring_model_output
    
    regex = regex_model_output

    types = {
        'parallel_io':YesNo,
        'log_file_level':str,
        'restart_time_interval':str,
        'output_intermediate_files':YesNo,
        'output_time_start':str,
        'output_time_interval':float,
        'output_3d_concentration':YesNo,
        'output_3d_concentration_bins':YesNo,
        'output_surface_concentration':YesNo,
        'output_column_load':YesNo,
        'output_cloud_top':YesNo,
        'output_ground_load':YesNo,
        'output_ground_load_bins':YesNo,
        'output_wet_deposition':YesNo,
        'output_track_points':YesNo,
        'output_track_points_file':str,
        'output_concentrations_at_xcuts':YesNo,
        'x_values':str,
        'output_concentrations_at_ycuts':YesNo,
        'y_values':str,
        'output_concentrations_at_zcuts':YesNo,
        'z_values':str,
        'output_concentrations_at_fl':YesNo,
        'fl_values':str,
    }
        
    def __init__(self,
            parallel_io:YesNo,
            log_file_level:str,
            restart_time_interval:str,
            output_intermediate_files:YesNo,
            output_time_start:str,
            output_time_interval:float,
            output_3d_concentration:YesNo,
            output_3d_concentration_bins:YesNo,
            output_surface_concentration:YesNo,
            output_column_load:YesNo,
            output_cloud_top:YesNo,
            output_ground_load:YesNo,
            output_ground_load_bins:YesNo,
            output_wet_deposition:YesNo,
            output_track_points:YesNo,
            output_track_points_file:str,
            output_concentrations_at_xcuts:YesNo,
            x_values:str,
            output_concentrations_at_ycuts:YesNo,
            y_values:str,
            output_concentrations_at_zcuts:YesNo,
            z_values:str,
            output_concentrations_at_fl:YesNo,
            fl_values:str,
            ):


        assert(type(parallel_io)==YesNo)
        assert(type(log_file_level)==str)
        assert(type(restart_time_interval)==str)
        assert(type(output_intermediate_files)==YesNo)
        assert(type(output_time_start)==str)
        assert(type(output_time_interval)==float)
        assert(type(output_3d_concentration)==YesNo)
        assert(type(output_3d_concentration_bins)==YesNo)
        assert(type(output_surface_concentration)==YesNo)
        assert(type(output_column_load)==YesNo)
        assert(type(output_cloud_top)==YesNo)
        assert(type(output_ground_load)==YesNo)
        assert(type(output_ground_load_bins)==YesNo)
        assert(type(output_wet_deposition)==YesNo)
        assert(type(output_track_points)==YesNo)
        assert(type(output_track_points_file)==str)
        assert(type(output_concentrations_at_xcuts)==YesNo)
        assert(type(x_values)==str)
        assert(type(output_concentrations_at_ycuts)==YesNo)
        assert(type(y_values)==str)
        assert(type(output_concentrations_at_zcuts)==YesNo)
        assert(type(z_values)==str)
        assert(type(output_concentrations_at_fl)==YesNo)
        assert(type(fl_values)==str)


        self.parallel_io = parallel_io
        self.log_file_level = log_file_level
        self.restart_time_interval = restart_time_interval
        self.output_intermediate_files = output_intermediate_files
        self.output_time_start = output_time_start
        self.output_time_interval = output_time_interval
        self.output_3d_concentration = YesNo(output_3d_concentration)
        self.output_3d_concentration_bins = YesNo(output_3d_concentration_bins)
        self.output_surface_concentration = YesNo(output_surface_concentration)
        self.output_column_load = YesNo(output_column_load)
        self.output_cloud_top = YesNo(output_cloud_top)
        self.output_ground_load = YesNo(output_ground_load)
        self.output_ground_load_bins = YesNo(output_ground_load_bins)
        self.output_wet_deposition = YesNo(output_wet_deposition)
        self.output_track_points = YesNo(output_track_points)
        self.output_track_points_file = output_track_points_file
        self.output_concentrations_at_xcuts = YesNo(output_concentrations_at_xcuts)
        self.x_values = x_values
        self.output_concentrations_at_ycuts = YesNo(output_concentrations_at_ycuts)
        self.y_values = y_values
        self.output_concentrations_at_zcuts = YesNo(output_concentrations_at_zcuts)
        self.z_values = z_values
        self.output_concentrations_at_fl = YesNo(output_concentrations_at_fl)
        self.fl_values = fl_values




class ModelValidation(Section):

    fstring = fstring_model_validation
    
    regex = regex_model_validation

    types = {
        'observations_type':str,
        'observations_file':str,
        'observations_dictionary_file':str,
        'results_file':str,
        'column_mass_observations_threshold_gm2':float,
        'column_mass_observations_threshold_du':int,
        'ground_load_observation_threshold_kgm2':float 
    }

    def __init__(self,
                observations_type:str,
                observations_file:str,
                observations_dictionary_file:str,
                results_file:str,
                column_mass_observations_threshold_gm2:float,
                column_mass_observations_threshold_du:int,
                ground_load_observation_threshold_kgm2:float, 
                ):
        
        assert(type(observations_type)==str)
        assert(type(observations_file)==str)
        assert(type(observations_dictionary_file)==str)
        assert(type(results_file)==str)
        assert(type(column_mass_observations_threshold_gm2)==float)
        assert(type(column_mass_observations_threshold_du)==int)
        assert(type(ground_load_observation_threshold_kgm2)==float)

        self.observations_type = observations_type
        self.observations_file = observations_file
        self.observations_dictionary_file = observations_dictionary_file
        self.results_file = results_file
        self.column_mass_observations_threshold_gm2 = column_mass_observations_threshold_gm2
        self.column_mass_observations_threshold_du = column_mass_observations_threshold_du
        self.ground_load_observation_threshold_kgm2 = ground_load_observation_threshold_kgm2 

class Fall3DInputFile:

    types = {
                'time_utc': TimeUTC,
                'insertion_data': InsertionData,
                'meteo_data': MeteoData,
                'grid':Grid,
                'species':Species,
                'tephra_tgsd':TephraTgsd,
                'radionucleides_tgsd':RadionucleidesTgsd,
                'particle_aggregation':ParticleAggregation,
                'source':Source,
                'ensemble':Ensemble,
                'emsemble_postprocess': EnsemblePostprocess,
                'model_physics':ModelPhysics,
                'model_output':ModelOutput,
                'model_validation':ModelValidation
        
    }

    def __init__(self,
                time_utc: TimeUTC,
                insertion_data: InsertionData,
                meteo_data: MeteoData,
                grid:Grid,
                species:Species,
                tephra_tgsd:TephraTgsd,
                radionucleides_tgsd:RadionucleidesTgsd,
                particle_aggregation:ParticleAggregation,
                source:Source,
                ensemble:Ensemble,
                emsemble_postprocess: EnsemblePostprocess,
                model_physics:ModelPhysics,
                model_output:ModelOutput,
                model_validation:ModelValidation,
                file=None
                ):
        
        assert(type(time_utc)==TimeUTC)
        assert(type(insertion_data)==InsertionData)
        assert(type(meteo_data)==MeteoData)
        assert(type(grid)==Grid)
        assert(type(species)==Species)
        assert(type(tephra_tgsd)==TephraTgsd)
        assert(type(radionucleides_tgsd)==RadionucleidesTgsd)
        assert(type(particle_aggregation)==ParticleAggregation)
        assert(type(source)==Source)
        assert(type(ensemble)==Ensemble)
        assert(type(emsemble_postprocess)== EnsemblePostprocess)
        assert(type(model_physics)==ModelPhysics)
        assert(type(model_output)==ModelOutput)
        assert(type(model_validation)==ModelValidation)

        self.time_utc = time_utc
        self.insertion_data = insertion_data
        self.meteo_data = meteo_data
        self.grid = grid
        self.species = species
        self.tephra_tgsd = tephra_tgsd
        self.radionucleides_tgsd = radionucleides_tgsd
        self.particle_aggregation = particle_aggregation
        self.source = source
        self.ensemble = ensemble
        self.emsemble_postprocess = emsemble_postprocess
        self.model_physics = model_physics
        self.model_output = model_output
        self.model_validation = model_validation
        self.file = file


    @classmethod
    def from_string(cls, string:str):
        

        
        return(cls(
            time_utc= TimeUTC.from_string(),
            insertion_data = InsertionData.from_string(),
            meteo_data = MeteoData.from_string(),
            grid = Grid.from_string(),
            species = Species.from_string(),
            tephra_tgsd = TephraTgsd.from_string(),
            radionucleides_tgsd = RadionucleidesTgsd.from_string(),
            particle_aggregation = ParticleAggregation.from_string(),
            source = Source.from_string(),
            ensemble = Ensemble.from_string(),
            emsemble_postprocess = EnsemblePostprocess.from_string(),
            model_physics = ModelPhysics.from_string(),
            model_output = ModelOutput.from_string(),
            model_validation = ModelValidation.from_string()
        ))

    @classmethod
    def from_file(cls, file:str):
        """New Fall3DInputFile from file
        """

        with open(file) as f:
            lines = f.readlines()

    
        string_timeutc = "".join(lines[16:38]) 
        string_insertiondata = "".join(lines[38:52]) 
        string_meteodata = "".join(lines[52:78]) 
        string_grid = "".join(lines[78:110]) 
        string_species = "".join(lines[110:152]) 
        string_tephratgsd = "".join(lines[152:186]) 
        string_radionucleidestgsd = "".join(lines[186:224]) 
        string_particleaggregation = "".join(lines[224:246]) 
        string_source = "".join(lines[246:338]) 
        string_ensemble = "".join(lines[338:440]) 
        string_ensemblepostprocess = "".join(lines[440:469]) 
        string_modelphysics = "".join(lines[469:513])
        string_modeloutput = "".join(lines[513:564]) 
        string_modelvalidation = "".join(lines[564:])
    


        time_utc = TimeUTC.from_string(
                            string_timeutc
                            )
        
        insertion_data = InsertionData.from_string(
                            string_insertiondata
                            )
        
        meteo_data = MeteoData.from_string(
                            string_meteodata
                            )
        
        grid = Grid.from_string(
                            string_grid
                            )
        
        species = Species.from_string(
                            string_species
                            )
        
        tephra_tgsd = TephraTgsd.from_string(
                            string_tephratgsd
                            )
        
        radionucleides_tgsd = RadionucleidesTgsd.from_string(
                            string_radionucleidestgsd
                            )
        
        particle_aggregation = ParticleAggregation.from_string(
                            string_particleaggregation
                            )
        
        source = Source.from_string(
                            string_source
                            )

        ensemble = Ensemble.from_string(
                            string_ensemble
                            )
        
        emsemble_postprocess =  EnsemblePostprocess.from_string(
                            string_ensemblepostprocess
                            )
        
        model_physics = ModelPhysics.from_string(
                            string_modelphysics
                            )
        
        model_output = ModelOutput.from_string(
                            string_modeloutput
                            )

        model_validation = ModelValidation.from_string(
                            string_modelvalidation
                            )

        f3dif = cls(
            time_utc = time_utc,
            insertion_data = insertion_data,
            meteo_data = meteo_data,
            grid = grid,
            species = species,
            tephra_tgsd = tephra_tgsd,
            radionucleides_tgsd = radionucleides_tgsd,
            particle_aggregation = particle_aggregation,
            source = source,
            ensemble = ensemble,
            emsemble_postprocess =  emsemble_postprocess,
            model_physics = model_physics,
            model_output = model_output,
            model_validation = model_validation,
            file=file
        )
        
        

        return(f3dif)

        
    def to_string(self):
        
        string = "".join([
            fstring_boilerplate,
            self.time_utc.to_string(),
            self.insertion_data.to_string(),
            self.meteo_data.to_string(),
            self.grid.to_string(),
            self.species.to_string(),
            self.tephra_tgsd.to_string(),
            self.radionucleides_tgsd.to_string(),
            self.particle_aggregation.to_string(),
            self.source.to_string(),
            self.ensemble.to_string(),
            self.emsemble_postprocess.to_string(),
            self.model_physics.to_string(),
            self.model_output.to_string(),
            self.model_validation.to_string()
        ])

        return(string)

    def to_file(self, file):

        string = self.to_string()

        with open(file,"w+") as f:
            f.writelines(string)
            
        self.output_file = file

    def update(self,new_values:dict):
        """Update with values in dict
        """

        sections = {}

        # we iterate through each section of the file ....
        for section_name, section_type in self.types.items():
            
            section = getattr(self,section_name)

            
            # ... for each section of the file we get the current params ...
        
            params = {}
    
            for attribute_name, attribute_type in section.types.items():
                
                params[attribute_name] = getattr(section,attribute_name)

            # ...we then check to see if there are any values in that section to update ...
            if section_name in new_values.keys():

                # ... and if there are, we update them
                for name, value in new_values[section_name].items():
                    params[name] = value

            # ... finally we reinitialise the section with the new params ..
            section.__init__(**params)

            # ... and append it to our sections dict.
            sections[section_name] = section

        
        # ... finally, we reinitialise with the update parameters
        self.__init__(**sections)
        
    def get_meteodata(self):
        """Fetches Meteo data based on specification in time_utc, grid and meteo_data
        """

        source = get_MeteoSource(self)

        source.get_fall3d_input(self)


    def plot_on_map(self):
        """Plots meteodata extent, grid extent and source location on a map with a high resolution
        coastlinbe.
        """


        
        plt.figure("Test Map")
        ccrs.PlateCarree()
        crs = ccrs.PlateCarree()
        ax = plt.subplot(111, projection=crs)        
        #ax.add_feature(NaturalEarthFeature('physical', 'ocean', '50m'))
        ax.coastlines(resolution='10m',color='blue')

        ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)

        extent = self.meteo_data.plot_on_map(ax=ax)
        self.grid.plot_on_map(ax=ax)
        self.source.plot_on_map(ax=ax)

        ax.set_extent(extent, crs=crs)
        






                             

def get_MeteoSource(file:Fall3DInputFile): #time_utc:TimeUTC, grid:Grid, 
    """General class to fetch
    """

    data_format = file.meteo_data.meteo_data_format

    switch = {
        'WRF':WRFSource,
        'ERA5':ERA5Source,
        'GFS':GFSSource,
        'IFS':IFSSource,
        'CARRA':CARRASource,
        'ERA5ML':ERA5MLSource
    }

    source = switch[data_format]() #time_utc, grid, meteo_data

    return(source)
        
class ERA5MLSource:

    def __init__(self):#,time_utc:TimeUTC, grid:Grid, meteo_data:MeteoData):
        raise NotImplementedError

class ERA5Source:

    def __init__(self):#,time_utc:TimeUTC, grid:Grid, meteo_data:MeteoData):
        raise NotImplementedError

class WRFSource:

    def __init__(self):#,time_utc:TimeUTC, grid:Grid, meteo_data:MeteoData):
        raise NotImplementedError

class GFSSource:

    def __init__(self):#,time_utc:TimeUTC, grid:Grid, meteo_data:MeteoData):
        raise NotImplementedError


class IFSSource:

    def __init__(self):#,time_utc:TimeUTC, grid:Grid, meteo_data:MeteoData):
        raise NotImplementedError

class CARRASource:
    
    # we need a file on the full native CARRA West domain grid
    # so that we have the latitude and longitude grids, as well
    # as the projection information in the attributes to hand
    #file_native = "/media/talfan/Maxtor/test_interpolated.nc"
    file_native = "mnt/aux/CARRA_orography_west.nc"

    
    def __init__(self, local_storage="mnt/archive/"):
        
        self.ds_native = xr.open_mfdataset(self.file_native)['orog'].drop(['time','step','surface','valid_time'])
        self.local_storage = local_storage

    def get_fall3d_input(self, file:Fall3DInputFile):
        """
        """
        
        time_utc = file.time_utc
        grid = file.grid
        meteo_data = file.meteo_data
        
        if os.path.exists(meteo_data.meteo_data_file):
        	raise ValueError("File already exists: "+meteo_data.meteo_data_file)

        # get extent from Grid object
        latmax = grid.latmax
        latmin = grid.latmin
        # remember! Longitudes are specified -180 to 180 in a grid object
        lonmax = Longitude(grid.lonmax,180)
        lonmin = Longitude(grid.lonmin,180)

        # get time info from TimeUTC object ...
        year = time_utc.year
        month = time_utc.month
        day = time_utc.day
        run_start = time_utc.run_start
        run_end = time_utc.run_end

        
        # ... convert to start and end datetimes ...
        start = datetime.datetime(year=year,month=month,day=day,hour=run_start)
        
        duration = run_end-run_start
        
        end = start + datetime.timedelta(hours=duration)

        # ... and get the months we need to order.
        # Fuirst, we get the difference between the two dates in seconds 
        # by differencing the unix timestamps ...
        seconds = (end.timestamp()-start.timestamp())

        # .. which we convert to decimal days ...
        days = seconds/(24*60*60)

        # ... round up to the nearest whole day ...
        days = np.ceil(days)

        # .. then convert to int for iterating over ...
        days = int(days)
        
        # ... to get a list of all daily dates between the start and end date ...
        months = []
        for day in range(days):
            date = start + datetime.timedelta(days=day)

            # ... and for each date we save the year and month ...
            months.append((date.year,date.month))

        # .. so that by performiung a set operation we can get the uniuque
        # year month pairs we need to order using the cdsapi ...
        months = set(months)

        # ... we order the data to get a list of datasets ...
        ds_months = [self.get_month(lonmin, lonmax, latmin, latmax, year, month) for year, month in months]

        # .. which we concatenate into a single one
        ds = xr.concat(ds_months,dim='time')

        # SUBSET BY DATE HERE!!!!

        # write to the file specified in the meteo_data section
        ds.to_netcdf(meteo_data.meteo_data_file)

        return(ds)
    
 
        
    def get_month(self,lonmin:Longitude, lonmax:Longitude, latmin:float, latmax:float, year:int, month:int):
        """NOTE! We order by the month due to CDSAPI putting orders longer than a month
         to the back of the queue!!!
         """
        # first we check local storage to see if we already have a file for that
        # year and month that covers our spatial extent
        
        file = self.check_local_storage(lonmin, lonmax, latmin, latmax, year, month)
        
        if file is None:
            print("No appropriate file found in local storage for ",year, month, lonmin._180, lonmax._180, latmin,latmax)
            print("Submitting order to CDSAPI")
            file = self.make_request(lonmin, lonmax, latmin, latmax, year, month)
            
        else:
            print("Local file found:",file)
        ds = xr.open_dataset(file)
        
        return(ds)
            

    def get_local_storage(self):
        
        # search for netcdfs in the local archive
        files = glob.glob("mnt/archive/*.nc")

        # if there are no .nc files, we return an empty dataframe ..
        if len(files) ==0:
            
            df = pd.DataFrame(columns=['file','year','month','lonmin','lonmax','latmin','latmax'])

        # ... otherwise we populate the dataframe with dates and extents
        else:
            
            data = [re.findall(r"(.*(\d{4})(\d{2})_(.+)_(.+)_(.+)_(.+)__CARRA\.nc)",file)[0]  for file in files]    

            df = pd.DataFrame(
                data, 
                columns=['file','year','month','lonmin','lonmax','latmin','latmax']
            )

            df['year'] = df['year'].astype(int)
            df['month'] = df['month'].astype(int)
            df['lonmin'] = df['lonmin'].astype(float)
            df['lonmax'] = df['lonmax'].astype(float)
            df['latmin'] = df['latmin'].astype(float)
            df['latmax'] = df['latmax'].astype(float)

        return(df)

    
            
    
    
    def check_local_storage(self,lonmin:Longitude, lonmax:Longitude, latmin:float, latmax:float, year:int, month:int):
        """Check to see if we already have a file that matches the request
        """
       
        df = self.get_local_storage()

        # if there are no files we wil return None, otherwise ...
        if len(df)==0:
            
            file = None

        # we check to see if we have a file that meets our requirements:
        else:
        
            check = (
                    lambda r: 
                        (r['year']== year)&
                        (r['month']== month)&
                        (r['lonmin']<= lonmin._360)&
                        (r['lonmax']>= lonmax._360)&
                        (r['latmin']<= latmin)&
                        (r['latmax']>= latmax)
                        )
    
            df['match'] = df.apply(check, axis=1)

            df_match = df[df['match']]

            # if there is 1 or more files that meets our requirements we will return the first one ...
            if len(df_match)>=1:
                file = df_match['file'].item()

            # ... otherwise we return None as before
            else:
                file = None
        
            
        return(file)

      
    
    def basic_request(self, lonmin:Longitude, lonmax:Longitude, latmin:float, latmax:float, year:int, month:int)->dict:
        """Gets the part of the cdsapi request dict that is common to both single and pressure level
        requests.
        """
        #print(lonmin._180)
        #print(lonmin._360)
        
        self.lonmin = lonmin
        self.lonmax = lonmax
        self.latmin = latmin
        self.latmax = latmax
        self.year = year
        self.month = month
        
        # then we get the filenames we will use - NOTE we use 0-360 long here
        self.filename = "_".join([
                str(year) + str(month).zfill(2),
                str(self.lonmin._360),
                str(self.lonmax._360),
                str(self.latmin),
                str(self.latmax),
                '_CARRA.nc'
            ])
        
        
        
        # a visualisation of the problem:
        # we need to order a CARRA subset 
        # specified in lat lon (outer box) 
        # big enough to contain the extent
        # in CARRA x-y coordinates (middle box)
        # that encompasse our domain of interest
        # (inner box)
        #
        #    *-----------------------*
        #    |       *               |
        #    |      +   +            |        
        #    |     +       +         |
        #    |    +*----------*      |    
        #    |   + |          |  +   |
        #    |  +  |          |    * |   
        #    | +   |          |   +  |  
        #    |*    |          |  +   |   
        #    |  +  *----------* +    |     
        #    |     +           +     |             
        #    |        +       +      |          
        #    |           +   +       |      
        #    |              *        |         
        #    *-----------------------*
        #
        
        
        # In other words, we need the range of lat lon to order from cdsapi, that
        # encompases a box in native grid (x,y) that in turn
        # encmopases the range of lat lon we want for Fall3D. It's complicated!
        # (native lon is 0-360)
        domain_xy=(
                    (self.ds_native.latitude>latmin)&
                    (self.ds_native.latitude<latmax)&
                    (self.ds_native.longitude>lonmin._360)&
                    (self.ds_native.longitude<lonmax._360)
                    )
        
        # we get the min and max x and y values
        # we will use to clip the data we download after
        # we order it
        xs = np.where(domain_xy.values.any(axis=0))[0]
        xmin = xs.min()
        xmax = xs.max()

        ys = np.where(domain_xy.values.any(axis=1))[0]
        ymin = ys.min()
        ymax = ys.max()
        
        # we buffer these values to be on the safe side
        xmin -= 2
        xmax += 2
        ymin -= 2
        ymax += 2

        # Now, we need to find a range of (lat, lon) values
        # that completely enclose this range of xy values

        domain_xy = domain_xy.sel(x=slice(xmin, xmax),y=slice(ymin, ymax))
        
        # we need to save domain_xy for resampling later
        self.domain_xy = domain_xy

        latmax_for_cdsapi = domain_xy.latitude.max().values.item()
        latmin_for_cdsapi = domain_xy.latitude.min().values.item()
        # longitude goes 0-360 here for some reason
        lonmax_for_cdsapi = Longitude( domain_xy.longitude.max().values.item(), 360)
        lonmin_for_cdsapi = Longitude(domain_xy.longitude.min().values.item(),360)

        # buffer the cdsapi bounds by 0.1 of a degree just in case
        latmax_for_cdsapi += 0.1 
        latmin_for_cdsapi -= 0.1 
        lonmax_for_cdsapi = Longitude( lonmax_for_cdsapi._360 + 0.1, 360)
        lonmin_for_cdsapi = Longitude( lonmin_for_cdsapi._360 - 0.1, 360)

        # finally, remember that longitude for cdsapi is -180-180, not 0-360
        #lon360_to_180 = lambda lon:( lon + 180) % 360 - 180
        #lonmax_for_cdsapi = lon360_to_180(lonmax_for_cdsapi)
        #lonmin_for_cdsapi = lon360_to_180(lonmin_for_cdsapi)

        # the subset area for cdsapi is specified like this
        # NOTE! lon is specified -180 - 180 here
        area = [latmax_for_cdsapi, lonmin_for_cdsapi._180, latmin_for_cdsapi, lonmax_for_cdsapi._180 ]
        
        # now we need to specify the resolution
        # this needs to be a number in degrees that is the result of dividing
        # 90 by an integer:
        num_samples_in_90_degrees = 3000

        resolution = 90/num_samples_in_90_degrees

        # we need the number of days in the month, which we find by subtracting
        # the 1st of the next month from the 1st of the current month
        start = datetime.datetime(year=year,month=month,day=1)

        if month==12:
            end = datetime.datetime(year=year+1,month=1,day=1)
        else:
            end = datetime.datetime(year=year,month=month+1,day=1)

        days_in_month = (end-start).days

        days = [day for day in range(1,days_in_month+1)]


        # now we turn all the date integers to strings
        day = [str(day).zfill(2) for day in days]

        year = str(year)

        month = str(month).zfill(2)

        # we always get all three hourly data for every day
        # these will never change
        time =[
                    '00:00', '03:00', '06:00',
                    '09:00', '12:00', '15:00',
                    '18:00', '21:00',
                ]


        # ... we will make two requests, one for 
        basic_request= {
            'format': 'grib',
            'domain': 'west_domain',
            'format': 'grib',
            'product_type': 'analysis',
            'grid':[resolution,resolution],
            'area': area,
            'time': time,
            'year': year,
            'month': month,
            'day': day
        }
        
        return(basic_request)
    
    def levels_request(self, lonmin:Longitude, lonmax:Longitude, latmin:float, latmax:float, year:int, month:int)->dict:
        
        levels_request = self.basic_request(lonmin, lonmax, latmin, latmax, year, month)
        
        levels_request['variable'] = [
                    'geometric_vertical_velocity','geopotential', 
                    'relative_humidity', 'temperature',
                    'u_component_of_wind', 'v_component_of_wind',
                ]
        
        levels_request['pressure_level'] = [
                                                '800', '825', '850',
                                                '875', '900', '925',
                                                '950', '1000',
                                            ]
        
        return(levels_request)
        

    def single_request(self, lonmin:Longitude, lonmax:Longitude, latmin:float, latmax:float, year:int, month:int)->dict:
    
        single_request = self.basic_request(lonmin, lonmax, latmin, latmax, year, month)
        
        single_request['variable']= [
                    '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_relative_humidity',
                    '2m_temperature', 'land_sea_mask',
                     'orography', 'surface_pressure', 'surface_roughness'
                ]
        
        single_request['level_type']= 'surface_or_atmosphere'
        
        return(single_request)
    
    def make_request(self, lonmin:Longitude, lonmax:Longitude, latmin:float, latmax:float, year:int, month:int):
        
        # first we check to see if we already have the data
        #file = local_storage(self,lonmin, lonmax, latmin, latmax, year, month)
        
        #if file is None:
        #    return(file)
                             
        # first we get the requests dicts
        levels_request = self.levels_request(lonmin, lonmax, latmin, latmax, year, month)
        
        single_request = self.single_request(lonmin, lonmax, latmin, latmax, year, month)
        
        # then we get the path to store the data
        path = os.path.join(self.local_storage, self.filename)
        
        # now we order the data. First we initialise the cdsapi client ...

        c = cdsapi.Client()

        # and we make the two requests
        c.retrieve(
            'reanalysis-carra-pressure-levels',
            levels_request,
            'pressure_levels.grib')

        c.retrieve(
            'reanalysis-carra-single-levels',
            single_request,
            'single_levels.grib')
        
        # Once we have the two gribs, we open them as xarray datasets ...
        gribs_as_ds = cfgrib.open_datasets('single_levels.grib')
        
        gribs_as_ds += cfgrib.open_datasets('pressure_levels.grib')

        # ... and iterate over them, interpolating them back to the CARRA grid
        # from the lat lon subset we we have ordered, so Fall3D can use it ...
        dss = []

        for i,ds in enumerate(gribs_as_ds):

            # resample
            ds_resampled = (
                            ds
                            .interp({
                                'latitude':self.domain_xy.latitude,
                                'longitude':self.domain_xy.longitude
                                    })
                            .drop(['surface','heightAboveGround'],errors='ignore')
                            )


            dss.append(ds_resampled)
            
        # ... merge the resampled datasets ...
        ds = xr.merge(dss)
        
        # ... and add the projection information back
        for name in ds:
            #ds[name].attrs = ds_works[name].attrs
            ds[name].attrs['GRIB_gridType']= 'lambert'
            ds[name].attrs['GRIB_gridDefinitionDescription']= 'Lambert conformal '
            ds[name].attrs['GRIB_LaDInDegrees'] = 72.0
            ds[name].attrs['GRIB_LoVInDegrees'] = 324.0
            ds[name].attrs['GRIB_DyInMetres']= 2500.0
            ds[name].attrs['GRIB_DxInMetres'] = 2500.0
            ds[name].attrs['GRIB_Latin2InDegrees'] = 72.0
            ds[name].attrs['GRIB_Latin1InDegrees']= 72.0
            #ds[name].attrs['GRIB_latitudeOfSouthernPoleInDegrees'] = 0.0
            #ds[name].attrs['GRIB_longitudeOfSouthernPoleInDegrees'] = 0.0
            
        ds.to_netcdf(path)
        
        return(path)


class Fall3DBatch:

    def __init__(self, name:str, basefile:str|Fall3DInputFile, df:pd.DataFrame, basedir="mnt/runs",n_parallel = 5):
        """Initialise a new batch run object
        """

        # name has to be a valid directory name
        # https://stackoverflow.com/questions/59672062/elegant-way-in-python-to-make-sure-a-string-is-suitable-as-a-filename
        ok = ".-_0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        
        if not all([ c in ok for c in name]):
            raise TypeError("Name variable contains characters illegal in a path")
        
        
        self.name = name
        self.basefile = basefile
        self.df = df
        self.basedir = basedir
        self.n_parallel = n_parallel

    def initialise(self):
        """Initialise a new batch run - create directories
        """

        # make the path that will itself contain a directory for each run
        self.path = os.path.join(self.basedir,self.name)
        
        # we don't want to risk overwriting previous batch runs by accident
        if os.path.isdir(self.path):
        	raise ValueError("Path already exists for batch name: "+self.name)
        
        os.mkdir(self.path)
        

        # save the dataframe containing the specification for each run
        self.df_file = os.path.join(self.path, self.name+".csv")
        
        if os.path.exists(self.df_file):
        	raise ValueError("File already exists for batch name: "+self.name)
        
        self.df.to_csv(self.df_file)

        # create directory and input file for each run
        self.input_files = []
        
        # each run is a row in the dataframe, so we iterate over them...
        for i, r in self.df.iterrows():
            
            # ... convert each row to a dict ...
            rdict =r.to_dict()

            # ... now the update function for a Fall3DinputFile object
            # takes a nexted dict, with one subdict for each sectiom.
            # So we restructure the row dict into a nested dict:
            update = {}
            
            for key, value in rdict.items():
            
                section, attribute = key.split('.')
                
                if section in update.keys():
            
                    update[section][attribute] = value
            
                else:
                    update[section] = {attribute:value}

            # now we have our nested dict, we make a copy of our base file
            # and update it with the new values. This reinitialises the object
            # so any update that is of wrong type / invalid value / clashes with
            # another setting in the inputfile should be caught
        
            newf = copy.deepcopy(self.basefile)
        
            newf.update(update)

            # now we need the folder to store this run in.
            # This is just the basepath we defined earlier
            # plus the uuid, which is the indes of our row, i:
        
            new_dir = os.path.join(self.path,i)
            
            if os.path.isdir(new_dir):
            	raise ValueError("Path already exists: "+new_dir)
            
            os.mkdir(new_dir)
            
            # and then save the Fall3D input file to that folder:
            input_file = os.path.join(new_dir,i+".inp")
            
            if os.path.exists(input_file):
            	raise ValueError("File already exists: "+input_file)
            
            newf.to_file(input_file)

            # save input file
            self.input_files.append(input_file)

    def get_meteo_data(self):
        """
        """

        for file in tqdm(self.input_files):

            f3if = Fall3DInputFile.from_file(file)
            
            meteosource = get_MeteoSource(f3if)
        
            meteosource.get_fall3d_input(
                    #time_utc=f3if.time_utc, 
                    #grid=f3if.grid, 
                    #meteo_data=f3if.meteo_data
                    f3if
                )

    def oldrun(self):
    


        for file in tqdm(self.input_files):

            subprocess.run([path_fall3d, "ALL",file]) 
            
    def run(self):
        

        processes = { }

        files = copy.deepcopy(self.input_files)

        # initialise first n_parallel runs
        for i in range(self.n_parallel):
            file = files.pop()
            processes[i] = subprocess.Popen(["/fall3d/bin/Fall3d.r8.x","All",file])

        # while there is more than 1 file left ....
        while len(files)>0:
        
            # ... iterate over each parallel slot ...
            for i, p in processes.items():
            
                # ... and when that run has finished ...
                if p.poll() is not None:
                
                    # ... get the next run in the list ...
                    file = files.pop()
                    
                    #print(len(files))
                    # ... put it in that slot and start it ...
                    processes[i] = subprocess.Popen(["/fall3d/bin/Fall3d.r8.x","All",file])
            
        # ... and once we've used up all the files we just need to wait until the last one is finished
        finished = False
       
        while not finished:
       
            # the batch is finished when all of the last Popen objects stop returning None
            finished  = all([(p.poll() is not None) for i, p in processes.items()])
           
        print("Batch completed")
       
            
            
    def get_meteo_and_run(self): 
    
    	for file in tqdm(self.input_files):
            f3if = Fall3DInputFile.from_file(file)
            
            print("********************************************************************")
            print("FETCHING METEO DATA")
            print("********************************************************************")
            meteosource = get_MeteoSource(f3if)
            
            print("********************************************************************")
            print("RUN FALL3D")
            print("********************************************************************")
        
            meteosource.get_fall3d_input(
                    #time_utc=f3if.time_utc, 
                    #grid=f3if.grid, 
                    #meteo_data=f3if.meteo_data
                    f3if
                )
                
            subprocess.run(["/fall3d/bin/Fall3d.r8.x", "All",file]) 

            


