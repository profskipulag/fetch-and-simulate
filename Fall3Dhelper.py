import re
from regex_strings import *
from fstrings import *




# awesome function by 'Anon' : https://stackoverflow.com/a/51981596
def is_valid_date(year, month, day):
    day_count_for_month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if year%4==0 and (year%100 != 0 or year%400==0):
        day_count_for_month[2] = 29
    return (1 <= month <= 12 and 1 <= day <= day_count_for_month[month])


class YesNo:
    """Class to associate various spellings of "yes" and "no" with True and False
    """
    
    def __init__(self,val):
        
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
    
    def __init__(self,val):
        
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
        # the value of the OnOff object that is inserted into the fstring
        return(self.val)
    
    def __str__(self):
        # the value of the YesNo object that is inserted into the fstring
        return(self.val)
    
    def __bool__(self):
        # the value of the OnOff object that is used for logical operations
        return(self.bool)
    




class Section:
    """Generic class for a Fall3D input file section.  
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
        
        variables = {}
        
        for key in self.types.keys():
            
            variables[key] = getattr(self,key)
        
        string = self.fstring.format(**variables)
        
        return(string)




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
                model_validation:ModelValidation
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
            model_validation = model_validation
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
        


class Vent:

    def __init__(self, lat, lon, Z):

        self.lat = lat
        self.lon = lon
        self.Z = Z


class RandomSourceGenerator(Source):

    def __init__(self):
        pass

