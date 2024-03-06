from pyf3d import TimeUTC, Grid, MeteoData, Source, Fall3DInputFile, CARRASource, Fall3DBatch
import os
import copy
import uuid
import time
import numpy as np
import datetime
import pandas as pd

print("Running script", flush=True)
print("Also running script", flush=True)
print(os.listdir())
# initialise random number generator

# seed for test2
#seed = 12345

# seed for test3
#seed = 5676

# seed for test 4
#seed = 475979

# seed fpor test 5
seed = 39864
rng = np.random.default_rng(seed)

# we load a default SO2 input file
f3if = Fall3DInputFile.from_file("mnt/examples/default_so2_reykjanes.inp")

# update what we need to update
f3if.update({
		'model_output':{'output_track_points_file':'mnt/examples/stations.pts'},
		'meteo_data':{'meteo_data_dictionary_file':'mnt/examples/CARRA.tbl'}
		})

# decide how many training samples we want
size= 20

# Get random starting dates
start_date = datetime.datetime(year=1991,month=1,day=1)

end_date = datetime.datetime(year=2022,month=1,day=1)

interval_in_hours = (end_date - start_date).days * 24

random_hours = rng.integers( low=0, high=interval_in_hours, size=size)

random_dates  = [start_date + datetime.timedelta(hours=int(h)) for h in random_hours]

years = [d.year for d in random_dates]
months = [d.month for d in random_dates]
days = [d.day for d in random_dates]
run_starts = [d.hour for d in random_dates]

# Get random plume heights
min_height_above_vent = 500.0
max_height_above_vent = 1500.0
height_above_vents = rng.uniform(
                                size=size, 
                                low = min_height_above_vent, 
                                high = max_height_above_vent
                            )


# Get random gas fluxes
min_mass_flow_rate = 0.0
max_mass_flow_rate = 500.0
mass_flow_rates = rng.uniform(
                                size=size,
                                low = min_mass_flow_rate,
                                high = max_mass_flow_rate
                            )
#  get random identifier for run
uids = [str(uuid.uuid4()) for i in range(size)]

#name = "test2"#str(uuid.uuid4())
#name = "test4"
name = "test5"
meteo_data_files = [os.path.join("mnt","runs",name,u,u+"_meteo.nc") for u in uids]

df = pd.DataFrame({
        'time_utc.year':years,
        'time_utc.month':months,
        'time_utc.day':days,
        'time_utc.run_start':run_starts,
        'source.height_above_vent':height_above_vents,
        'source.mass_flow_rate':mass_flow_rates,
        'meteo_data.meteo_data_file':meteo_data_files
    },index=uids)
    
    


batch = Fall3DBatch(name=name, basefile=f3if, df=df, basedir="mnt/runs")

batch.initialise()

batch.get_meteo_data()

batch.run()

 
 
 
