from Fall3Dhelper import *

#############################################################################
#
# Custom data type tests
#
#############################################################################


def test_yesno():

    assert(YesNo('YES'))
    assert(YesNo('Yes'))
    assert(YesNo('yes'))

    assert(not YesNo('NO'))
    assert(not YesNo('No'))
    assert(not YesNo('no'))

    assert(YesNo(True).__str__() == "Yes") 
    assert(YesNo(False).__str__() == "No")

def test_onoff():

    assert(OnOff('ON'))
    assert(OnOff('On'))
    assert(OnOff('on'))

    assert(not OnOff('OFF'))
    assert(not OnOff('Off'))
    assert(not OnOff('off'))

    assert(OnOff(True).__str__() == "On") 
    assert(OnOff(False).__str__() == "Off")

#############################################################################
#
# test Longitude type
#
#############################################################################

def test_lon360():
    lon = Lon360(350)

    assert(lon._180 == -10)
    assert(lon._360 == 350)

def test_lon180():
    lon = Lon180(-30)

    assert(lon._180 == -30)
    assert(lon._360 == 330)

def test_longitude():
    lon = Longitude(-10, 180)

    assert(lon._180==-10)
    assert(lon._360==350)

#############################################################################
#
# SECTION TO/FROM STRING TESTS
#
#############################################################################


def test_section_to_from_string():

    with open("Example.inp") as f:
        lines = f.readlines()

    sections = {
        'TimeUTC':TimeUTC, 
        'InsertionData':InsertionData, 
        'MeteoData':MeteoData, 
        'Grid':Grid, 
        'Species':Species, 
        'TephraTgsd':TephraTgsd, 
        'RadionucleidesTgsd':RadionucleidesTgsd, 
        'ParticleAggregation':ParticleAggregation, 
        'Source':Source, 
        'Ensemble':Ensemble, 
        'EnsemblePostprocess':EnsemblePostprocess, 
        'ModelPhysics':ModelPhysics,
        'ModelOutput':ModelOutput, 
        'ModelValidation':ModelValidation
        }

    test_strings ={
        'TimeUTC':"".join(lines[16:38]), 
        'InsertionData':"".join(lines[38:52]), 
        'MeteoData':"".join(lines[52:78]), 
        'Grid':"".join(lines[78:110]), 
        'Species':"".join(lines[110:152]), 
        'TephraTgsd':"".join(lines[152:186]), 
        'RadionucleidesTgsd':"".join(lines[186:224]), 
        'ParticleAggregation':"".join(lines[224:246]), 
        'Source':"".join(lines[246:338]), 
        'Ensemble':"".join(lines[338:440]), 
        'EnsemblePostprocess':"".join(lines[440:469]), 
        'ModelPhysics':"".join(lines[469:513]),
        'ModelOutput':"".join(lines[513:564]), 
        'ModelValidation':"".join(lines[564:])
    }

    for name, section in sections.items():
        print(name)

        input_string = test_strings[name]

        output_string = section.from_string(input_string).to_string()

        assert(input_string == output_string)

def test_update_section():
    with open("Example.inp") as f:
        lines = f.readlines()

    sections = {
        'TimeUTC':TimeUTC, 
        'InsertionData':InsertionData, 
        'MeteoData':MeteoData, 
        'Grid':Grid, 
        'Species':Species, 
        'TephraTgsd':TephraTgsd, 
        'RadionucleidesTgsd':RadionucleidesTgsd, 
        'ParticleAggregation':ParticleAggregation, 
        'Source':Source, 
        'Ensemble':Ensemble, 
        'EnsemblePostprocess':EnsemblePostprocess, 
        'ModelPhysics':ModelPhysics,
        'ModelOutput':ModelOutput, 
        'ModelValidation':ModelValidation
        }

    test_strings ={
        'TimeUTC':"".join(lines[16:38]), 
        'InsertionData':"".join(lines[38:52]), 
        'MeteoData':"".join(lines[52:78]), 
        'Grid':"".join(lines[78:110]), 
        'Species':"".join(lines[110:152]), 
        'TephraTgsd':"".join(lines[152:186]), 
        'RadionucleidesTgsd':"".join(lines[186:224]), 
        'ParticleAggregation':"".join(lines[224:246]), 
        'Source':"".join(lines[246:338]), 
        'Ensemble':"".join(lines[338:440]), 
        'EnsemblePostprocess':"".join(lines[440:469]), 
        'ModelPhysics':"".join(lines[469:513]),
        'ModelOutput':"".join(lines[513:564]), 
        'ModelValidation':"".join(lines[564:])
    }

    test_string = test_strings['TimeUTC']

    time_utc = TimeUTC.from_string(test_string)

    time_utc.update({'year':1992})

    assert(time_utc.year==1992)
    
    


#############################################################################
#
# Fall3DInputFile tests
#
#############################################################################

def test_fall3dinputfile_to_from_file():

    in_file = "Example.inp"
    
    out_file = "test.inp"


    Fall3DInputFile.from_file(in_file).to_file(out_file)

    with open(in_file) as f:
        in_string = f.read()

    with open(out_file) as f:
        out_string = f.read()

    assert(in_string==out_string) 

    
def test_fall3dinputfile_update():

    in_file = "Example.inp"

    f3if = Fall3DInputFile.from_file(in_file)

    update = {
            'time_utc': {'year': 2010, 'month': 6, 'day': 27, 'run_start': 16},
            'source': {'height_above_vent': 4970.768719125381, 'mass_flow_rate': 343.217891263361},
            'meteo_data': {'meteo_data_file': 'RUNS/fabb6b6f-ea69-4e78-864d-f236d84e0398/fabb6b6f-ea69-4e78-864d-f236d84e0398_meteo.nc'}
    }

    f3if.update(update)

    assert(f3if.time_utc.year == 2010)
    assert(f3if.time_utc.month == 6)
    assert(f3if.time_utc.day == 27)
    assert(f3if.time_utc.run_start == 16)
    
    assert(f3if.source.height_above_vent == 4970.768719125381)
    assert(f3if.source.mass_flow_rate == 343.217891263361)
    
    assert(f3if.meteo_data.meteo_data_file == 'RUNS/fabb6b6f-ea69-4e78-864d-f236d84e0398/fabb6b6f-ea69-4e78-864d-f236d84e0398_meteo.nc')





