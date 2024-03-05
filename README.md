# SS5401 fetch_and_simulate.py
Pacakge for managing batch runs of Fall3D, in order to generate training data for DT-GEO WP5 DTC4.



    SS5401/
    ├── docker-compose.yaml                     - orchestrates containers, ports, volumes
    ├── Dockerfile                              - specifies the container that runs the API
    ├── .dockerignore                           - files to be excluded from the build context
    ├── environment.yaml                        - the required python packages
    ├── .github                                 - github workflows (coverage, CI, etc.)
    │   └── workflows
    │       └── python-package.yml
    ├── .gitignore                              - files to be ignored by git (e.g. .cdsapirc)
    ├── mnt
    │   ├── archive                             - local storage for reanalysis data from CDS
    │   │   └── README.md
    │   ├── aux                                 - extra files needed
    │   │   ├── CARRA_orography_west.nc         - USER SUPPLIED any CARRA dataset with full lat lon grids for reprojection    
    │   │   └── README.md
    │   ├── examples                            - example input files
    │   │   ├── default_so2_reykjanes.inp       - example input file used for SO2 dispersion
    │   │   ├── README.md
    │   │   └── stations.pts                    - example inpui file specifying station locationsa
    │   ├── runs                                - where the training data will be stored
    │   │   └── README.md
    │   └── secrets
    │       ├── .cdsapirc                       - USER SUPPLIED cdsapi login credentials 
    │       └── README.md
    ├── notebook.ipynb                          - notebook illustrating example usage (note requires manual setup of environment with conda)
    ├── pyf3d                                   - package for manipulating Fall3D input files and generating random batch runs
    │   ├── __init__.py
    │   ├── fstrings.py                         - formatting strings for writing a Fall3D input file
    │   ├── regex_strings.py                    - regex strings for reading a Fall3D input file
    │   └── source.py				- class definitions for file sections, the input file, meteo data sources, batch runs
    ├── README.md                               - this file
    ├── fetch_and_simulate.py                   - the script that calls pyf3d to generate input files and runsa Fall3D
    └── test_pyf3d.py                           - unit tests - run with pytest



## To download the repository
Clone the repository to your machine

    git clone https://github.com/profskipulag/SS5401.git

You will be asked for your username and password. For the password github now requires a token:
- on github, click yur user icon in the top right corner
- settings -> developer settings -> personal access tokens -> Tokens (classic) -> Generate new token -> Generate new token (classic) 
- enter you authentifcation code
- under note give it a name, click "repo" to select al check boxes, then click generate token
- copy result enter it as password

## To run the juputer notebook
Create a new conda environment from the environment.yaml file:

    conda env create -f environment.yml

Activate the enmironment

    conda activate ss5401
    
Launch the notebook server

    jupyter notebook
    
Navigate to the SS5401 directory and click the file notebook.ipynb to launch it.

## To run the docker container
The user needs to provide two files, a netcdf file of the whole of the west CARRA domain (e.g. orography)

    mnt/aux/CARRA_orography_west.nc
    
which can be easily downloaded from Copernicus Climate Data Center (CDS) here https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-carra-single-levels?tab=form. The software also needs the users CDS API key, stored as .cdsapirc - see https://cds.climate.copernicus.eu/api-how-to for more details .

    mnt/secrets/.cdsapirc
    
then in the SS5401 directory edit fetch_and_simulate.py to generate random Fall3D runs of your choice, and run:

    sudo docker compose up --build 
    
  
