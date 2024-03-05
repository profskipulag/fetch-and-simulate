FROM continuumio/miniconda3
RUN apt-get --allow-releaseinfo-change update
RUN apt-get update && apt-get install -y libgl1-mesa-glx ffmpeg libsm6 libxext6 libarchive-dev gawk build-essential gfortran libnetcdf-dev libnetcdff-dev mpich makedepf90 automake netcdf-bin vim
RUN git clone https://gitlab.com/fall3d-suite/fall3d.git
RUN /bin/bash -c "cd fall3d; git pull"
RUN /bin/bash -c "cd fall3d; ./configure; make; make install"
RUN /bin/bash -c "cd fall3d/Example; wget https://gitlab.com/fall3d-distribution/testsuite/-/raw/master/example-8.0/InputFiles/example-8.0.wrf.nc"
RUN conda install -c conda-forge mamba
COPY environment.yaml /app/environment.yaml
RUN /bin/bash -c "mamba env create -f /app/environment.yaml"
COPY pyf3d/ /app/pyf3d/
COPY fetch_and_simulate.py /app/fetch_and_simulate.py
ENTRYPOINT ["conda", "run", "--no-capture-output", "--cwd","/app","-n", "ss5401","python","fetch_and_simulate.py"]

#ENTRYPOINT ["ls","app"]
#CMD cd fall3d/bin; ./Fall3d.r8.x All /RUNS/Example.inp  > /RUNS/log.txt
# TO LOGIN! sudo docker run --rm -it --entrypoint /bin/bash fall3d-fall3d
