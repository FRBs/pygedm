# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster
RUN apt-get update 
RUN apt-get install build-essential git gfortran f2c pkg-config libhdf5-dev wget -y
RUN pip3 install numpy ipython scipy astropy 
RUN pip3 install pyne2001
RUN pip3 install git+https://github.com/FRBs/ne2001
RUN pip3 install git+https://github.com/FRBs/pygedm

# DASH APP
RUN mkdir /app
RUN mkdir /app/data
COPY ./app.py /app/app.py
COPY ./assets/* /app/assets
WORKDIR /app
RUN wget 'https://zenodo.org/record/4708481/files/skymap.h5?download=1' -O /app/data/skymap.h5
RUN pip3 install dash dash-bootstrap-components pandas xarray spectra tinycss2 gunicorn Pillow requests
RUN pip3 install h5py --no-binary=h5py
RUN pip3 install git+https://github.com/plotly/dash-labs
EXPOSE 8050/tcp
CMD gunicorn -w 8 -b 0.0.0.0:8050 app:server


