# syntax=docker/dockerfile:1

FROM python:3.14-slim-trixie
RUN apt-get update
RUN apt-get install build-essential git gfortran f2c pkg-config libhdf5-dev wget -y
RUN mkdir /app
RUN mkdir /app/data
COPY ./app.py /app/app.py
COPY ./requirements.txt /app/requirements.txt
COPY ./assets/* /app/assets/
WORKDIR /app

RUN pip3 install -r requirements.txt

# DASH APP
# Copy local gedm_dist_maps.hkl if it exists in app/data/
RUN --mount=type=bind,source=.,target=/tmp/build \
    if [ -f /tmp/build/data/gedm_dist_maps.hkl ]; then cp /tmp/build/data/gedm_dist_maps.hkl /app/data/gedm_dist_maps.hkl; else wget -O /app/data/gedm_dist_maps.hkl "https://zenodo.org/records/18779007/files/gedm_dist_maps.hkl?download=1"; fi || true

EXPOSE 8050/tcp
CMD ["gunicorn", "-w", "8", "-b", "0.0.0.0:8050", "app:server"]
