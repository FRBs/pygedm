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
# --preload: load the app once BEFORE forking workers so skymap data
#   lives in shared copy-on-write memory instead of being duplicated.
# -w 2: two workers is plenty for a Dash app; 8 workers each loading
#   ~60 MB of data was consuming ~500 MB+ of RAM unnecessarily.
# --timeout 120: data loading + first callback can take a while on cold start.
# --access-logfile -: emit HTTP request logs to stdout for observability.
CMD ["gunicorn", "--preload", "-w", "2", "--timeout", "120", "--access-logfile", "-", "-b", "0.0.0.0:8050", "app:server"]
