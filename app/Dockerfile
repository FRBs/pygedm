# syntax=docker/dockerfile:1

FROM python:3.14-slim-trixie

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential git gfortran f2c pkg-config libhdf5-dev wget && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /app/data /app/assets
COPY ./app.py /app/app.py
COPY ./requirements.txt /app/requirements.txt
COPY ./download_data.py /app/download_data.py
COPY ./assets/ /app/assets/
WORKDIR /app

RUN pip3 install --no-cache-dir -r requirements.txt


# Download gedm_dist_maps.hkl from Zenodo (or use a local copy if in ./assets)
RUN --mount=type=bind,target=/tmp/ctx,source=. \
    python3 /app/download_data.py

EXPOSE 8050/tcp
CMD ["gunicorn", "--preload", "-w", "4", "--timeout", "120", "--access-logfile", "-", "-b", "0.0.0.0:8050", "app:server"]
