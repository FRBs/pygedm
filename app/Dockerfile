# syntax=docker/dockerfile:1

FROM python:3.14-slim-trixie

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential git gfortran f2c pkg-config libhdf5-dev wget && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /app/data /app/assets
COPY ./app.py /app/app.py
COPY ./requirements.txt /app/requirements.txt
COPY ./assets/ /app/assets/
WORKDIR /app

RUN pip3 install --no-cache-dir -r requirements.txt

# Try to copy a local gedm_dist_maps.hkl if the builder has one or download from Zenodo
ENV DATA_URL="https://zenodo.org/records/18779007/files/gedm_dist_maps.hkl?download=1"

# Copy data dir contents if present
COPY data/gedm_dist_maps.hk[l] /app/data/

RUN if [ ! -f /app/data/gedm_dist_maps.hkl ]; then \
        echo "gedm_dist_maps.hkl not found locally — downloading..." && \
        wget --progress=dot:giga -O /app/data/gedm_dist_maps.hkl "$DATA_URL"; \
    else \
        echo "gedm_dist_maps.hkl found locally — skipping download."; \
    fi && \
    # Sanity check: file must exist and not be empty
    test -s /app/data/gedm_dist_maps.hkl || \
        { echo "FATAL: gedm_dist_maps.hkl is missing or empty"; exit 1; }

EXPOSE 8050/tcp
CMD ["gunicorn", "--preload", "-w", "4", "--timeout", "120", "--access-logfile", "-", "-b", "0.0.0.0:8050", "app:server"]
