# syntax=docker/dockerfile:1
FROM python:3.8-slim-buster
RUN apt-get update 
RUN apt-get install build-essential git f2c pkg-config  -y
COPY . /app
WORKDIR /app
RUN pip3 install -r requirements.txt
RUN pip3 install -r requirements_test.txt
RUN pip3 install .




