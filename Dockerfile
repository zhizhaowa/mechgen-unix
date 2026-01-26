FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    bison \
    flex \
    autoconf \
    wget \
    tar \
    vim \
    telnet \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

#COPY MechGen_online_v5_unix.db /app/mechgen/in.db
COPY src /app
COPY unix /app
#COPY MOO-1.8.1 /app/mechgen/MOO

# Work directory for runtime
WORKDIR /app
