FROM --platform=linux/arm/v7 debian:bookworm

RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    git \
    ca-certificates \
    wget \
    libasound2-dev \
    puredata-dev \
    && rm -rf /var/lib/apt/lists/*
