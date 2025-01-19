# particlemap

Data structure for particle-based data which coordinates of particles are defined as sites of Worley noise.

## How to use

Enable the feature `particlemap` in your `Cargo.toml`.

```toml
[dependencies]
worley-particle = { ..., features = ["particlemap"] }
```

## Format Specification (.particlemap)

This library also provides a data format to store particle-based data.

This format uses [Protocol Buffers](https://protobuf.dev/) to define the data structure. See [particlemap.proto](proto/particlemap.proto) for the format specification.

## Features

This library provides `ParticleMap` class to parse and handle the data.

- Read data from file
- Interpolation (Nearest, IDW)
- Rastarisation
- Vectorisation (using [contour-rs](https://crates.io/crates/contour))

## Preview

![map-rasterise](data/output/map-rasterise.png)

![map-isobands](data/output/map-isobands.png)
