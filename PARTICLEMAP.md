# particlemap

Data structure for particle-based data which particles are defined as sites of Worley noise.

## How to use

Enable the feature `particlemap` in your `Cargo.toml`.

```toml
[dependencies]
worley-particle = { ..., features = ["particlemap"] }
```

## Features

This library provides `ParticleMap` class to parse and handle the data.

- Serialization and Deserialization
- Interpolation (Nearest, IDW)
- Rastarisation
- Vectorisation (using [contour-rs](https://crates.io/crates/contour))

## Format Specification

This library uses a custom binary format based on [Protocol Buffers](https://protobuf.dev/) to serialize and deserialize the data.

See [particlemap.proto](proto/particlemap.proto) for the format specification.


## Preview

![map-rasterise](data/output/map-rasterise.png)

![map-isobands](data/output/map-isobands.png)
