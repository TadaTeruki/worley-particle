# worleyMap

Data format parser and handler for particle-based data which coordinates of particles are defined as sites of Worley noise.

## Format Specification

`x` and `y` is the sites of each particle.

```plaintext:example.worleymap
seed,min_randomness,max_randomness,scale
[{x},{y}]{other data}
...
```

### Example

```
seed:150,min_randomness:0.5,max_randomness:0.8,scale:0.05
[23.8,45.2]355.7,mountain,rocky
[18.5,32.1]2.3,plain,grass
[3.2,7.8]-0.5,sea,sand
```

## Features

This library provides `WorleyMap` class to parse and handle the data.
- Read data from file
- Interpolation
- Rastarisation
- Create Contour (using [contour-rs](https://github.com/mthh/contour-rs))