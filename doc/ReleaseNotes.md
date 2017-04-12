# v01-02

* 2017-03-29 Daniel Hynds ([PR#3](https://github.com/iLCSoft/ConformalTracking/pull/3))
  - Collections are added sequentially to the list of hits being used to make tracks, including unused hits from the previous collections. Previously hits were considered individually before this step, leading to missed hits in the interaction region.

* 2017-03-28 Andre Sailer ([PR#2](https://github.com/iLCSoft/ConformalTracking/pull/2))
  - Removed "CellIDDecoderString" parameter, not necessary by using encoding_string()

# v01-00
- **General**:The package contains a first implementation of a fully hermitic pattern recognition using all tracking detectors (developed for the CLIC detector model 2015 containing only silicon tracking). Hits from all detectors are considered together, without  knowledge of geometry, using a two-dimensional conformal mapping to reduce the track search to a straight line search. Cellular automaton is used in this conformal space to produce tracks. 

- **ConformalTracking**: The main pattern recognition processor. It will try to produce tracks by sequentially adding the input hit collections, which should be provided from inner to outer radius.

- **KDTree and kdtree2**: Implementation of a fast nearest neighbour search algorithm (from BOOST) and a wrapper which interfaces the KDHit class.

- **KDCluster**: A simple class representing a hit in conformal space. It stores the uv co-ordinates in cartesian and polar notation, to allow for angular searches (during seeding) and spatial nearest neighbour searches.

- **Cell**: The cell class used for cellular automaton. Holds some simple information concerning the cell weight, other cells to which it is connected, and the hits from which it is made.

