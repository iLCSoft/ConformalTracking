# v01-05

* 2017-10-20 Daniel Hynds ([PR#25](https://github.com/iLCSoft/ConformalTracking/pull/25))
  - Require 5 hits/track for displaced track finding, to reduce fake rate

* 2017-10-13 Daniel Hynds ([PR#24](https://github.com/iLCSoft/ConformalTracking/pull/24))
  - Functions updated to allow tracker to vertex tracking, to allow displaced track reconstruction

* 2017-12-04 Andre Sailer ([PR#28](https://github.com/iLCSoft/ConformalTracking/pull/28))
  - Performance optimisations: avoiding temporaries, divisions, sqrts
  - protect against too many tracks

* 2017-12-12 Andre Sailer ([PR#29](https://github.com/iLCSoft/ConformalTracking/pull/29))
  - Performance optimisations: 
     - avoiding copies of `shared_ptrs`, avoiding temporary objects (use `emplace_back`, direct construction etc.)
     - move filtering of kdtree search results to before sorting
     - use `vdt::fast_atan` in `Cell:getAngle[RZ]` (enabled if ROOT provides `vdt`)
     - add option to turn off sorting of `kdtree` search results, this changes the outcome of the reconstruction, thus sorting is enabled by default

* 2017-11-27 Emilia Leogrande ([PR#26](https://github.com/iLCSoft/ConformalTracking/pull/26))
  - Parameter struct to set pattern recognition parameters
    - `_maxCellAngle`, `_maxCellAngleRZ`, `_chi2cut`, `_minClustersOnTrack`, `_maxDistance`, `_highPTfit`, `_onlyZSchi2cut`

* 2018-03-13 Marko Petric ([PR#30](https://github.com/iLCSoft/ConformalTracking/pull/30))
  -  Fix for iLCSoft/LCIO#35

# v01-04

* 2017-09-20 Daniel Hynds ([PR#20](https://github.com/iLCSoft/ConformalTracking/pull/20))
  - Distance Of Closest Approach to the origin added, used as a cut on new seed cells in order to reduce combinatorics

* 2017-07-11 Daniel Hynds ([PR#15](https://github.com/iLCSoft/ConformalTracking/pull/15))
  - Rewrote track fit in SZ, proper phi treatment 
  - Updated track strategy including extension of hits throughout the detector
  - Utility tool added for debugging

* 2017-09-21 Andre Sailer ([PR#21](https://github.com/iLCSoft/ConformalTracking/pull/21))
  - Fix gcc compiler warnings

* 2017-08-18 Andre Sailer ([PR#19](https://github.com/iLCSoft/ConformalTracking/pull/19))
  - errors for endcaps updated and forward flag added
  - high pt track extension now uses CA (over 10 GeV/c pt)
  - strategy slightly changed for displaced tracks to stop overloading in case of high occupancy events, fixes #16

* 2017-10-02 Andre Sailer ([PR#22](https://github.com/iLCSoft/ConformalTracking/pull/22))
  - Use smart pointers to wrap objects and simplify memory management.

* 2017-07-27 Andre Sailer ([PR#17](https://github.com/iLCSoft/ConformalTracking/pull/17))
  - Added temporary workaround when a way too large number of tracks would be created. Add skipEvent flag to allow for proper cleanup before exception is thrown.

* 2017-10-06 Andre Sailer ([PR#23](https://github.com/iLCSoft/ConformalTracking/pull/23))
  - Drop unused and no longer existing header includes AidaSoft/DD4hep#241

# v01-03

* 2017-06-13 Daniel Hynds ([PR#11](https://github.com/iLCSoft/ConformalTracking/pull/11))
  - Added clang-format to cmake lists and added extra target. "make format" will now format the code, and will be prepared automatically if clang-format exists.
  - Updated tracking strategy and added some improvements in terms of refitting tracks each time a hit is added, and included a maximum chi2 criteria for adding hits

* 2017-06-13 Andre Sailer ([PR#10](https://github.com/iLCSoft/ConformalTracking/pull/10))
  - ConformalTracking: Fix some memory leaks

* 2017-06-28 Andre Sailer ([PR#14](https://github.com/iLCSoft/ConformalTracking/pull/14))
  -  Add format checker as a CI check

* 2017-06-28 Andre Sailer ([PR#13](https://github.com/iLCSoft/ConformalTracking/pull/13))
  - ConformalTracking: cleanup final conformal tracks

* 2017-05-31 Andre Sailer ([PR#7](https://github.com/iLCSoft/ConformalTracking/pull/7))
  - Make all output DEBUG for now, because there is too much output for the grid to handle

* 2017-06-20 Andre Sailer ([PR#12](https://github.com/iLCSoft/ConformalTracking/pull/12))
  - Adapt to DD4hep namespace change, get BField from utility function

* 2017-06-03 Daniel Hynds ([PR#8](https://github.com/iLCSoft/ConformalTracking/pull/8))
  - Put track fit back in, accidentally missing during debugging

* 2017-05-24 Daniel Hynds ([PR#6](https://github.com/iLCSoft/ConformalTracking/pull/6))
  - Chi2 calculations completely reworked
  - Significant code cleanup
  - Code modularised to allow different tracking strategies 
  - Default state currently still vertex-only track construction

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

