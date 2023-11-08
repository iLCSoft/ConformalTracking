# v01-12

* 2023-11-08 Leonhard Reichenbach ([PR#62](https://github.com/ilcsoft/ConformalTracking/pull/62))
  - Emit an error message and skip further steps when too many tracks are created instead of throwing a SkipEventException

# v01-11-01

* 2023-05-25 Andre Sailer ([PR#60](https://github.com/iLCSoft/ConformalTracking/pull/60))
  - ParameterParser: fix linking issue with recent versions of boost
  - ParameterParser: move to non-deprecated header file location

# v01-11

* 2021-09-27 Andre Sailer ([PR#59](https://github.com/iLCSoft/ConformalTracking/pull/59))
  - CI: Move to github actions

# v01-10

* 2019-10-08 Remi Ete ([PR#55](https://github.com/iLCSoft/ConformalTracking/pull/55))
  - ConformalTracking processor: Use streamlog macros instead of raw calls

* 2019-09-09 Andre Sailer ([PR#54](https://github.com/iLCSoft/ConformalTracking/pull/54))
  - Add Parameter TooManyTracks to configure what is considered too many, defaults to 5e5
  - Fix the abort logic so that after ten tries the event is actually skipped, or if retryTooManyTracks is False immediately, instead of never or until the number of created tracks is small enough

* 2019-08-22 Erica Brondolin ([PR#53](https://github.com/iLCSoft/ConformalTracking/pull/53))
  - Include cut in slope in z for CA building part to reduce number of combinatorics
    - In ConformalTrackingV2: configurable via per step parameter slopeZRange
    - Reduces CPU time
  - Change sum of chi2 from square-root of the squared sum to linear sum
  - No big changes expected in eff/fakes

# v01-08

* 2019-02-07 Emilia Leogrande ([PR#49](https://github.com/iLCSoft/ConformalTracking/pull/49))
  - Clone checker: if equal length, keep the one with the better chi2
  - Debug9 prints clone checker info
  - Fix address of the track in the debug printouts among functions

* 2018-10-30 Erica Brondolin ([PR#48](https://github.com/iLCSoft/ConformalTracking/pull/48))
  - Solving bug introduced with the PR #47 : the `fitError` variable was overwritten even when the inverted fit was not better than the normal one.
  - The maximum number of track hit to try the inverted fit is now a parameter.

* 2018-10-26 Erica Brondolin ([PR#47](https://github.com/iLCSoft/ConformalTracking/pull/47))
  - If the track is too short (usually less than 7 hits correspond to a vertex track) the fit in the inverse direction is tried. If the track is longer, then it is kept as the final one.
  - This is important especially for tracks reconstructed in the very forward region where the extrapolation VTX-TRK is longer.

* 2018-10-11 Erica Brondolin ([PR#46](https://github.com/iLCSoft/ConformalTracking/pull/46))
  - Introduce `extendTracksPerLayer` function to find all the hits layer by layer
  - Filter the neighbours that are in the wrong direction of the track and that are in layers too far
  - The best candidates are inserted in the track if their chi2 is within a certain range, otherwise only the hit with the best chi2 is kept
  - Results in the bbar show improved number of hits for tracks with pT > 1 GeV and basically unchanged for tracks with pT < 1 GeV

* 2018-09-27 Emilia Leogrande ([PR#45](https://github.com/iLCSoft/ConformalTracking/pull/45))
  - ConfomalTrackingV2
      - Kalman fit direction can now be set per reconstruction step, to allow for the foreseen possibility to fit prompt tracks forward and non-prompt backward
      - added flag "KalmanFitBackward" and "KalmanFitForward" that can be parsed from steering file; if nothing is set, the default is forward
   - Kalman fit direction added as member variable of the KDTrack class, with getter/setter

* 2018-09-26 Emilia Leogrande ([PR#44](https://github.com/iLCSoft/ConformalTracking/pull/44))
  - in ConformalTracking::runStep, marking hits as used is now repositioned at the end of every step
  - implemented debug streamlog printouts to follow pattern recognition
  -- streamlog_out( DEBUG9 ) for 'main' functions: processEvent, buildNewTracks [to do: extendTracks, but function is intended to be upgraded]
  -- streamlog_out( DEBUG8 ) for 'secondary' functions: extendSeedCells, createTracksNew, getFittedTracks, getLowestChi2

* 2018-08-09 Emilia Leogrande ([PR#42](https://github.com/iLCSoft/ConformalTracking/pull/42))
  - The pt threshold to run extendHighPt is now a parameter, to be set in the steering file
    - It has effect only on the steps when extendTracks is run (buildNewTracks does not use extendHighPt)
    - Default values at the moment: 10 GeV for extending in the vertex endcap, 1 GeV for extending in the trackers

* 2018-07-24 Emilia Leogrande ([PR#41](https://github.com/iLCSoft/ConformalTracking/pull/41))
  - ConformalTrackingV2: removed the call to init()
    - init() is already done in ConformalTracking.cc
    - since in init() the step parameters are parsed, they were parsed twice, resulting in twice the number of steps in the pattern recognition chain

* 2018-06-28 Andre Sailer ([PR#40](https://github.com/iLCSoft/ConformalTracking/pull/40))
  - New processor ConformalTrackingV2, which allows complete freedom to configure the steps of the pattern recognition
    - See the README.md File
    - The ConformalTracking processor is kept for backward compatibility but should be considered deprecated
    - Refactored ConformalTracking for inheritance, loop over vector of parameters
    - Cleaned up header includes and CMakeLists

* 2018-06-25 Emilia Leogrande ([PR#39](https://github.com/iLCSoft/ConformalTracking/pull/39))
  - Conformal Tracking restructured to be run in steps with function runStep
  - runStep allows to decide which functions to run per step for the reconstruction, i.e. combineCollections, buildNewTracks, extendTracks
  - prior to runStep, the set of parameters for the reconstruction step is initialized

* 2018-05-18 Emilia Leogrande ([PR#37](https://github.com/iLCSoft/ConformalTracking/pull/37))
  - ConformalTracking: fixed two issues, one in the pattern recognition, one in the choice among clones
  - when building tracks in the combined vertex barrel + endcap, the step with looser cut must be done all at one (open the angles and increase the chi2 cut at the same time, not in two subsequent steps), otherwise tracks are made with less hits than desired
  - once increased, keep the chi2 cut the same for subsequent steps
  - in presence of clones, it is not good to prefer shorter tracks with better chi2 (the comparison is also progressive, so one can easily go from a 8 hits track to a 4 hits track). To avoid in order not to have double tracks per particle

* 2018-05-16 Emilia Leogrande ([PR#36](https://github.com/iLCSoft/ConformalTracking/pull/36))
  - new strategy: hits on the same subdetector layer are accepted (otherwise the leftover of the pairs can make a second track)
  - if hits on the same *sensor* of the same subdetector layer, only one (smaller radius) is taken
  - flag to enable tight cuts in the combined vertex endcap + barrel added as processor parameter (Temporary, will be removed in the future)

* 2018-05-09 Andre Sailer ([PR#35](https://github.com/iLCSoft/ConformalTracking/pull/35))
  - ConformalTracking.cc: [extendTracks]fix for hits from the same subdetector layers not to be accepted in the same track, when tracks are extended to vertex endcap and to trackers


# v01-07

* 2018-03-29 Andre Sailer ([PR#33](https://github.com/ilcsoft/ConformalTracking/pull/33))
  - combineCollections: Properly protect against empty collections, fixes #31

* 2018-04-17 Andre Sailer ([PR#34](https://github.com/ilcsoft/ConformalTracking/pull/34))
  - ConformalTracking: change setting of input collections, use collection names instead of indices, check for mistakes in processor::init
  
     * Drop Processorparameters: AllCollectionIndices, TrackerHitCollectionIndices
     * Add ProcessorParameters: MainTrackerHitCollectionNames, VertexBarrelHitCollectionNames, VertexEndcapHitCollectionNames
     * Automatically obtain collection indices, throw exception if there is a collection name mismatch

# v01-06

* 2018-03-28 Frank Gaede ([PR#31](https://github.com/iLCSoft/ConformalTracking/pull/31))
  - make compatible for other detectors, e.g. ILD
       - introduce new parameters with indices of tracker hit collections:
            - AllHitCollectionIndices
            - TrackerHitCollectionIndices
        - default values as used for CLIC
        - protect against empty collections

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

