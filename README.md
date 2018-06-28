# ConformalTracking
[![Build Status](https://travis-ci.org/iLCSoft/ConformalTracking.svg?branch=master)](https://travis-ci.org/iLCSoft/ConformalTracking)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12348/badge.svg)](https://scan.coverity.com/projects/ilcsoft-conformaltracking)

Package for running pattern recognition based on conformal mapping and cellular automaton. This is not tied to a given geometry, but has been developed for the CLIC detector model 2015.

ConformalTracking is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)


## License and Copyright
Copyright (C), ConformalTracking Authors

ConformalTracking is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.



## ConformalTrackingV2

The **ConformalTrackingV2** processor is fully configurable via the steering
file. The parameter *Steps* takes a "Config" like String.  For example:

    [VXDBarrel]
    @Collections = VXDTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, BuildNewTracks
    [VXDEncap]
    @Collections = VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, ExtendTracks
    [TCVC]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits
    @Parameters = MaxCellAngle = 0.035; MaxCellAngleRZ = 0.035;  Chi2Cut = 300;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker
    @Functions = CombineCollections, BuildNewTracks
    [LowerCellAngle1]
    @Collections =
    @Parameters = MaxCellAngle = 0.175; MaxCellAngleRZ = 0.175;  Chi2Cut = 300;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = BuildNewTracks
    [LowerCellAngle2]
    @Collections =
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 6000;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = BuildNewTracks
    [LowerHitNumbers]
    @Collections =
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 6000;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = BuildNewTracks, SortTracks
    [Tracker]
    @Collections = ITrackerHits,OTrackerHits,ITrackerEndcapHits,OTrackerEndcapHits
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 6000;  MinClustersOnTrack = 4;  MaxDistance = 0.015;
    @Flags = HighPTFit, VertexToTracker, RadialSearch
    @Functions = CombineCollections, ExtendTracks
    [Displaced]
    @Collections = VXDTrackerHits, VXDEndcapTrackerHits, ITrackerHits, OTrackerHits, ITrackerEndcapHits, OTrackerEndcapHits
    @Parameters = MaxCellAngle = 0.35; MaxCellAngleRZ = 0.35;  Chi2Cut = 3000;  MinClustersOnTrack = 5;  MaxDistance = 0.015;
    @Flags = OnlyZSchi2cut, RadialSearch
    @Functions = CombineCollections, BuildNewTracks

Each steps starts with a *[Name]* in brackets. Following are the *@Collections*,
*@Parameters*, *@Flags*, and *@Functions* entries. All *@* have to be
present. The *@Collections* and *@Flags* can be empty. It does not make sense to
leave the *@Functions* empty. The *@Parameters* needs to be present and all
parameters need to be set. Parameters and values need to be separated by an
equal sign, separation between parameters can be a comma, semicolon, or colon.

The *@Parameters* are

  * MaxCellAngle
  * MaxCellAngleRZ
  * Chi2Cut
  * MinClustersOnTrack
  * MaxDistance

Possible *@Flags* are

  * OnlyZSchi2cut
  * VertexToTracker
  * RadialSearch
  * HighPTFit

Possible *@Functions* are

  * CombineCollections
  * BuildNewTracks
  * ExtendTracks
  * SortTracks

