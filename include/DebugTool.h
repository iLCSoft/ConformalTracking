#ifndef DEBUGTOOL_H
#define DEBUGTOOL_H 1

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <math.h>
#include "KDCluster.h"
#include "KDTrack.h"

// ------------------------------------------------------------------------------------
// The ConformalDebugger is a helper class to find associations between MCTruth and the
// objects used in the ConformalTracking processor. It will create links from
// MCParticles to KDClusters and then give purity values etc. for given tracks.
// ------------------------------------------------------------------------------------

using SLCRelationNavigator = std::shared_ptr<LCRelationNavigator>;

class ConformalDebugger {
public:
  //--- Constructor and destructor
  ConformalDebugger(){};
  ~ConformalDebugger(){};

  ConformalDebugger(const ConformalDebugger&) = delete;
  ConformalDebugger& operator=(const ConformalDebugger&) = delete;

  //--- Member variables
  LCCollection*                     m_particleCollection = nullptr;  // MC particles
  std::vector<SLCRelationNavigator> m_relations{};                   // relations from tracker hits to MCP
  std::vector<LCCollection*>        m_trackerHitCollections{};       // tracker hits

  std::map<int, SharedKDClusters>         m_kdClusters{};    // list of kd clusters
  std::map<SKDCluster, TrackerHitPlane*>  m_kdClusterMap{};  // their link to "real" hits
  std::map<SKDCluster, MCParticle*>       m_kdParticles{};   // link from conformal hit to MC particle
  std::map<MCParticle*, SharedKDClusters> m_particleHits{};  // list of conformal hits on each MC particle

  //--- Set member variables
  void setMCParticles(LCCollection* particleCollection) { m_particleCollection = particleCollection; }
  void setRelations(std::vector<SLCRelationNavigator> relations) { m_relations = relations; }
  void setTrackerHits(std::vector<LCCollection*> trackerHitCollections) { m_trackerHitCollections = trackerHitCollections; }
  void setKDClusters(int collection, SharedKDClusters clusters) { m_kdClusters[collection] = clusters; }

  //--- Clear all member variables
  void clear() {
    m_particleCollection = NULL;
    m_relations.clear();
    m_trackerHitCollections.clear();
    m_kdClusters.clear();
    m_kdClusterMap.clear();
    m_kdParticles.clear();
    m_particleHits.clear();
  }

  //--- Register a tracker hit and its kd hit
  void registerHit(int collection, SKDCluster kdHit, TrackerHitPlane* hit) {
    // Store the link between the kdhit and tracker hit
    m_kdClusterMap[kdHit] = hit;

    // Store the link to the associated MCParticle
    const LCObjectVec& simHitVector = m_relations[collection]->getRelatedToObjects(hit);

    // Take the first hit only (TODO: this should be changed?)
    SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));

    MCParticle* particle = simHit->getMCParticle();
    m_kdParticles[kdHit] = particle;

    // Don't save hits from secondaries to the expected particle hits
    if (simHit->isProducedBySecondary())
      return;

    m_particleHits[particle].push_back(kdHit);
  }

  //--- Get the particle that corresponds to the largest number of hits on this track
  MCParticle* getAssociatedParticle(UKDTrack& track) {
    // Temporarily store all mcparticles associated to this track
    std::vector<MCParticle*> particles;
    std::map<MCParticle*, double> particleHits;
    double nHits = 0.;

    // Get the clusters from this track
    SharedKDClusters clusters = track->m_clusters;

    // Loop over all hits and see which particle they are associated to
    for (size_t itCluster = 0; itCluster < clusters.size(); itCluster++) {
      // Get the hit
      SKDCluster cluster = clusters[itCluster];
      nHits++;

      // If we already have hits on this particle, then just increment the counter
      if (particleHits.count(m_kdParticles[cluster]))
        particleHits[m_kdParticles[cluster]]++;
      else {
        // Otherwise register the new particle and set the number of hits to 1
        particles.push_back(m_kdParticles[cluster]);
        particleHits[m_kdParticles[cluster]] = 1;
      }
    }

    // Now look how many hits are on each particle and calculate the purity
    double      bestHits           = 0.;
    MCParticle* associatedParticle = NULL;
    for (size_t iPart = 0; iPart < particles.size(); iPart++) {
      if (particleHits[particles[iPart]] > bestHits) {
        bestHits           = particleHits[particles[iPart]];
        associatedParticle = particles[iPart];
      }
    }

    // Return the associated MC particle
    return associatedParticle;
  }

  //--- Check if a hit is associated to the given MC particle
  bool isAssociated(SKDCluster const& kdHit, MCParticle* particle) {
    if (std::find(m_particleHits[particle].begin(), m_particleHits[particle].end(), kdHit) != m_particleHits[particle].end())
      return true;
    else
      return false;
  }

  SharedKDClusters getAssociatedHits(MCParticle* particle) { return m_particleHits[particle]; }
};

#endif
