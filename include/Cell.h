#ifndef CELL_H
#define CELL_H 1

#include <math.h>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <vector>
#include "KDCluster.h"

#ifdef CF_USE_VDT
#include <vdt/atan.h>
#endif

// ------------------------------------------------------------------------------------
// The Cell class is a simple object which connects two points in 2D space. It is
// used in Cellular Automaton tracking to create tracks, by connecting all plausable
// hit points with cells, and linking these cells together in a chain. Each cell
// needs to know what cells it is connected to, and hold a weight determined by
// its position in the chain.
// ------------------------------------------------------------------------------------

class Cell {
public:
  typedef std::vector<std::weak_ptr<Cell>> WeakCells;
  typedef std::shared_ptr<Cell>            SCell;
  typedef std::weak_ptr<Cell>              WCell;

public:
  // Constructors, main initialisation is with two kd hits
  Cell() { m_weight = 0; }
  Cell(const Cell&) = delete;
  Cell(Cell&&)      = default;
  Cell& operator=(const Cell&) = delete;

  Cell(SKDCluster const& hit1, SKDCluster const& hit2)
      : m_weight(0),
        m_gradient((hit2->getV() - hit1->getV()) / (hit2->getU() - hit1->getU())),
        m_gradientRZ((hit2->getRadius() - hit1->getRadius()) / (hit2->getZ() - hit1->getZ())),
        m_start(hit1),
        m_end(hit2) {}

  // Destructor
  ~Cell() {}

  // Weight of the cell (first cell in a chain has weight 0, and each subsequent link has weight +1)
  int  getWeight() const { return m_weight; }
  void setWeight(int weight) { m_weight = weight; }

  // Gradient of the cell connecting two hits
  double getGradient() const { return m_gradient; }
  void setGradient(double gradient) { m_gradient = gradient; }
  double                  getGradientRZ() const { return m_gradientRZ; }

  // Angle between two cells. This is assumed to be less than 90 degrees
  double getAngle(SCell const& cell2) const { return getAngle(*(cell2.get())); }
  double getAngleRZ(SCell const& cell2) const { return getAngleRZ(*(cell2.get())); }

  double getAngle(Cell const& cell2) const {
#ifdef CF_USE_VDT
    return fabs(vdt::fast_atan((cell2.m_gradient - m_gradient) / (1 + m_gradient * cell2.m_gradient)));
#else
    return fabs(std::atan((cell2.m_gradient - m_gradient) / (1 + m_gradient * cell2.m_gradient)));
#endif
  }

  double getAngleRZ(Cell const& cell2) const {
#ifdef CF_USE_VDT
    return fabs(vdt::fast_atan((cell2.m_gradientRZ - m_gradientRZ) / (1 + m_gradientRZ * cell2.m_gradientRZ)));
#else
    return fabs(std::atan((cell2.m_gradientRZ - m_gradientRZ) / (1 + m_gradientRZ * cell2.m_gradientRZ)));
#endif
  }

  // Start and end points of the cell
  SKDCluster const& getStart() const { return m_start; }
  SKDCluster const& getEnd() const { return m_end; }

  // Increment the cell weight (usually if the chain length is extended upstream of this cell)
  void update(SCell const& cell2) {
    if ((cell2->getWeight() + 1) > m_weight)
      m_weight = cell2->getWeight() + 1;
  }

  // The cell has a memory of all cells that connect to it, and all cells that it connects to. If several cells point to this
  // cell, then the weight taken from the highest weighted of those (longest chain)
  void setFrom(SCell const& cell2) {
    m_from.push_back(WCell(cell2));
    m_weights.push_back(cell2->getWeight());
    if ((cell2->getWeight() + 1) > m_weight)
      m_weight = cell2->getWeight() + 1;
  }
  void setTo(SCell const& cell2) { m_to.push_back(WCell(cell2)); }
  WeakCells&              getFrom() { return m_from; }
  WeakCells&              getTo() { return m_to; }

  double doca() const {
    double intercept = m_start->getV() - m_start->getU() * m_gradient;
    double doca      = fabs(intercept) / sqrt(m_gradient * m_gradient + 1.);
    return doca;
  }

private:
  // Each cell contains a weight, a gradient, two hits which it connects
  // and a list of cells that it connects to or from
  int              m_weight     = 0;
  double           m_gradient   = 0.0;
  double           m_gradientRZ = 0.0;
  SKDCluster       m_start      = nullptr;
  SKDCluster       m_end        = nullptr;
  WeakCells        m_from{};
  std::vector<int> m_weights{};
  WeakCells        m_to{};
};

using SCell                = Cell::SCell;
using cellularTrack        = std::vector<SCell>;
using SharedCells          = std::vector<std::shared_ptr<Cell>>;
using UcellularTrack       = std::unique_ptr<cellularTrack>;
using UniqueCellularTracks = std::vector<std::unique_ptr<cellularTrack>>;
#endif
