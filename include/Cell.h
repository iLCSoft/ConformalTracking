#ifndef CELL_H
#define CELL_H 1

#include <iosfwd>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "KDCluster.h"
#include <math.h>

// ------------------------------------------------------------------------------------
// The Cell class is a simple object which connects two points in 2D space. It is
// used in Cellular Automaton tracking to create tracks, by connecting all plausable
// hit points with cells, and linking these cells together in a chain. Each cell
// needs to know what cells it is connected to, and hold a weight determined by
// its position in the chain.
// ------------------------------------------------------------------------------------

class Cell
{
	public:
  
  	// Constructors, main initialisation is with two kd hits
		Cell(){m_weight=0;}
	  Cell(KDCluster* hit1, KDCluster* hit2){
			m_start = hit1;
			m_end = hit2;
			m_gradient = (hit2->getV() - hit1->getV())/(hit2->getU() - hit1->getU());
			m_weight=0;
		}
  
  	// Destructor
		virtual ~Cell(){
      m_from.clear();
      m_to.clear();
    }
	
  	// Weight of the cell (first cell in a chain has weight 0, and each subsequent link has weight +1)
		int getWeight(){return m_weight;}
		void setWeight(int weight){m_weight = weight;}
		
  	// Gradient of the cell connecting two hits
  	double getGradient(){return m_gradient;}
		void setGradient(double gradient){m_gradient=gradient;}

  	// Angle between two cells. This is assumed to be less than 90 degrees
  	double getAngle(Cell* cell2){
			return fabs(std::atan( (cell2->getGradient()-m_gradient)/(1+m_gradient*cell2->getGradient()) ));
		}
  
  	// Start and end points of the cell
		KDCluster* getStart(){return m_start;}
		KDCluster* getEnd(){return m_end;}
  
  	// Increment the cell weight (usually if the chain length is extended upstream of this cell)
  	void update(Cell* cell2){
    	if((cell2->getWeight()+1)>m_weight) m_weight = cell2->getWeight()+1;
  	}
  
  	// The cell has a memory of all cells that connect to it, and all cells that it connects to. If several cells point to this
		// cell, then the weight taken from the highest weighted of those (longest chain)
  	void setFrom(Cell* cell2){m_from.push_back(cell2); m_weights.push_back(cell2->getWeight()); if((cell2->getWeight()+1)>m_weight) m_weight = cell2->getWeight()+1;}
    void setTo(Cell* cell2){m_to.push_back(cell2);}
		std::vector<Cell*>* getFrom(){return &m_from;}
		std::vector<Cell*>* getTo(){return &m_to;}
	
	private:

  	// Each cell contains a weight, a gradient, two hits which it connects
  	// and a list of cells that it connects to or from
	  int m_weight;
		double m_gradient;
	  KDCluster* m_start;
  	KDCluster* m_end;
		std::vector<Cell*> m_from;
		std::vector<int> m_weights;
		std::vector<Cell*> m_to;
	
};

// Vector of cells
typedef std::vector<Cell*> cellularTrack;

#endif
