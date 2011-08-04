/*



*/


#include "Definitions.h"
#include "ComponentState.h"
#include "PropertiesModel.h"


class TwoPhaseDynamic {
	
	Public:
	int TwoPhaseDynamic(ComponentState* initialconditions, PropertiesModel* props);  //Constructor - Fills in blanks in initialconditions
	int Timestep(ComponentState* state, double interval, double timestep);  //Run component with with given flow variables for timestep

	Private:
	// internal parameters and conditions	
	double V, Fill, heatrate, mass_liquid, mass_vapor, h_liquid, h_vapor, x_liquid, x_vapor, den_liquid, den_vapor,
	       P, T, E;
	// flow conditions
	double m_l_in, x_l_in, h_l_in; 
	double m_v_in, x_v_in, h_v_in;
	double m_l_out, x_l_out, h_l_out; 
	double m_v_out, x_v_out, h_v_out;

};
