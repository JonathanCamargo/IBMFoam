#ifndef IBM_H
	#define IBM_H

	#include "mesh.h"
	#include "fluid.h"

	using namespace std;

namespace Tera{
	void interpolation(fluid fluido, mesh membrana, int x, int y, int z);
	void spread(fluid fluido, mesh membrana, int x, int y, int z);
	double dirac_2(double *x);
	double dirac_3(double *x);
	double dirac_4(double *x);
	double phi_2(double r);
	double phi_3(double r);
	double phi_4(double r);
} //end namespace
#endif
