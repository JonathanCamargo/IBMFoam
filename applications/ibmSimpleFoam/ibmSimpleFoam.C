 /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow modified to use immerse boundary
    method with Terablood.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

#include "interpolation.H"

#include "mesh.h"
#include "ibm.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createAdditionalFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);
    dimensionedScalar iteraciones=runTime.endTime();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        int N=1;
        Tera::mesh membrana[N];
        Tera::mesh referencia[N];
        double dt = 1.0;
        double dx = 1.0;
        int X = 40;
        int Y = 100;
        int Z = 40;
        int VTK = 10;

        // Parametros adimensionales
	dimensionedScalar rho("rho",dimensionSet(1,-3,0,0,0,0,0),1.0); 
        double nu = 1./6.;
        double Re = 50.0;
        double G = 0.1;
        double R = 10;
	double scale=0.1;
	double gamma_dot = (Re*nu)/(rho.value()*Foam::pow(R,2));
        double ks = (gamma_dot*nu*R)/(G);
        double kb = ks*1.0e-6;
        double kp = (gamma_dot)/(G);
        double STEPS = 100000/kp;


	Info<<"Generando Geometría de la membrana"<<endl;
	 // Membrana
        for(int i = 0; i < N; i++)
        {
                membrana[i].setID(i);
                membrana[i].mesh_refine_tri4();
                membrana[i].mesh_refine_tri4();
                membrana[i].mesh_refine_tri4();
                membrana[i].proyectarEsfera(R);
                membrana[i].proyectarRBC(R);
                membrana[i].rotarEstructura(0.0,0.0,0.0);
                membrana[i].moverCentro(0,0,10.0*(i+1));
                membrana[i].iniciarGeometria();
                referencia[i].setID(i);
                referencia[i].mesh_refine_tri4();
                referencia[i].mesh_refine_tri4();
                referencia[i].mesh_refine_tri4();
                referencia[i].proyectarEsfera(R);
                referencia[i].proyectarRBC(R);
                referencia[i].rotarEstructura(0.0,0.0,0.0);
                referencia[i].moverCentro(0,0,10.0*(i+1));
                referencia[i].iniciarGeometria();
                referencia[i].actualizarGeometria();

        }


              
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;
    
    for(int ts = 0 ; ts < STEPS ; ts++)
    { 


    dictionary interpolationDict=mesh.schemesDict().subDict("interpolationSchemes");

    Info<<"\nStarging convergece loop for time: "<<ts<<endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
    
	// --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End of fluid iteration\n" << endl;
    runTime.setEndTime(runTime.endTime()+iteraciones);
    Info<<runTime.endTime()<<endl;

    Info<< "Interpolating velocities on the nodes of the membrane"<<endl;
    // ------------------------------------------------------------$
    // 1. Interpolation
    // ------------------------------------------------------------$
           
	autoPtr< interpolation<vector> > Uinterpolation=interpolation<vector>::New(interpolationDict,U);
	double postmp[3];
	for (int j=0;j<membrana[0].darNumeroNodos();j++)
	{
		membrana[0].darPosNodo(j,postmp);
		//Info<<postmp[0]<<","<<postmp[1]<<","<<postmp[2]<<tab ;
		point pos(postmp[0],postmp[1],postmp[2]);
		label celli=mesh.findCell(pos);
		vector Uj=Uinterpolation->interpolate(pos,celli);
		//Info<<"U:"<<Uj<<endl;
		membrana[0].setVelocidad(j,Uj.x(),Uj.y(),Uj.z());
	}

	Info<<"Moving nodes to new positions"<<endl;
	 // ------------------------------------------------------------$
         // 2. Encontrar nuevas posiciones de la membrana
         // ------------------------------------------------------------$
         for(int i = 0; i < N; i++)
          {
          membrana[i].moverNodos(dt, dx, Y);
          }


	Info<<"Running membrane FEM model"<<endl;
          // ------------------------------------------------------------$
          // 3. Calcular fuerzas en los nodos de la membrana
          // ------------------------------------------------------------$
          for(int i = 0; i < N; i++)
           {
           //membrana[i].calcularFuerzasHelfrich(kb);
           membrana[i].calcularFuerzasFEM(referencia[i], ks);
           }

	   Info<<"Calculating force excerted to fluid for the next fluid iteration loop"<<endl;
           // ------------------------------------------------------------$
           // 4. Propagar la densidad de fuerza hacia el fluido
           // ------------------------------------------------------------$
	   F=F*0;
	   volVectorField X=mesh.C(); //Cell center distance field
	   scalarField V=mesh.V().field();
	   forAll (X,i){
		//Info<<X[i]<<endl;
		  for (int j=0;j<membrana[0].darNumeroNodos();j++)
		        {
			double fNodo[3]={0,0,0};
			double postmp[3]={0,0,0};
			double delta;
                	membrana[0].darPosNodo(j,postmp);
			membrana[0].darFuerzaNodo(j,fNodo);
        	        point pos(postmp[0],postmp[1],postmp[2]);
			vector r=X[i]-pos; //distance to RBC node
			double distance[3]={r.x(),r.y(),r.z()};
			delta=Tera::dirac_4(distance);
			vector df(fNodo[0],fNodo[1],fNodo[2]);
			df=df*delta;
			F[i]+=df;
			}
		F[i]=F[i]/V[i];
//		if (F[i]!=vector::zero){
//			Info<<"Hola soy "<<i<<"no soy cero:"<<F[i]<<endl;
//		}
	   }


   	   // -----------------------------------------------------------------------//
           // 7. Calcular propiedades macro de la membrana
           // -----------------------------------------------------------------------//
           for(int i = 0; i<N;i++)
                {
                        membrana[i].calcularCambioArea(referencia[i]);
                        membrana[i].actualizarGeometria();
                }



		// -----------------------------------------------------------------------//
                // 9. Visualización
                // -----------------------------------------------------------------------//
                        for(int i=0;i<N;i++)
                        {
                                membrana[i].guardarVTU(ts);
                        }

	}








    return 0;
}




// ************************************************************************* //
