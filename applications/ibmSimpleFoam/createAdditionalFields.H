    Info<< "Membrane force source term p\n" << endl;
    volVectorField F
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedVector("F",dimensionSet(1,-2,-2,0,0,0,0),vector::zero)
    );

