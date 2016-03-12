    if (runTime.outputTime())
    {
        volVectorField gradPhi(fvc::grad(phi));	

	
        volScalarField gradPhix
        (
            IOobject
            (
                "gradPhix",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradPhi.component(vector::X)
        );

        volScalarField gradPhiy
        (
            IOobject
            (
                "gradPhiy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradPhi.component(vector::Y)
        );

        volScalarField gradPhiz
        (
            IOobject
            (
                "gradPhiz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradPhi.component(vector::Z)
        );

	S = mag(gradPhi.component(vector::X)) + mag(gradPhi.component(vector::Y)) + mag(gradPhi.component(vector::Z));
	
        runTime.write();
    }
