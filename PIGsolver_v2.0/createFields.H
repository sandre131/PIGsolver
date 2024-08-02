    Info<< "Reading plasma properties\n" << endl;

    IOdictionary plasmaProperties
    (
        IOobject
        (
            "plasmaProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
/*
    dimensionedScalar Da
    (
        plasmaProperties.lookup("Da")
    );
*/

//    dimensionedScalar NUi
//    (
//        plasmaProperties.lookup("NUi")
//    );
/*
    dimensionedScalar LAMBDAe
    (
        plasmaProperties.lookup("LAMBDAe")
    );
*/

// me/2ma
    dimensionedScalar DELTA
    (
        plasmaProperties.lookup("DELTA")
    );

/*
    dimensionedScalar Esquared
    (
        plasmaProperties.lookup("Esquared")
    );
*/

/*
    dimensionedScalar SIGMA
    (
        plasmaProperties.lookup("SIGMA")
    );
*/
    dimensionedScalar Ei
    (
        plasmaProperties.lookup("Ei")
    );

    dimensionedScalar NUc
    (
        plasmaProperties.lookup("NUc")
    );

    dimensionedScalar podgonsec
    (
        plasmaProperties.lookup("podgonsec")
    );

    dimensionedScalar podgonkelvin
    (
        plasmaProperties.lookup("podgonkelvin")
    );

    dimensionedScalar podgonkg
    (
        plasmaProperties.lookup("podgonkg")
    );

    dimensionedScalar podgonmetr
    (
        plasmaProperties.lookup("podgonmetr")
    );


    dimensionedScalar podgonforNUi
    (
        plasmaProperties.lookup("podgonforNUi")
    );
       
    dimensionedScalar podgonDa
    (
        plasmaProperties.lookup("podgonDa")
    );
    
    dimensionedScalar podgonLAMBDAe
    (
        plasmaProperties.lookup("podgonLAMBDAe")
    );       
/*
    dimensionedScalar sigma
    (
        financialProperties.lookup("sigma")
    );

    dimensionedScalar sigmaSqr = sqr(sigma);
*/

    Info<< nl << "Reading field Ne" << endl;

    volScalarField Ne
    (
        IOobject
        (
            "Ne",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field Te" << endl;

    volScalarField Te
    (
        IOobject
        (
            "Te",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volVectorField Va
    (
        IOobject
        (
            "Va",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
        //для volfield - mesh.C()
        //mesh.Cf()
    );   

    volScalarField Pa
    (
        IOobject
        (
            "Pa",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
    );


    volScalarField Ta
    (
        IOobject
        (
            "Ta",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
    );


/*
    surfaceVectorField Pf
    (
        IOobject
        (
            "Pf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.Cf()
    );
*/
    
    //surfaceVectorField Va
/*
    volVectorField Va
    (
        IOobject
        (
            "Va",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.C()
        //для volfield - mesh.C()
        //mesh.Cf()
    );
*/

/*    
    volVectorField UMean
    (
        IOobject
        (
            "UMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
        //для volfield - mesh.C()
        //mesh.Cf()
    );    
*/

/*
//надо делать так и присваивать в проге с помощью == чтобы глюков не было
    volVectorField UMean
    (
        IOobject
        (
            "UMean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
        //для volfield - mesh.C()
        //mesh.Cf()
    );
    
    volVectorField Va
    (
        IOobject
        (
            "Va",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
        //для volfield - mesh.C()
        //mesh.Cf()
    );    

    volScalarField P
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
    );
    
    volScalarField Pa
    (
        IOobject
        (
            "Pa",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
    );
    
//читаю Ta из рассчитанного
//для этого в кейсе 2 файла с температурой - bondaryT и overallT вначале одинаковые

    volScalarField overallT
    (
        IOobject
        (
            "overallT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField Ta
    (
        IOobject
        (
            "Ta",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        //!!надо обязательно mesh тогда размерность сохраняется!
        mesh
    );    
*/

//читаю Va из рассчитанного
/*
    volVectorField Va
    (
        IOobject
        (
            "Va",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C()
        //для volfield - mesh.C()
        //mesh.Cf()
    );


    volScalarField Pa
    (
        IOobject
        (
            "Pa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
	    //IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/

//читаю давление Pa из рассчитанного
/*
    volScalarField Pa
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
	    //IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
        mesh
    );
*/

//Ve только для записи
/*
    volVectorField Ve
    (
        IOobject
        (
            "Ve",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C()
    ); 
*/
  

//!!похоже этот файл только для чтения из папки 0, после рассчетов переопределять переменные из новых временных папок о не будет!!
  
 
/*
    volScalarField Ta
    (
        IOobject
        (
            "Ta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
*/


/*
    volScalarField Esquared
    (
        IOobject
        (
            "Esquared",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
*/


    	dimensionedScalar Esqrcoeff
	(
		"Esqrcoeff",
		dimensionSet(0, 0, -2, 2, 0, 0 ,0),
		1E6
	);
/*
    	dimensionedScalar podgonmetr
	(
		"podgonmetr",
		dimensionSet(0, -1, 0, 0, 0, 0 ,0),
		1.0
	);
*/
    volScalarField Esquared
    (
        IOobject
        (
            "Esquared",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
//E^2=E_r^2+E_fi^2+E_z^2
	Esqrcoeff*(Foam::pow(5*Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2))*(-16*Foam::pow(podgonmetr*mesh.C().component(2),2)+8*podgonmetr*mesh.C().component(2)),2)     +     Foam::pow((-44.4*Foam::pow(Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2)),2)+4.44*Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2))+0.888)*(-16*Foam::pow(podgonmetr*mesh.C().component(2),2)+8*podgonmetr*mesh.C().component(2)),2)     +      Foam::pow((1-Foam::pow(Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2)),2))*2*podgonmetr*mesh.C().component(2),2))
//      Делаем по z в 4 раза больше чтобы сгладить провал в центре
	//Esqrcoeff*(Foam::pow(5*Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2))*(-16*Foam::pow(podgonmetr*mesh.C().component(2),2)+8*podgonmetr*mesh.C().component(2)),2)     +     Foam::pow((-44.4*Foam::pow(Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2)),2)+4.44*Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2))+0.888)*(-16*Foam::pow(podgonmetr*mesh.C().component(2),2)+8*podgonmetr*mesh.C().component(2)),2)     +      4*Foam::pow((1-Foam::pow(Foam::sqrt(Foam::pow(podgonmetr*mesh.C().component(0),2)+Foam::pow(podgonmetr*mesh.C().component(1),2)),2))*2*podgonmetr*mesh.C().component(2),2))

//	Esqrcoeff*Foam::exp(-1000*podgonmetr*podgonmetr*mesh.C().component(0)*mesh.C().component(0) - 1000*podgonmetr*podgonmetr*mesh.C().component(1)*mesh.C().component(1)  -(podgonmetr*mesh.C().component(2)-0.1)*(podgonmetr*mesh.C().component(2)-0.1)+0.00001)
    );


//частота ионизации

    	dimensionedScalar NUicoeff
	(
		"NUicoeff",
		dimensionSet(0, 0, -1, 0, 0, 0 ,0),
		1E7
	);
/*
    volScalarField NUi
    (
        IOobject
        (
            "NUi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	NUicoeff*Foam::exp(-1000*podgonmetr*podgonmetr*mesh.C().component(0)*mesh.C().component(0) - 1000*podgonmetr*podgonmetr*mesh.C().component(1)*mesh.C().component(1)  -(podgonmetr*mesh.C().component(2)-0.1)*(podgonmetr*mesh.C().component(2)-0.1)+0.00001) 
    );
*/
/*
    volScalarField NUi
    (
        IOobject
        (
            "NUi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/
/*
    volVectorField Ve
    (
        IOobject
        (
            "Ve",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	Va-(Da/Ne)*fvc::grad(Ne)
        //mesh.C()
        //для volfield - mesh.C()
        //mesh.Cf()
    );
*/

//-    V == max(strike - P.x(), dimensionedScalar("0", V.dimensions(), 0.0));
/*
    V == max
    (
        P.component(Foam::vector::X) - strike,
        dimensionedScalar("0", V.dimensions(), 0.0)
    );
*/

/*
    volScalarField delta
    (
        IOobject
        (
            "delta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(V)().component(Foam::vector::X)
    );
*/
