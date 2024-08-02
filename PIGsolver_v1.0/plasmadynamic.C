/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Created by Aleksandr Shemakhin (c)

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
    plasmadynamic

Description
    Direct simulation Monte Carlo (DSMC) solver for 3D, transient, multi-
    species flows. Applied F force.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//концентрация электронов
double n_e(double x, double y, double z, double n_e0)
{
  //d_sh - diametr strui
  double d_sh = 0.024+z/10.0;
  //radius plasmatrona
  double xmax = 0.06;
   //  Info<< z << "d_sh" << d_sh << endl;
  return n_e0*(xmax-Foam::sqrt(x*x+y*y))/d_sh;
};

//отдает плотность ионов в завис от плотности атомов
double n_ecur(double currho, double cur_ne_deg)
{
    return currho*cur_ne_deg;
};


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

//    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Constructing dsmcCloud " << endl;

    dsmcCloud dsmc("dsmc", mesh);

//    IOdictionary
/*
    fileName dictName("dsmcInitialisePlasma");
    dictionary* dsmcInitialisePlasma = new dictionary(dictName);
    const scalar n_e0
    (
        readScalar(dsmcInitialisePlasma->lookup("n_e0"))
    );
*/


    IOdictionary* dsmcInitialisePlasma = new IOdictionary( IOobject
        (
            "dsmcInitialisePlasma",
            dsmc.mesh().time().constant(),
            dsmc.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ));

    const scalar n_e_deg
    (
        readScalar(dsmcInitialisePlasma->lookup("n_e_deg"))
    );

    const scalar T_e
    (
        readScalar(dsmcInitialisePlasma->lookup("T_e"))
    );

    const scalar C_p
    (
        readScalar(dsmcInitialisePlasma->lookup("C_p"))
    );



    Info<< "\nStarting time loop\n" << endl;

    Info<< "\nn_e0=\n" << endl;
    Info<< n_e_deg << endl;
    Info<< "\nT_e=\n" << endl;
    Info<< T_e << endl;
    Info<< "\nC_p=\n" << endl;
    Info<< C_p << endl;


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

//блок изменения скорости за счет подогрева электронами(создание плазмы)
	volScalarField rhoMcur = dsmc.rhoM();

	//volScalarField A = dsmc.internalE();
	scalar deltaT = dsmc.mesh().time().deltaTValue();
	scalar k_Bolz = 1.3806488E-23;
	//scalar k_Bolz = Foam::constant::physicoChemical::k;
	//scalar n_e0 = 1.0E+18; //nado e+16-e+18
	//scalar T_e = 40000; //nado 10000 - 40000
	//teploemkost argona Dj/kg*gradus
	//scalar C_p = 530;
	scalar m_e = 9.10938291E-31;
	//scalar m_e = Foam::constant::atomic::me;



	const List<DynamicList<dsmcParcel*> > &CurcellOccupancy = dsmc.cellOccupancy();

	forAll(CurcellOccupancy, curcellI)
	{
	//ParcelType = "dsmcParcel";
          const DynamicList<dsmcParcel*>& cellParcels(CurcellOccupancy[curcellI]);
/*
    	    Info<< nl << "curcellI " << curcellI << "cellParcels.size()" << cellParcels.size()
            << nl << endl;
*/
        //label nC(cellParcels.size());

     //   if (nC > 1)
     //   {

          forAll(cellParcels, i)
          {

           dsmcParcel& p = *cellParcels[i];

	   scalar delta_m = m_e/(2*dsmc.constProps(p.typeId()).mass());

	    //physicoChemical::k.value()
	   vector& Ui = p.U();
	   vector Ui_new = p.U() ;

	    //   Ui_new.x() = 1;

	   //какая-то проблемма с rhoMcur postavil>e-8
	   if(rhoMcur[p.cell()]>0.00000001)
	   {

	   //  double n_e_cur = n_e(p.position().x(),p.position().y(),p.position().z(),n_e0);
	     double n_e_cur = n_ecur(rhoMcur[p.cell()], n_e_deg);
/*
    Info<< "\nn_e_cur=\n" << endl;
    Info<< n_e_cur << endl;
    Info<< "\nrhocur\n" << endl;
    Info<< rhoMcur[p.cell()] << endl;
*/
	   //  Info<< "n_e_cur" << n_e_cur << endl;
	
	     scalar Voldsquare = p.U().x()*p.U().x() + p.U().y()*p.U().y() +p.U().z()*p.U().z();	     

	     Ui_new.x() = Foam::sqrt(Voldsquare+4.5*(k_Bolz*k_Bolz*delta_m*n_e_cur/(C_p*rhoMcur[p.cell()]))*(T_e-dsmc.constProps(p.typeId()).mass()*(Voldsquare)/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass())*(p.U().x()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));

	     Ui_new.y() = Foam::sqrt(Voldsquare+4.5*(k_Bolz*k_Bolz*delta_m*n_e_cur/(C_p*rhoMcur[p.cell()]))*(T_e-dsmc.constProps(p.typeId()).mass()*(Voldsquare)/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass())*(p.U().y()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));

	     Ui_new.z() = Foam::sqrt(Voldsquare+4.5*(k_Bolz*k_Bolz*delta_m*n_e_cur/(C_p*rhoMcur[p.cell()]))*(T_e-dsmc.constProps(p.typeId()).mass()*(Voldsquare)/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass())*(p.U().z()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));
/*
	    Info<< nl << "nz=  " << p.U().z()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z())  << "ny=  " << p.U().y() << "nx=  " << p.U().x()
            << nl << endl;
*/

      //	    Info<< nl << "p.U().z() " << p.U().z() << " Ui_new.z " << Ui_new.z() 
      //      << nl << endl;

	    Ui= Ui_new;

      //	    Info<< nl << "p.U().z() new " << p.U().z()
      //      << nl << endl;


//	   vector Ui_new=sqrt(p.U() & p.U()+4.5*k_Bolz*k_Bolz*delta_m*n_e*(T_e-dsmc.constProps(p.typeId()).mass()*(p.U() & p.U())/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass());

	    //  p.U().x() = Ui_new.x();

//
//           vector relPos = p.position();
//      	    Info<< nl << "Particle p position " << pos(relPos.x()) << " "
//            << nl << endl;
//

	    }; //if
	  }; // last forall

	//};
      };



        dsmc.evolve();

        dsmc.info();

        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
