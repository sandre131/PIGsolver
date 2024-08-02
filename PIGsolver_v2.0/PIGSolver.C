/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    PIGsolver v1.0

Description
    Direct simulation Monte Carlo (DSMC) solver for 3D, transient, multi-
    species flows +
    Solves the N_e T_e equations to plasma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"

#include "IOdictionary.H"

//к N_e T_e
#include "writeCellGraph.H"
#include "OSspecific.H"

//добавил для чтения Umean
//#include "dsmcFields.H"

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

    //к N_e T_e
    #include "createFields.H"

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
    
//чтобы считать writeInterval
/*
    IOdictionary* controlDict = new IOdictionary( IOobject
        (
            "../system/controlDict",
            dsmc.mesh().time().constant(),
            dsmc.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ));
*/    
    
/*
    const scalar n_e_deg
    (
        readScalar(dsmcInitialisePlasma->lookup("n_e_deg"))
    );

    const scalar T_e
    (
        readScalar(dsmcInitialisePlasma->lookup("T_e"))
    );
*/
    const scalar C_p
    (
        readScalar(dsmcInitialisePlasma->lookup("C_p"))
    );


    scalar k_Bolz = 1.3806488E-23;
    scalar m_e = 9.10938291E-31;
    //scalar k_Bolz = Foam::constant::physicoChemical::k;
    //scalar n_e0 = 1.0E+18; //nado e+16-e+18
    //scalar T_e = 40000; //nado 10000 - 40000
    //teploemkost argona Dj/kg*gradus
    //scalar C_p = 530;
    //scalar m_e = Foam::constant::atomic::me;
    
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
    
/*
volVectorField Va("Va",mesh.C());

              for(int i=0;i<Va.size();i++)
              {
                if(Va[i].x()<0.0001)
                {
                  Va[i].x()=0.0001;
                  //Info << "Va!!!!=" << Va[i] << "Va!!!!";
                };
                if(Va[i].y()<0.0001)
                {
                  Va[i].y()=0.0001;
                };
                if(Va[i].z()<0.0001)
                {
                  Va[i].z()=0.0001;
                };               
              };
*/              
            
             
//volScalarField Pa("Pa",Va&Va,dimensionSet(-2,2,4,0,0,0,0));
 /*   volScalarField Pa2
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

 
volScalarField Pa = Pa2;

              for(int i=0;i<Pa.size();i++)
              {
                if(Pa[i]<0.0001)
                {
                  Pa[i]=0.0001;
                };               
              };

 */

   
    dimensionedScalar Kb("Kb",constant::physicoChemical::k.value());


//к N_e T_e    
    Info<< nl << "Calculating value(price of comodities)" << endl;    
    volVectorField NeVa("NeVa", Ne*Va);
    

    //частота ионизации порядка 1E7
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
        //2.275*podgonforNUi*pow(10,6.)*(Pa/133.0)*pow(1.5*Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,5.34)
	2.275*podgonforNUi*pow(10,4.)*(Pa/133.0)*pow(1.5*Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,5.34)
    );       
    //volScalarField NUi("NUi", 2.275*podgonforNUi*pow(10,6.)*Pa*pow(1.5*Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,5.34));
    
    //Da-амбиполярная диффузия
    volScalarField Da
    (
        IOobject
        (
            "Da",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //0.1*podgonDa*Kb*Te/Pa/constant::electromagnetic::e.value()
	10.0*podgonDa*Kb*Te/Pa/constant::electromagnetic::e.value()
    );
    
    volScalarField Pa_old("Pa_old", Pa);    
    //volScalarField Da("Da", 0.1*podgonDa*Kb*Te/Pa/constant::electromagnetic::e.value());
    
    //volScalarField NUi("NUi", Pa);
   
    volVectorField DIVCoeff("DIVCoeff", 2.5*Kb*(Ne*Va - Da*fvc::grad(Ne)));
    volScalarField Source2("Source2", 1.5*Kb*DELTA*NUc*Ne);
   //oldold volScalarField SIGMA("SIGMA", Ne*pow(constant::electromagnetic::e.value(),2.0)*NUc/(constant::atomic::me.value()*(pow(NUc,2.0)+pow(2*Foam::constant::mathematical::pi*1.76*1E6,2.0))));
 
    volScalarField SIGMA
    (
        IOobject
        (
            "SIGMA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ne*constant::electromagnetic::e.value()*constant::electromagnetic::e.value()*NUc/(constant::atomic::me.value()*(pow(podgonsec*2.0*constant::mathematical::pi*1.76*pow(10,6.),2.)+pow(NUc,2.)))
    );      
    //volScalarField SIGMA("SIGMA", Ne*pow(constant::electromagnetic::e.value(),2.)*NUc/(constant::atomic::me.value()*(pow(podgonsec*2*constant::mathematical::pi*1.76*1E6,2.)+pow(NUc,2.))));

    //LAMBDAe-электронная теплопроводность
    volScalarField LAMBDAe
    (
        IOobject
        (
            "LAMBDAe",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	//5.65*podgonLAMBDAe*pow(10,-14.)*Ne
        5.65*podgonLAMBDAe*pow(10,-14.)*Ne*pow(Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,0.5)
    );

//тепловыделение??
    volScalarField Qj
    (
        IOobject
        (
            "Qj",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	//5.65*podgonLAMBDAe*pow(10,-14.)*Ne
        SIGMA*Esquared
    );

    //volScalarField LAMBDAe("LAMBDAe", 5.65*podgonLAMBDAe*pow(10,-8)*Ne*pow(Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,0.5));

    //для заполенения всех значений массива
//    Ve == Va;

//
    
    Info<< "\nStarting time loop\n" << endl;
/*
    Info<< "\nn_e0=\n" << endl;
    Info<< n_e_deg << endl;
    Info<< "\nT_e=\n" << endl;
    Info<< T_e << endl;
*/
    Info<< "\nC_p=\n" << endl;
    Info<< C_p << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

//блок изменения скорости за счет подогрева электронами(создание плазмы)
	volScalarField rhoMcur = dsmc.rhoM();
        
        volScalarField rhoNcur = dsmc.rhoN();

	//volScalarField A = dsmc.internalE();
	scalar deltaT = dsmc.mesh().time().deltaTValue();


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
//Info << "delta_m!!!!=" << delta_m << "delta_m!!!!";
	    //physicoChemical::k.value()
	   vector& Ui = p.U();
	   vector Ui_new = p.U() ;

	    //   Ui_new.x() = 1;

	   double n_e_cur = Ne[p.cell()];
	   double T_e = Te[p.cell()];
           //мб Ta будет редко обновляться
           double T_a = Ta[p.cell()];
	   //n_ecur(rhoMcur[p.cell()], n_e_deg);U.mesh().C()[cellI]
	   //какая-то проблемма с rhoMcur postavil>e-8
	   if(rhoMcur[p.cell()]>0.00000001)
	   {

	   //  double n_e_cur = n_e(p.position().x(),p.position().y(),p.position().z(),n_e0);
	    // double n_e_cur = n_ecur(rhoMcur[p.cell()], n_e_deg);
/*
    Info<< "\nn_e_cur=\n" << endl;
    Info<< n_e_cur << endl;
    Info<< "\nrhocur\n" << endl;
    Info<< rhoMcur[p.cell()] << endl;
*/
	   //  Info<< "n_e_cur" << n_e_cur << endl;
	
          //считаем как dE=E/dt=(m_a*V_new^2/2dt-m_a*V_old^2/2dt)=(3/2)k*delta*NUc*Ne*(Te-Ta) отсюда V_new^2 = ..       
        
        
	     scalar Voldsquare = p.U().x()*p.U().x() + p.U().y()*p.U().y() +p.U().z()*p.U().z();

// разобраться с формулой!!
//V_new^2 = V_old^2+(3/m_a)*k_Bolz*delta*NUc*Ne*(Te-Ta)*deltaV*dt*napravlenie_V ; deltaV - objem jachejki мб m_a-масса всехчастиц в тек ячейке? m_a = m_a*chislo_chastic
	     Ui_new.x() = Foam::sqrt(Voldsquare+3.0*(k_Bolz*delta_m*NUc.value()*n_e_cur/(dsmc.constProps(p.typeId()).mass()*rhoNcur[p.cell()]))*(T_e-T_a)*deltaT*mesh.V()[p.cell()])*(p.U().x()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));

	     Ui_new.y() = Foam::sqrt(Voldsquare+3.0*(k_Bolz*delta_m*NUc.value()*n_e_cur/(dsmc.constProps(p.typeId()).mass()*rhoNcur[p.cell()]))*(T_e-T_a)*deltaT*mesh.V()[p.cell()])*(p.U().y()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));

	     Ui_new.z() = Foam::sqrt(Voldsquare+3.0*(k_Bolz*delta_m*NUc.value()*n_e_cur/(dsmc.constProps(p.typeId()).mass()*rhoNcur[p.cell()]))*(T_e-T_a)*deltaT*mesh.V()[p.cell()])*(p.U().z()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));


//Info << "Voldsquare!!!!+=" << Voldsquare <<"+" << 3.0*(k_Bolz*delta_m*NUc.value()*n_e_cur/(dsmc.constProps(p.typeId()).mass()*rhoNcur[p.cell()]))*(T_e-T_a)*deltaT*mesh.V()[p.cell()]<< "UIOLD!!!!=";
             

 //старый вариант почему так не помню
/*
	     Ui_new.x() = Foam::sqrt(Voldsquare+4.5*(k_Bolz*k_Bolz*delta_m*n_e_cur/(C_p*rhoMcur[p.cell()]))*(T_e-dsmc.constProps(p.typeId()).mass()*(Voldsquare)/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass())*(p.U().x()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));

	     Ui_new.y() = Foam::sqrt(Voldsquare+4.5*(k_Bolz*k_Bolz*delta_m*n_e_cur/(C_p*rhoMcur[p.cell()]))*(T_e-dsmc.constProps(p.typeId()).mass()*(Voldsquare)/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass())*(p.U().y()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));

	     Ui_new.z() = Foam::sqrt(Voldsquare+4.5*(k_Bolz*k_Bolz*delta_m*n_e_cur/(C_p*rhoMcur[p.cell()]))*(T_e-dsmc.constProps(p.typeId()).mass()*(Voldsquare)/(3*k_Bolz))*deltaT/dsmc.constProps(p.typeId()).mass())*(p.U().z()/Foam::sqrt(p.U().x()*p.U().x()+p.U().y()*p.U().y()+p.U().z()*p.U().z()));
*/
	     
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

              //т.к. проблема с вычислением градиента или Ne слишком большое считаю Ve вручную надо проверять мб неверно!
              //Ve == Va-(Da/Ne)*fvc::grad(Ne);
/*              
              if(i>0)
              {
                dsmcParcel& p2 = *cellParcels[i-1];
                if(Ne[p.cell()]<1E19)
                {
                  
                  Ve[p.cell()].x() = Va[p.cell()].x() - Da[p.cell()]*(Ne[p.cell()]-Ne[p2.cell()])/(p.position().x()-p2.position().x())/Ne[p.cell()];
                  Ve[p.cell()].y() = Va[p.cell()].y() - Da[p.cell()]*(Ne[p.cell()]-Ne[p2.cell()])/(p.position().y()-p2.position().y())/Ne[p.cell()];
                  Ve[p.cell()].z() = Va[p.cell()].z() - Da[p.cell()]*(Ne[p.cell()]-Ne[p2.cell()])/(p.position().z()-p2.position().z())/Ne[p.cell()];
                
                };
              }
*/

	    }; //if
	  }; // last forall

	//};
      };



        dsmc.evolve();

        dsmc.info();

        //вставил из конца
        runTime.write();
        
//    help dimensions        
//1 	Mass 	kilogram (kg) 	
//2 	Length 	metre (m) 	
//3 	Time 	second (s) 
//4 	Temperature 	Kelvin (K) 	
//5 	Quantity 	kilogram-mole (kgmol) 	
//6 	Current 	ampere (A)
//7 	Luminous intensity  candela (cd)

       
//расчет N_e T_e
        NeVa == Ne*Va;
//        Info << "NeVa!!!!=" << NeVa << "NeVa!!!!";
//проблемма с размерностью в дивергенции Ne*Va (на NUi там домножать не надо!)




if(runTime.value()>0.001)
{
  //readScalar(Foam::debug::controlDict::lookup("writeInterval"))
 //doubleScalar writetime1 = readScalar(runTime.controlDict().lookup("writeInterval"));
 //doubleScalar writetime1 = runTime.outputTime() -0.0001;
 //word& previoustimedir = runTime.findClosestTime(runTime.value()).name() ;
 //runTime
 
 //Foam::Time::setWriteInterval(0.000002);
 word previoustimedir = runTime.timeName(runTime.findClosestTime(runTime.value()).value());
 //Info << "runTime.timeName(-1) =" <<  runTime.outputTime() - runTime.controlDict().lookup("writeInterval") <<endl;
 //Info << "runTime.timeName(-1) =" <<  runTime.timeName(runTime.findClosestTime(runTime.value()).value()) << "runTime.timeName(-1)!!!" <<endl;
 //Info << "runTime.timeName(-1) =" <<  runTime.outputTime() << "runTime.timeName(-1)!!!" <<endl;
 
 //чтобы по одному разу вокруг точек записи считались Ne Te надо 0.000002
//  if (((runTime.value() - runTime.findClosestTime(runTime.value()).value())<0.005)&&(!runTime.outputTime()))
  if (!runTime.outputTime())
  {
      

  //читаем U P T здесь чтобы обновить из последних Va Ta Pa    
    volVectorField UMean
    (
        IOobject
        (
            "UMean",
            previoustimedir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
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
            previoustimedir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
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
            previoustimedir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

  
//по нормальному надо доставть поля из dsmc
//блок объявления переменных из dsmc для использования при расчете Ne Te
            if (min(mag(dsmc.rhoN())).value() > VSMALL)
            {
              //if(runTime.timeName()==0.002)
              //{
              //переопределяю поля из dsmc как это там сделано т.к. читать из следующей по времени директории не получается новопосчитанные поля
              //Va = Foam::dsmcFields::read("UMean");
              /*
        if (runTime.outputTime())
        {              
              volScalarField rhoM(IOobject("rhoM",
runTime.timeName(), mesh,
IOobject::MUST_READ,
IOobject::NO_WRITE), mesh);
              
              volVectorField momentum(IOobject("momentum",
runTime.timeName(), mesh,
IOobject::MUST_READ,
IOobject::NO_WRITE), mesh);
              
                      Va = momentum/rhoM;
        };
        */
              
              //Info << "Umean!!!!=" << UMean << "Umean!!!!";
          //надо ставить == !!! и границы не менять и меньше нуля тоже - иначе глюки! т.к. там покомпонентно и там мб отриц знач
          //можеть быть надо ==
              Va == UMean;
              forAll( Va.boundaryField(), facej)
              {
               // Pa.boundaryField()[facej]=(max(mag(Pa)).value()+min(mag(Pa)).value())/2.0;
                Va.boundaryField()[facej]=UMean.boundaryField()[facej];
              };              
              //Info << "Va!!!!=" << Va << "Va!!!!";       

     /* 
              const volScalarField& rhoMMean = Foam::dictionary::lookup
              (
                  "rhoMMean",1,1
              );
      
              const volVectorField& momentumMean = Foam::dictionary::lookup
              (
                  "momentumMeanName",1,1
              );
       */       
              //Va = (1.0)*dsmc.momentum()/dsmc.rhoM();
              
              //forAll(Va, i)
              //{
/*
              for(int i=0;i<Va.size();i++)
              {
                if(Va[i].x()<0.0001)
                {
                  Va[i].x()=0.0001;
                  //Info << "Va!!!!=" << Va[i] << "Va!!!!";
                };
                if(Va[i].y()<0.0001)
                {
                  Va[i].y()=0.0001;
                };
                if(Va[i].z()<0.0001)
                {
                  Va[i].z()=0.0001;
                };               
              };
              
              //ставим на границах типа везде скорость =1 чтобы уе для Ne считалось
              forAll( Va.boundaryField(), facej)
              {
                Va.boundaryField()[facej]=vector(0.0001,0.0001,0.0001);
              }
*/              
              //Info << "Va!!!!=" << Va << "Va!!!!";
              //Info << min(mag(Va & Va)).value() << " minmag!!!";
              
              //Info << nl << "!!!!VAminmag= " << min(mag(Va)).value();
              //if (min(mag(Va)).value() > VSMALL)
              //{


//Info << "Boundary =" <<  Pa.boundaryField().patchInternalField() <<endl;              
      
              //You are right about the temperature units. I would guess that the units of the Boltzmann constant are not passed through the .value() function (even though it is correctly set up in etc/controlDict). If k was dimensionless, the units would match those given in the temperature files. The same thing would happen with pressure p but in this case k is both in the numerator and denominator and therefore it's ok. I wouldn't worry too much about that anyway. 
              //overallT и p размерность не K и Па т.к. Кb почему-то безразмерна см. etc/controlDict
//!!!тут надо сделать нормальные размерности у Ta и Pa!!! - вроде бы сделал!
              Ta ==  overallT*podgonkelvin*podgonkg*podgonmetr*podgonmetr/podgonsec/podgonsec;
              //volScalarField trT = 2.0/(3.0*Kb*dsmc.rhoN())*(dsmc.linearKE() - 0.5*dsmc.rhoM()*(Va & Va)*podgonsec*podgonsec);
              Pa ==  P/podgonmetr;

              
              //Ta.boundaryField() == overallT.boundaryField();
              //Pa.boundaryField() == P.boundaryField();
            //  Ta ==  podgonkg*podgonmetr*podgonmetr*podgonkelvin/(podgonsec*podgonsec)*2.0/(Kb*(3.0*dsmc.rhoN() + dsmc.iDof()))*(dsmc.linearKE() - 0.5*dsmc.rhoM()*(Va & Va) + dsmc.internalE());
              //volScalarField trT = 2.0/(3.0*Kb*dsmc.rhoN())*(dsmc.linearKE() - 0.5*dsmc.rhoM()*(Va & Va)*podgonsec*podgonsec);
            //  Pa ==  Kb*dsmc.rhoN()*(2.0/(3.0*Kb*dsmc.rhoN())*(dsmc.linearKE() - 0.5*dsmc.rhoM()*(Va & Va)));

           
              forAll(Pa, i)
              {              
              //for(int i=0;i<Pa2.size();i++)
              //{
                if(Pa[i]<0.0001)
                {
                  Pa[i]=0.1;
                //Info << Pa[i] << "!!!!";
                };
              };
              
//!!!!!problem!!!!!
              //!!!!поставил на границе для Pa 1 Па надо как-то это обойти т.е. заставить считать только по внутренней сетке или поставить calculated в типе гр поля!!!!
              forAll( Pa.boundaryField(), facej)
              {
               // Pa.boundaryField()[facej]=(max(mag(Pa)).value()+min(mag(Pa)).value())/2.0;
                Pa.boundaryField()[facej]=(1.0);
              };

              //Info << "Ta!!!!=" << Ta << "Ta!!!!";
/*              
              forAll(Ta, i)
              {              
              //for(int i=0;i<Pa2.size();i++)
              //{
                if(Ta[i]<0.0001)
                {
                  Ta[i]=0.0001;
                //Info << Pa[i] << "!!!!";
                };
              };
                
              //!!!!поставил на границе для Ta 300 К надо как-то это обойти!!!!
              forAll( Ta.boundaryField(), facej)
              {
               // Pa.boundaryField()[facej]=(max(mag(Pa)).value()+min(mag(Pa)).value())/2.0;
                Ta.boundaryField()[facej]=(300.0);
              };
*/
/*
              forAll( Pa.boundaryField(), facej)
              {
               // Pa.boundaryField()[facej]=(max(mag(Pa)).value()+min(mag(Pa)).value())/2.0;
                Pa.boundaryField()[facej]==P.boundaryField()[facej];;
              };
*/
              
              forAll( Ta.boundaryField(), facej)
              {
               // Pa.boundaryField()[facej]=(max(mag(Pa)).value()+min(mag(Pa)).value())/2.0;
                Ta.boundaryField()[facej]==overallT.boundaryField()[facej];
              };

              //Pa.write();
            //writeCellGraph(Pa, runTime.graphFormat());
              //};
              
            };
//конец расчета блока объявления переменных из dsmc для использования при расчете Ne Te





//!!то что под pow надо чтобы было double иначе проблема - неправильно вычисляет!!
//Info << Pa/133.0 << "!!!!";
        //nui esli bolshe v 10 raz uze rashoditsa !!!проверить эту формулу!!!мб 0.01*Pa заменить на Pa/133
        
//!!! менять сверху эти функции в определениях файлов для записи!!! 
        //NUi == 2.275*podgonforNUi*pow(10,6.)*(Pa/133.0)*pow(1.5*Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,5.34);
	//разница между Da и NUi: NUi=(10^5 или 10^6)*Da приблизительно
        NUi == 2.275*podgonforNUi*pow(10,4.)*(Pa/133.0)*pow(1.5*Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,5.34);
        //Da == 0.1*podgonDa*Kb*Te/Pa/constant::electromagnetic::e.value();
	Da == 10.0*podgonDa*Kb*Te/Pa/constant::electromagnetic::e.value();
        SIGMA == Ne*constant::electromagnetic::e.value()*constant::electromagnetic::e.value()*NUc/(constant::atomic::me.value()*(pow(podgonsec*2.0*constant::mathematical::pi*1.76*pow(10,6.),2.)+pow(NUc,2.)));
        //Info << "SIGMA!!!!=" << SIGMA << "SIGMA!!!!";
        //Info << "Ne!!!!=" << Pa.internalField() << "Pa!!!!";
       // Info << "Ne!!!!=" << Ne << "Ne!!!!";
        //NUi == Pa;

	Qj == SIGMA*Esquared;
        
//Info << "fvm::laplacian(Da, Ne)!!!!=" << (fvm::laplacian(Da, Ne))<< "fvm::laplacian(Da, Ne)!!!!";

//if(runTime.value()>0.01 + readScalar(runTime.controlDict().lookup("writeInterval")))
if(runTime.value()>0.26 + 0.00005)
{
  
      Info << "!!!!!ZapuskNeTe" << runTime.value() - runTime.findClosestTime(runTime.value()).value() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl;
      
        solve
        (
//явная схема

            fvm::ddt(Ne)
          - fvm::laplacian(Da, Ne)  
	  - NUi*Ne
         ==
          - fvc::div(NeVa)

//неявная схема

//            fvm::ddt(Ne)
//          - fvm::laplacian(Da, Ne)  
//         ==
//            fvm::Sp(NUi, Ne)
//          - podgonsec*fvc::div(NeVa)

        );
//Info << "NeVa!!!!=" << NeVa[100] << "NeVa!!!!";


        //DIVCoeff == 2.5*Kb*Ne*Ve;
        DIVCoeff == 2.5*Kb*(Ne*Va - Da*fvc::grad(Ne)); //делить на podgonsec не надо - размерность:(
        //Ve находится так:	
	//Ve == Va-(Da/Ne)*fvc::grad(Ne);
        Source2 == 1.5*Kb*DELTA*NUc*Ne;
        //NUiEi == NUi*Ei;
        //SIGMA == Ne*pow(constant::electromagnetic::e.value(),2)*NUc/(constant::atomic::me.value()*(pow(NUc,2.0)+pow(2*Foam::constant::mathematical::pi*1.76*1E6,2.)));
        //SIGMA = Ne*(e^2)*NUc/(me(NUc^2+omega^2))
        
        
        LAMBDAe == 5.65*podgonLAMBDAe*pow(10,-14.)*Ne*pow(Kb*Te/constant::electromagnetic::e.value()/podgonkelvin,0.5);
	//LAMBDAe == 5.65*podgonLAMBDAe*pow(10,-14.)*Ne;

        // Sp - источник для неявной схемы       
        solve
        (
/*          
            2.5*Kb*Te*constant::atomic::me.value()*Ne*fvm::ddt(Te)
          + fvm::laplacian(LAMBDAe, Te) 
	  - podgonkelvin*fvc::div(DIVCoeff*Te) //podgonsec*podgonkelvin - пробл с размерностью
          + Source2*Te
         ==
	    Source2*Ta
	  + SIGMA*Esquared
	  - NUi*Ei*Ne
*/

            2.5*Kb*Te*constant::atomic::me.value()*Ne*fvm::ddt(Te)
          - fvm::laplacian(LAMBDAe, Te) 
	  + podgonkelvin*fvc::div(DIVCoeff*Te) //podgonsec*podgonkelvin - пробл с размерностью
          + Source2*Te
         ==
	    Source2*Ta
	  + SIGMA*Esquared
	  - NUi*Ei*Ne
          
	   //пробл с разм - NUi*Ei*Ne
	  //размерность неправ - fvc::div(DIVCoeff*Te)
        );

        //потом записываю скорость электронов из расчета
        //  Ve == Va-(Da/Ne)*fvc::grad(Ne);
};
        
   //if  в районе записываемой точки     
     };

//end if runtime>..
};

      //  runTime.write();

      if (!runTime.outputTime())
      {
        volScalarField Pa_old("Pa_old", Pa);
      };
      
        if (runTime.outputTime())
        {
          
          
	 // raznost mezhdu tochnum i pribl
         volScalarField Pa_new_minus_Pa_old("Pa_new_minus_Pa_old", Pa - Pa_old);
         //Ne_new_minus_Ne_old == Ne - Ne_old;
         Pa_old == Pa;
        


        runTime.write();

        scalar maxRes = (max(Pa_new_minus_Pa_old)).value();
        scalar minRes = (min(Pa_new_minus_Pa_old)).value();
        if(abs(minRes)>abs(maxRes))
        {
          maxRes == minRes;
        };
        
        scalar AbsPogr = maxRes/(max(Pa)).value();
        Info<< "||C|| Max(Pa_new_minus_Pa_old)/max(Pa_new) = " << AbsPogr << " !!" << endl;
        
 


            writeCellGraph(Ne, runTime.graphFormat());
            writeCellGraph(Te, runTime.graphFormat());
            
//            Info << nl << "!!!!VSMALL = " << VSMALL ;
//            Info << nl << "!!!!RHONminmag= "<< min(mag(dsmc.rhoN())).value() ;

            
//Pa == Foam::dsmcFields::overallT;
//Pa.read();
        }

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
