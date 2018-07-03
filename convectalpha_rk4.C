#include "fvCFD.H"
//#include "simpleControl.H"
//#include "surfaceInterpolationScheme.H"
//#include "linear.H"
#include "labelList.H"
#include "fixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensionedScalar SMALL_NUMBER("small", dimless, 5e-16);


template<typename T, typename Rhs>
void rk3(const T& x0, double dt, Rhs & rhs, T& x)
{
    T& u1 = rhs.u1Ref();  // pobierz referencje do pola u1 pracuje na polach zdefiniowanych w dYdt
    u1 = x0 + dt * rhs.eval(x0);

    T& u2 = rhs.u2Ref();
    u2 = 3./4 * x0 + 1./4* u1 + 1./4 * dt * rhs.eval(u1);

    x = 1./3 * x0 + 2./3* u2 + 2./3 * dt * rhs.eval(u2);
}


struct dYdt
{
    double u1, u2;

    double& u1Ref()
    {
        return u1;
    }
    double& u2Ref()
    {
        return u2;
    }

    double eval(const double &y0)
    {
        return y0;
    }
};


void test_rk3()
{
    double y0 = 1;

    double dt = 0.0001;

    dYdt rhs;

    double y = 0;

    for(int i = 0; i < 20000; ++i)
    {
        rk3(y0, dt, rhs, y);

        Info <<"time "<<(i+1)*dt <<"\t error=" << Foam::mag( y - Foam::exp( (i+1)*dt) )<< endl;

        y0 = y;
    }
    Info << "time " << endl;

}


//class FVM_Helper
//{
//    const fvMesh& mesh;

//    volScalarField k1, k2;
//    volScalarField mGradPsi;
//    volVectorField gradPsi;
//    surfaceVectorField gradPsiFace, snGradPsi;
//    surfaceScalarField mGradPsiFace, mSnGradPsi, alphaFace, phiR;

//    int Ntau, PsiCritic;


//    void LimitGradPsi2(const volScalarField & Psi, volScalarField & mGradPsi, double epsh, double PsiCritic)
//    {
//        forAll( mesh.cellCentres(), cellId)
//        {
//            if( mag(Psi[cellId]) > (PsiCritic+3.0)*epsh )
//            {
//                 mGradPsi[cellId] = scalar(1.0);
//            }
//        }

//        forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
//        {
//            const fvPatch& patch = mesh.boundary()[patchi];
//            fvPatchScalarField& gradPsiPatch = mGradPsi.boundaryField()[patchi];
//            const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

//            forAll(patch, faceId) // petla po centrach objetosci danego patcha
//            {
//                if ( mag(PsiPatch[faceId]) > (PsiCritic+3.0)*epsh )
//                {
//                    gradPsiPatch[faceId] = scalar(1.0);
//                }
//            }
//        }
//    }


//public:
//    FVM_Helper(Foam::Time & runTime, Foam::fvMesh & mesh_): mesh(mesh_),
//        k1
//        (
//            IOobject
//            (
//                "k1",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("k1", dimless, scalar(0))
//        ),
//        k2
//        (
//            IOobject
//            (
//                "k2",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("k2", dimless, scalar(0))
//        ),

//        phiR
//        (
//            IOobject
//            (
//                "phiR",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::AUTO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("phiR", (dimLength*dimLength*dimLength)/dimTime, scalar(0.0))
//        ),

//        Ntau( runTime.controlDict().lookupOrDefault("Ntau", 0) ),
//        PsiCritic( runTime.controlDict().lookupOrDefault("PsiCritc", 37) )
//    {
//    }


////==================================================================================================

//    volScalarField& u1Ref(){ return k1;}
//    volScalarField& u2Ref(){ return k2;}


//    tmp<volScalarField> rhs(const volScalarField& alpha){
//        const fvMesh & mesh = alpha.mesh();

//            alphaFace == linearInterpolate(alpha);
//            AlphaToPsi(k1, Psi, epsh);
//        //        evaluateBoundaryConditions(Psi);

//            mGradPsi == Foam::mag(fvc::grad(Psi));
//            gradPsi == fvc::grad(Psi);
//        //        setBoundaryValues(gradPsi, vector(1, 0, 0));
//        //        gradPsi.boundaryField().evaluate();

//            gradPsiFace == linearInterpolate( gradPsi );
//            LimitGradPsi2(mesh, Psi, lv.mGradPsi, epsh.value(), PsiCritic);

//        //        mGradPsiFace == linearInterpolate( mGradPsi );
//            mGradPsiFace == Foam::mag(gradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

//            snGradPsiFun(Psi, snGradPsi);           //  liczy grad dyfuzyjny miedzy dwoma sr komorek i przypisuje do pola snGradpsi

//            mSnGradPsi = Foam::mag(snGradPsi);
//        //        mGradPsiFaceCondition.boundaryField().evaluate();
//            phiR = Cf * alphaFace * ( scalar(1.0) - alphaFace ) * ( (mSnGradPsi - scalar(1.0) ) * gradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();

//            return fvc::div(phiR);
//    }
//};



int main(int argc, char *argv[])
{

    test_rk3();


    return 0;
}

















/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

//#include "fvCFD.H"
////#include "simpleControl.H"
////#include "surfaceInterpolationScheme.H"
////#include "linear.H"
//#include "labelList.H"
//#include "fixedGradientFvPatchField.H"
//#include "zeroGradientFvPatchField.H"

//// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//dimensionedScalar SMALL_NUMBER("small", dimless, 1e-16);

////!!!!!!!!!!!!!!!!!!!
//void psiInit(volScalarField& Psi)
//{
//    forAll(Psi.mesh().cellCentres(), cellI) // loop over cell centers
//    {
//        Psi[cellI] = Psi.mesh().cellCentres()[cellI].x();           //-0.5;
//    }
////    forAll(Psi.mesh().boundary(), patchi)  // loop over all boundary patches
////    {
////        const fvPatch& patch = Psi.mesh().boundary()[patchi];
////        fvPatchScalarField& psiPatch = Psi.boundaryField()[patchi];
////        forAll(patch, faceId)       //  loop over c.v. centers of the given patch
////        {
////            psiPatch[faceId] = patch.Cf()[faceId].x();          //-0.5;
////        }
////    }
//}


//void LimitGradPsi2(const fvMesh& mesh, const volScalarField & Psi, volScalarField & mGradPsi, double epsh, double PsiCritic)
//{
//    forAll( mesh.cellCentres(), cellId)
//    {
//        if( mag(Psi[cellId]) > (PsiCritic+3.0)*epsh )
//        {
//             mGradPsi[cellId] = scalar(1.0);
//        }
//    }

//    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
//    {
//        const fvPatch& patch = mesh.boundary()[patchi];
//        fvPatchScalarField& gradPsiPatch = mGradPsi.boundaryField()[patchi];
//        const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

//        forAll(patch, faceId) // petla po centrach objetosci danego patcha
//        {
//            if ( mag(PsiPatch[faceId]) > (PsiCritic+3.0)*epsh )
//            {
//                gradPsiPatch[faceId] = scalar(1.0);
//            }
//        }
//    }
//}

// first interpolate to faces then compute Psi at f from Alpha
//void AlphaToPsiF(const surfaceScalarField & AlphaF, surfaceScalarField & PsiF, dimensionedScalar epsh )
//{
//    double arg = 0;
//    forAll(AlphaF, CellId)
//    {
//        arg = (AlphaF[CellId] + SMALL_NUMBER.value())/(1.0 - AlphaF[CellId] + SMALL_NUMBER.value());
//        arg = Foam::max(0.0, arg);
//        PsiF[CellId] = epsh.value()*( Foam::log( arg ) );
//    }

//    forAll(AlphaF.mesh().boundary(), patchi)
//    {
//        const fvPatch& patch = AlphaF.mesh().boundary()[patchi];
//        fvsPatchScalarField& PsiPatchF = PsiF.boundaryField()[patchi];
//        const fvsPatchScalarField& AlphaPatchF = AlphaF.boundaryField()[patchi];
//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            arg = (AlphaPatchF[faceId] + SMALL_NUMBER.value())/(1.0 - AlphaPatchF[faceId] + SMALL_NUMBER.value());
//            arg = Foam::max(0.0, arg);
//            PsiPatchF[faceId] = epsh.value()*( Foam::log( arg ) );
//        }
//    }
//}

// // map Alpha to Psi in centeres of CV's
//void AlphaToPsi(const volScalarField & Alpha, volScalarField & Psi, dimensionedScalar epsh )
//{
//    double arg = 0;
//    forAll(Alpha, CellId)
//    {
//        arg = (Alpha[CellId] + SMALL_NUMBER.value())/(1.0 - Alpha[CellId] + SMALL_NUMBER.value());
//        arg = Foam::max(0.0, arg) + 1.0e-32;
//        Psi[CellId] = epsh.value()*( Foam::log( arg ) );
//    }

//    forAll(Alpha.mesh().boundary(), patchi)
//    {
//        // get patch reference from (Alpha).mesh volume field
//        const fvPatch& patch = Alpha.mesh().boundary()[patchi];

//        // get reference to PsiPatch from Psi boundary field
//        fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

//        // get reference to AlphaPatach from Alpha.boundary field
//        const fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
//        forAll(patch, faceId) // lopp over faces centers belonging to a given patch
//        {
//            arg = (AlphaPatch[faceId] + SMALL_NUMBER.value())/(1.0 - AlphaPatch[faceId] + SMALL_NUMBER.value());
//            arg = Foam::max(0.0, arg) + 1e-32;
//            PsiPatch[faceId] = epsh.value()*( Foam::log( arg ) );
//        }
//    }
//}

//void snGradPsiFun(const volScalarField & Psi, surfaceScalarField & snGradPsi)
//{
////    const unallocLabelList& owner = Psi.mesh().owner();
////    const unallocLabelList& neighbour = Psi.mesh().neighbour();
////    const vectorField& Sf = Psi.mesh().Sf();

////    forAll(owner, facei)
////    {
////        double d = Foam::sqrt( Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].x() - Psi.mesh().cellCentres()[ neighbour[facei] ].x()) + Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].y() - Psi.mesh().cellCentres()[ neighbour[facei] ].y()) + Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].z() - Psi.mesh().cellCentres()[ neighbour[facei] ].z()));
////        snGradPsi[facei] = (Psi[neighbour[facei]] - Psi[owner[facei]])/ d;               // owner zwraca numer komorki wlasciciela jezeli chce dostac wartosc np. wsp X tej komurki to musze ja wyciagnac z meszu
////    }

//    snGradPsi = fvc::snGrad(Psi);
//}

//void limitPhiR(surfaceScalarField & phiR, surfaceScalarField & mGradPsiFace)
//{
//        forAll(phiR, ip)
//        {
//            if (Foam::mag(mGradPsiFace[ip]) == 0.)
//            {
//                phiR[ip] = 0.;
//            }
//        }
//}


//void setBoundaryValues(volVectorField & f, const vector & value){
//    forAll(f.boundaryField(), bid)
//    {
//        f.boundaryField()[bid] = value;
//    }
//}

//template<typename T>
//void evaluateBoundaryConditions(GeometricField<T, fvPatchField, volMesh> & field){

//    forAll(field.boundaryField(), fid)
//    {

//        if( field.boundaryField()[fid].type() == "fixedGradient" )
//        {
//            fixedGradientFvPatchField<T>& bf = (fixedGradientFvPatchField<T>&)field.boundaryField()[fid];
//            bf == (bf.patchInternalField() + bf.gradient() / bf.patch().deltaCoeffs())();
//        }
//        else if( field.boundaryField()[fid].type() == "zeroGradient" )
//        {
//            zeroGradientFvPatchField<T>& bf = (zeroGradientFvPatchField<T>&)field.boundaryField()[fid];
//            bf == bf.patchInternalField();
//        }
//        else
//            field.boundaryField()[fid].evaluate();
//    }
//}



//int main(int argc, char *argv[])
//{
//#   include "setRootCase.H"
//#   include "createTime.H"
//#   include "createMesh.H"
////    simpleControl simple(mesh);
//#   include "createFields.H"

//    std::fstream f, fconv, fB, f0, fB0;
//    f.open("./fFields.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
//    fconv.open("./fConv.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
//    fB.open("./valuesBond.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
//    f0.open("./values0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
//    fB0.open("./valuesBond0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);

//// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USTAWIENIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    int Ntau = runTime.controlDict().lookupOrDefault("Ntau", 0);
//    int PsiCritic = runTime.controlDict().lookupOrDefault("PsiCritc", 37);

//    scalar dx = scalar(1.0)/scalar(128.0);
//    dimensionedScalar epsh( "epsh", dimLength, dx / scalar(4.0) );
//    dimensionedScalar dtau( "dtau", dimTime, epsh.value() / scalar(1.0) );
//    dimensionedScalar dtauDimmless( "dtauDimmless", dimless, epsh.value() / scalar(1.0) );
//    Info << "End\n" << endl;
//    psiInit(PsiZero);
////    PsiZero.boundaryField().evaluate();
//    evaluateBoundaryConditions(PsiZero);

//    Psi == PsiZero;


////    PsiZero.correctBoundaryConditions();
//    AlphaAnalit == scalar(1.0)/(scalar(1.0)+Foam::exp(-PsiZero/epsh));
////    AlphaAnalit =  scalar(0.5)*(scalar(1.0)+Foam::tanh(PsiZero/(scalar(2.)*epsh)));


////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INICJALIZACJA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////    forAll(mesh.cellCentres(), cellI )
////    {
////        Alpha[cellI] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero[cellI]/(scalar(2.0)*epsh.value())) );  // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//////        Alpha[cellI] = 1.0/(1.0+Foam::exp(-PsiZero[cellI]/(2.0*epsh.value())));  // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////    }
////    forAll(mesh.boundary(), patchi )
////    {
////         const fvPatch& patch = mesh.boundary()[patchi];
////         fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
////         fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];
////         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
////         {
////            AlphaPatch[faceId] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZeroPatch[faceId]/(scalar(2.0)*epsh.value())) );   // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//////            Alpha[faceId] = 1.0/(1.0+Foam::exp(-PsiZero[faceId]/(2.0*epsh.value())));   // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////         }
////    }
//    Alpha == scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero/(scalar(2.0)*epsh)) );
//    AlphaZero == Alpha;
//    snGradPsiFun(Psi, snGradPsi);
////    mGradPsiFace = Foam::mag(snGradPsi);
////    Info << mGradPsiFace << endl;
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ZAPIS  PO INICJALIZACJI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//    forAll(mesh.cellCentres(), cellI )
//    {
//    f0 << mesh.cellCentres()[cellI].x()
//      << " " << mesh.cellCentres()[cellI].y()
//      << " " << std::setprecision(15) << Alpha[cellI]
//      << " " << std::setprecision(15) << PsiZero[cellI]
//      << " " << std::setprecision(15) << Psi[cellI]
//      << " " << std::setprecision(15) << AlphaAnalit[cellI]
//      << std::endl;
//    }

//    forAll(mesh.boundary(), patchi )
//    {
//        const fvPatch& patch = mesh.boundary()[patchi];

//        if( patch.type() == "empty" )
//        {
//            continue;
//        }

//        fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
//        const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
//        fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];

//        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//        {
//            fB0 << patch.Cf()[faceId].x()
//                << " " << patch.Cf()[faceId].y()
//                << " " << std::setprecision(15) << AlphaPatch[faceId]
//                << " " << std::setprecision(15) << PZPatch[faceId]
//                << " " << std::setprecision(15) << PPatch[faceId]
//                << " " << std::setprecision(15) << AlphaAnalit[faceId]
//                << std::endl;
//        }
//    }
//    fB0.close();

////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALKOWANE W PSEUDO CZASIE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    for ( int i = 0; i < Ntau; i++ )
//    {
//        Info<< "reinitialization 3thO iter: " << i << endl;

//        AlphaOld == Alpha;
//        /*--------------------------------------------------- 1-RK-step --------------------------------------------------------*/
//        AlphaFace == linearInterpolate(Alpha);
//        AlphaToPsi(Alpha, Psi, epsh);
////        evaluateBoundaryConditions(Psi);

//        gradPsi == fvc::grad(Psi);
//        mGradPsi == Foam::mag(gradPsi);

////        setBoundaryValues(gradPsi, vector(1, 0, 0));

//        //gradPsi.boundaryField().evaluate();
//        evaluateBoundaryConditions(gradPsi);
//        GradPsiFace == linearInterpolate( gradPsi );
//        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);


////        mGradPsiFace == linearInterpolate( mGradPsi );
//        mGradPsiFace == Foam::mag(GradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

//        snGradPsiFun(Psi, snGradPsi);
//        mGradPsiFaceCondition = Foam::mag(snGradPsi);
//        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
//        limitPhiR(phiR, mGradPsiFace);  // opcjonalnie
//        k1 == Alpha + fvc::div(phiR)*dtau;
//        k1 == Foam::min(Foam::max(k1, scalar(0.0)), scalar(1.0) );
////        evaluateBoundaryConditions(k1);
//        k1f == linearInterpolate(k1);

//        /*--------------------------------------------------- 2-RK-step --------------------------------------------------------*/
//        AlphaToPsi(k1, Psi, epsh);
////        evaluateBoundaryConditions(Psi);

//        mGradPsi == Foam::mag(fvc::grad(Psi));
//        gradPsi == fvc::grad(Psi);
////        setBoundaryValues(gradPsi, vector(1, 0, 0));
////        gradPsi.boundaryField().evaluate();
//        evaluateBoundaryConditions(gradPsi);
//        GradPsiFace == linearInterpolate( gradPsi );
//        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);

////        mGradPsiFace == linearInterpolate( mGradPsi );
//        mGradPsiFace == Foam::mag(GradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

//        snGradPsiFun(Psi, snGradPsi);           //  liczy grad dyfuzyjny miedzy dwoma sr komorek i przypisuje do pola snGradpsi

//        mGradPsiFaceCondition = Foam::mag(snGradPsi);
////        mGradPsiFaceCondition.boundaryField().evaluate();
//        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
//        limitPhiR(phiR, mGradPsiFace);   // opcjonalnie
//        k2 == scalar(0.75)* Alpha + scalar(0.25)* k1 + scalar(0.25)* fvc::div(phiR)*dtau;
//        k2 == Foam::min(Foam::max(k2, scalar(0.0)), scalar(1.0) );
////        evaluateBoundaryConditions(k2);
//        k2f == linearInterpolate(k2);

//        /*--------------------------------------------------- 3-RK-step --------------------------------------------------------*/
//        AlphaToPsi(k2, Psi, epsh);
////        evaluateBoundaryConditions(Psi);

//        mGradPsi == Foam::mag(fvc::grad(Psi));
//        gradPsi == fvc::grad(Psi);
////        setBoundaryValues(gradPsi, vector(1, 0, 0));
////        gradPsi.boundaryField().evaluate();
//        evaluateBoundaryConditions(gradPsi);
//        GradPsiFace == linearInterpolate( gradPsi );
//        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);

//        //mGradPsiFace == linearInterpolate( mGradPsi );
//        mGradPsiFace == Foam::mag(GradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

//        snGradPsiFun(Psi, snGradPsi);
//        mGradPsiFaceCondition = Foam::mag(snGradPsi);
//        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
//        limitPhiR(phiR, mGradPsiFace);    // opcjonalnie
//        Alpha ==  scalar(100.0)/scalar(300.0)* Alpha + scalar(200.0)/scalar(300.0)* k2 + scalar(200.0)/scalar(300.0)*fvc::div(phiR)*dtau;
//        Alpha == Foam::min(Foam::max(Alpha, scalar(0.0)), scalar(1.0) );
////        evaluateBoundaryConditions(Alpha);

//        double norm1c = Foam::sum(Foam::mag(Alpha-AlphaOld)).value() / Alpha.size();
//        Info << "Norma conv = " << norm1c << endl;
//        double norm1 = Foam::sum(Foam::mag(AlphaAnalit-Alpha)).value() / Alpha.size();
//        Info << "Norma 1 analit = " << norm1 << endl;
//        fconv << i << "   " << norm1 << std::endl;
//    }

////    runTime.write();
////    Alpha.write();
//    std::cout << std::endl;
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ZAPIS PO OBLICZENIACH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//     //   AlphaToPsi(Alpha,Psi,epsh);
//     //  AlphaToPsi(AlphaZero,PsiZero,epsh);
//        mGradPsiLimitedPlot = Foam::mag(fvc::grad(Psi));
//        forAll(mesh.cellCentres(), cellI )
//        {
//        f << mesh.cellCentres()[cellI].x()
//          << " " << mesh.cellCentres()[cellI].y()
//          << " " << std::setprecision(32) << AlphaZero[cellI]
//          << " " << std::setprecision(32) << Alpha[cellI]
//          << " " << std::setprecision(32) << PsiZero[cellI]
//          << " " << std::setprecision(32) << Psi[cellI]
////          << " " << std::setprecision(32) << gradPsi.x[cellI] << std::endl;
////          << " " << std::setprecision(32) << (gradPsi[cellI].x()) << std::endl;
//          << " " << std::setprecision(32) << mGradPsi[cellI]
//          << " " << std::setprecision(32) << mGradPsiLimitedPlot[cellI]
//          << " " << std::setprecision(32) << phiR[cellI]
//          << std::endl;
//        }

//        forAll(mesh.boundary(), patchi )
//        {
//            const fvPatch& patch = mesh.boundary()[patchi];
//            if( patch.type() == "empty" )
//            { continue; }
//            const fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
//            const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
//            const fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];
//            const fvPatchScalarField& mGradPP = mGradPsi.boundaryField()[patchi];

//            forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//            {
//                fB << patch.Cf()[faceId].x()
//                    << " " << patch.Cf()[faceId].y()
//                    << " " << std::setprecision(32) << AlphaPatch[faceId]
//                    << " " << std::setprecision(32) << PZPatch[faceId]
//                    << " " << std::setprecision(32) << PPatch[faceId]
//                    << " " << std::setprecision(32) << mGradPP[faceId]
//                    << std::endl;
//            }
//        }
//    f.close();
//    fB.close();
//    f0.close();
//    fconv.close();
//    return 0;
//}

// ************************************************************************* //

