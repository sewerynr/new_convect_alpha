/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "simpleControl.H"
//#include "surfaceInterpolationScheme.H"
//#include "linear.H"
#include "labelList.H"
#include "fixedGradientFvPatchField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensionedScalar SMALL_NUMBER("small", dimless, 5e-16);

//!!!!!!!!!!!!!!!!!!!
void psiInit(volScalarField& Psi)
{
    forAll(Psi.mesh().cellCentres(), cellI) // loop over cell centers
    {
        Psi[cellI] = -Psi.mesh().cellCentres()[cellI].x() + 0.5;
    }
//    forAll(Psi.mesh().boundary(), patchi)  // loop over all boundary patches
//    {
//        const fvPatch& patch = Psi.mesh().boundary()[patchi];
//        fvPatchScalarField& psiPatch = Psi.boundaryField()[patchi];
//        forAll(patch, faceId)       //  loop over c.v. centers of the given patch
//        {
//            psiPatch[faceId] = patch.Cf()[faceId].x();          //-0.5;
//        }
//    }
}

//        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);

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

void LimitMGradPsiFaceCondition(const fvMesh& mesh, const surfaceScalarField & PsiF, surfaceScalarField&  mGradPsiFaceCondition, double epsh, double PsiCritic)
{
    forAll( PsiF, FaceId)
    {
        if( mag(PsiF[FaceId]) > (PsiCritic+3.0)*epsh )
        {
             mGradPsiFaceCondition[FaceId] = scalar(1.0);
        }
    }
    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvsPatchScalarField& gradPsiPatch = mGradPsiFaceCondition.boundaryField()[patchi];
        const fvsPatchScalarField& PsiPatch = PsiF.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(PsiPatch[faceId]) > (PsiCritic+3.0)*epsh )
            {
                gradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }
}

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

// map Alpha to Psi in centeres of CV's
void AlphaToPsi(const volScalarField & Alpha, volScalarField & Psi, dimensionedScalar epsh )
{
    double arg = 0;
    forAll(Alpha, CellId)
    {
        arg = (Alpha[CellId] + SMALL_NUMBER.value())/(1.0 - Alpha[CellId] + SMALL_NUMBER.value());
        arg = Foam::max(0.0, arg) + 1.0e-32;
        Psi[CellId] = epsh.value()*( Foam::log( arg ) );
    }

    forAll(Alpha.mesh().boundary(), patchi)
    {
        // get patch reference from (Alpha).mesh volume field
        const fvPatch& patch = Alpha.mesh().boundary()[patchi];

        // get reference to PsiPatch from Psi boundary field
        fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

        // get reference to AlphaPatach from Alpha.boundary field
        const fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
        forAll(patch, faceId) // lopp over faces centers belonging to a given patch
        {
            arg = (AlphaPatch[faceId] + SMALL_NUMBER.value())/(1.0 - AlphaPatch[faceId] + SMALL_NUMBER.value());
            arg = Foam::max(0.0, arg) + 1e-32;
            PsiPatch[faceId] = epsh.value()*( Foam::log( arg ) );
        }
    }
}

void snGradPsiFun(const volScalarField & Psi, surfaceScalarField & snGradPsi)
{
//    const unallocLabelList& owner = Psi.mesh().owner();
//    const unallocLabelList& neighbour = Psi.mesh().neighbour();
//    const vectorField& Sf = Psi.mesh().Sf();

//    forAll(owner, facei)
//    {
//        double d = Foam::sqrt( Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].x() - Psi.mesh().cellCentres()[ neighbour[facei] ].x()) + Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].y() - Psi.mesh().cellCentres()[ neighbour[facei] ].y()) + Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].z() - Psi.mesh().cellCentres()[ neighbour[facei] ].z()));
//        snGradPsi[facei] = (Psi[neighbour[facei]] - Psi[owner[facei]])/ d;               // owner zwraca numer komorki wlasciciela jezeli chce dostac wartosc np. wsp X tej komurki to musze ja wyciagnac z meszu
//    }

    snGradPsi = fvc::snGrad(Psi);
}

void limitPhiR(surfaceScalarField & phiR, surfaceScalarField & mGradPsiFace)
{
        forAll(phiR, ip)
        {
            if (Foam::mag(mGradPsiFace[ip]) == 0.)
            {
                phiR[ip] = 0.;
            }
        }
}


void setBoundaryValues(volVectorField & f, const vector & value)
{
    forAll(f.boundaryField(), bid)
    {
        f.boundaryField()[bid] = value;
    }
}

template<typename T>
void evaluateBoundaryConditions(GeometricField<T, fvPatchField, volMesh> & field)
{
    forAll(field.boundaryField(), fid)
    {

        if( field.boundaryField()[fid].type() == "fixedGradient" )
        {
            fixedGradientFvPatchField<T>& bf = (fixedGradientFvPatchField<T>&)field.boundaryField()[fid];
            bf == (bf.patchInternalField() + bf.gradient() / bf.patch().deltaCoeffs())();
        }
        else if( field.boundaryField()[fid].type() == "zeroGradient" )
        {
            zeroGradientFvPatchField<T>& bf = (zeroGradientFvPatchField<T>&)field.boundaryField()[fid];
            bf == bf.patchInternalField();
        }
        else
            field.boundaryField()[fid].evaluate();
    }
}

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
//    simpleControl simple(mesh);
#   include "createFields.H"

    std::fstream f, fconv, fB, f0, fB0;
    f.open("./fFields.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fconv.open("./fConv.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB.open("./valuesB.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f0.open("./values0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB0.open("./valuesB0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USTAWIENIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int Ntau = runTime.controlDict().lookupOrDefault("Ntau", 0);
    int PsiCritic = runTime.controlDict().lookupOrDefault("PsiCritic", 30);

    scalar dx = scalar(1.0)/scalar(129.0);
    dimensionedScalar epsh( "epsh", dimLength, dx / scalar(4.0) );
    dimensionedScalar dtau( "dtau", dimTime, epsh.value() / scalar(1.0) );
    dimensionedScalar dtauDimmless( "dtauDimmless", dimless, epsh.value() / scalar(1.0) );
    Info << "End\n" << endl;
    psiInit(PsiZero);
//    PsiZero.boundaryField().evaluate();
    evaluateBoundaryConditions(PsiZero);

    Psi == PsiZero;

//    PsiZero.correctBoundaryConditions();
//    AlphaAnalit == scalar(1.0)/(scalar(1.0)+Foam::exp(-PsiZero/epsh));
    AlphaAnalit == scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero/(scalar(2.0)*epsh)) );

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INICJALIZACJA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    forAll(mesh.cellCentres(), cellI )
//    {
//        Alpha[cellI] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero[cellI]/(scalar(2.0)*epsh.value())) );  // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////        Alpha[cellI] = 1.0/(1.0+Foam::exp(-PsiZero[cellI]/(2.0*epsh.value())));  // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    }
//    forAll(mesh.boundary(), patchi )
//    {
//         const fvPatch& patch = mesh.boundary()[patchi];
//         fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
//         fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];
//         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//         {
//            AlphaPatch[faceId] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZeroPatch[faceId]/(scalar(2.0)*epsh.value())) );   // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////            Alpha[faceId] = 1.0/(1.0+Foam::exp(-PsiZero[faceId]/(2.0*epsh.value())));   // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//         }
//    }
//    Alpha == AlphaAnalit;
    Alpha == scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero/(scalar(2.0)*epsh)) );
    AlphaZero == Alpha;
    snGradPsiFun(Psi, snGradPsi);
//    mGradPsiFace = Foam::mag(snGradPsi);
//    Info << mGradPsiFace << endl;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ZAPIS  PO INICJALIZACJI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    forAll(mesh.cellCentres(), cellI )
    {
    f0 << mesh.cellCentres()[cellI].x()
      << " " << mesh.cellCentres()[cellI].y()
      << " " << std::setprecision(15) << Alpha[cellI]
      << " " << std::setprecision(15) << PsiZero[cellI]
      << " " << std::setprecision(15) << Psi[cellI]
      << " " << std::setprecision(15) << AlphaAnalit[cellI]
      << std::endl;
    }

    forAll(mesh.boundary(), patchi )
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if( patch.type() == "empty" )
        {
            continue;
        }

        const fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
        const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
        const fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];

        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
            fB0 << patch.Cf()[faceId].x()
                << " " << patch.Cf()[faceId].y()
                << " " << std::setprecision(15) << AlphaPatch[faceId]
                << " " << std::setprecision(15) << PZPatch[faceId]
                << " " << std::setprecision(15) << PPatch[faceId]
                << " " << std::setprecision(15) << AlphaAnalit[faceId]
                << std::endl;
        }
    }
    fB0.close();

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALKOWANE W PSEUDO CZASIE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for ( int i = 0; i < Ntau; i++ )
    {
        Info<< "reinitialization 3thO iter: " << i << endl;

        AlphaOld == Alpha;
        /*--------------------------------------------------- 1-RK-step --------------------------------------------------------*/
        AlphaFace == linearInterpolate(Alpha);
        AlphaToPsi(Alpha, Psi, epsh);
//        evaluateBoundaryConditions(Psi);
                Psi.boundaryField() == PsiZero.boundaryField();

        gradPsi == fvc::grad(Psi);
        mGradPsi == Foam::mag(gradPsi);
//        setBoundaryValues(gradPsi, vector(1, 0, 0));
//        gradPsi.boundaryField().evaluate();
//        evaluateBoundaryConditions(gradPsi);
        GradPsiFace == linearInterpolate( gradPsi );

//        mGradPsiFace == linearInterpolate( mGradPsi );
        mGradPsiFace == Foam::mag(GradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

        PsiF == linearInterpolate(Psi);
        snGradPsiFun(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        LimitMGradPsiFaceCondition(mesh, PsiF, mGradPsiFaceCondition, epsh.value(), PsiCritic);

        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
        limitPhiR(phiR, mGradPsiFace);  // opcjonalnie
        k1 == Alpha + fvc::div(phiR)*dtau;
        k1 == Foam::min(Foam::max(k1, scalar(0.0)), scalar(1.0) );
//        evaluateBoundaryConditions(k1);
                k1.boundaryField() == AlphaAnalit.boundaryField();
        k1f == linearInterpolate(k1);

        /*--------------------------------------------------- 2-RK-step --------------------------------------------------------*/
        AlphaToPsi(k1, Psi, epsh);
//        evaluateBoundaryConditions(Psi);
                Psi.boundaryField() == PsiZero.boundaryField();

        mGradPsi == Foam::mag(fvc::grad(Psi));
        gradPsi == fvc::grad(Psi);
//        setBoundaryValues(gradPsi, vector(1, 0, 0));
//        gradPsi.boundaryField().evaluate();
//        evaluateBoundaryConditions(gradPsi);
//        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);
        GradPsiFace == linearInterpolate( gradPsi );
//        mGradPsiFace == linearInterpolate( mGradPsi );
        mGradPsiFace == Foam::mag(GradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

        PsiF == linearInterpolate(Psi);
        snGradPsiFun(Psi, snGradPsi);           //  liczy grad dyfuzyjny miedzy dwoma sr komorek i przypisuje do pola snGradpsi

        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        LimitMGradPsiFaceCondition(mesh, PsiF, mGradPsiFaceCondition, epsh.value(), PsiCritic);

//        mGradPsiFaceCondition.boundaryField().evaluate();
        phiR = Cf * k1f*( scalar(1.0) - k1f ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
        limitPhiR(phiR, mGradPsiFace);   // opcjonalnie
        k2 == scalar(0.75)* Alpha + scalar(0.25)* k1 + scalar(0.25)* fvc::div(phiR)*dtau;
        k2 == Foam::min(Foam::max(k2, scalar(0.0)), scalar(1.0) );
//        evaluateBoundaryConditions(k2);
                k2.boundaryField() == AlphaAnalit.boundaryField();

        k2f == linearInterpolate(k2);

        /*--------------------------------------------------- 3-RK-step --------------------------------------------------------*/
        AlphaToPsi(k2, Psi, epsh);
//        evaluateBoundaryConditions(Psi);
                Psi.boundaryField() == PsiZero.boundaryField();

        mGradPsi == Foam::mag(fvc::grad(Psi));
        gradPsi == fvc::grad(Psi);
//        setBoundaryValues(gradPsi, vector(1, 0, 0));
//        gradPsi.boundaryField().evaluate();
//        evaluateBoundaryConditions(gradPsi);
        GradPsiFace == linearInterpolate( gradPsi );

        //mGradPsiFace == linearInterpolate( mGradPsi );
        mGradPsiFace == Foam::mag(GradPsiFace); //mag z interpolowanego, a nie interpolacja z mag

        PsiF == linearInterpolate(Psi);
        snGradPsiFun(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        LimitMGradPsiFaceCondition(mesh, PsiF, mGradPsiFaceCondition, epsh.value(), PsiCritic);

        phiR = Cf * k2f*( scalar(1.0) - k2f ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
        vector(1.,0.,0.);
        limitPhiR(phiR, mGradPsiFace);    // opcjonalnie
        Alpha ==  scalar(100.0)/scalar(300.0)* Alpha + scalar(200.0)/scalar(300.0)* k2 + scalar(200.0)/scalar(300.0)*fvc::div(phiR)*dtau;
        Alpha == Foam::min(Foam::max(Alpha, scalar(0.0)), scalar(1.0) );
//        evaluateBoundaryConditions(Alpha);

        double norm1c = Foam::sum(Foam::mag(Alpha-AlphaOld)).value() / Alpha.size();
        Info << "Norma conv = " << norm1c << endl;
        double norm1 = Foam::sum(Foam::mag(AlphaAnalit-Alpha)).value() / Alpha.size();
        Info << "Norma 1 analit = " << norm1 << endl;
        fconv << i << "   " << norm1 << std::endl;

        Alpha.boundaryField() == AlphaAnalit.boundaryField();
        Psi.boundaryField() == PsiZero.boundaryField();
    }

//    runTime.write();
//    Alpha.write();
    std::cout << std::endl;
//    Info <<  Psi.mesh().size() << endl;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ZAPIS PO OBLICZENIACH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     //   AlphaToPsi(Alpha,Psi,epsh);
     //  AlphaToPsi(AlphaZero,PsiZero,epsh);
        mGradPsiLimitedPlot = Foam::mag(fvc::grad(Psi));
        forAll(mesh.cellCentres(), cellI )
        {
        f << mesh.cellCentres()[cellI].x()
          << " " << mesh.cellCentres()[cellI].y()
          << " " << std::setprecision(32) << AlphaZero[cellI]
          << " " << std::setprecision(32) << AlphaAnalit[cellI]
          << " " << std::setprecision(32) << Alpha[cellI]
          << " " << std::setprecision(32) << PsiZero[cellI]
          << " " << std::setprecision(32) << Psi[cellI]
//          << " " << std::setprecision(32) << gradPsi.x[cellI] << std::endl;
//          << " " << std::setprecision(32) << (gradPsi[cellI].x()) << std::endl;
          << " " << std::setprecision(32) << mGradPsi[cellI]
          << " " << std::setprecision(32) << mGradPsiLimitedPlot[cellI]
          << " " << std::setprecision(32) << phiR[cellI]
          << std::endl;
        }

        forAll(mesh.boundary(), patchi )
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if( patch.type() == "empty" )
            { continue; }
            const fvPatchScalarField& AlphaPatchZ = AlphaZero.boundaryField()[patchi];
            const fvPatchScalarField& AlphaPatchAn = AlphaAnalit.boundaryField()[patchi];
            const fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
            const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
            const fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];
            const fvPatchScalarField& mGradPP = mGradPsi.boundaryField()[patchi];

            forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
            {
                fB << patch.Cf()[faceId].x()
                    << " " << patch.Cf()[faceId].y()
                    << " " << std::setprecision(32) << AlphaPatchZ[faceId]
                    << " " << std::setprecision(32) << AlphaPatchAn[faceId]
                    << " " << std::setprecision(32) << AlphaPatch[faceId]
                    << " " << std::setprecision(32) << PZPatch[faceId]
                    << " " << std::setprecision(32) << PPatch[faceId]
                    << " " << std::setprecision(32) << mGradPP[faceId]
                    << std::endl;
            }
        }
    f.close();
    fB.close();
    f0.close();
    fconv.close();
    return 0;
}

// ************************************************************************* //
