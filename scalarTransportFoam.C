/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "simpleControl.H"
//#include "surfaceInterpolationScheme.H"
//#include "linear.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensionedScalar SMALL_NUMBER("small", dimless, 5e-16);

//!!!!!!!!!!!!!!!!!!!
void psiInit(volScalarField& Psi)
{
    forAll(Psi.mesh().cellCentres(), cellI) // loop over cell centers
    {
        Psi[cellI] = Psi.mesh().cellCentres()[cellI].x();           //-0.5;

    //    Info << cellI << endl;
    }
    forAll(Psi.mesh().boundary(), patchi)  // loop over all boundary patches
    {
        const fvPatch& patch = Psi.mesh().boundary()[patchi];
        fvPatchScalarField& psiPatch = Psi.boundaryField()[patchi];
        forAll(patch, faceId)       //  loop over c.v. centers of the given patch
        {
            psiPatch[faceId] = patch.Cf()[faceId].x();          //-0.5;
        }
    }
          //  Info << Psi << endl;
}


void LimitGradPsi2(const fvMesh& mesh, const volScalarField & Psi, volScalarField & mGradPsi, double epsh, double PsiCritic)
{
    forAll( mesh.cellCentres(), cellId)
    {
        if( mag(Psi[cellId]) > (PsiCritic+3.0)*epsh )
        {
             mGradPsi[cellId] = scalar(1.0);
        }
    }

    forAll(mesh.boundary(), patchi)  // mesh.boundary() returns addresses to b.c.
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        fvPatchScalarField& gradPsiPatch = mGradPsi.boundaryField()[patchi];
        const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

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
void AlphaToPsiF(const surfaceScalarField & AlphaF, surfaceScalarField & PsiF, dimensionedScalar epsh )
{
    double dz = 1.;
    double arg = 0;
    forAll(AlphaF, CellId)
    {
        arg = ((AlphaF[CellId] + SMALL_NUMBER.value()/dz)/(1.0 - AlphaF[CellId] + SMALL_NUMBER.value()/dz));
        arg = Foam::max(0.0, arg);
        PsiF[CellId] = epsh.value()*( Foam::log( arg ) );
    }

    forAll(AlphaF.mesh().boundary(), patchi)
    {
        const fvPatch& patch = AlphaF.mesh().boundary()[patchi];
        fvsPatchScalarField& PsiPatchF = PsiF.boundaryField()[patchi];
        const fvsPatchScalarField& AlphaPatchF = AlphaF.boundaryField()[patchi];
        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
            arg = ((AlphaPatchF[faceId] + SMALL_NUMBER.value()/dz)/(1.0 - AlphaPatchF[faceId] + SMALL_NUMBER.value()/dz));
            arg = Foam::max(0.0, arg);
            PsiPatchF[faceId] = epsh.value()*( Foam::log( arg ) );
        }
    }
}


// map Alpha to Psi in centeres of CV's
void AlphaToPsi(const volScalarField & Alpha, volScalarField & Psi, dimensionedScalar epsh )
{
    double dz = 1.;
    double arg = 0;
    forAll(Alpha, CellId)
    {
        arg = ((Alpha[CellId] + SMALL_NUMBER.value()/dz)/(1.0 - Alpha[CellId] + SMALL_NUMBER.value()/dz));
        arg = Foam::max(0.0, arg);// + 1.0e-32;
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
            arg = ((AlphaPatch[faceId] + SMALL_NUMBER.value()/dz)/(1.0 - AlphaPatch[faceId] + SMALL_NUMBER.value()/dz));
            arg = Foam::max(0.0, arg);//+ 1e-32;
            PsiPatch[faceId] = epsh.value()*( Foam::log( arg ) );
        }
    }
}

//void AlphaToPsi3(const volScalarField& T, volScalarField& Psi, const dimensionedScalar& epsH)
//{
//    double small = 1.;
//    double mult = 10.;
//     forAll(T, CellId)
//     {
//         if(T[CellId] > SMALL_NUMBER.value()*small && T[CellId] < (1.0 - SMALL_NUMBER.value()*small  ) )
//         {
//            Psi[CellId] = epsH.value()*( Foam::log(T[CellId]) - Foam::log(1.0 - T[CellId])) ;
//         }
//         else
//         {
//            Psi[CellId] = epsH.value()*( Foam::log(T[CellId] +SMALL_NUMBER.value()*mult) - Foam::log(1.0 - T[CellId] +SMALL_NUMBER.value()*mult)) ;
//         }
//     }
//     forAll(T.mesh().boundary(), patchi)
//     {
//         const fvPatch& patch = T.mesh().boundary()[patchi];
//         fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];
//         const fvPatchScalarField& TPatch = T.boundaryField()[patchi];
//         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
//         {
//             if(TPatch[faceId] > SMALL_NUMBER.value()*small  && TPatch[faceId] < (1.0 - SMALL_NUMBER.value()*small ) )
//             {
//                 PsiPatch[faceId] = epsH.value()*( Foam::log(TPatch[faceId]) - Foam::log(1.0 - TPatch[faceId]));
//             }
//             else
//             {
//                 PsiPatch[faceId] = epsH.value()*( Foam::log(TPatch[faceId]+SMALL_NUMBER.value()*mult) - Foam::log(1.0 - TPatch[faceId]+SMALL_NUMBER.value()*mult));
//             }
//         }
//     }
//}

void snGradPsiFun(const volScalarField & Psi, surfaceScalarField & snGradPsi)
{
    const unallocLabelList& owner = Psi.mesh().owner();
    const unallocLabelList& neighbour = Psi.mesh().neighbour();
    const vectorField& Sf = Psi.mesh().Sf();

    forAll(owner, facei)
    {
        double d = Foam::sqrt( Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].x() - Psi.mesh().cellCentres()[ neighbour[facei] ].x()) + Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].y() - Psi.mesh().cellCentres()[ neighbour[facei] ].y()) + Foam::sqr(Psi.mesh().cellCentres()[ owner[facei] ].z() - Psi.mesh().cellCentres()[ neighbour[facei] ].z()));
        snGradPsi[facei] = (Psi[neighbour[facei]] - Psi[owner[facei]])/ d;               // owner zwraca numer komorki wlasciciela jezeli chce dostac wartosc np. wsp X tej komurki to musze ja wyciagnac z meszu
    }
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
    fB.open("./valuesBond.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f0.open("./values0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB0.open("./valuesBond0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USTAWIENIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int Ntau = runTime.controlDict().lookupOrDefault("Ntau", 0);
    int limit = runTime.controlDict().lookupOrDefault("limit", 0.0);

    int PsiCritic = runTime.controlDict().lookupOrDefault("PsiCritc", 37);

    scalar dx = scalar(1.0)/scalar(128.0);
    dimensionedScalar epsh( "epsh", dimLength, dx / scalar(4.0) );
    dimensionedScalar dtau( "dtau", dimTime, epsh.value() / scalar(1.0) );
    dimensionedScalar dtauDimmless( "dtauDimmless", dimless, epsh.value() / scalar(1.0) );
    Info << "End\n" << endl;
    psiInit(PsiZero);
    Psi == PsiZero;
    AlphaAnalit = scalar(1.0)/(scalar(1.0)+Foam::exp(-PsiZero/epsh));
//    AlphaAnalit =  scalar(0.5)*(scalar(1.0)+Foam::tanh(PsiZero/(scalar(2.)*epsh)));

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INICJALIZACJA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    forAll(mesh.cellCentres(), cellI )
    {
        Alpha[cellI] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero[cellI]/(scalar(2.0)*epsh.value())) );  // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//        Alpha[cellI] = 1.0/(1.0+Foam::exp(-PsiZero[cellI]/(2.0*epsh.value())));  // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    forAll(mesh.boundary(), patchi )
    {
         const fvPatch& patch = mesh.boundary()[patchi];
         fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
         fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];
         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
         {
            AlphaPatch[faceId] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZeroPatch[faceId]/(scalar(2.0)*epsh.value())) );   // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//            Alpha[faceId] = 1.0/(1.0+Foam::exp(-PsiZero[faceId]/(2.0*epsh.value())));   // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         }
    }
    AlphaZero == Alpha;
    snGradPsiFun(Psi, snGradPsi);
    mGradPsiFace = Foam::mag(snGradPsi);
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
      << " " << std::setprecision(15) << mGradPsiFace[cellI]
      << std::endl;
    }

    forAll(mesh.boundary(), patchi )
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        if( patch.type() == "empty" )
        { continue; }
        fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
        const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
        fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];

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
    double psiLimiter = Foam::max(Psi).value() -1.5*dx;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALKOWANE W PSEUDO CZASIE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for ( int i = 0; i < Ntau; i++ )
    {
        Info<< "reinitialization 3thO iter: " << i << endl;

        /*--------------------------------------------------- 1-RK-step --------------------------------------------------------*/
        AlphaOld == Alpha;
        AlphaFace == linearInterpolate(Alpha);

        AlphaToPsi(Alpha, Psi, epsh);
        mGradPsi == Foam::mag(fvc::grad(Psi));

        gradPsi == fvc::grad(Psi);
        GradPsiFace == linearInterpolate( gradPsi );
        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);
        mGradPsiFace == linearInterpolate( mGradPsi );
//        mGradPsiFace == Foam::mag(fvc::snGrad(Psi));
        snGradPsiFun(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();

        limitPhiR(phiR, mGradPsiFace);  // opcjonalnie

        k1 == Alpha + fvc::div(phiR)*dtau;
        k1 == Foam::min(Foam::max(k1, scalar(0.0)), scalar(1.0) );
        k1f == linearInterpolate(k1);

        /*--------------------------------------------------- 2-RK-step --------------------------------------------------------*/
        AlphaToPsi(k1, Psi, epsh);
        mGradPsi == Foam::mag(fvc::grad(Psi));

        gradPsi == fvc::grad(Psi);
        GradPsiFace == linearInterpolate( gradPsi );
        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);
        mGradPsiFace == linearInterpolate( mGradPsi );
//        mGradPsiFace == Foam::mag(fvc::snGrad(Psi));
        snGradPsiFun(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();

        limitPhiR(phiR, mGradPsiFace);

        k2 == scalar(0.75)* Alpha + scalar(0.25)* k1 + scalar(0.25)* fvc::div(phiR)*dtau;
        k2 == Foam::min(Foam::max(k2, scalar(0.0)), scalar(1.0) );
        k2f == linearInterpolate(k2);

        /*--------------------------------------------------- 3-RK-step --------------------------------------------------------*/
        AlphaToPsi(k2, Psi, epsh);
        mGradPsi == Foam::mag(fvc::grad(Psi));

        gradPsi == fvc::grad(Psi);
        GradPsiFace == linearInterpolate( gradPsi );
        LimitGradPsi2(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);
        mGradPsiFace == linearInterpolate( mGradPsi );
//        mGradPsiFace == Foam::mag(fvc::snGrad(Psi));
        snGradPsiFun(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();

        limitPhiR(phiR, mGradPsiFace);

        Alpha ==  scalar(100.0)/scalar(300.0)* Alpha + scalar(200.0)/scalar(300.0)* k2 + scalar(200.0)/scalar(300.0)*fvc::div(phiR)*dtau;
        Alpha == Foam::min(Foam::max(Alpha, scalar(0.0)), scalar(1.0) );

        double norm1c = Foam::sum(Foam::mag(Alpha-AlphaOld)).value() / Alpha.size();
        Info << "Norma conv = " << norm1c << endl;

        double norm1 = Foam::sum(Foam::mag(AlphaAnalit-Alpha)).value() / Alpha.size();
        Info << "Norma 1 analit = " << norm1 << endl;
        fconv << i << "   " << norm1 << std::endl;

    }

    std::cout << std::endl;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ZAPIS PO OBLICZENIACH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        AlphaToPsi(Alpha,Psi,epsh);
        AlphaToPsi(AlphaZero,PsiZero,epsh);
        mGradPsiLimitedPlot = Foam::mag(fvc::grad(Psi));
        forAll(mesh.cellCentres(), cellI )
        {
        f << mesh.cellCentres()[cellI].x()
          << " " << mesh.cellCentres()[cellI].y()
          << " " << std::setprecision(32) << AlphaZero[cellI]
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
            const fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
            const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
            const fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];
            const fvPatchScalarField& mGradPP = mGradPsi.boundaryField()[patchi];

            forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
            {
                fB << patch.Cf()[faceId].x()
                    << " " << patch.Cf()[faceId].y()
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
