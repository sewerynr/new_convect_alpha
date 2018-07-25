/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensionedScalar SMALL_NUMBER("small", dimless, SMALL);

void psiInit(volScalarField& Psi);

// TW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
//            PsiPatch[faceId] = -100;
        }
    }
//    psiInit(Psi);
}

scalar funInitT4( double PsiZero, const dimensionedScalar& epsH, double par, const double x,const double y, const double z )
{
    const static scalar R = 0.15;
    scalar eps = epsH.value();
    scalar ret = 0.0;
    scalar r = Foam::sqrt(x*x + (y-0.25)*(y-0.25));
    return 1 - 0.5*(1.0 + Foam::tanh((r-R)/2/eps));
}


//!!!!!!!!!!!!!!!!!!!
void psiInit(volScalarField& Psi)
{
    forAll(Psi.mesh().cellCentres(), cellI) // loop over cell centers
    {
        Psi[cellI] = - Psi.mesh().cellCentres()[cellI].x() + 0.5;
    }
    forAll(Psi.mesh().boundary(), patchi)  // loop over all boundary patches
    {
        const fvPatch& patch = Psi.mesh().boundary()[patchi];
        fvPatchScalarField& psiPatch = Psi.boundaryField()[patchi];
        forAll(patch, faceId)       //  loop over c.v. centers of the given patch
        {
            psiPatch[faceId] = - patch.Cf()[faceId].x() + 0.5;
        }
    }
}

void LimitGradPsi(const fvMesh& mesh, const volScalarField & Psi, volScalarField & mGradPsi, double epsh, double PsiCritic)
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
        fvPatchScalarField& mgradPsiPatch = mGradPsi.boundaryField()[patchi];
        const fvPatchScalarField& PsiPatch = Psi.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(PsiPatch[faceId]) > (PsiCritic+3.0)*epsh )
            {
                mgradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }
}

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
        fvsPatchScalarField& mGradPsiPatch = mGradPsiFaceCondition.boundaryField()[patchi];
        const fvsPatchScalarField& PsiPatch = PsiF.boundaryField()[patchi];

        forAll(patch, faceId) // petla po centrach objetosci danego patcha
        {
            if ( mag(PsiPatch[faceId]) > (PsiCritic+3.0)*epsh )
            {
                mGradPsiPatch[faceId] = scalar(1.0);
            }
        }
    }
}

void setInternalToBoundary(volVectorField & field)
{

//    field = value;

    forAll(field.boundaryField(), fid)
    {
            fvPatchVectorField& bf = field.boundaryField()[fid];
            forAll(bf.patch(), faceId)
            {
                bf[faceId] = bf.patchInternalField()()[faceId];
            }
    }
}


void snGradFunPsi(const volScalarField & field, surfaceScalarField & result)
{
    const fvMesh & mesh = field.mesh();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        label ocell = owner[facei];
        label ncell = neighbour[facei];
        scalar Vdiff = - field[ocell] + field[ncell];
        vector Cdiff = mesh.cellCentres()[ocell] - mesh.cellCentres()[ncell];
        vector nf = mesh.Sf()[facei];
        nf = nf / mag(nf);
        result[facei] = Vdiff * (Cdiff & nf) / (Cdiff & Cdiff);
    }
    forAll(field.boundaryField(), patchId)
    {
            const fvPatchScalarField& PsiPach = field.boundaryField()[patchId];
            fvsPatchScalarField& snGradPsiPach = result.boundaryField()[patchId];
            forAll(snGradPsiPach.patch(), faceId)
            {
                snGradPsiPach[faceId] = 0.;
            }
    }
// to samo ! co wyzej
//    result == fvc::snGrad(field);
}

tmp< surfaceScalarField > linearScalarInterpolate(const volScalarField & vf)
{
    const fvMesh & mesh = vf.mesh();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    tmp< surfaceScalarField > tsf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    surfaceScalarField & result = tsf();


    forAll(result, facei)
    {
        label ocell = owner[facei];
        label ncell = neighbour[facei];

        vector ofc = mesh.cellCentres()[ocell] - mesh.Cf()[facei];
        vector nfc = mesh.cellCentres()[ncell] - mesh.Cf()[facei];

        scalar oWeight = mag(ofc) / ( mag(ofc) + mag(nfc));
        scalar nWeight = 1.0 - oWeight;

        result[facei] = oWeight * vf[ocell] + nWeight * vf[ncell];
    }


    //Wymuszamy wiedząc że na brzegu jest zeroGrad
    forAll(vf.boundaryField(), patchId)
    {
        //const fvPatchScalarField
        const typename volScalarField::PatchFieldType & vfPatch = vf.boundaryField()[patchId];

        //const fvsPatchScalarField
        typename surfaceScalarField::PatchFieldType& sfPatch = result.boundaryField()[patchId];

        sfPatch = vfPatch.patchInternalField();
    }

    return tsf;
}

tmp< surfaceVectorField > linearVectorInterpolate(const volVectorField & vf)
{
    const fvMesh & mesh = vf.mesh();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    tmp< surfaceVectorField > tsf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    surfaceVectorField & result = tsf();


    forAll(result, facei)
    {
        label ocell = owner[facei];
        label ncell = neighbour[facei];

        vector ofc = mesh.cellCentres()[ocell] - mesh.Cf()[facei];
        vector nfc = mesh.cellCentres()[ncell] - mesh.Cf()[facei];

        scalar oWeight = mag(ofc) / ( mag(ofc) + mag(nfc));
        scalar nWeight = 1.0 - oWeight;

        result[facei] = oWeight * vf[ocell] + nWeight * vf[ncell];
    }


    //Wymuszamy wiedząc że na brzegu jest zeroGrad
    forAll(vf.boundaryField(), patchId)
    {

        const typename volVectorField::PatchFieldType & vfPatch = vf.boundaryField()[patchId];


        typename surfaceVectorField::PatchFieldType& sfPatch = result.boundaryField()[patchId];

        sfPatch = vfPatch.patchInternalField();
    }

    return tsf;
}

tmp<volVectorField> gaussGrad(const volScalarField & vf)
{
    const fvMesh & mesh = vf.mesh();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    tmp<surfaceScalarField> surfTmp = linearScalarInterpolate(vf);
    surfaceScalarField & surfaceInterpolation = surfTmp();

    tmp<volVectorField> tmpResult
    (
        new volVectorField
        (
            IOobject
            (
                "gaussGrad("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions() / dimLength
        )
    );

    volVectorField & result = tmpResult();
    result = vector(0.0, 0.0, 0.0);

    forAll(surfaceInterpolation, facei)
    {
        label ocell = owner[facei];
        label ncell = neighbour[facei];

        vector sf = mesh.Sf()[facei];
        result[ocell] += surfaceInterpolation[facei] * sf;
        result[ncell] -= surfaceInterpolation[facei] * sf;
    }

    forAll(result, celli)
    {
        result[celli] /= mesh.V()[celli];
    }
//    result = result / mesh.V().field();

    // Nic nie robimy na brzegach - gradient nie zdefiniowany
    // dla green-gaussa. Chamsko damy 0 value dla gradientu, ale
    // nie powinniśmy nigdzie po drodze korzystać z tych wartości
    forAll(result.boundaryField(), patchId)
    {
            fvPatchVectorField& resultPatch = result.boundaryField()[patchId];

            //zakładamy że pochodna gradientu też jest 0 - bierzemy wartości z centrów
            resultPatch = resultPatch.patchInternalField()();

//            forAll(resultPatch.patch(), faceId)
//            {
//                resultPatch[faceId] = vector(0.0, 0.0, 0.0);
//            }
    }

    return tmpResult;
}


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
//    simpleControl simple(mesh);
#   include "createFields.H"

    std::fstream f, fB, f0, fB0, fconv;
    f.open("./values.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB.open("./valuesB.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    f0.open("./values0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fB0.open("./valuesB0.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);
    fconv.open("./fconv.vtk", std::fstream::in | std::fstream::out | std::fstream::trunc);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USTAWIENIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int iloscKrCz = runTime.controlDict().lookupOrDefault("iloscKrCz", 0);
//    int limitGradPsi = runTime.controlDict().lookupOrDefault("limitGradPsi", 64);
    int PsiCritic = runTime.controlDict().lookupOrDefault("PsiCritic", 30);

    scalar dx = scalar(1.0)/scalar(128.0);
    dimensionedScalar epsh( "epsh", dimLength, dx / scalar(4.0) );
    dimensionedScalar dtau( "dtau", dimTime, epsh.value() / scalar(2.0) );
//    dimensionedScalar dtauDimmless( "dtauDimmless", dimless, epsh.value() / scalar(1.0) );
    Info << "End\n" << endl;
    psiInit(PsiZero);
    Psi == PsiZero;
    AlphaAnalit = scalar(1.0)/(scalar(1.0)+Foam::exp(-PsiZero/epsh));
//    AlphaAnalit =  scalar(0.5)*(scalar(1.0)+Foam::tanh(PsiZero/(scalar(2.)*epsh)));
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INICJALIZACJA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    forAll(mesh.cellCentres(), cellI )
    {
//        Alpha[cellI] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZero[cellI]/(scalar(2.0)*epsh.value())) );  // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Alpha[cellI] = 1.0/(1.0+Foam::exp(-PsiZero[cellI]/(epsh.value())));  // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
//    Alpha.boundaryField().evaluate();
    forAll(mesh.boundary(), patchi )
    {
         const fvPatch& patch = mesh.boundary()[patchi];
         fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
         fvPatchScalarField& AlphaPatchAnalit = AlphaAnalit.boundaryField()[patchi];

         fvPatchScalarField& PsiZeroPatch = PsiZero.boundaryField()[patchi];
         forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
         {
//            AlphaPatch[faceId] = scalar(0.5)*( scalar(1.0) + Foam::tanh(PsiZeroPatch[faceId]/(scalar(2.0)*epsh.value())) );   // 2 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            AlphaPatch[faceId] = 1.0/(1.0+Foam::exp(-PsiZeroPatch[faceId]/(epsh.value())));   // 1 analit ~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            AlphaPatchAnalit[faceId] = AlphaPatch[faceId];
         }
    }
    AlphaZero == Alpha;
//    Alpha = 0.0;
//    Alpha.boundaryField() = 0.0;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ZAPIS  PO INICJALIZACJI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    forAll(mesh.cellCentres(), cellI )
    {
    f0 << mesh.cellCentres()[cellI].x()
      << " " << mesh.cellCentres()[cellI].y()
      << " " << std::setprecision(20) << Alpha[cellI]
      << " " << std::setprecision(20) << AlphaAnalit[cellI]
      << " " << std::setprecision(20) << PsiZero[cellI]
      << " " << std::setprecision(20) << Psi[cellI]
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
        const fvPatchScalarField& AlphaPatchAnalit = AlphaAnalit.boundaryField()[patchi];
        const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
        const fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];
        forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
        {
//            if(AlphaPatch[faceId] == 1)
//                Info << patch.name() << endl;

            fB0 << patch.Cf()[faceId].x()
                << " " << patch.Cf()[faceId].y()
                << " " << std::setprecision(20) << AlphaPatch[faceId]
                << " " << std::setprecision(20) << AlphaPatchAnalit[faceId]
                << " " << std::setprecision(20) << PZPatch[faceId]
                << " " << std::setprecision(20) << PPatch[faceId]
                << std::endl;
        }
    }

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALKOWANE W PSEUDO CZASIE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for ( int i = 0; i < iloscKrCz; i++ )
    {
        Info<< "reinitialization 3thO iter: " << i << endl;

        AlphaOld == Alpha;
//        AlphaOld.boundaryField() = Alpha.boundaryField();

        //AlphaFace == linearInterpolate(Alpha);
        AlphaFace == linearScalarInterpolate(Alpha);


        AlphaToPsi(Alpha, Psi, epsh);
        //PsiF == linearInterpolate(Psi);
        //Psi.boundaryField() == PsiZero.boundaryField();


//        gradPsi = fvc::grad(Psi);
//        setInternalToBoundary(gradPsi);
        gradPsi = gaussGrad(Psi);

//        gradPsi.boundaryField().evaluate();
//        extrapolateToBoundaries(gradPsi);

        //GradPsiFace == linearInterpolate( gradPsi );
        GradPsiFace == linearVectorInterpolate(gradPsi );

        mGradPsi = Foam::mag(gradPsi);
        //LimitGradPsi(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);

        //mGradPsiFace == linearInterpolate( mGradPsi );
        mGradPsiFace == linearScalarInterpolate( mGradPsi );

//!!!!!!!!!!!!!!!!!!!!!!
        snGradFunPsi(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        //LimitMGradPsiFaceCondition(mesh, PsiF, mGradPsiFaceCondition, epsh.value(), PsiCritic);
//!!!!!!!!!!!!!!!!!!!!!!
        phiR = Cf * AlphaFace*( scalar(1.0) - AlphaFace ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
        k1 == Alpha + fvc::div(phiR)*dtau;
        k1 == Foam::min(Foam::max(k1, scalar(0.0)), scalar(1.0) );
//        k1.boundaryField() == AlphaZero.boundaryField();
//        k1f == linearInterpolate(k1);
        k1f == linearScalarInterpolate(k1);


        AlphaToPsi(k1, Psi, epsh);
        //PsiF == linearInterpolate(Psi);
        //Psi.boundaryField() == PsiZero.boundaryField();

//        gradPsi = fvc::grad(Psi);
//        setInternalToBoundary(gradPsi);
        gradPsi = gaussGrad(Psi);


//        gradPsi.boundaryField().evaluate();
//        extrapolateToBoundaries(gradPsi);

        //GradPsiFace == linearInterpolate( gradPsi );
        GradPsiFace == linearVectorInterpolate( gradPsi );


        mGradPsi = Foam::mag(gradPsi);
        //LimitGradPsi(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);

        //mGradPsiFace == linearInterpolate( mGradPsi );
        mGradPsiFace == linearScalarInterpolate( mGradPsi );

        snGradFunPsi(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        //LimitMGradPsiFaceCondition(mesh, PsiF, mGradPsiFaceCondition, epsh.value(), PsiCritic);

        phiR = Cf * k1f*( scalar(1.0) - k1f ) * ( (mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / (mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
        k2 == scalar(0.75)* Alpha + scalar(0.25)* k1 + scalar(0.25)* fvc::div(phiR)*dtau;
        k2 == Foam::min(Foam::max(k2, scalar(0.0)), scalar(1.0) );
//        k2.boundaryField() == AlphaZero.boundaryField();
        //k2f == linearInterpolate(k2);
        k2f == linearScalarInterpolate(k2);


        AlphaToPsi(k2, Psi, epsh);
        //PsiF == linearInterpolate(Psi);
        //Psi.boundaryField() == PsiZero.boundaryField();
//        gradPsi = fvc::grad(Psi);
//        setInternalToBoundary(gradPsi);
        gradPsi = gaussGrad(Psi);

//        gradPsi.boundaryField().evaluate();
//        extrapolateToBoundaries(gradPsi);

        //GradPsiFace == linearInterpolate( gradPsi );
        GradPsiFace == linearVectorInterpolate( gradPsi );

        mGradPsi = Foam::mag(gradPsi);
        //LimitGradPsi(mesh, Psi, mGradPsi, epsh.value(), PsiCritic);

//        mGradPsiFace == linearInterpolate( mGradPsi );
        mGradPsiFace == linearScalarInterpolate( mGradPsi );

        snGradFunPsi(Psi, snGradPsi);
        mGradPsiFaceCondition = Foam::mag(snGradPsi);
        //LimitMGradPsiFaceCondition(mesh, PsiF, mGradPsiFaceCondition, epsh.value(), PsiCritic);

        phiR = Cf * k2f*( scalar(1.0) - k2f ) * ( ( mGradPsiFaceCondition - scalar(1.0) ) * GradPsiFace / ( mGradPsiFace + SMALL_NUMBER) ) & mesh.Sf();
        Alpha ==  scalar(1.0)/scalar(3.0)* Alpha + scalar(2.0)/scalar(3.0)* k2 + scalar(2.0)/scalar(3.0)*fvc::div(phiR)*dtau;
        Alpha == Foam::min(Foam::max(Alpha, scalar(0.0)), scalar(1.0) );

//        Alpha.boundaryField() == AlphaZero.boundaryField();

        double norm1c = Foam::sum(Foam::mag(Alpha-AlphaOld)).value() / Alpha.size();
        Info << "Norma conv = " << norm1c << endl;

        double norm1 = Foam::sum(Foam::mag(AlphaAnalit-Alpha)).value() / Alpha.size();
//        double norm2 = Foam::pow(Foam::sum(Foam::pow(TAnalit-T, 2)).value(), 0.5 ) / T.size();
//        double norm3 = Foam::max(Foam::mag(TAnalit-T)).value();
        Info << "Norma 1 analit = " << norm1 << endl;
//            Info << gradPsi;
            std::cout << std::endl;
        fconv << i << "   " << norm1 << std::endl;
    }
    Info << PsiCritic*epsh << " !!!!!!!!!!!!!      " << 9.*dx;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ZAPIS PO OBLICZENIACH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     snGradFunPsi(Psi, snGradPsi);
//     Info << snGradPsi << endl;
//gradPsi == fvc::grad(Psi);
//Info << gradPsi << endl;



        forAll(mesh.cellCentres(), cellI )
        {
        f << mesh.cellCentres()[cellI].x()
          << " " << mesh.cellCentres()[cellI].y()
          << " " << std::setprecision(26) << Alpha[cellI]
          << " " << std::setprecision(26) << AlphaAnalit[cellI]
          << " " << std::setprecision(26) << PsiZero[cellI]
          << " " << std::setprecision(26) << Psi[cellI]
          << " " << std::setprecision(26) << gradPsi[cellI].y()
          << " " << std::setprecision(26) << AlphaZero[cellI]
          << std::endl;
        }

        forAll(mesh.boundary(), patchi )
        {
            const fvPatch& patch = mesh.boundary()[patchi];
            if( patch.type() == "empty" )
            { continue; }
            fvPatchScalarField& AlphaPatch = Alpha.boundaryField()[patchi];
            const fvPatchScalarField& PZPatch = PsiZero.boundaryField()[patchi];
            const fvPatchScalarField& AlphaPatchAnalit = AlphaAnalit.boundaryField()[patchi];
            fvPatchScalarField& PPatch = Psi.boundaryField()[patchi];

            forAll(patch, faceId)                       // lopp over faces centers belonging to a given patch
            {
                fB << patch.Cf()[faceId].x()
                    << " " << patch.Cf()[faceId].y()
                    << " " << std::setprecision(26) << AlphaPatch[faceId]
                    << " " << std::setprecision(26) << AlphaPatchAnalit[faceId]
                    << " " << std::setprecision(26) << PZPatch[faceId]
                    << " " << std::setprecision(26) << PPatch[faceId]
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
// NNNNNNNNNNNNNNNNNNNNNNNNNNEEEEEEEEEEEEEEEEEEEEEEEEWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
