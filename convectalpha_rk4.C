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
    T& u1 = rhs.u1Ref();
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

    double eval(const double &y0){
        return y0;
    }
};


void test_rk3()
{
    double y0 = 1;

    double dt = 0.01;

    dYdt rhs;

    double y = 0;

    for(int i=0; i<100; ++i)
    {
        rk3(y0, dt, rhs, y);

        Info <<"time "<<(i+1)*dt <<"\t error=" << Foam::mag( y - Foam::exp( (i+1)*dt) )<< endl;

        y0 = y;
    }
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
