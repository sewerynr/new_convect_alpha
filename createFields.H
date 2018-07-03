    volScalarField Alpha
    (
        IOobject
        (
            "Alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField AlphaOld
    (
        IOobject
        (
            "AlphaOld",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Alpha
    );

    volScalarField AlphaZero
    (
        IOobject
        (
            "AlphaZero",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Alpha
    );

//    volScalarField Psi
//    (
//        IOobject
//        (
//            "Psi",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("Psi", dimLength, scalar(0.0))
//    );

    volScalarField Psi
    (
        IOobject
        (
            "Psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
        //dimensionedScalar("Psi", dimLength, scalar(0.0))
    );

    surfaceScalarField PsiF
    (
        IOobject
        (
            "Psi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Psi)
    );

    volScalarField PsiZero
    (
        IOobject
        (
            "PsiZero",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
//        dimensionedScalar("PsiZero", dimLength, scalar(0.0))
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField C
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("C", dimLength/dimTime, scalar(1.0))
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );


    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    surfaceScalarField Cf
    (
        IOobject
        (
            "Cf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(C)
    );

    surfaceScalarField phiR
    (
        IOobject
        (
            "phiR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiR", (dimLength*dimLength*dimLength)/dimTime, scalar(0.0))
    );

    surfaceScalarField AlphaFace
    (
        IOobject
        (
            "AlphaFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Alpha)
    );

    surfaceScalarField PsiFace
    (
        IOobject
        (
            "PsiFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(Psi)
    );

    volScalarField mGradPsi
    (
        IOobject
        (
            "mGradPsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mGradPsi", dimless, scalar(-100000000.0))
    );

    volScalarField mGradPsiLimitedPlot
    (
        IOobject
        (
            "mGradPsiLimitedPlot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mGradPsiLimitedPlot", dimless, scalar(-100000000.0))
    );

    volVectorField gradPsi
    (
        IOobject
        (
            "gradPsi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField gradPsiTemp
    (
        IOobject
        (
            "gradPsiTemp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(Psi)
    );

    surfaceVectorField GradPsiFace
    (
        IOobject
        (
            "GradPsiFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(fvc::grad(Psi))
    );



    surfaceScalarField  snGradPsi
    (
        IOobject
        (
            "snGradPsi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::snGrad(Psi)
    );

    surfaceScalarField  mGradPsiFace
    (
        IOobject
        (
            "mGradPsiFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Foam::mag(linearInterpolate(fvc::grad(Psi)))
    );

    surfaceScalarField  mGradPsiFaceCondition
    (
        IOobject
        (
            "mGradPsiFaceCondition",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Foam::mag(linearInterpolate(fvc::grad(Psi)))
    );

    volScalarField k1
    (
        IOobject
        (
            "k1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("k1", dimless, scalar(0))
    );

    volScalarField k2
    (
        IOobject
        (
            "k2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("k2", dimless, scalar(0))
    );

    surfaceScalarField k1f
    (
        IOobject
        (
            "k1f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("k1f", dimless, scalar(0))
    );

    surfaceScalarField k2f
    (
        IOobject
        (
            "k2f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("k2f", dimless, scalar(0))
    );

    volScalarField AlphaAnalit
    (
        IOobject
        (
            "AlphaAnalit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("AlphaAnalit", dimless, scalar(-100000000.0))
        //0.5*(1.+Foam::tanh(PsiZero/(2.*epsH)))
    );
