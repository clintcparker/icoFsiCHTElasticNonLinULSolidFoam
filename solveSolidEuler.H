    {
#       include "readStressedFoamControls.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 0;

        lduMatrix::debug = 0;

#       include "EulerCoeffs.H"

        do
        {
            DU.storePrevIter();

            fvVectorMatrix DUEqn
            (
                Cn*rho*fvm::ddt(DU)
              - Co*rho*DV.oldTime()
             ==
                fvm::laplacian(2*mu + lambda, DU, "laplacian(DDU,DU)")
              - fvc::laplacian(mu + lambda, DU, "laplacian(DDU,DU)")
              + fvc::div
                (
                    mu*gradDU.T()
                  + lambda*(I*tr(gradDU))
                  + mu*(gradDU&gradDU.T())
                  + 0.5*lambda*(I*tr(gradDU & gradDU.T()))
                  + (sigma & DF.T())
                  + (DSigma & DF.T()),
                    "div(sigma)"
                )
            );

            solverPerf = DUEqn.solve();

            DU.relax();

            if(iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            gradDU = fvc::grad(DU);

            DF = gradDU.T();

#           include "calculateDSigma.H"
        } 
        while
        (
            solverPerf.initialResidual() > convergenceTolerance 
         && ++iCorr < nCorr
        );

        Info << "Solving for " << DU.name() 
            << ", Initial residual = " << initialResidual 
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr << endl;

        DV = fvc::ddt(DU);

        lduMatrix::debug = 1;

//Code addition which will hopefully add the CHT portion to the solid (following block till end is taken straight from Cht/solid/solvesolid.H

        if (finalIter)
        {
            mesh.data::add("finalIteration", true);
        }

        {
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                tmp<fvScalarMatrix> TEqn
                (
                    fvm::ddt(rho*cp, T)
                  - fvm::laplacian(K, T)
                );
                TEqn().relax();
                TEqn().solve(mesh.solver(T.select(finalIter)));
            }

            Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
        }

        thermo.correct();

        if (finalIter)
        {
            mesh.data::remove("finalIteration");
        }


    }

