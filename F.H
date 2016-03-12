    List<scalar> nonDiag; // sums of non-diagonal elements

    nonDiag.resize(psi.size());
    forAll(nonDiag, i)
    {
	nonDiag[i] = 0.0;
    }
	
    const lduAddressing& addr = psiEqn.lduAddr();
    const labelList& lowerAddr = addr.lowerAddr();
    const labelList& upperAddr = addr.upperAddr();
	
    forAll(lowerAddr, i)
    {
// 	    A[lowerAddr[i]][upperAddr[i]] = pEqn.upper()[i];
// 	    A[upperAddr[i]][lowerAddr[i]] = pEqn.lower()[i];
	nonDiag[lowerAddr[i]] += psiEqn.upper()[i];
	nonDiag[upperAddr[i]] += psiEqn.lower()[i];
    }
	
    forAll(psi.boundaryField(),I)
    {
	const fvPatch &ptch=psi.boundaryField()[I].patch();
	forAll(ptch,J)
	{
	    int w=ptch.faceCells()[J];
	    nonDiag[w]+=abs(-psiEqn.internalCoeffs()[I][J]);
// 	    Info << "nonDiag" << "[" << w << "]" << " = " << nonDiag[w]  << nl << endl;
	}
    
    }
	
// F = alpha1;
// 	forAll(p_rghEqn.diag(),i)
    forAll(psiEqn.diag(),i)
    {

	      F[i] = (abs(1-psiEqn.diag()[i]))+nonDiag[i];
// 	F[i] = pEqn.diag()[i];
// 	      F[i] = abs(nonDiag[i]);
// 	      F[i] = abs(p_rghEqn.diag()[i])+abs(nonDiag[i]);
// 	      Info << "p_rghEqn.diag()[i]=" << p_rghEqn.diag()[i] << nl << endl;
// 	      Info << "UEqn.diag()=" << UEqn.diag() << nl << endl; 
    }
