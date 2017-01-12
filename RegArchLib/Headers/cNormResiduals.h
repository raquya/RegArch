#pragma once
#ifndef _CNORMRESIDUALS_H_
#define _CNORMRESIDUALS_H_
#include "cAbstResiduals.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
 \file cNormResiduals.h
 \brief Definition of the class for N(0, 1) conditional distribution.

 \author Jean-Baptiste DURAND, Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {

	/*! 
	 * \class cNormResiduals
	 * \brief  Class to implement the N(0, 1) residuals
	 */

	class _DLLEXPORT_ cNormResiduals: public cAbstResiduals
	{
	public :
		cNormResiduals(const cDVector* theDistrParameter=NULL, bool theSimulFlag=true) ; ///< a simple constructor
		virtual ~cNormResiduals() ; ///< a simple destructor
		virtual cAbstResiduals* PtrCopy() const ; /// < Return a copy of *this
		void Print(ostream& theOut) const ; ///< print the distribution type
		void Generate(uint theNSample, cDVector& theEpst) const ; ///< Draw a sample from residual distribution 
		double LogDensity(double theX) const ;
		/** Return the number of parameters in distribution */
		uint GetNParam(void) const ;
		/** Compute the derivative of log density with respect to the random variable (theGradData[0]) \e and the gradient
            of log density with respect to the model parameters (other components in theGradData) */
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData) const;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex) const ;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
	} ;

}

#endif //_CNORMRESIDUALS_H_

