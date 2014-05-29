#include <Rcpp.h>
using namespace Rcpp;

double L2DistCPP(double x, NumericVector ys) 
{
    int n = ys.size();
    double out = 0;
    for(int i = 0; i < n; ++i) {
      out += (ys[i] - x)*(ys[i] - x);
    }
    return out;
}

double L1DistCPP(double x, NumericVector ys) 
{
    int n = ys.size();
    double out = 0;
    for(int i = 0; i < n; ++i) {
      out += fabs(ys[i] - x);
    }
    return out;
}

double crossEntropyCPP(NumericVector pX) 
{
    int nX  = pX.size();
	double totalX = 0;
    for(int i = 0; i < nX; ++i) 
	{
      if (pX[i] != 0) { totalX +=  pX[i]*log(pX[i]); }	  
    }
	return -(totalX);
}

double giniCPP(NumericVector pX) 
{
    int nX  = pX.size();
	double totalX = 0;
    for(int i = 0; i < nX; ++i) 
	{  totalX +=  pX[i]*(1 - pX[i]);  }
	return totalX;
}

NumericMatrix XMinMaxCPP(NumericMatrix X)
{
	int ncolX = X.ncol();
	NumericMatrix Y(2, ncolX);
		
	for(int j = 0; j < ncolX; ++j) 	
	{ Y(0,j) =  min(X(_,j)); Y(1,j)= max(X(_,j));	}
	
	return Y;
}

NumericVector checkUniqueObs(NumericMatrix X)
{
	int ncolX = X.ncol();
	NumericVector uniqueObs(ncolX);
		
	for(int j = 0; j < ncolX; ++j) 	
	{ 
		uniqueObs[j] = 2;
		if (min(X(_,j)) == max(X(_,j)))
		{ uniqueObs[j] = 1; }
	}
	
	return uniqueObs;
}

NumericVector classifyCPP( NumericMatrix treeObject, NumericVector X)
{
	int k = 0; 
	    	
	while (treeObject(k,4) >  0)
	{
		if ( X[treeObject(k,2)-1] >  treeObject(k,3) ) 
		{	k = treeObject(k,1)-1; }
		else		
		{  k = treeObject(k,0)-1; }
	}
			
	NumericVector ZZ(3);
	ZZ[0] = treeObject(k,5);
    ZZ[1] =	treeObject(k,6);
	ZZ[2] = k + 1;
	
	return ZZ;
}

NumericMatrix classifyMatrixCPP( NumericMatrix treeObject, NumericMatrix X)
{
	int k = 0; 
	int nrowX = X.nrow();
	
	NumericMatrix ZZ(nrowX,3);
	NumericVector Xtmp;
	    	
	for(int i = 0; i < nrowX; ++i) 	
	{
		Xtmp = X(i,_);
		while (treeObject(k,4) >  0)
		{
			if ( Xtmp[treeObject(k,2)-1] >  treeObject(k,3) ) 
			{	k = treeObject(k,1)-1; }
			else		
			{  k = treeObject(k,0)-1; }
		}
		
		ZZ(i,0) = treeObject(k,5);
		ZZ(i,1) = treeObject(k,6);
		ZZ(i,2) = k + 1;
		
		k = 0;
	}
		
	return ZZ;
}

NumericVector runifMatrixCPP(NumericMatrix X)
{
	int ncolX = X.ncol();
	NumericVector res(ncolX);
	RNGScope scope; 
	
	for(int j = 0; j < ncolX; ++j) 	{	res[j] = runif(10,  min(X(_,j)),  max(X(_,j)))[1];	}
	
	return res;
}

NumericVector L2InformationGainCPP(NumericVector Y, NumericMatrix X, NumericVector splitPoint)
{
	int ncolX = X.ncol();
	int nrowX = X.nrow();
	NumericVector err(ncolX);
	double mYLow;		
	double mYHigh;
	int a;
	int b;
	
	for (int j = 0; j < ncolX; ++j)
	{
		a = 0; b = 0;
		NumericVector YLow(nrowX); 
		NumericVector YHigh(nrowX);
	    YLow[0] = 0; YHigh[0] = 0;
		mYLow = 0; mYHigh = 0;
		for (int i = 0; i < nrowX; ++i)
		{
			if (X(i,j) <= splitPoint[j])
			{	YLow[a] = Y[i];  mYLow += Y[i];	a = a + 1;}
			else
			{	YHigh[b] = Y[i]; mYHigh += Y[i];  b = b + 1; }
		}
	
		if (b == 0)
		{	err[j]= L2DistCPP(mYLow/a, YLow[Range(0,a-1)]);	}
		else
		{ 	err[j]= L2DistCPP(mYLow/a, YLow[Range(0,a-1)]) + L2DistCPP(mYHigh/b, YHigh[Range(0,b-1)]);	}
	}
		
	return err;
}

NumericVector L1InformationGainCPP(NumericVector Y, NumericMatrix X, NumericVector splitPoint)
{
	int ncolX = X.ncol();
	int nrowX = X.nrow();
	NumericVector err(ncolX);
	double mYLow;		
	double mYHigh;
	int a;
	int b;
	
	for (int j = 0; j < ncolX; ++j)
	{
		a = 0; b = 0;
		NumericVector YLow(nrowX); 
		NumericVector YHigh(nrowX);
	    YLow[0] = 0; YHigh[0] = 0;
		mYLow = 0; mYHigh = 0;
		for (int i = 0; i < nrowX; ++i)
		{
			if (X(i,j) <= splitPoint[j])
			{	YLow[a] = Y[i];  mYLow += Y[i];	a = a + 1;}
			else
			{	YHigh[b] = Y[i]; mYHigh += Y[i];  b = b + 1; }
		}
	
		if (b == 0)
		{	err[j]= L1DistCPP(mYLow/a, YLow[Range(0,a-1)]);	}
		else
		{ 	err[j]= L1DistCPP(mYLow/a, YLow[Range(0,a-1)]) + L1DistCPP(mYHigh/b, YHigh[Range(0,b-1)]);	}
	}
		
	return err;
}

double conditionalCrossEntropyCPP(NumericVector pkLow, NumericVector pkHigh,double pLow, double pHigh) 
{
	return (crossEntropyCPP(pkLow)*pLow + crossEntropyCPP(pkHigh)*pHigh);
}

double conditionalGiniCPP(NumericVector pkLow, NumericVector pkHigh, double pLow, double pHigh) 
{
	return (giniCPP(pkLow)*pLow + giniCPP(pkHigh)*pHigh);
}

NumericVector entropyInformationGainCPP(NumericVector Y, NumericMatrix X, NumericVector splitPoint, NumericVector classes, int nClasses, 
	NumericVector eClasswt, int entropy)
{
	double ncolX = X.ncol();
	double nrowX = X.nrow();
	NumericVector err(ncolX);
	NumericVector pLow(nClasses);
	NumericVector pHigh(nClasses);
	double nLow;
	double nHigh;
	
	for (int j = 0; j < ncolX; ++j)
	{
	    nLow = 0;
		nHigh = 0;
		
		for (int k = 0; k < nClasses; ++k)
		{   pLow[k] = 0; pHigh[k] = 0; }
		
		for (int i = 0; i < nrowX; ++i)
		{
			if (X(i,j) <= splitPoint[j])
			{	
				nLow = nLow + 1;
				for (int k = 0; k < nClasses; ++k)
				{
					if (Y[i] == classes[k])
					{  	pLow[k] = pLow[k] + 1; }
				}
			}
			else
			{	
				nHigh = nHigh + 1;
				for (int k = 0; k < nClasses; ++k)
				{
					if (Y[i] == classes[k])
					{  	pHigh[k] = pHigh[k]  + 1; }
				}
			}
		}
			
		if (eClasswt[0] > 0)
		{
			for (int ec = 0; ec < nClasses; ++ec)
			{
				pLow[ec] = pLow[ec] * eClasswt[ec];
				pHigh[ec] = pHigh[ec] * eClasswt[ec];
			}
		}
	
		if (nHigh == 0)
		{	nHigh = nHigh  + 1; }
		
		for (int k = 0; k < nClasses; ++k)
		{
			pLow[k] = pLow[k]/nLow; 
			pHigh[k] = pHigh[k]/nHigh; 
		}	
		
		double pPrior1 = nLow/nrowX;
		double pPrior2 = 1 - pPrior1;
					
		if (entropy == 1)
		{	err[j] = conditionalCrossEntropyCPP(pLow, pHigh,  pPrior1, pPrior2); }
		else
		{ err[j] = conditionalGiniCPP(pLow, pHigh,  pPrior1, pPrior2); }
	}
	
	return err;
}

NumericVector sortCPP(NumericVector x) 
{
   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

RCPP_MODULE(rUniformForestCppClass)
{
	using namespace Rcpp; 
	
	function("L2DistCPP", &L2DistCPP);
	function("L1DistCPP", &L1DistCPP);
	function("crossEntropyCPP", &crossEntropyCPP);
	function("giniCPP", &giniCPP);
	function("XMinMaxCPP", &XMinMaxCPP);
	function("checkUniqueObsCPP", &checkUniqueObs);
	function("classifyCPP", &classifyCPP);
	function("classifyMatrixCPP", &classifyMatrixCPP);
	function("runifMatrixCPP", &runifMatrixCPP);
	function("L2InformationGainCPP", &L2InformationGainCPP);
	function("L1InformationGainCPP", &L1InformationGainCPP);
	function("conditionalCrossEntropyCPP", &conditionalCrossEntropyCPP);
	function("conditionalGiniCPP", &conditionalGiniCPP);
	function("entropyInformationGainCPP", &entropyInformationGainCPP);
	function("sortCPP", &sortCPP);
}
