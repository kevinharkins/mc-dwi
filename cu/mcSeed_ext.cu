#define MAG(a,b) (sqrt((a)*(a)+(b)*(b)))
#define TWO (2147483648.0)

/* 
 * linear congruent generator for device side random number calculations
*/ 
__device__ int lcg(int seed)
{
    return seed*1103515245+12345 % 2147483648;
}

/* seeds individual spins for monte carlo simulation into onlyl extra-axonal 
 * space
 *
 *   xi - output x,y locations
 *   rn1 - seeds for computation of x,y
 *   rn2 - seeds for determining M0r
 *   roi - outer,inner radii
 *   xo - x0,y0 for axons
 *   nA - number of axons
 *   Lx - boundary of geometry
 *   M0r - ratio of myelin to non-myelin density
 *   N - number of spins
 */
__global__ void mcSeed_ext(float *xi, float *rn1, float *rn2, float* roi, 
        float* xo, int nA, float Lx, float M0r, int N) 
{
	int p = threadIdx.x + blockDim.x*blockIdx.x;
    
    if (p<N)
    {
        float r1, r2, r3, xt, yt;
        int m;
        
        r1 = rn1[p];
        r2 = rn1[p+N];
        r3 = rn2[p];
        
        xt = r1*Lx;
        yt = r2*Lx;
        
        int state=1;

        // loop until the spin is placed
        while (state)
        {
            state = 0;
            
            // Is this spin within an axon?
            for (m=0;m<nA;m++) 
            {
                if ( MAG(xt-xo[m],yt-xo[m+nA])<roi[m] )     
                {
                    // in side axon or myelin, reseeding...
                    state = 1;
                    r1 = abs(lcg(round(r1*TWO)))/TWO;
                    r2 = abs(lcg(round(r2*TWO)))/TWO;
                    r3 = abs(lcg(round(r3*TWO)))/TWO;
                    
                    xt = r1*Lx;
                    yt = r2*Lx;
                }
            }
        }
    
        xi[p] = xt;
        xi[p+N] = yt;
    }
}
