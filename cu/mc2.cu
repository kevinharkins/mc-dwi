#define MAG(a,b) (sqrt((a)*(a)+(b)*(b)))
#define ANG(a,b) (atan2((b),(a)))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define PI (3.14159265359)
#define TOL (1e-4)
#define INF __int_as_float(0x7f800000)
#define TWO (2147483648.0)

// how far do we jump forward before we hit this circle?
// x1 and y1 are relative to x0 and y0
__device__ float jumpLin(float x1,float y1,float ang,float ro)
{ 
    float b = 2*(x1*cos(ang) + y1*sin(ang));
    float c = x1*x1 + y1*y1 - ro*ro;
    if (b*b-4*c < 0) return INF;
    float a = (-b+sqrt(b*b-4*c))/2;
    a = a<=0? INF:a;
    b = (-b-sqrt(b*b-4*c))/2;
    b = b<=0? INF:b;
    return MIN(a,b);
}

// how far do we jump forward before we hit this circle?
// but for myelin (in radial coordinates)
__device__ float jumpRadFrac(float r1,float jnorm,float ro,float ri)
{ 
    if (jnorm > 0)
        return MIN((ro-r1)/jnorm,1);
    else
        return MIN((ri-r1)/jnorm,1);
}

/* 
 * linear congruent generator for device side random number calculations
*/ 
__device__ int lcg(int seed)
{
    return seed*1103515245+12345 % 2147483648;
}

/*
 * mc2 calculates 1 timestep of a Monte Carlo diffusion simulation
 *
 *  xi - spin x,y locations (input & output)
 *  sig - signal magntidue & phase (input & output) 
 *  rn1 - seeds for computation of x,y (input)
 *  Dn - diffusion coefficient in the direction normal to myelin
 *  Dp - diffusion coefficient in the direction perpendicular to myelin
 *  Do - diffusion coefficient in non-myelin
 *  G - the amplitude of the gradient waveform at this timestep
 *  A - the area of the gradient waveform at this timestep
 *  roi - inner, outer radius of the axon
 *  xo - x0,y0 for axons (input)
 *  nA - number of axons (input)
 *  Lx - boundary of geometry (input())
 *  M0r - ratio of myelin to non-myelin density (input)
 *  R2m - 1/T1 of myelin (input)
 *  R2o - 1/T2 of intra & extra-axonal space (input)
 *  N - number of spins (input)
*/
__global__ void mc2(float *xi, float *sig, float* rn1, float Dn, float Dp, float Do,
                    float G, float A, float* roi, float* xo, int nA, float Lx, 
                    float M0r, float R2m, float R2o, int N) 
{
	int p = threadIdx.x + blockDim.x*blockIdx.x;
    
    if (p<N)
    {
        float dr, da, tmp1, rn, ang, xt, yt, jump, j1;
        int n=0,m,thA=-1,tmp;
        int state=-1;
        ang = rn1[p]*2*PI;
        rn = rn1[p+N];
        
        xt = xi[p];
        yt = xi[p+N];
        jump = 0;
        
        // which axon are we interacting with?
        for (m=0;m<nA;m++) 
        {
            if ( MAG(xt-xo[m],yt-xo[m+nA]) < roi[m] )
            { // this spin is inside an axon
                thA = m;
                if (MAG(xt-xo[m],yt-xo[m+nA])>roi[m+nA]) state = 3; // inside myelin
                else state = 0; // inside intra-axonal space
                break;
            }
        }
        
        //////////////////////////////////////////////////////
        // main loop for multiple interactions  
        //////////////////////////////////////////////////////
        while ((jump<1) && n<20)
        {
            n++; // count to get rid of race conditions and precision issues
            
            if (thA==-1)
            { 
                // This spin is outside of an axon
                
                tmp = -1; // index to the the potential new axon

                // loop over all axons, looking for interactions with this spin
                for (m=0;m<nA;m++)
                {
                    j1 = 1-jump;

                    // does this interact with the outside of an axon?
                    if (MAG(xt-xo[m],yt-xo[m+nA])<(roi[m]+Do*j1)) 
                        j1 = MIN(jumpLin(xt-xo[m],yt-xo[m+nA],ang,roi[m])/Do,j1);
                    if (j1 < (1-jump))
                    {
                        tmp = m;
                        break;
                    }
                }

                // to fix floating point precision issues (?), do the same  
                // calculation as before for a radius + some tolderance
                if (tmp == -1)
                {
                    for (m=0;m<nA;m++)
                    {
                        j1 = 1-jump;
                        // does this interact with the outside of an axon?
                        if (MAG(xt-xo[m],yt-xo[m+nA])<(roi[m]+Do*j1+TOL)) 
                            j1 = MIN(jumpLin(xt-xo[m],yt-xo[m+nA],ang,roi[m]+TOL)/Do,j1);
                        if (j1 < (1-jump))
                        {
                            tmp = m; // tmp is the potential new axon
                            break;
                        }
                    }
                }
                
                if (tmp>=0)
                    j1 = MIN(jumpLin(xt-xo[tmp],yt-xo[tmp+nA],ang,roi[tmp])/Do,1-jump);
                else
                    j1 = 1-jump;
                yt = yt+sin(ang)*j1*Do;
                xt = xt+cos(ang)*j1*Do;
                sig[p] = sig[p]*(1-(1-R2o)*j1);
                jump = j1+jump;
                
                // if it hit a boundary, does it reflect or go through?
                if (tmp >= 0)
                {
                    if (jump<1)
                    {
                        if (rn>(Dn/Do*M0r)) 
                            ang = PI-ang+2*ANG(xt-xo[tmp],yt-xo[tmp+nA]); // reflects
                        else // goes into myelin
                        {
                            state = 3;
                            ang -= ANG(xt-xo[tmp],yt-xo[tmp+nA]);
                            thA = tmp;
                        }
                        rn = abs(lcg(round(rn*TWO)))/TWO;
                    }
                    if (abs(MAG(xt-xo[tmp],yt-xo[tmp+nA])-roi[tmp])<=TOL)
                    {
                        if (state==3) j1 = TOL;
                        else j1 = TOL;
                        yt = yt+j1*sin(ANG(xt-xo[tmp],yt-xo[tmp+nA]));
                        xt = xt+j1*cos(ANG(xt-xo[tmp],yt-xo[tmp+nA]));
                    }
                }
            }
            else if (state==3)
            {
                // This spin is inside myelin

                // put into relative coordinates for the jump calculation
                xt = xt-xo[thA];
                yt = yt-xo[thA+nA];
                
                // calculate the jump
                if (rn>((MAG(xt,yt)+Dn/2*(1-jump)*cos(ang))/MAG(xt,yt))) ang += PI-2*ang;
                rn = abs(lcg(round(rn*TWO)))/TWO;
                dr = Dn*cos(ang);
                da = Dp*sin(ang)/MAG(xt,yt); 
                j1 = MIN(jumpRadFrac(MAG(xt,yt),Dn*cos(ang),roi[thA],roi[thA+nA]),(1-jump));
                tmp1 = (MAG(xt,yt)+dr*j1)*cos(ANG(xt,yt)+da*j1);
                yt  = (MAG(xt,yt)+dr*j1)*sin(ANG(xt,yt)+da*j1);
                xt = tmp1;
                sig[p] = sig[p]*(1-(1-R2m)*j1);
                jump = j1+jump;
                
                // back into absolute coordinates
                xt = xt+xo[thA];
                yt = yt+xo[thA+nA];
                
                // if it hit a boundary, does it reflect or go through?
                tmp = thA;
                if (jump<1)
                {
                    if (rn>(Do/Dn/M0r)) 
                    {
                        // reflects back
                        ang += PI-2*ang;
                    } 
                    else 
                    {
                        // leaves myelin
                        ang += ANG(xt-xo[thA],yt-xo[thA+nA]);
                        state = 0;
                        if (abs(MAG(xt-xo[thA],yt-xo[thA+nA])-roi[thA])<TOL) 
                        {
                            // and leaves the axon
                            thA=-1; 
                        }
                    }

                    rn = abs(lcg(round(rn*TWO)))/TWO;
                }

                // fix fp precision issues near the inner & outer boundary
                if (abs(MAG(xt-xo[tmp],yt-xo[tmp+nA])-roi[tmp])<TOL)
                {
                    j1=TOL;
                    if (state==3) j1 = -TOL;
                    yt = yt+j1*sin(ANG(xt-xo[tmp],yt-xo[tmp+nA]));
                    xt = xt+j1*cos(ANG(xt-xo[tmp],yt-xo[tmp+nA]));
                }

                if (abs(MAG(xt-xo[tmp],yt-xo[tmp+nA])-roi[tmp+nA])<TOL)
                {
                    j1=TOL;
                    if (state==0) j1 = -TOL;
                    yt = yt+j1*sin(ANG(xt-xo[tmp],yt-xo[tmp+nA]));
                    xt = xt+j1*cos(ANG(xt-xo[tmp],yt-xo[tmp+nA]));
                }
            }
            else if (state==0) 
            {
                // This spin is in intra-axonal space

                // put into relative coordinates
                xt = xt-xo[thA];
                yt = yt-xo[thA+nA];
                
                // calculate the jump
                j1 = MIN(jumpLin(xt,yt,ang,roi[m+nA])/Do,1-jump);
                yt = yt+sin(ang)*j1*Do;
                xt = xt+cos(ang)*j1*Do;
                sig[p] = sig[p]*(1-(1-R2o)*j1);
                jump = j1+jump;
                
                // if it hit a boundary, does it reflect or go through?
                if (jump<1)
                {
                    if (rn>(Dn/Do*M0r)) 
                    {
                        // reflects back
                        ang = PI-ang+2*ANG(xt,yt); 
                    }
                    else 
                    {
                        // goes into myelin
                        state = 3;
                        ang -= ANG(xt,yt);
                    }
                    rn = abs(lcg(round(rn*TWO)))/TWO;
                }

                // fix FP preciions issues near the boundary 
                if (abs(MAG(xt,yt)-roi[thA+nA])<TOL)
                {
                    if (state==0) j1 = -1;
                    else j1 = 1;
                    yt = yt+j1*sin(ANG(xt,yt))*TOL;
                    xt = xt+j1*cos(ANG(xt,yt))*TOL;
                }
                // back into absolute coordinates
                xt = xt+xo[thA];
                yt = yt+xo[thA+nA];
            }
        }

        // increment the signal phase
        sig[p+N] += xt*G;

        // enforce periodic boundary conditions
        if (xt>Lx)
        {
            xt -= Lx;
            sig[p+N] -= Lx*A;
        }
        if (xt<0)
        {
            xt += Lx;
            sig[p+N] += Lx*A;
        }
        if (yt<0) yt+=Lx;
        if (yt>Lx) yt-=Lx;
    
        // output new x,y coordinates
        xi[p] = xt;
        xi[p+N] = yt;
        rn1[p] = ang/2/PI;
        rn1[p+N] = rn;
    }
}
