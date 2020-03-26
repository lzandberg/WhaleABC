package org.ChaffinchABC;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

public class Priors {

	private static final long DOUBLE_MASK = (1L << 53) - 1;
        private static final double NORM_53 = 1. / (1L << 53);
        private long state0, state1;
	
	//public double[] priorMax, priorMin, parameters, variables;
	//int[] priorType;
	int nParams;
        //MAINLAND EUROPE
        //public double[] variables=new double[] {2,1500,101,500,0.4,6, 2, 1.1, 0.05};  
        
        //MAINLAND TENTSMUIR
        //public double[] variables=new double[] {2,3000,23,100,0.4,6, 2, 1.1, 0.05};  
        
        //public double[] variables=new double[] {2,500, 533, 101,0.4,6, 2, 1.1, 0.05};
        
        //Azores Faial
        //public double[] variables=new double[] {2,500,120,80,0.4,6, 2, 1.1, 0.05};  
        
        //TENERIFE COMMON
        //public double[] variables=new double[] {2,1000, 540, 30,0.4,6, 2, 1.1, 0.05};
		
        
       //public double[] variables=new double[] {2,2000,120,80,0.4,6, 2, 1.1, 0.05};
         
         //BLUE CHAFFINCH PRIORS:
        public double[] variables=new double[] {2,500,120,80,0.4,4, 2, 1.1, 0.05}; 
        
        
	int[] priorType= new int[]{1, 1, 1, 1, 1};
	double[] priorMax=new double[] {0.05, 1, 5, 40, 2.5};
	double[] priorMin=new double[] {0.0001, 0.0001, 0.5, 0.5, 0.8};
        
        //double[] priorMin=new double[] {0.00723, 0.00527, 1.109, 3.510, 0.893};
        //double[] priorMax=new double[] {0.00724, 0.00528, 1.11, 3.511, 0.894};

        //double[] priorMin=new double[] {0.0231, 0.054, 1.82, 0.83, 2.07};
        //double[] priorMax=new double[] {0.0232, 0.055, 1.83, 0.84, 2.08};
        
        //double[] priorMin=new double[] {0.005, 0.32, 1.64, 0.98, 2.27};
	//double[] priorMax=new double[] {0.006, 0.33, 1.65, 0.99, 2.28};
        
        //double[] priorMin=new double[] {0.0029, 0.010, 1.99, 1.98, 1.00};
	//double[] priorMax=new double[] {0.003, 0.011, 2, 1.99, 1.01};
        
       
        
        
        public Priors(){
            nParams=priorType.length;
        }
       
	public Priors(long seed) {
		setSeed(seed);
		
		
		
		//priorMax=new double[] {2, 1000, 533, 100, 0.4, 6, 0.1, 0, 3, 10, 1.1, 0.06, 2, 0.02, 198};
		//priorMin=new double[] {2, 1000, 533, 100, 0.4, 6, 0.0002, 0, 0.5, 1, 1.1, 0.06, 0.5, 0.02, 198};
		
		//priorMax=new double[] {2, 1000, 500, 100, 0.4, 6, 0.1, 0, 3, 10, 1.1, 0.06, 2, 0.02, 51};
		//priorMin=new double[] {2, 1000, 500, 100, 0.4, 6, 0.0002, 0, 0.5, 1, 1.1, 0.06, 0.5, 0.02, 51};
		
		nParams=priorType.length;
		
		//modelType=(int)p[0];
		//nYears=(int)p[1];
		//nx=(int)p[2];
		//ny=(int)p[3];
		//popSize=nx*ny;
		//mortalityRate=p[4];
		//sylsPerSong=(int)p[5];
		//mutationVar=p[6];
		//recomRate=p[7];
		//betaa=p[8];
		//betab=p[9];
		//neighthresh=p[10];
		//typeThresh=p[11];
		//confBias=p[12];
		//minDist=p[13];
		//numTypes=(int)p[14];
			
	}
	
	public double[] sampleFromPriors() {
		double[]parameters=new double[nParams];
		
		for (int i=0; i<nParams; i++) {
			if (priorMax[i]==priorMin[i]) {
				parameters[i]=priorMax[i];
			}
			else {
				if (priorType[i]==0) {
					parameters[i]=sampleUniformInteger(i);
				}
				else if (priorType[i]==1) {
					parameters[i]=sampleLogUniform(i);
				}
				else if (priorType[i]==2) {
					parameters[i]=sampleUniform(i);
				}
			}
			
		}
		return parameters;
	}
	
	public double sampleUniformInteger(int i) {
		int q=(int)(priorMax[i]-priorMin[i]);
		double p=nextInt(q)+priorMin[i];
		return p;
		
	}
	
	public double sampleLogUniform(int i) {
		double q=Math.log(priorMax[i])-Math.log(priorMin[i]);
		double p=nextDouble()*q+Math.log(priorMin[i]);
		return Math.exp(p);
	}
	
	public double sampleUniform(int i) {
		double q=priorMax[i]-priorMin[i];
		double p=nextDouble()*q+priorMin[i];
		return p;
	}
	
	public double[] transform(double[] v) {
		double[] v3=new double[nParams];
		for (int i=0; i<nParams; i++) {
			if (priorType[i]==1) {
				v3[i]=Math.log(v[i]);
			}
			else {
				v3[i]=v[i];
			}
		}
		return v3;
	}
	
	public double[] drawFromProposal(double[][] cov, double[] v){	
		double[] v2=null;
		boolean allok=false;
		
		MultivariateNormalDistribution mnd=new MultivariateNormalDistribution(v, cov);
		mnd.reseedRandomGenerator(nextLong());
		
		while(!allok){	
			v2= mnd.sample();
			allok=true;
			for (int i=0; i<nParams; i++) {
				if (priorType[i]==1) {
					v2[i]=Math.exp(v2[i]);
				}
				else {
					v2[i]=v2[i];
				}
				
				if (v2[i]>priorMax[i]) {allok=false;}
				if (v2[i]<priorMin[i]) {allok=false;}
			}
		}
		return v2;
	}
	
	public long nextLong() {
	       final long s0 = state0;
	       long s1 = state1;
	       final long result = s0 + s1;

	       s1 ^= s0;
		   state0 = Long.rotateLeft(s0, 55) ^ s1 ^ (s1 << 14); // a, b
		   state1 = Long.rotateLeft(s1, 36); // c

		   return result;
		}

		    /**
		     * Can return any int, positive or negative, of any size permissible in a 32-bit signed integer.
		     * @return any int, all 32 bits are random
		     */
		    public int nextInt() {
		        return (int)nextLong();
		    }

		    /**
		     * Exclusive on the upper bound.  The lower bound is 0.
		     * @param bound the upper bound; should be positive
		     * @return a random int less than n and at least equal to 0
		     */
		    public int nextInt( final int bound ) {
		        if ( bound <= 0 ) return 0;
		        int threshold = (0x7fffffff - bound + 1) % bound;
		        for (;;) {
		            int bits = (int)(nextLong() & 0x7fffffff);
		            if (bits >= threshold)
		                return bits % bound;
		        }
		    }

		    public double nextDouble() {
		        return (nextLong() & DOUBLE_MASK) * NORM_53;
		    }


		    /**
		     * Sets the seed of this generator using one long, running that through LightRNG's algorithm twice to get the state.
		     * @param seed the number to use as the seed
		     */
		    public void setSeed(final long seed) {

		        long state = seed + 0x9E3779B97F4A7C15L,
		                z = state;
		        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
		        state0 = z ^ (z >>> 31);
		        state += 0x9E3779B97F4A7C15L;
		        z = state;
		        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
		        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
		        state1 = z ^ (z >>> 31);
		    }		
		
			public double nextGaussian(){
				double x1, x2, w, y1;
				 
		         do {
		        	 x1 = 2.0 * nextDouble() - 1.0;
		        	 x2 = 2.0 * nextDouble() - 1.0;
		        	 w = x1 * x1 + x2 * x2;
		         } while ( w >= 1.0 );

		         //w = StrictMath.sqrt( (-2.0 * StrictMath.log( w ) ) / w );
		         w = Math.sqrt( (-2.0 * Math.log( w ) ) / w );
		         y1 = x1 * w;
		         return y1;
			}
			

	/**
	 * A port of Blackman and Vigna's xoroshiro 128+ generator; should be very fast and produce high-quality output.
	 * Testing shows it is within 5% the speed of LightRNG, sometimes faster and sometimes slower, and has a larger period.
	 * It's called XoRo because it involves Xor as well as Rotate operations on the 128-bit pseudo-random state.
	 * <br>
	 * Machines without access to efficient bitwise rotation (such as all desktop JREs and some JDKs run specifying the
	 * {@code -client} flag or that default to the client VM, which includes practically all 32-bit Windows JREs but almost
	 * no 64-bit JREs or JDKs) may benefit from using XorRNG over XoRoRNG. LightRNG should continue to be very fast, but has
	 * a significantly shorter period (the amount of random numbers it will go through before repeating), at
	 * {@code pow(2, 64)} as opposed to XorRNG and XoRoRNG's {@code pow(2, 128)}, but LightRNG also allows the current RNG
	 * state to be retrieved and altered with {@code getState()} and {@code setState()}. For most cases, you should decide
	 * between LightRNG and XoRoRNG based on your needs for period length and state manipulation (LightRNG is also used
	 * internally by almost all StatefulRNG objects).
	 * <br>
	 * Original version at http://xoroshiro.di.unimi.it/xoroshiro128plus.c
	 * Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
	 *
	 * @author Sebastiano Vigna
	 * @author David Blackman
	 * @author Tommy Ettinger
	 */
	
	
}
