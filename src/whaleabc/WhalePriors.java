/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package whaleabc;

/**
 *
 * @author robertlachlan
 */

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

public class WhalePriors {
    
        int mode=1; //mode=0 means that parameters are sampled from previous round and error added - as per normal ABC procedure
                    //mode=1 means that parameters are sampled but no error added. If you want to re-run models.
    
    
        int hemisphere=0; //hemisphere: 0=SH, 1=NH
        double[] means={14.8012164243119,8.2719274768121,0.341104619113366,0.022167776346596,0.18404937278838,0.174967906166869};
        double[] sds={8.28742378754896,1.13837582411023,0.187320031596461,0.052091076961061,0.126785416689095,0.100758794774356};
        //double[] sds={8.28742378754896,1.13837582411023,0.187320031596461,0.026045538480531,0.126785416689095,0.100758794774356};

        double[]empstats={0.128205128,0.925579976,1.666666667,0.934126984,6.576846453,1.211447786,1,2,0.593434343,0.648584299,
            0.377083333,0.971367521,0.212365591,0.212365591,0.443357868,0.450983809,0.233644444,0.786648079,0.118510965,0.413398407,0.324846903,0.553003635,0.443357868};

        int[] statindices={2, 4,17,20,21,22};
    

	private static final long DOUBLE_MASK = (1L << 53) - 1;
        private static final double NORM_53 = 1. / (1L << 53);
        private long state0, state1;

	int nParams;

        //public double[] variables=new double[] {1000, 8, 10, 0.1, 0, 0.06};
        //public double[] variables=new double[] {500, 11, 10, 2, 0, 0.06}; //FOR NORMAL SIMULATION
        
        public double[] variables=new double[] {500, 11, 10, 0.002, 0, 0.06, hemisphere}; //FOR SIMPLE SIMULATION!!!
        //public double[] variables=new double[] {5000, 11, 10, 0.002, 0, 0.06, hemisphere}; //FOR LONG SIMULATION!!!

        // SH pops are ordered W to E starting with S American Atlantic population
        public int[] popsizessh={6000,6000, 240, 4000,4000,10000,10000,1125,750,375,1000,5000};     //RIGHT WAY ROUND!!
        public int[] minpopssh={250, 750, 35, 350, 950, 400, 100, 23, 16, 7, 20, 350};
        public int[] kpopssh={12000, 9000, 2000, 4500, 4500, 11000, 13000, 2500, 1600, 800, 2150, 6000};
              
        
        //NH pops are ordered: Okinawa, Hawaii, Cen Am, Mexico, Caribbean, Cape Verde
        
        public int[] popsizesnh={750,1750, 5500, 325,7370, 130};     //RIGHT WAY ROUND!!
        public int[] minpopsnh={20, 30, 400, 10, 1000, 10};
        public int[] kpopsnh={5308, 12385, 38926, 2300, 52160, 920};
        
       
       
        
        //public int[] popsizes={5000,1000,375,750,1125,10000,10000,4000,4000,240,6000,6000, 6000,6000};
        //public int[] minpops={2350,20,7,16,23,100,400,950,350,35,750,250};
        //public int[] kpops={6000,2150,800,1600,2500,13000,11000,4500,4000,2000,9000,12000};
        
        
/* Should be these numbers (with population BSB2 included, but not the smaller subpops of Oceania)
        public int[] popsizes={6000,6000,240,4000,4000,10000,10000,3000,5000};
        public int[] minpops={250, 750, 35, 350, 950, 400, 100, 65, 350};
        public int[] kpops={12000, 9000, 2000, 4000, 4500, 11000, 13000, 7000, 6000};
*/        
        
	int[] priorType= new int[]{1, 2, 1, 1, 2, 2};
	//double[] priorMax=new double[] {0.1, 2, 0.05, 20, -1, 0.001};//Normal priors
	//double[] priorMin=new double[] {0.002, -0.5, 0.001, 3,  8, 0.1};
        
        //double[] priorMax=new double[] {0.001, 2, 0.05, 20, 2, 0.02};   //simple priors
	//double[] priorMin=new double[] {0.00001, -0.5, 0.0001, 3,  12, 0.0001};
        
        double[] priorMax=new double[] {0.00001, 2, 0.05, 20, 14, 1.5};   //simple priors
	double[] priorMin=new double[] {0.000000001, -0.5, 0.0001, 3,  6, 0.75};
        
        //double[] priorMax=new double[] {0.000000069, 1.7, 0.00032, 7.7, 10.07, 1.25};   //best fit priors
	//double[] priorMin=new double[] {0.000000068, 1.6, 0.00031, 7.69,  10.06, 1.245};

        //double[] priorMax={1.3846435725525E-07,1.36735620892202,0.004020576690663,3.41342995831741,12.18714946223,1.18347948327436};//best fit priors 2
        //double[] priorMin={1.38386435725525E-07,1.36535620892202,0.004018076690663,3.41330995831741,12.18604946223,1.18330948327436};
        
        
        //double[] priorMax=new double[] {0.0000069, 1.7, 0.00032, 7.7, 10.07, 1.25};   //simple priors
	//double[] priorMin=new double[] {0.0000068, 1.6, 0.00031, 7.69,  10.06, 1.245};
        
        //double[] priorMax=new double[] {0.00334, 1.614,0.00175,13, -5, 0.001};
	//double[] priorMin=new double[] {0.00333, 1.613,0.00174, 12, 5, 0.1};
        
        
        public WhalePriors(){
            nParams=priorType.length;
        }
       
	public WhalePriors(long seed) {
		setSeed(seed);
		nParams=priorType.length;
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
			//System.out.println(parameters[i]);
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
		double p=nextDouble()*q;
		return Math.exp(p+Math.log(priorMin[i]));
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
            
            if (mode==1){
                boolean found=false;
                double[] v2=new double[nParams];
                for (int i=0; i<nParams; i++) {
                    if (priorType[i]==1) {
                        v2[i]=Math.exp(v[i]);
                    }
                    else {
                        v2[i]=v[i];
                    }
                }
                return v2;
            }
            
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