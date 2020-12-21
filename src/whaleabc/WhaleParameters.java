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
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.Random;
//import org.whale.ZigguratNormalizedGaussianSampler;

public class WhaleParameters implements Serializable {


    
    private long DOUBLE_MASK = (1L << 53) - 1;
    private double NORM_53 = 1. / (1L << 53);
    private long state0, state1;
    private double nextNextGaussian;
    private boolean haveNextNextGaussian = false;
    
    
    private ZigguratNormalizedGaussianSampler zng=new ZigguratNormalizedGaussianSampler(this);
    
        int hemisphere=0;
    
        int minThemes=2;
        int maxThemes=11;
    
    
	int modelType=1;
	int nYears=10;
        int runinperiod=0;
	int epochsperyear=10;
        int loglength=11;
        
	double mortalityRate=0.4;
	
	double betaa=1;
	double betab=3;
	double neighthresh=1.1;
        
        int[] repSizes= {1,2,2,2,2,3,3,3,4,4,5};
        //int maxRep=1; //Whale 
        
	double confBias=1.5;
        double novbias=1;
	public int[] popsizes={1000,1000,1000,1000,1000,1000,1000,1000,1000,1000};
        public int[] minpops, kpops;
        double popgrowth=0.05;
	public int popSize;
	int sylsPerSong=11;
        int numdims=2;
        int memorylength=50;
        int memorysize=numdims*sylsPerSong*memorylength;
        int ntutors=10;
	double mutationVar=0.08;
	int maxlearn=10;
        int yearsThreshold=0;
	double typeThresh=0.1;
	double problearn1=0;
        double problearn2=0;
        
        double dropparam=0;
        double probadd=0.01;
	
	float[] songBuffer;
	int[] songFreq;
	double[] cumFreq, freq;
	double[] powerLookUp, powlookUp;
	double[] tutsongsim;
	
	
	public WhaleParameters(){
		long seed=System.currentTimeMillis();
		setUp(seed);
	}
	
	public WhaleParameters(double[] p, double[] q, int[] r, int[] s, int[]t, long seed){

                popsizes=r;
                minpops=s;
                kpops=t;
                popSize=0;
                for (int i=0; i<t.length; i++){
                    popSize+=kpops[i];
                }
                //System.out.println(popSize);
		nYears=(int)q[0];
                
                yearsThreshold=(nYears/maxlearn)-11;
                
                
		sylsPerSong=(int)q[1];
                //System.out.println(sylsPerSong + " = sylsPerSong parameters");
                numdims=(int)q[2];
		typeThresh=q[3];
                problearn2=q[4];
                popgrowth=q[5];
                hemisphere=(int)q[6];
                
		mutationVar=p[0];
		novbias=p[1];
                problearn1=p[2];
                ntutors=(int)p[3];
                dropparam=p[4];
                
                probadd=Math.exp(-1*p[5]*dropparam);         //THIS IS A STATISTICAL SHORTCUT!!! I FOUND THAT THIS IS A WAY TO GENERATE SONGS OF THE RIGHT LENGTH
                
                //probadd=p[5];
                //System.out.println(mutationVar+"" +novbias+" "+problearn1+" "+ntutors);
                
                
		setUp(seed);
		
	}
	
	
	public WhaleParameters(WhaleParameters p){
		
		setUp(p.nextLong());
	}
	
	public void setUp(long seed) {
		//System.out.println("WhaleParameters.setup");
		songBuffer=new float[maxThemes*20*sylsPerSong*numdims];
                tutsongsim=new double[100];
		cumFreq=new double[maxThemes*20];
                powerLookUp=new double[30];
                //System.out.println(powerLookUp.length);
		for (int i=0; i<powerLookUp.length; i++){
			powerLookUp[i]=Math.pow(i, novbias);
		}
		
                //powerLookUp2=new double[100000];
                //for (int i=0; i<powerLookUp2.length; i++){
                //    powerLookUp2[i]=Math.pow(20*i/100000.0, novbias);
                //}
                
                
                makePowerLookUp();
		setSeed(seed);

	}
        
        
        public void makePowerLookUp(){
            powlookUp=new double[500];
            for (int i=0; i<500; i++){
                powlookUp[i]=(1/(1+ Math.pow(0.5,(i+dropparam))));
            }   
        }
        
        public void setPopulationSizes(int[][] pop){
            /*
            nx=pop[0][0];
            ny=pop[0][1];
            int rl=(pop.length-1)/2;
            removes=new int[rl][2][2];
            for (int i=0; i<rl; i++){
                removes[i][0][0]=pop[1+(i*2)][0];
                removes[i][0][1]=pop[1+(i*2)][1];
                removes[i][1][0]=pop[2+(i*2)][0];
                removes[i][1][1]=pop[2+(i*2)][1];
            }
            
            popSize=nx*ny;
            */
        }
	
	
	public void reportParameters(){
		
		DecimalFormat myFormat=new DecimalFormat("0.000");
		
		//double[] v=getPVec(modelType);
		//for (int i=0; i<v.length; i++){
		//	System.out.print(myFormat.format(v[i])+" ");
		//}
		//System.out.println();
	}
	
	public long rotateLeft(long x, int amt) {
		return (x << amt) | (x >>> 64-amt);		
	}
	
	public long nextLong() {
	       final long s0 = state0;
	       long s1 = state1;
	       final long result = s0 + s1;

	       s1 ^= s0;
		   
		   state0 = (s0 << 55 | s0 >>> 9)^ s1 ^ (s1 << 14); // a, b
		   state1 = (s1 << 36 | s1 >>> 28); // c

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
                    
                    public float nextFloat(){
                        return (float)(nextDouble());
                    }
                    
                    public float nextFloat(float min, float max) {
                    Random rand = new Random();
                    return rand.nextFloat() * (max - min) + min;

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
                        return zng.sample();
                    }
		    
		    public double nextGaussian2() {
		    	if (haveNextNextGaussian) {
		    		haveNextNextGaussian = false;
		    		return nextNextGaussian;
		    	} else {
		    	     double v1, v2, s;
		    	     do {
		    	       v1 = 2 * nextDouble() - 1;   // between -1.0 and 1.0
		    	       v2 = 2 * nextDouble() - 1;   // between -1.0 and 1.0
		    	       s = v1 * v1 + v2 * v2;
		    	     } while (s >= 1 || s == 0);
		    	     //double multiplier = StrictMath.sqrt(-2 * StrictMath.log(s)/s);
		    	     double multiplier = Math.sqrt(-2 * Math.log(s)/s);
		    	     nextNextGaussian = v2 * multiplier;
		    	     haveNextNextGaussian = true;
		    	     return v1 * multiplier;
		    	}
		    }
		    
			/*public double nextGaussian(){
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
			}*/
			

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