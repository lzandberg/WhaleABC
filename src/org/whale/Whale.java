/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.whale;

/**
 *
 * @author lieszandberg
 */
import java.util.Random;

public class Whale {
	
  
        boolean verbose=true;
	WhaleParameters param;
	WhaleIndividual[] inds;
	WhalePopulation population;
        WhaleMeasureStatistics stats;
	
	public int[][]locs;
	public int[] repSizes;
        public int[][] popSizes;
	public int iter;
        public int counter;
        public int maxlearn = 10; //total number of learning instances per year
	public int[][][][] out;
        public double[] means, sds;
        public double[][]loadings;
       
	
	public String fileLocation="/data/home/btw774/";
	
	public String[] paramNames={ "MutationVar", "RecombRate", "Gamma1", "Gamma2", "ConfBias"};
                
	public String[] statNames= {"ModelType", "NYrs", "nx", "ny", "MortalityRate", "Sylls", 
			"MutationVar", "RecombRate", "Gamma1", "Gamma2", "NeighbThresh", "TypeThresh", 
			"ConfBias", "MinDist", "NTypes"};
	public String[] measNames= {"AverageShare", "Sing", "Rare", "Int", "Common", "MaxF", "h1", 
			"Alpha1", "FitA", "FitB", "Share1", "Share2", "LAverageShare", "LSing", "LRare", "LInt", 
                        "LCommon", "LMaxF", "Lh1", "LAlpha1", "LFitA", "LFitB",  "LShare1", "LShare2", 
			"JaccAll", "Jacc1", "Jacc2", "Jacc3", "Jacc4", "Jacc5", 
			"WTDAll", "WTD1", "WTD2", "WTD3", "WTD4", "WTD5", "SylAverageSim", "SylSingletons", "SylRares", 
			"SylInts", "SylCommons", "SylMaxFreq", "Sylh1", "SylShZero", "SylShSingleton", "SylShMults",
			"SylShAll", "SylShMean", "SylRatio"};

	

	//long t1, t2, t3,t4, t5, t6, it1, it2, it3, it4, it5;
	
        public Whale (double[] x, double[] y, int[] z, int[]w, long seed){ 
                System.out.println("Whale");
                //Random random=new Random(System.currentTimeMillis());
            	
            	param=new WhaleParameters(x, y, z, w, seed);
		int[] rs={1};
		param.setRepertoireSizes(rs);
                int[][]ps={{1000, 10}};
		param.setPopulationSizes(ps);
		int [][] locs={{0,0}};
                runSimulation();
                population.makeEmpPop();
                //System.out.println("completed");
                WhaleMeasureStatistics ms=new WhaleMeasureStatistics(population, param);
                //out=ms.output;
                ms=null;
                population=null;
                inds=null;
                
        }
	
	public Whale(String fileLocation, long seed, int[][] locs, int[] rep, int[][] popSizes, double[] p, double[] q, int[] r, int[]s){
		this.fileLocation=fileLocation;
		param=new WhaleParameters(p,q, r, s, seed);
		param.setRepertoireSizes(rep);
                param.setPopulationSizes(popSizes);
                this.inds=inds;
		this.locs=locs;
		//this.p=p;
		//ed=new EmpData(param);
		
		runSimulation();
		//System.out.println("Measuring stats");
		//MeasureStatistics ms=new MeasureStatistics(population, param, ed);
                //WhaleMeasureStatistics ms=new WhaleMeasureStatistics(population, param);
                //stats.calculateThresholdSpectrum(inds);
	}
        	
	public Whale(String fileLocation, double[] p, double[] q, int[]r, int[]s, long seed) {
		this.fileLocation=fileLocation;
		EmpData ed=new EmpData(fileLocation);
		locs=ed.locs;
		repSizes=ed.repSizeD;
                popSizes=ed.popSizes;
		param=new WhaleParameters(p, q, r, s, seed);
		param.setRepertoireSizes(repSizes);
                param.setPopulationSizes(popSizes);
		WhaleMeasureStatistics ms=new WhaleMeasureStatistics(population, ed, param);
		//stats=ms.out;
		//for (int i=0; i<stats.length; i++) {System.out.println(measNames[i]+" "+stats[i]);}
                
                //System.out.println(ed.empDiss.length+" "+ed.empSylDiss.length);
                
	}
        
	
	
	public void runSimulation(){
                //System.out.println("runSimulation");
		long a=System.currentTimeMillis();

		//System.out.println("Initiating...");
		initiateSimulation();
		//long b=System.currentTimeMillis();
		int nYears=param.nYears;
		//System.out.println("Initiation took: "+(b-a));
		//System.out.println("Simulation running for "+nYears);
		int j=0;
                counter=0;
		/*
		t1=0;
		t2=0;
		t3=0;
		t4=0;
		it1=0;
		it2=0;
		it3=0;
		it4=0;
		it5=0;
		*/
		for (int i=0; i<nYears; i++){
                    //if(i%100==0){
                    if (verbose){
                        System.out.println("iterateSimulation = " + i + "/" + nYears);
                    }
                    //}
			/*
			if (j==100) {
				j=0;
                                MeasureStatistics ms=new MeasureStatistics(population, param);
                                
                                for (int k=0; k<ms.out.length; k++){
                                    System.out.print(ms.out[k]+" ");
                                }
                                System.out.println();
                                
                                
                                
				//System.out.println("Year: "+i);
			}
			j++;
                        */
                        
                        iterateSimulation(counter);
                        counter++;
                        if(counter>(maxlearn-1)){
                            counter=0;
                        }
			
		}
               
                
                //population.printrec();
		/*
		double t1x=(t3-t2)*0.000000001;
		double t2x=(t4-t3)*0.000000001;
		double it1x=it1*0.000000001;
		double it2x=it2*0.000000001;
		double it3x=it3*0.000000001;
		double it4x=it4*0.000000001;
		double it5x=it5*0.000000001;
		
		
		System.out.println(t1x+" "+t2x+" "+it1x+" "+it2x+" "+it3x+" "+it4x+" "+it5x);
		*/
		//long c=System.currentTimeMillis();
		
		//long d=System.currentTimeMillis();
		//System.out.println((b-a)+" "+(c-b)+" "+(d-c));

	}
        
	
	public void initiateSimulation(){
		
		//System.out.println("initiateSimulation");
		inds=new WhaleIndividual[param.popSize];
		
		//System.out.println("Making individuals");

		for (int i=0; i<param.popSize; i++){
			inds[i]=new WhaleIndividual(i, param, -1);
		}
		
		//System.out.println("Making population structure");
		population=new WhalePopulation(inds, param, locs);
		
		for (int i=0; i<param.popSize; i++){
			inds[i].setPopulation(population); 
                        inds[i].initiate(); //initiate repertoire (newRep) per subpop, updateRepert, and initiate memory
		}
		//int[] x=population.calculateEmpIDs();
                
                for (int i=0; i<inds.length; i++){ 
                    inds[i].mortalityRate=1;
		}
                
                for (int i=0;  i<param.runinperiod; i++){
                    for (int j=0; j<inds.length; j++){
			inds[j].mortality();
                    }
                    for (int j=0; j<inds.length; j++){
			inds[j].learnSongs();
                    }
		}
                    
                for (int i=0; i<inds.length; i++){
                    inds[i].mortalityRate=param.mortalityRate;
                    inds[i].mutationVariance=param.mutationVar;
                    //inds[i].recombinationRate=param.recomRate; //might want to change this?
		}
                
		//System.out.println("CHECK1: "+x.length);
		//System.out.println("Population initiated");
	}
	
	public void iterateSimulation(int counter){

		/* No mortality in the Whales..?
                
		for (int i=0; i<inds.length; i++){
			inds[i].mortality();
		}*/ 
                //System.out.println("counter = " + counter);
                population.setSeason(counter);
		//System.out.println("Song learning...");
                
                population.learnSongs();
		
                
                //for (int i=0; i<inds.length; i++){
		//	inds[i].learnSongs();
			/*double[] dt=
			
			if (dt!=null) {
				it1+=dt[0];
				it2+=dt[1];
				it3+=dt[2];
				it4+=dt[3];
				it5+=dt[4];
			}
			*/
		//}
		
		//t3+=System.nanoTime()-t1;
		//System.out.println("Repertoire update...");
		//for (int i=0; i<inds.length; i++){
		//	inds[i].updateRepertoire();
		//}
                
                
                if (verbose){
                for (int i=0; i<population.subpopsize.length; i++){
                    int ii=population.subpopstarts[i];
                    System.out.print(i+" "+ii+" "+counter+" "+population.currentpopsizes[i]+" "+population.breeding+" "+" "+population.problearn1+" "+param.mutationVar+" "+param.novbias+" "+param.ntutors+" ");
                    for (int j=0; j<inds[ii].ns; j++){
                        System.out.print(inds[ii].newRepertoire[j]+" ");
                    }
                    System.out.println();
                    ii+=1;
                    System.out.print(i+" "+ii+" "+counter+" "+population.currentpopsizes[i]+" "+population.breeding+" "+" "+population.problearn1+" "+param.mutationVar+" "+param.novbias+" "+param.ntutors+" ");
                    for (int j=0; j<inds[ii].ns; j++){
                        System.out.print(inds[ii].newRepertoire[j]+" ");
                    }
                    System.out.println();
                    ii+=1;
                    System.out.print(i+" "+ii+" "+counter+" "+population.currentpopsizes[i]+" "+population.breeding+" "+" "+population.problearn1+" "+param.mutationVar+" "+param.novbias+" "+param.ntutors+" ");
                    for (int j=0; j<inds[ii].ns; j++){
                        System.out.print(inds[ii].newRepertoire[j]+" ");
                    }
                    System.out.println();
                    
                }
                }
		//t4+=System.nanoTime()-t1;

	}
		
		
		

		public static void main (String args[]) {
                    //System.out.println("mainWhale");
                    //new Whale();
                   
                    /*
			String fileLocation="/home/rflachlan/Dropbox/ChaffMainlandN/";
                        Priors p=new Priors(System.currentTimeMillis());
                        double[] x=p.sampleFromPriors();
			new ChaffinchABC(fileLocation, x, p.variables, System.currentTimeMillis());
			*/
		}
		

	}