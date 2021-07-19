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


/**
 *
 * @author lieszandberg
 */
import java.util.LinkedList;
import java.util.Random;

public class Whale {
	
    //boolean breakSimulation=false;
        boolean verbose=false;
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
       
	//LinkedList<double[]> resultstore;
        
	public String fileLocation="/data/home/btw774/";
	
	public String[] paramNames={ "MutationRate", "NovBias", "PLearnNeighbour", "NTutors", "DropParam", "AddParam"};
                
	public String[] measNames= { "SongSingletons",  "SongMaxFreq",     "ThemeSingletons", "ThemeMaxFreq",    "MeanSongLength",  "PSDSongLength",   "NumMinLength",   
 "NumMaxLength",    "SharePL",        "SharePA",         "SharePR",         "ShareCL",        "ShareCA",         "ShareCR",        
"TsharePL",        "TsharePA",        "TsharePR",        "TshareCL",        "TshareCA",        "TshareCR",       "TShareDP",        "TTurnover",        "TMaxShareP" };

        WhaleMeasureStatistics wms;
	int sim=0;
        SimulationRunner sr;
        
        public Whale (double[] x, double[] y, int[] z, int[]w, int[] t, long seed, int sim, SimulationRunner sr){ 
                //System.out.println("Whale");
                //Random random=new Random(System.currentTimeMillis());
            	this.sim=sim;
                this.sr=sr;
                //if (verbose){resultstore=new LinkedList<double[]>();}
            	param=new WhaleParameters(x, y, z, w, t, seed);
                int[][]ps={{1000, 10}};
		param.setPopulationSizes(ps);
		int [][] locs={{0,0}};
                runSimulation();
                //population.makeEmpPop();
                //System.out.println("completed");
                wms=new WhaleMeasureStatistics(population, param);
                //out=ms.output;
                population=null;
                inds=null;
                
        }
	
	public Whale(String fileLocation, long seed, int[][] locs, int[] rep, int[][] popSizes, double[] p, double[] q, int[] r, int[]s, int[] t){
		this.fileLocation=fileLocation;
		param=new WhaleParameters(p,q, r, s, t, seed);
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
                        iterateSimulation(counter, i);
                        counter++;
                        if(counter>(maxlearn-1)){
                            counter=0;
                        }
                    //if (breakSimulation){
                    //    i=nYears;
                        
                    //}
			
		}
                
                
                //if (!breakSimulation){
                //    for (double[] x : resultstore){
                //        for (int i=0; i<x.length; i++){
                //            System.out.print(x[i]+" ");
                //        }
                //        System.out.println();  
                //    }
                //}
                
                //breakSimulation=false;
               
               
	}
        
	
	public void initiateSimulation(){
		
		//System.out.println("initiateSimulation");
		inds=new WhaleIndividual[param.popSize];
		
		for (int i=0; i<param.popSize; i++){
			//inds[i]=new WhaleIndividual(i, param);
                        inds[i]=new SimpleIndividual(i, param);
		}
		
		population=new WhalePopulation(inds, param, locs);
		
		for (int i=0; i<param.popSize; i++){
			inds[i].setPopulation(population); 
                        inds[i].initiate();
		}

                for (int i=0;  i<param.runinperiod; i++){
                    for (int j=0; j<inds.length; j++){
			inds[j].learnSongs(true);
                    }
		}
                    
                for (int i=0; i<inds.length; i++){
                    inds[i].mutationVariance=param.mutationVar;
		}
	}
	
	public void iterateSimulation(int counter, int p){

                population.setSeason(counter);                
                population.learnSongs();
                
                if (verbose&&(counter==9)){
                    if (p>20){
                        WhaleMeasureStatistics ms=new WhaleMeasureStatistics(population, param, true);
                        System.out.println(sr.id+" "+sim+" "+((p-9)/10)+" "+ms.out4[0]+" "+ms.out4[1]+" "+ms.out4[2]+" "+ms.out4[3]);
                        //double[]x={sr.id,sim,((p-9)/10),ms.out4[0],ms.out4[1],ms.out4[2],ms.out4[3]};
                        //resultstore.add(x);
                        //for (int i=0; i<ms.out4.length; i++){
                          //  System.out.print(ms.out4[i]+" ");
                        //}
                        /*
                        for (int i=0; i<ms.out1.length; i++){
                            System.out.print(ms.out1[i]+" ");
                        }
                        for (int i=0; i<ms.out2.length; i++){
                            System.out.print(ms.out2[i]+" ");
                        }
                        for (int i=0; i<ms.out3.length; i++){
                            System.out.print(ms.out3[i]+" ");
                        }
                        */
                        //System.out.println();
                        //System.out.println("Sharing: "+p+" "+ms.out[0]+" "+ms.out[1]+" "+ms.out[2]);
                    }
                    //System.out.println(p);
                    /*if (p==499){
                        WhaleMeasureStatistics ms=new WhaleMeasureStatistics(population, param);
                        double[] empstats=sr.empstats;
                        double[] simstats=new double[sr.p.statindices.length];
                        for (int j=0; j<sr.p.statindices.length; j++){
                            simstats[j]=ms.out[sr.p.statindices[j]];
                        }
                        double score=sr.calculateDifference(empstats, simstats);
                        if (score>2){
                            //breakSimulation=true;
                            System.out.println(sr.id+" Reject run "+score);
                        } 
                        else{
                            System.out.println(sr.id+" Accept run "+score);
                            ms.traceThemes(sr.id);
                        }
                    }*/
                    /*
                for (int i=0; i<population.subpopsize.length; i++){
                    int ii=population.subpopstarts[i];
                    System.out.print(p+" "+i+" "+ii+" "+counter+" "+population.currentpopsizes[i]+" "+population.breeding+" "+" "+population.problearn1+" "+param.mutationVar+" "+param.novbias+" "+param.ntutors+" ");
                    for (int j=0; j<inds[ii].ns; j++){
                        System.out.print(inds[ii].newRepertoire[j]+" ");
                    }
                    System.out.println();
                    
                    /*
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
                    */
                    
                //}
                
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
