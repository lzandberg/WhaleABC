package ABCRunner;
import SimulationRunner.SimulationRunner;
import org.whale.WhalePriors;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

import org.whale.*;
import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

public class ABCRunner {
	
    LinkedList<long[]> records=new LinkedList<long[]>();
	
	double[][] pvals1, pvals2;
	double[] weights, cumweights;
	double[][] cov;

        String[] paramNames, measNames;
        WhalePriors priors=new WhalePriors();
	
	//int numReps=5000;
	int nsamps=1280; //1000
        
        //double[] thresholds={15, 11, 10, 9, 8, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3, 2.5, 2};
        
        //double[] thresholds={7, 6, 5.5, 5, 4.5,4.25, 4,3.75, 3.5, 3, 2.5, 2};
	
        double[] thresholds={8, 7, 6, 5.5, 5, 4.5, 4, 3.5, 3.25, 3,2.75, 2.5, 2.25, 2.125, 2, 1.875, 1.75, 1.625, 1.5, 1.375, 1.25, 1.125, 1, 0.875, 0.75};
        
        int[] selected;
	//public String fileLocation="/data/home/btw774/Chaff/";
	//public String fileLocation="/users/Rob/Desktop/";
	//public String fileLocation="/home/rflachlan/Dropbox/ChaffTenerife/";
        //public String fileLocation="/home/rflachlan/Dropbox/ChaffMainlandN/";
        //public String fileLocation="/home/rflachlan/Dropbox/ChaffTentsmuir/";
        /*public String[] fileLocation={"/home/rflachlan/Dropbox/ChaffinchFaial/",
            "/home/rflachlan/Dropbox/ChaffFlores/",
            "/home/rflachlan/Dropbox/ChaffinchPico/",
            "/home/rflachlan/Dropbox/ChaffMiguel/"
        };*/
       // public String[] fileLocation={"/home/rflachlan/Dropbox/ChaffinchFaial/"};
        //public String[] fileLocation={"/home/rflachlan/Dropbox/ChaffTenerife/"};
        //public String[] fileLocation={"/home/rflachlan/Dropbox/ChaffTentsmuir/"};
        //public String[] fileLocation={"/home/rflachlan/Dropbox/ChaffGranCan/"};
        //public String[] fileLocation={"/home/rflachlan/Dropbox/ChaffBlue/"};
        public String[] fileLocation={"/Users/lieszandberg/STACK/PostDoc QMUL/Phase III - Comparative study/Whales/Results"};
	
        int ncores=64;
	int np=0;
	int npc=0;
	
	LinkedList<threadRunner[]> trl=new LinkedList<threadRunner[]>();
	
	int currentRound=1; //10
	
	int nRounds=thresholds.length;
	
	double[][] stats, params;
	double[] scores;
	
	Random random=new Random(System.currentTimeMillis());
	
	double threshold=0;
        
        Whale[] cemp;
        int numpops;
        
        
	public ABCRunner() {

		np=ncores;
		//npc=(int)(Math.ceil(numReps/(np+0.0)));
		System.out.println("HERE WE ARE! "+npc);
		//numReps=npc*np;
                
                npc=50; 
                numpops=fileLocation.length;
                
                //ChaffinchABC ct=new ChaffinchABC();
                
                
                
                WhalePriors p=new WhalePriors(random.nextLong());
                double[] x=p.sampleFromPriors();
                cemp=new Whale[numpops];
                for (int i=0; i<numpops; i++){
                    //cemp[i]=new Whale(fileLocation[i], x, p.variables, System.currentTimeMillis()); Applicable for ChaffinchABC, not Whale?
                }
                this.paramNames=cemp[0].paramNames;
                this.measNames=cemp[0].measNames;
                
                if (currentRound>0){
                    readRound(fileLocation[0]);
                }
                /*
                tr=new threadRunner[ncores];
                for (int i=0; i<32; i++){
                    tr[i]=new threadRunner(npc, i, 2);
                    tr[i].readRound();
                }                
                readRound(fileLocation);
                collateResults();
		selectBestRuns();
		weights=calculate_Weights(weights, cov, pvals1, pvals2);
                cumweights=calculateCumWeights(weights);
		cov=calculate_CovMat(weights, pvals2);
		updatePVals();	    
                writeRound();
                writeResults(3);
                */
               
                
		for (int i=0; i<nRounds; i++) {
			iterateRound();
		}		
      
	}
        
        
        public ABCRunner(int a){
            System.out.println("ABCRUNNER");
            np=ncores;
            npc=200; 
            iterateRound(0);
        }
        
        public void writeRound(String fl){
            	DocumentSave ds=new DocumentSave(fl+"/cvmat.csv", ", ");
		
		for (int i=0; i<cov.length; i++){
                    for (int j=0; j<cov[i].length; j++){
                        ds.writeDouble(cov[i][j]);
                    }
                    ds.writeLine();
		}
		ds.finishWriting();
                
                ds=new DocumentSave(fl+"/pvals.csv", ", ");
		
		for (int i=0; i<pvals2.length; i++){
                    for (int j=0; j<pvals2[i].length; j++){
                        ds.writeDouble(pvals2[i][j]);
                    }
                    ds.writeDouble(cumweights[i]);
                    ds.writeLine();
		}
		ds.finishWriting();       
        }
        
        public void writeNews(String fl, int a, int b, int c){
            /*
            LinkedList<String[]> locs=new LinkedList<String[]>();
            try{
                File file=new File(fl+"/newsupdate.csv");
                
                String cvsSplitBy = ",";
                BufferedReader reader=new BufferedReader(new FileReader(file));
                String line=null;
                while((line=reader.readLine())!=null){
				
                    String[] s=line.split(cvsSplitBy);
                    locs.add(s);
                }
			
                reader.close();
            }
            catch(Exception e){}
            */
            
            long [] newresults=new long[4];
            newresults[0]=a;
            newresults[1]=b;
            newresults[2]=c;
            newresults[3]=System.currentTimeMillis();
            records.add(newresults);
            
            try{
            	DocumentSave ds=new DocumentSave(fl+"/newsupdate.csv", ", ");
		
                for (long[] s: records){
                    for (int i=0; i<s.length; i++){
                        ds.writeLong(s[i]);
                    }
                    ds.writeLine();
                }
                ds.finishWriting(); 
                /*
		ds.writeInt(a);
                ds.writeInt(b);
                ds.writeInt(c);
                ds.writeLong(System.currentTimeMillis());
                ds.writeLine();
		    
                */
            }
            catch(Exception e){e.printStackTrace();}
        }
	
	
	public void iterateRound(){

		iterateRound(currentRound);
	
		collateResults();
		selectBestRuns();
		
		if (currentRound==0){
			weights=calculateWeightsRoundOne(nsamps);
			
		}
		else{
			weights=calculate_Weights(weights, cov, pvals1, pvals2);
			
		}
		cumweights=calculateCumWeights(weights);
                
                
		cov=calculate_CovMat(weights, pvals2);
		updatePVals();	
                
                writeRound(fileLocation[0]);
		
		currentRound++;
                writeResults(fileLocation[0], currentRound);
		System.out.println("Finished Round: "+currentRound+" "+threshold);
		
	}
	
	
	public void iterateRound(int round) {
            System.out.println("IterateRound");
            int numpassed=0;
            
            trl=new LinkedList<threadRunner[]>();
                        
            while (numpassed<nsamps){
            
            
		threadRunner[] tr=new threadRunner[np];
		for (int cores=0; cores<np; cores++){
			tr[cores]=new threadRunner(npc, cores, round);
			tr[cores].setPriority(Thread.MAX_PRIORITY);
			tr[cores].start();
		}
		
		try{
                    for (int cores=0; cores<np; cores++){
			tr[cores].join();
                    }
		}
		catch (Exception f){
			f.printStackTrace();
		}
                
                
                for (int i=0; i<np; i++){
                    numpassed+=tr[i].numPassThresh;
                }
                trl.add(tr);
                
                writeNews(fileLocation[0], (currentRound+1), trl.size(), numpassed);
                
                System.out.println("NUM REPS: "+trl.size()+" "+numpassed);
                
            }      
	}
        
        public void readRound(String fileloc){
            System.out.println("Reading round...");
            File file=new File(fileloc+"/cvmat.csv");
            cov=null;
            try{
                String cvsSplitBy = ",";
                BufferedReader reader=new BufferedReader(new FileReader(file));
                String line=null;
			
                LinkedList<double[]> locs=new LinkedList<double[]>();
                while((line=reader.readLine())!=null){
				
                    String[] s=line.split(cvsSplitBy);
				
                    double[] t=new double[s.length-1];
                    for (int i=0; i<s.length-1; i++){
                        t[i]=Double.parseDouble(s[i]);
                    }
                    locs.add(t);
                }
			
                reader.close();
			
                cov=new double[locs.size()][];
                for (int i=0; i<cov.length; i++) {
                    cov[i]=locs.get(i);	
                }
            }
		
            catch(Exception e){
                e.printStackTrace();
            }
        
            file=new File(fileloc+"/pvals.csv");
            pvals1=null;
            weights=null;
            try{
                String cvsSplitBy = ",";
                BufferedReader reader=new BufferedReader(new FileReader(file));
                String line=null;
			
                LinkedList<double[]> locs=new LinkedList<double[]>();
                while((line=reader.readLine())!=null){
				
                    String[] s=line.split(cvsSplitBy);
				
                    double[] t=new double[s.length-1];
                    for (int i=0; i<s.length-1; i++){
                        t[i]=Double.parseDouble(s[i]);
                    }
                    locs.add(t);
                }
			
                reader.close();
			
                pvals1=new double[locs.size()][];
                weights=new double[locs.size()];
                cumweights=new double[locs.size()];
                for (int i=0; i<weights.length; i++) {
                    double[] x=locs.get(i);
                    pvals1[i]=new double[x.length-1];
                    System.arraycopy(x, 0, pvals1[i], 0, x.length-1);
                    cumweights[i]=x[x.length-1];
                    if (i==0){
                        weights[i]=cumweights[i];
                    }
                    else{
                        weights[i]=cumweights[i]-cumweights[i-1];
                    }
                }
            }
		
            catch(Exception e){
                e.printStackTrace();
            }
        
        
            for (int i=0; i<cov.length; i++){
                System.out.println("COV: "+i+" "+cov[i][0]+" "+cov[i].length);
            }
            for (int i=0; i<weights.length; i++){
                System.out.println("WEI: "+i+" "+weights[i]);
            }
        
        }
	
	public void collateResults() {
                int tn=trl.size();
		stats=new double[np*npc*tn][];
		params=new double[np*npc*tn][];
		scores=new double[np*npc*tn];
		
		int k=0;
                for (threadRunner[] tr : trl){
                    
                    for (int cores=0; cores<np; cores++){	
			for (int i=0; i<tr[cores].params.length; i++) {
				params[k]=tr[cores].params[i];
                                //System.out.println(cores+" "+i+" "+params[k].length);
				stats[k]=tr[cores].stats[i];
				scores[k]=tr[cores].score[i];
				k++;
			}
                    }
		}
	}
	
	public void selectBestRuns() {
		int n=scores.length;
		double[] x=new double[n];
                selected=new int[n];
		System.arraycopy(scores, 0, x, 0, n);
		Arrays.sort(x);
				
		threshold=x[nsamps];
		
		//System.out.println("SELECT RUNS: "+n+" "+threshold+" "+nsamps+" "+x[0]);
		int ns=0;
                for (int i=0; i<scores.length; i++){
                    if (scores[i]<threshold){
                        ns++;
                    }
                }
                
		pvals2=new double[ns][];
		
		int k=0;
		for (int i=0; i<n; i++) {		
			if (scores[i]<threshold) {		
				pvals2[k]=priors.transform(params[i]);
                                //System.out.println("TRANSFORM: "+pvals2[k][0]+" "+pvals2[k][1]);
                                selected[i]=1;
				k++;
			}
		}	
	}
	
	
	public void writeResults(String fl, int cr) {
		DocumentSave ds=new DocumentSave(fl+"results"+cr+".csv", ", ");
		
		for (int i=0; i<paramNames.length; i++) {
			ds.writeString(paramNames[i]);
		}
                for (int i=0; i<numpops; i++){
                    for (int j=0; j<measNames.length; j++) {
			ds.writeString(measNames[j]+"_"+(i+1));
                    }
		}
                ds.writeString("Score");
                ds.writeString("Selected");
                ds.writeString("Weight");
		ds.writeLine();
		
		//double[][] stats=new double[numReps][];
		//int a=0;
                int k=0;
                int kk=0;
		for (int i=0; i<params.length; i++){
                    for (int j=0; j<params[i].length; j++) {
			ds.writeDouble(params[i][j]);
                    }
                    for (int j=0; j<stats[i].length; j++) {
			ds.writeDouble(stats[i][j]);
                    }
                    ds.writeDouble(scores[i]);
                    ds.writeDouble(selected[i]);
                    if (selected[i]==1){
                        ds.writeDouble(weights[kk]);
                        kk++;
                    }
                    ds.writeLine();
		}
		ds.finishWriting();
		//long btime=System.currentTimeMillis();
	}
	
	
	class threadRunner extends Thread{
		int n, id, round;
		double[][] stats, params;
                double[] score;
                int numPassThresh=0;
		
		public threadRunner(int n, int id, int round) {
			this.n=n;
			this.id=id;
			this.round=round;
		}
		
		public synchronized void run(){
                    
                    //SimulationRunner sr=new SimulationRunner(cemp, numpops, n, id, round, random.nextLong(), fileLocation, cumweights, pvals1, cov);
                    
                    SimulationRunner sr=new SimulationRunner(n, id, round, random.nextLong());
                    sr.runSimulation();
                    
                    stats=sr.stats;
                    score=sr.score;
                    params=sr.params;
                    
                    //readRound();
                    numPassThresh=sr.passThresh(thresholds[round]);
			
		}
                
                public void readRound(){

                    File file=new File(fileLocation[0]+"/outputstats"+id+".csv");
                    stats=new double[n][];
                    score=new double[n];
                    try{
                        String cvsSplitBy = ",";
                        BufferedReader reader=new BufferedReader(new FileReader(file));
                        String line=null;
			int a=0;
                        while((line=reader.readLine())!=null){
				
                            String[] s=line.split(cvsSplitBy);
                            double[] t=new double[s.length-1];
                            for (int i=0; i<s.length-1; i++){
                                t[i]=Double.parseDouble(s[i]);
                            }
                            stats[a]=new double[t.length-1];
                            System.arraycopy(t, 0, stats[a], 0, t.length-1);
                            score[a]=t[t.length-1];
                            a++;
                        }
			
                        reader.close();
			
            
                    }
		
                    catch(Exception e){
                        e.printStackTrace();
                    }
        
                    file=new File(fileLocation[0]+"/outputparams"+id+".csv");
                    params=new double[n][];
                    try{
                        String cvsSplitBy = ",";
                        BufferedReader reader=new BufferedReader(new FileReader(file));
                        String line=null;
			int a=0;
                        while((line=reader.readLine())!=null){
				
                            String[] s=line.split(cvsSplitBy);
				
                            double[] t=new double[s.length-1];
                            for (int i=0; i<s.length-1; i++){
                                t[i]=Double.parseDouble(s[i]);
                            }
                            params[a]=new double[t.length];
                            System.arraycopy(t, 0, params[a], 0, t.length);
                            a++;
                        }
			
                        reader.close();
			
           
                    }
		
                    catch(Exception e){
                        e.printStackTrace();
                    }
        
                }	
	}
        
        public double[] calculateCumWeights(double[] w){
            int n=w.length;
            double[] out=new double[n];
            out[0]=w[0];
            for (int i=1; i<n; i++){
                out[i]=out[i-1]+w[i];
            }
            return out;
        }
	
	public double[] calculateWeightsRoundOne(int n){
		double out[]=new double[n];
		
		for (int i=0; i<n; i++){
			out[i]=1/(n+0.0);
		}	
		return out;
	}
	
	
	//weights=calculate_Weights(weights, cov, pvals1, pvals2);
	public double[] calculate_Weights(double[] pw, double[][] pc, double[][] thetp, double[][] thetc){
		
		int n=thetp.length;
	
		double out[]=new double[n];
		
		for (int i=0; i<n; i++){
			MultivariateNormalDistribution mnd=new MultivariateNormalDistribution(thetc[i], pc);					
			for (int j=0; j<n; j++){				
				double d=mnd.density(thetp[j]);
				d*=pw[j];
				out[i]+=d;			
			}
			out[i]=1/out[i];
		}
		
		double t=0;
		for (int i=0; i<n; i++){
			t+=out[i];
		}
		for (int i=0; i<n; i++){
			out[i]/=t;
		}
		
		return out;
	}
	
	public double[][] calculate_CovMat(double[] pw, double[][] thetas){
		
		int n=pw.length;
		int m=thetas[0].length;
		
		double s1=0;
		double s2=0;
		for (int i=0; i<n; i++){
			s1+=pw[i];
			s2+=pw[i]*pw[i];
		}
		double s3=s1*s1;
		
		double s4=s1/(s3-s2);
		
		//System.out.println("COVMAT S CALC: "+s1+" "+s2+" "+s3+" "+s4);
		
		
		double[] wm=new double[m];
		
		for (int i=0; i<m; i++){
			double t=0;
			for (int j=0; j<n; j++){
				wm[i]+=thetas[j][i]*pw[j];
				t+=pw[j];
			}
			wm[i]/=t;
		}
		
		
		double[][] out=new double[m][m];
		
		for (int i=0; i<m; i++){
			for (int j=0; j<m; j++){
				for (int k=0; k<n; k++){
					out[i][j]+=pw[k]*(thetas[k][i]-wm[i])*(thetas[k][j]-wm[j]);
				}
				out[i][j]*=2*s4;
				
				//System.out.print(out[i][j]+" ");
			}
			//System.out.println();
		}
	
		return out;		
	}
	
	public double[][] doubleCovMat(double[][] cov){
		int n=cov.length;
		double[][] cov2=new double[n][n];
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				cov2[i][j]=2*cov[i][j];
			}
		}
		
		return cov2;
	}
	
	public void updatePVals(){
		pvals1=new double[nsamps][];
		for (int i=0; i<pvals1.length; i++){
			pvals1[i]=new double[pvals2[i].length];
			System.arraycopy(pvals2[i], 0, pvals1[i], 0, pvals2[i].length);
		}
	}
	
	
	
	
	
	
	
	public static void main (String args[]) {
		new ABCRunner(1);
		
	}
	
	
}
