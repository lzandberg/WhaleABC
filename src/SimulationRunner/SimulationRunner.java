/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SimulationRunner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedList;
import org.whale.Whale;
import org.whale.WhalePriors;

/**
 *
 * @author rflachlan
 */
public class SimulationRunner {
    int type=0;
    int n, id, round;
    int numpops=1;
    long seed;
    String[] fileloc;
    
    public double[][] results;
    public double[][] stats;
    public double[][] params;
    public double[] score;
    
    double[][] pvals1, cov;
    double[] weights;
    
    WhalePriors p;
    Whale[] cemp;
    int[][][] locs;
    int[][] repSizes;
    int[][][] popSizes;
    String currentpop="";
    
    double[] means, sds;
    double[][] loadings;

    
    public SimulationRunner(Whale[] cemp, int numpops, int n, int id, int round, long seed, String[] fileloc, double[] weights, double[][] pvals1, double[][]cov){
    
	this.n=n;
        this.numpops=numpops;
	this.id=id;
	this.round=round;
	this.results=new double[n][];
        this.seed=seed;
        this.fileloc=fileloc;
        p=new WhalePriors(seed);
        
        this.cemp=cemp;
        
        means=cemp[0].means;
        sds=cemp[0].sds;
        loadings=cemp[0].loadings;
        
        //double[] x=p.sampleFromPriors();
        //cemp=new ChaffinchABC(fileloc, x, p.variables, System.currentTimeMillis());
	locs=new int[numpops][][];
        repSizes=new int[numpops][];
        popSizes=new int[numpops][][];
        for (int i=0; i<numpops; i++){
            locs[i]=cemp[i].locs;
            repSizes[i]=cemp[i].repSizes;
            popSizes[i]=cemp[i].popSizes;
        }
        
        if (weights!=null){
            this.weights=new double[weights.length];
            System.arraycopy(weights, 0, this.weights, 0, weights.length);
        }
        
        if (pvals1!=null){
            this.pvals1=new double[pvals1.length][];
            for (int i=0; i<pvals1.length; i++){
                this.pvals1[i]=new double[pvals1[1].length];
                System.arraycopy(pvals1[i], 0, this.pvals1[i], 0, pvals1[i].length);
            }
        }
        
        if (cov!=null){
            this.cov=new double[cov.length][];
            for (int i=0; i<cov.length; i++){
                this.cov[i]=new double[cov[1].length];
                System.arraycopy(cov[i], 0, this.cov[i], 0, cov[i].length);
            }
        }
        
    }
    
    public SimulationRunner(int n, int id, int round, long seed){
        this.n=n;
	this.id=id;
	this.round=round;
	this.results=new double[n][];
        this.seed=seed;
        this.fileloc=fileloc;
        p=new WhalePriors(seed);
        type=1;
    }
    
    public void readRound(){
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
            for (int i=0; i<weights.length; i++) {
                double[] x=locs.get(i);
                pvals1[i]=new double[x.length-1];
                System.arraycopy(x, 0, pvals1[i], 0, x.length-1);
                weights[i]=x[x.length-1];
            }
        }
		
	catch(Exception e){
            e.printStackTrace();
	}
        
        
       // for (int i=0; i<cov.length; i++){
        //    System.out.println("COV: "+i+" "+cov[i][0]+" "+cov[i].length);
        //}
        //for (int i=0; i<weights.length; i++){
        //    System.out.println("WEI: "+i+" "+weights[i]);
        //}
        
    }
		
    public  void runSimulation(){
        
        if (type==0){
            //runSimulationC();
        }
        else{
            runSimulationW();
        }
    }
     
/*    
    public void runSimulationC(){
	stats=new double[n][];
	params=new double[n][];
	score=new double[n];
	Whale cabc;
	
			
        //if (round>0){
          //  readRound();
        //}
        
	for (int i=0; i<n; i++) {
            //System.out.println("Start: "+id+" "+round+" "+i);
				
            double[] x;
            if (round==0) {
		x=p.sampleFromPriors();		
            }
            else {
		double[] y=pvals1[pickParams()];
		x=p.drawFromProposal(cov, y);
            }
				
            params[i]=x;
            int xx=cemp[0].stats.length;
            stats[i]=new double[xx*numpops];
            score[i]=0;
            for (int j=0; j<numpops; j++){
                currentpop=fileloc[j];
                cabc=new Whale(fileloc[j], p.nextLong(), locs[j], repSizes[j], popSizes[j], x, p.variables);		
                System.arraycopy(cabc.stats, 0, stats[i], j*xx, xx);
                score[i]+=calculateDifference(cabc.stats, cemp[j].stats);
            }		
            score[i]/=numpops+0.0;
            System.out.println("End: "+id+" "+round+" "+i+" "+score[i]+" "+x[0]+" "+x[1]+" "+x[2]+" "+x[3]+" "+x[4]);
            cabc=null;
            System.gc();
        }
       // writeRound();
     
    }
    */
    
    public void runSimulationW(){
	stats=new double[n][];
	params=new double[n][];
	score=new double[n];
	Whale whale;
        p=(WhalePriors)p;
	for (int i=0; i<n; i++) {				
            double[] x;
            if (round==0) {
		x=p.sampleFromPriors();		
            }
            else {
		double[] y=pvals1[pickParams()];
		x=p.drawFromProposal(cov, y);
            }
            params[i]=x;
            score[i]=0;
            for (int j=0; j<numpops; j++){
                //whale=new Whale(x, p.variables, p.popsizes, p.minpops, p.nextLong());
                whale=new Whale(x, p.variables, p.kpops, p.minpops, p.nextLong());
                score[i]=0;
            }		
            score[i]/=numpops+0.0;
            //System.out.println("End: "+id+" "+round+" "+i+" "+score[i]+" "+x[0]+" "+x[1]+" "+x[2]+" "+x[3]+" "+x[4]);
            whale=null;
            System.gc();
        }
       // writeRound();
     
    }
                
    public int pickParams() {
        int n=weights.length;
        double x=p.nextDouble()*weights[n-1];
        int loc=0;
        for (int i=0; i<n; i++){
            if (weights[i]>x){
		loc=i;
		i=n;
            }
        }	
        return loc;
		
    }
                
    public double calculateDifference(double[] x1, double[] x2) {
		
        double[] compsx=calculatePLSComponents(x1);
        double[] compsy=calculatePLSComponents(x2);
        for (int i=0; i<x1.length; i++){
            if (Double.isNaN(x1[i])){System.out.println("NaN error: xs: "+i);}
            if (Double.isNaN(x2[i])){System.out.println("NaN error: ys: "+i);}
            if (Double.isInfinite(x1[i])){System.out.println("Inf error: xs: "+i+" "+currentpop);}
            if (Double.isInfinite(x2[i])){System.out.println("Inf error: ys: "+i+" "+currentpop);}
        }
		
        double d[]=new double[compsx.length];
        double r=0;
        for (int i=0; i<compsx.length; i++){
			//d[i]=Math.abs(compsx[i]-compsy[i]);
            d[i]=compsx[i]-compsy[i];
			//System.out.println(i+" "+compsx[i]+" "+compsy[i]);
            if (Double.isNaN(d[i])){
		System.out.println("NaN error: "+i+" "+currentpop);
		d[i]=1000000000;
            }
            r+=d[i]*d[i];
        }
		//System.out.println(r);
        return Math.sqrt(r);	
    }
	
    public double[] calculatePLSComponents(double[]x) {
        
        
	double[] out=new double[loadings[0].length];
		
		//System.out.println(loadings.length+" "+loadings[0].length);
		//System.out.println(x.length);
		
		//for (int j=0; j<loadings.length; j++){
		//	System.out.println(j+" "+x[j]+" "+means[j]+" "+sd[j]+" "+(x[j]-means[j])/sd[j]);
		//}
		
		
		
	for (int i=0; i<loadings[0].length; i++){
            for (int j=0; j<loadings.length; j++){
		out[i]+=((x[j]-means[j])/sds[j])*loadings[j][i];
                
                if(Double.isNaN(x[j])){System.out.println(currentpop+" "+j);}
				//System.out.println(i+" "+j+" "+out[i]+" "+x[j]+" "+means[j]+" "+sd[j]+" "+loadings[j][i]);
            }
        }
                
                /*
                for (int i=0; i<x.length; i++){
                    System.out.println("IN: "+x[i]);
                }
                for (int i=0; i<out.length; i++){
                    System.out.println("OUT: "+out[i]);
                }
                */
	return out;
    }
    
    
    public int passThresh(double threshold){
        int count=0;
        for (int i=0; i<score.length; i++){
            if (score[i]<threshold){count++;}
        }
        
        return count;
    }
    
   /* 
     public void writeRound(){
            	DocumentSave ds=new DocumentSave(fileloc+"/outputstats"+id+".csv", ", ");

                for (int i=0; i<stats.length; i++){
                    for (int j=0; j<stats[i].length; j++){
                        ds.writeDouble(stats[i][j]);
                    }
                    ds.writeDouble(score[i]);
                    ds.writeLine();      
                }
		
		ds.finishWriting();
                
                ds=new DocumentSave(fileloc+"/outputparams"+id+".csv", ", ");
		
                for (int i=0; i<params.length; i++){
                    for (int j=0; j<params[i].length; j++){
                        ds.writeDouble(params[i][j]);
                    }
                    ds.writeLine();   
                }	
		ds.finishWriting();       
        }
  */  
    
    public static void main (String args[]) {
       // System.out.println(args[0] + " "+ args[1]+ " "+args[2]+" "+args[3]);
        int a=Integer.parseInt(args[0]);
        int b=Integer.parseInt(args[1]);
        int c=Integer.parseInt(args[2]);
        long d=Long.parseLong(args[3]);
        String[] e=new String[1];
        e[0]=args[4];
        int f=1;
        WhalePriors p=new WhalePriors(System.currentTimeMillis());
        double[] x=p.sampleFromPriors();
        Whale[] cemp=new Whale[1];
        //cemp[0]=new Whale(e[0], x, p.variables, System.currentTimeMillis());
        
        
	SimulationRunner sr=new SimulationRunner(cemp, 1, a,b,c,d,e, null, null, null);
	sr.runSimulation();
        
        
    }
		
	
}
        
        
        
    

