/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package whaleabc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedList;


/**
 *
 * @author rflachlan
 */
public class SimulationRunner {
    boolean printBestTrajectories=false; 
   int type=0;
    int n, id, round;
    long seed;
    
    public double[][] results;
    public double[][] stats;
    public double[][] params;
    public double[] score;
    
    double[][] pvals1, cov;
    double[] weights;
    
    double[] means, sds;
    
    WhalePriors p;  
    double[] empstats;
    
    public SimulationRunner(int n, int id, int round, long seed, double[] weights, double[][] pvals1, double[][]cov){
    
	this.n=n;
	this.id=id;
	this.round=round;
	this.results=new double[n][];
        this.seed=seed;
        p=new WhalePriors(seed);
        
        this.means=p.means;
        this.sds=p.sds;

        
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
        p=new WhalePriors(seed);
        type=1;
    }
    
    
    public  void runSimulation(){

	stats=new double[n][];
	params=new double[n][];
	score=new double[n];
	Whale whale;
        p=(WhalePriors)p;
        
        empstats=new double[p.statindices.length];
        for (int i=0; i<p.statindices.length; i++){
            empstats[i]=p.empstats[p.statindices[i]];
        }
        
	for (int i=0; i<n; i++) {	
            
            double[] x;
            if (round==0) {
		x=p.sampleFromPriors();		
            }
            else {
		double[] y=pvals1[pickParams()];
                while (y[0]<Math.log(0.000001)){y=pvals1[pickParams()];}
		x=p.drawFromProposal(cov, y);
            }
            params[i]=x;
            if (p.hemisphere==0){
                whale=new Whale(x, p.variables, p.popsizessh, p.minpopssh, p.kpopssh, p.nextLong(), i, this);
            }
            else if (p.hemisphere==1){
                whale=new Whale(x, p.variables, p.popsizesnh, p.minpopsnh, p.kpopsnh, p.nextLong(), i, this);
            }
            else if (p.hemisphere==2){
                whale=new Whale(x, p.variables, p.popsizessh2, p.minpopssh2, p.kpopssh2, p.nextLong(), i, this);
            }
            else {
                whale=new Whale(x, p.variables, p.popsizessh3, p.minpopssh3, p.kpopssh3, p.nextLong(), i, this);
            }
            double[] simstats=new double[p.statindices.length];
                for (int j=0; j<p.statindices.length; j++){
                simstats[j]=whale.wms.out[p.statindices[j]];
            }
            //simstats[p.statindices.length]=whale.param.mutationCounter;
            //simstats[p.statindices.length+1]=whale.param.insertionCounter;
            //simstats[p.statindices.length+2]=whale.param.deletionCounter;
            int[] tr=p.transf;
            score[i]=calculateDifference(empstats, simstats, tr);
            //System.out.println(id+" "+i+" "+n+" "+score[i]+" "+simstats[0]+" "+simstats[1]+" "+simstats[2]+" "+simstats[3]+" "+simstats[4]+" "+simstats[5]);
            System.out.println(id+" "+i+" "+n+" "+score[i]+" "+x[0]+" "+x[1]+" "+x[2]+" "+x[3]+" "+x[4]+" "+x[5]+" "+whale.param.mutationCounter+" "+whale.param.insertionCounter+" "+whale.param.deletionCounter+" "+(whale.wms.out[14]-whale.wms.out[18])+" "+simstats[simstats.length-3]);
            if(printBestTrajectories){
                if (score[i]<2){
                    System.out.println(id+" Accept run "+score);
                    whale.wms.traceThemes(id);
                }
            }
            
            stats[i]=new double[whale.wms.out.length+3];
            System.arraycopy(whale.wms.out, 0, stats[i], 0, whale.wms.out.length);
            stats[i][whale.wms.out.length]=whale.param.mutationCounter;
            stats[i][whale.wms.out.length+1]=whale.param.insertionCounter;
            stats[i][whale.wms.out.length+2]=whale.param.deletionCounter;
            whale=null;
            System.gc();
        }
     
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
                
    public double calculateDifference(double[] x1, double[] x2, int[] t) {
		
        double[] compsx=calculatePLSComponents(x1, t);
        double[] compsy=calculatePLSComponents(x2, t);
        for (int i=0; i<x1.length; i++){
            if (Double.isNaN(x1[i])){System.out.println("NaN error: xs: "+i);}
            if (Double.isNaN(x2[i])){System.out.println("NaN error: ys: "+i);}
            if (Double.isInfinite(x1[i])){System.out.println("Inf error: xs: "+i);}
            if (Double.isInfinite(x2[i])){System.out.println("Inf error: ys: "+i);}
        }
		
        double d[]=new double[compsx.length];
        double r=0;
        for (int i=0; i<compsx.length; i++){
			//d[i]=Math.abs(compsx[i]-compsy[i]);
            d[i]=compsx[i]-compsy[i];
			//System.out.println(i+" "+compsx[i]+" "+compsy[i]);
            if (Double.isNaN(d[i])){
		System.out.println("NaN error: "+i);
		d[i]=1000000000;
            }
            //double q=Math.sqrt(d[i]*d[i]);
            //if (q<0.001){q=0.001;}
            //r+=Math.log(q);
            r+=d[i]*d[i];
        }
	//System.out.println(x1[compsx.length-3]+" "+x2[compsx.length-3]+" "+compsx[compsx.length-3]+" "+compsy[compsx.length-3]+" "+sds[compsx.length-3])	;//System.out.println(r);
        return Math.sqrt(r);	
        //return Math.exp(r/(compsx.length+0.0));
    }
	
    public double[] calculatePLSComponents(double[]x, int[]t) {
        
        
	double[] out=new double[x.length];
	
	for (int i=0; i<x.length; i++){
            double y=x[i];
            if (t[i]==1){y=Math.sqrt(y);}
            if (t[i]==2){y=Math.log(y+0.001);}
            
            out[i]+=((y-means[i])/sds[i]);
            if(Double.isNaN(x[i])){System.out.println(i);}
            //System.out.print(x[i]+" ");
        }
        //System.out.println();

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
        
        
	SimulationRunner sr=new SimulationRunner(a,b,c,d, null, null, null);
	sr.runSimulation();
        
        
    }
		
	
}
        
        
        
    
