/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.Whale;

/**
 *
 * @author lieszandberg
 */

import java.util.Arrays;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.function.Logistic;
import org.apache.commons.math3.fitting.SimpleCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

public class WhaleMeasureStatistics extends org.ChaffinchABC.MeasureStatistics {
        
    
    

	WhalePopulation population;
        WhalePopulation emppop;
	WhaleParameters param;
        WhaleIndividual[] inds;
        Whale wh;
        WhaleIndividual individual;
	double[] shareprofile;
	public double ntypes;
        int memorylength=100;
	
	//public double meanSylDiss, sylaveraget, sylsingletons, sylrares, sylintermediates, sylcommons, sylmaxfreq, sylh1, sylLD, sylshzero, sylshsing, sylshmult, sylshall, sylshmean, hdiff; 
	
        public double jaccardTot, withinTot, sylRatio;
	double[] geogThresh= {1.01, 3.5, 7, 70, 250};
        double geogthresh=5;
	
	double[] jaccard, withinTypeDiss;
	//double clusterDepth, clusterDepthEmp;
	double[][] diss, dissSyl;
	
	double[]out;
        int[][][][] output;
      
        int maxlearn=10;
        int[] ids;
        int[] syllab;
        int[] indbuffer;
        int[] membuffer;
        int[][][] output1;
        int[][][] output2;
        double[][][] output3;
        int[] subpopsize;
        int[][][] sharedsongs;
        double[][][] output4;

        


	
	double dsim=0.04;
	
	boolean typeEmpirical=true;
	
	public WhaleMeasureStatistics(org.Whale.WhalePopulation population, WhaleParameters param) {
            this.population=population;
            this.param=param;     
  //          this.individual=individual;
            typeEmpirical=true;
            dsim=param.typeThresh;
            this.subpopsize=population.subpopsize;

            
            //diss=population.calculateEmpDissimilarityMatrix(0);
            //ids=population.calculateEmpIDs();
            //dissSyl=population.calculateEmpSyllDissimilarityMatrix(0);
            //syllab=makeEmpSylLabel();
            //calculateStats();
            calculateThresholdSpectrum(population.emppop);
            calculateSongSharing();
            calculateSongSharingAvg();
            
        }
        
        public WhaleMeasureStatistics(WhalePopulation population, EmpData ed, WhaleParameters param) {
		this.population=population;
		this.param=param;
		dsim=param.typeThresh;
                memorylength=param.memorylength;
		ids=population.calculateEmpIDs();
		dissSyl=ed.empSylDiss;
		syllab=ed.syllab;
		diss=ed.empDiss;
                //calculateStats();

                
                
	}
        
            
        public void calculateStats(){
            double[] x2=calculateThresholdStats(ids, population.nemp, geogthresh);
            double[] x1=calculateThresholdStats(ids, population.nemp, 100000); 
            double[] x3=calculateGeogStats(ids, population.nemp, diss, param.repSizes);
            double[] x4=calculateSyllableStats(syllab);
            //hdiff=x1[6]-x4[6];
            //double avdiff=x1[0]-x4[0];
            //double singdiff=x1[1]-x4[1];
            out=new double[12+12+12+13];
            System.arraycopy(x1, 0, out, 0, 12);
            System.arraycopy(x2, 0, out, 12, 12);
            System.arraycopy(x3, 0, out, 24, 12);
            System.arraycopy(x4, 0, out, 36, 13);
            //out[49]=hdiff;
            //out[50]=avdiff;
            //out[51]=singdiff;
	}

                
        

	public void calculateThresholdSpectrum(WhaleIndividual[] pop){
           output1= new int[pop.length][param.memorylength][];
           output2= new int[pop.length][param.memorylength][];
           indbuffer= new int[param.memorylength*pop.length];
           membuffer= new int[param.memorysize*pop.length];
           int ml=param.memorylength;
           for (int i=0; i<pop.length; i++){ //For each individual i in population length
               if(i%100==0){
               System.out.println("Individual = " + i);
               }
               float[] mema=pop[i].getSongMemory();
               for (int a=0; a<ml; a++){ //memlength //for each song in total memory
                  int k=0;                    
                  int indexa=pop[i].getMemoryIndex(a); //get the index of song a    
                  for (int j=0; j<pop.length; j++){ //for each other individual in the population
                      float[] memb=pop[j].getSongMemory();
                      if(i!=j){ //as long as it is not individual i itself
                          for (int b=0; b<ml; b++){ 
                            int indexb=pop[j].getMemoryIndex(b);
                            if(pop[i].matchSongs(mema, memb, indexa, indexb)){
                              indbuffer[k]=j;
                              membuffer[k]=b;
                              k++;   
                            }                     
                          }
                      }
                  }
                  output1[i][a]=new int[k];
                  output2[i][a]=new int[k];
                  System.arraycopy(indbuffer, 0, output1[i][a], 0, k);
                  System.arraycopy(membuffer, 0, output2[i][a], 0, k); 
                  //System.out.println(k);
               }
           } 
          
}
        public void calculateSongSharing(){
        output3= new double[population.emppop.length][param.memorylength][param.memorylength];
        for (int i=0; i<population.emppop.length; i++){ //for each individual!!
          if(i%100==0){
              System.out.println("calculateSongSharing ID = " + i);
          }  
          int a=population.getEmpPop(i); //get the population number
          for (int j=0; j<param.memorylength; j++){ //for each timepoint
            //int count=0;  
            for (int k=0; k<output1[i][j].length; k++){ //go through songs of output1
              int b=population.getEmpPop(output1[i][j][k]); //b is population number of each k (other ID)
              if (a==b){             
                output3[a][j][output2[i][j][k]]++;          
                
              }
            }
          }
        }
        for (int i=0; i<subpopsize.length; i++){            //for every population
          double x=population.sampleperpop*population.sampleperpop; //x=samplesize pop i * samplesize pop i
          output4 = new double[subpopsize.length][param.memorylength][param.memorylength];
          for (int j=0; j<param.memorylength; j++){             //for every time point
            for (int k=0; k<param.memorylength; k++){           //for every other time point 
              //System.out.println("pop = " + i + " t = " + j + "sharing = " + output3[i][j][k]);
                output3[i][j][k]/=x;
              //System.out.println("pop = "+i+" timepoint = "+j+ " shared song time= " + k + " songsharing " +  output4[i][j][k]); 
              
              
            }
          }
        }
        
        
        }      
        
        public void calculateSongSharingAvg(){
            double[] outputavg= new double[10];

            for(int i=0; i<maxlearn; i++){ //timepoint i in yearly cycle of length maxlearn
                double sum=0;
                int x=0;    
                for(int j=0; j<population.subpopsize.length;j++){    //j=pop
                    for(int k=0; k<memorylength;k++){ //timepoint a=k
                        if(k%10==i){
                            for(int l=0; l<memorylength;l++){ //timepoint b=l
                                if(k==l){
                                sum+=output3[j][k][l];
                                //System.out.println(output3[j][k][l]);
                                x++;
                                } 
                            }
                        } 
                    }
                }
                outputavg[i]=sum/x;
                System.out.println("Time = " + i + "  & sharing avg = " + outputavg[i]);
            }
            
            
        }
        
/* 
        public void calculatePopSongSharing(){
        int[][][][] outputpop= new int[population.emppop.length][population.subpopsize.length][param.memorylength][param.memorylength];
        for (int i=0; i<population.emppop.length; i++){ //for each individual!!
          int a=population.getEmpPop(i); //get the population number
          for (int j=0; j<param.memorylength; j++){ //for each timepoint
            //int count=0;  
            for (int l=0; l<3; l++){
                if(l==0){
                    a=a-1;
                }
                else if(l==1){
                    a=a;
            }
                else if(l==1){
                    a=a+1;
            }
                for (int k=0; k<output1[i][j].length; k++){ //go through songs of output1
                int b=population.getEmpPop(output1[i][j][k]); //b is population number of each k (other ID)
                if (a==b){             
                output3[a][j][output2[i][j][k]]++;          //output3[a][j][]=x++
                }
              }
            }
          }
        }
        for (int i=0; i<subpopsize.length; i++){            //for every population
          double x=population.sampleperpop*population.sampleperpop; //x=samplesize pop i * samplesize pop i
          output4 = new double[subpopsize.length][param.memorylength][param.memorylength];
          for (int j=0; j<param.memorylength; j++){             //for every time point
            for (int k=0; k<param.memorylength; k++){           //for every other time point 
              output4[i][j][k]=output3[i][j][k]/x; 
              //System.out.println("pop = "+i+" timepoint = "+j+ " shared song time= " + k + " songsharing " +  output4[i][j][k]); 
              
              
            }
          }
        }
        
        
        }  
*/        
        
        
}       
        
        
   
 



    
    

    

            /*int n=d.length;
            int[][] out=new int[n][];
            int[] buffer=new int[n];
            //System.out.println("SONGS: "+n);
            for (int i=0; i<n; i++) {
                int a=0;
                int b=0;
                for (int j=0; j<n; j++) {
                    int s1=ids[i];
                    int s2=ids[j];
                    double dist=population.getEmpDistance(s1,s2);
                    //if(d[i][j]<t){
                    //    System.out.println((s1==s2)+" "+d[i][j] +" "+i+" "+j+" "+s1+" "+s2);
                    //}
                    
                    //System.out.println(s1+" "+s2+" "+dist);
                    if ((s1!=s2)&&(dist<geogthresh)){
                        b++;
                        if(d[i][j]<t){
                            buffer[a]=j;
                            a++;	
                        }
                    }
                }
		out[i]=new int[a+1];
		System.arraycopy(buffer, 0, out[i], 0, a);
                out[i][a]=b;
            }
            return out;
	}*/
        
 
 //       public int[][] calculateThresholdSpectrum(WhaleIndividual[] pop){
      
            
            
            
            
            /*
            allsong= new int[allsonglength];
            for (int i=0; i<ninds; i++){ //for all inds
                for (int j=0; j<memorysize,j++){
                    allsong[k]=songmemory[j];
                    k++;                      
                }
            }  
            
            for (int i=0; i<ninds; i++){ //for each id
                for (int j=0; j<memorylength,j++){ //each song in its memory
                     for (int k=0; k<allsongnr,k++){ //compare it so all songs in allsonglength
                         int a=k*ns;
                         int b=0
                         if(matchSongs(allsong, songmemory, a, j){
                          System.arraycopy(allsong, a, tresholdsong, b, ns);  
                          b++
                         }
                         
                         
                         int id=ninds[i];
                         int songid=j;
                         int 
                     }
                    
                    
            }
            }   
            
            */
            /*
            int n=d.length;
            int[][] out=new int[n][];
            int[] buffer=new int[n];
		
            for (int i=0; i<n; i++) {
                int a=0;
                for (int j=0; j<n; j++) {
                    if(d[i][j]<t){
                        buffer[a]=j;
                        a++;
                    }
                }
		out[i]=new int[a];
		System.arraycopy(buffer, 0, out[i], 0, a);
            }
            return out;
	}
	
	public double[] calculateThresholdStats(int[] ids, int nemp, double geogthresh) {
		int[][] x=calculateThresholdSpectrum(diss, dsim, ids, geogthresh);
                
		int n=x.length;
                shareprofile=new double[n];
		
                //System.out.println(intThresh+" "+n);
		double singletons=0;
		double rares=0;
		double intermediates=0;
		double commons=0;
		double maxfreq=-1000000;
		double averaget=0;
		
		double h1=0;
		double alpha1=0;
		double log2=1/Math.log(2);

		int mf=0;
		double ci=0;
                double mni=0;
		for (int i=0; i<n; i++) {
                    double ni=x[i][x[i].length-1];
                    if(ni>0){
                        ci++;
                        mni+=ni;
                        int intThresh=(int)Math.round(ni/20);
                        int rareThresh=(int)Math.round(4*ni/n);
                        if (rareThresh<2){rareThresh=2;}
                        if (intThresh<3){intThresh=3;}
			int p=x[i].length;
                        shareprofile[i]=p/ni;
                        
			averaget+=p/ni;

			if (p==1) {
				singletons++;
			}
			else if (p<=rareThresh) {
				rares++;
			}
			else if (p<=intThresh) {
				intermediates++;
			}
			else {
				commons++;
			}
                        //double mfr=logit(p, ni);
                        //System.out.println(i+" "+p+" "+ni+" "+maxfreq);
			if (p>maxfreq) {
                            maxfreq=p;
			}
                        
                        //if (p>n/2){
                        //    System.out.println(i+" "+p+" "+n);
                            
                        //}
                        
			if (p>mf){mf=p;}
			double x2=p/ni;
			
			h1+=Math.log(x2)*log2;
                    }
		}
		
                mni/=ci;
                
		//System.out.println(geogthresh+" "+n+" "+singletons+" "+rares+" "+intermediates+" "+commons+" "+maxfreq+" "+dsim);
		
		//System.out.println("S "+singletons+ " "+n);
		singletons=logit(singletons, n);
		commons=logit(commons, n);
		intermediates=logit(intermediates,n);
		rares= logit(rares, n);
                maxfreq=logit(maxfreq, n);
                averaget=logit(averaget, n);
	
		double[] counts=new double[mf+1];
		for (int i=0; i<n; i++) {
			int p=x[i].length-1;
			counts[p]++;
			
		}
		double s=0;
		for (int i=1; i<counts.length; i++) {
			counts[i]/=i+0.0;
			alpha1+=counts[i]*Math.log(i*2.0);
			s+=counts[i];
		}
		if (alpha1>0){
                    alpha1=1+(s/alpha1);
                }
                else{alpha1=0;}
		double[] fit1=calculateFit(alpha1, counts);
		double fita=fit1[0];
		double fitb=fit1[1];
		
		
		int[] y=calculateShareSpectrum2(ids, nemp, geogthresh);
		double st=0;
		for (int i=1; i<10; i++){
			st+=y[i];
		}
		
		double shareprop=logit(st, y[0]+st); //PROPORTION OF INDIVIDUALS SHARING AT LEAST ONE SONG
                double shareprop2=logit(st-y[1], st); //OF THOSE SHARING AT LEAST ONE SONG, HOW MANY SHREA MORE THAN 1...
                if (st==0){shareprop2=logit(0,100);}
                if (y[1]==0){shareprop2=-10;}
                if (Double.isInfinite(shareprop2)){shareprop2=-10;}
                double[] out={averaget, singletons, rares, intermediates, commons, maxfreq, h1, alpha1, fita, fitb, shareprop, shareprop2};
                return out;
	
	}
        
}
        
        
        
	
        /*
	public void calculateFrequencyStats(int[] classification, int[] ids, int nemp) {
		//System.out.println(classification.length+" "+ids.length);
		int[] f=calculateFrequencies(classification);
		int n=f.length;
		double m=classification.length;
		//double m=population.nemp;
		int intThresh=(int)Math.round(0.05*m);
		
		singletons=0;
		rares=0;
		intermediates=0;
		commons=0;
		maxfreq=0;
		ntypes=logit(n, m);
		
		h1=0;
		alpha1=0;
		double log2=1/Math.log(2);
		double[] counts1=new double[classification.length+1];
		int t=0;
		
		for (int i=0; i<n; i++) {
			t+=f[i];
			if (f[i]==1) {
				singletons++;
			}
			else if (f[i]<=4) {
				rares+=f[i];
			}
			else if (f[i]<=intThresh) {
				intermediates+=f[i];
			}
			else {
				commons+=f[i];
			}
			if (f[i]>maxfreq) {
				maxfreq=f[i];
			}
			
			double x2=f[i]/m;
			
			h1+=x2*Math.log(x2)*log2;
			counts1[f[i]]++;
			
			alpha1+=Math.log(f[i]*2.0);
			
			
			//System.out.println("FREQCHECK: "+f[i]+" "+m);
		}
		
		singletons=logit(singletons, m);
		commons=logit(commons, m);
		intermediates=logit(intermediates,m);
		rares= logit(rares, m);
		

		maxfreq=logit(maxfreq, m);
		
		alpha1=1+(f.length/alpha1);
		double[] fit1=calculateFit(alpha1, counts1);
		fita=fit1[0];
		fitb=fit1[1];
		
		//int[] p=new int[t];
		
		//t=0;
		//for (int i=0; i<f.length; i++){
		//	for (int j=0; j<f[i]; j++){
		//		p[t]+=f[i];
		//		t++;
		//	}
		//}
		
		//Arrays.sort(p);
		
		int[] y=calculateShareSpectrum(classification, ids, nemp);
		double st=0;
		for (int i=1; i<10; i++){
			st+=y[i];
		}
		
		shareprop=logit(st, y[0]+st);
		shareprop2=logit(st-y[1], st);
		
		
		
		
		//System.out.println("FREQSTATS: "+ntypes+" "+singletons+" "+intermediates+" "+commons+" "+maxfreq);
	}
	*/
        
/*	public int[] calculateShareSpectrum(int[] f, int[]g, int h){
		int[] shareSpectrum=new int[20];
		
		int[][]c=new int[h][h];
		
		for (int i=0; i<f.length; i++) {
			for (int j=0; j<i; j++) {
				if (f[i]==f[j]) {
					if (g[i]>=g[j]) {
						c[g[i]][g[j]]++;
					}
					else {
						c[g[j]][g[i]]++;
					}
				}
			}
		}
		
		
		for (int i=0; i<c.length; i++){
			for (int j=0; j<i; j++){
				if (c[i][j]>=20) {
					System.out.println("ALERT SHARE: "+i+" "+j+" "+c[i][j]);
				}
				else {
					shareSpectrum[c[i][j]]++;	
				}
			}
		}	
		
		return shareSpectrum;
	}
	
//	public int[] calculateShareSpectrum2(int[]g, int h, double gt){
		int[] shareSpectrum=new int[20];
		
		int[][]c=new int[h][h];
		
		for (int i=0; i<diss.length; i++) {
			for (int j=0; j<i; j++) {
                            int s1=g[i];
                            int s2=g[j];
                            double dist=population.getEmpDistance(s1,s2);
                            if (dist<gt){
				if (diss[i][j]<dsim) {
					if (g[i]>=g[j]) {
						c[g[i]][g[j]]++;
					}
					else {
						c[g[j]][g[i]]++;
					}
				}
                            }
			}
		}
		
		
		for (int i=0; i<c.length; i++){
			for (int j=0; j<i; j++){
                            double dist=population.getEmpDistance(i,j);
                            if (dist<gt){
				if (c[i][j]>=20) {
					//System.out.println(i+" "+j+" "+c[i][j]);
				}
				else {
					shareSpectrum[c[i][j]]++;	
				}
                            }
			}
		}	
		
		return shareSpectrum;
	}
	
	public double logit(double x, double y){
		if (x==0){x=1;}
		if (x==y){x=y-1;}
		
		double z=x/(y-x);
		
		return Math.log(z);
	}
	
//	public double[] calculateFit(double alpha, double[] x){
		double fmin=0;
		for (int i=0; i<1000; i++){
			fmin+=Math.pow(i+1, -1*alpha);
		}
		
		//for (int i=0; i<x.length; i++) {
		//	System.out.println(i+" "+x[i]);
		//}
		
		double[] cum=new double[x.length];
		
		for (int i=cum.length-2; i>0; i--){
			cum[i]=cum[i+1]+x[i];
		}
		double ts=cum[1];
		for (int i=0; i<cum.length; i++){
			cum[i]/=ts;
		}
		
		
		double[] fit={0,0};
		for (int i=0; i<cum.length; i++){		
			if ((cum[i]>0)&&(cum[i]<1)){
				
				double f=0;
				for (int j=0; j<50; j++){
					f+=Math.pow(j+i, -1*alpha);	
				}
				
				f/=fmin;
				
				
				
				double p=Math.abs(f - cum[i]);	
				//System.out.println(i+" "+cum[i]+" "+f+" "+p);
				if (p>fit[0]){fit[0]=p;}
				if (i==2){fit[1]=f - cum[i];}
			}	
		}
		return fit;
		
	}
	
//	public int[] calculateFrequencies(int[] d) {
		int x=0;
		for (int i=0; i<d.length; i++) {
			if (d[i]>x) {x=d[i];}
		}
		x++;
		int[] out=new int[x];
		for (int i=0; i<d.length; i++) {
			out[d[i]]++;
			//System.out.println("FREQUCALCCHK: "+d[i]+" "+out[d[i]]);
		}
		return out;
	}
	/*
	public double[] calculateGeogStats(int[] ids, int n, double[][] diss, int[] repSize) {
		
		double[] jaccard=new double[geogThresh.length+1];
		double[] withinTypeDiss=new double[geogThresh.length+1];
		double[] count=new double[geogThresh.length+1];
		double[] count2=new double[geogThresh.length+1];
		
		
		double[][] share=new double[n][n];
		double[][] wd=new double[n][n];
		
                /*
		for (int i=0; i<diss[0].length; i++) {
			boolean[] check=new boolean[n];
			for (int j=0; j<i; j++) {
				if ((!check[ids[j]])&&(diss[i][j]<param.typeThresh)) {
					share[ids[i]][ids[j]]++;
					wd[ids[i]][ids[j]]+=diss[i][j];
					check[ids[j]]=true;
				}
			}
		}
               
                for (int i=0; i<diss.length; i++){
                    boolean[] check=new boolean[n];
                    for (int j=0; j<i; j++){
                        if ((!check[ids[j]])&&(diss[i][j]<param.typeThresh)){
                            share[ids[i]][ids[j]]++;
                            wd[ids[i]][ids[j]]+=diss[i][j];
                            share[ids[j]][ids[i]]++;
                            wd[ids[j]][ids[i]]+=diss[i][j];
                            check[ids[j]]=true;
                        }
                    }
                }
		
	
		for (int i=0; i<n; i++) {
			for (int j=0; j<i; j++) {
				double dist=population.getEmpDistance(i,j);
				//System.out.println(i+" "+j+" "+dist);
				int c=geogThresh.length;;
				for (int k=0; k<geogThresh.length; k++) {
					if (dist<=geogThresh[k]) {
						c=k;
						k=geogThresh.length;
					}
				}
				
				double jacc=share[i][j]/(repSize[i]+repSize[j]-share[i][j]);
                                
				jaccard[c]+=jacc;
				count[c]++;
				
				if (share[i][j]>0) {
					withinTypeDiss[c]+=wd[i][j]/share[i][j];
					count2[c]++;
				}
				
				
			}
		}
		
		double jt=0;
                double wt=0;
		for (int i=0; i<geogThresh.length; i++) {
			jaccard[i]/=count[i];
                        if (count[i]==0){jaccard[i]=0;}
			if (count2[i]>0) {
				withinTypeDiss[i]/=count2[i];
			}
                        jt+=jaccard[i];
                        wt+=withinTypeDiss[i];
			//System.out.println("Geog stats: "+i+" "+jaccard[i]+" "+withinTypeDiss[i]+" "+count[i]);
		}
                
                //double jaccardTot=jt;
                //double withinTot=wt;
                for (int i=0; i<geogThresh.length; i++) {
                    if (jt>0){jaccard[i]/=jt;}
                    if (wt>0){withinTypeDiss[i]/=wt;}
                    else{withinTypeDiss[i]=dsim;}
                }
                
                double[] out=new double[12];
                out[0]=jt;
                System.arraycopy(jaccard, 0, out, 1, 5);
                out[6]=wt;
                System.arraycopy(withinTypeDiss,0, out, 7, 5);
                return out;
                
	}
*/
//	
/*	public int[] makeEmpSylLabel() {
		int n=dissSyl.length;
		int[] label=new int[n];
		int a=0;
		int b=0;
		for (int i=0; i<n; i++) {
			label[i]=a;
			b++;
			if (b==param.sylsPerSong) {
				b=0;
				a++;
			}
		}
		return label;
	}
	
	//public double[] calculateSyllableSharing(int[][] x, int[] y, int[] songlengths) {
		
		int nsong=songlengths.length;
                //System.out.println("NSONG: "+nsong);
		double[][] share=new double[nsong][nsong];
		
		
		for (int i=0; i<x.length; i++) {
                    for (int j=0; j<x[i].length; j++) {
                        share[y[i]][y[x[i][j]]]++;
                    }
		}
		
               
            
		double meanshare=0;
		double zeroes=0;
		double singles=0;
		double multiples=0;
		double alls=0;
		double count=0;
		for (int i=0; i<share.length; i++) {
			for (int j=0; j<i; j++) {
				double sc=(share[i][j]+share[j][i])*0.5;
                                double sp=shareprofile[i]*shareprofile[j];
				meanshare+=sc*sp;
				if (sc==0) {zeroes+=sp;}
				else if (sc<=1) {singles+=sp;}
				else {multiples+=sp;}
				if (sc>songlengths[i]/2) {alls+=sp;}
				count+=sp;
			}
		}
		
		
		double[] out= {meanshare/count, zeroes/count, singles/count, multiples/count, alls/count};
		
		return out;
	}
	
//	public int[] calculateSongLengths(int[] label) {
		int n=label.length;
		int lmax=0;
		for (int i=0; i<n; i++) {
			if (label[i]>lmax) {lmax=label[i];}
		}
		int[] out=new int[lmax+1];
		for (int i=0; i<n; i++) {
			out[label[i]]++;
		}
                
		return out;
	}
        
   //     public void calculateTransitions(int[] label){
            int n=dissSyl.length;
            int[][] x=calculateThresholdSpectrum(dissSyl, dsim);
            int[] slengths=calculateSongLengths(label);
            
        }
        

	//public double[] calculateSyllableStats(int[] label) {
		
		int n=dissSyl.length;
		
		int[][] x=calculateThresholdSpectrum(dissSyl, dsim);
		int[] slengths=calculateSongLengths(label);
		
                //double[] y=calculateTransThresholdSpectrum(dissSyl, dsim, label);
                
		/*
		meanSylDiss=0;
		double count=0;
		if (!typeEmpirical) {
			
			int a=0;
			for (int i=0; i<diss.length; i++) {
				for (int j=0; j<param.sylsPerSong; j++) {
					for (int k=0; j<i; j++) {
						int b=k*param.sylsPerSong+j;
						meanSylDiss+=dissSyl[a][b];
						count++;
					}
					a++;
				}
			}
		}
		//System.out.println(n+" "+meanSylDiss+" "+count);
		meanSylDiss/=count;
		
		//System.out.println(n+" "+meanSylDiss+" "+count);
		
		double[] ssh =calculateSyllableSharing(x, label, slengths);
		
		double sylshmean=ssh[0];
		double sylshzero=ssh[1];
		double sylshsing=ssh[2];
		double sylshmult=ssh[3];
		double sylshall=ssh[4];
                
                
          
                double sylRatio=0;
                if (sylshall==0){sylRatio=-10;}
                else if (sylshsing==0){sylRatio=-10;}
                else{
                    sylRatio=Math.log(sylshsing/sylshall);
                }
		
		//System.out.println("Done with LD");
		int intThresh=(int)Math.round(n/24);
		double sylsingletons=0;
		double sylrares=0;
		double sylintermediates=0;
		double sylcommons=0;
		double sylmaxfreq=0;
		double sylaveraget=0;
		
		double sylh1=0;
		
		double log2=1/Math.log(2);

		int mf=0;
		
		for (int i=0; i<n; i++) {
			int p=x[i].length;
			sylaveraget+=p/(n+0.0);
			
			if (p==1) {
				sylsingletons++;
			}
			else if (p<=4) {
				sylrares++;
			}
			else if (p<=intThresh) {
				sylintermediates++;
			}
			else {
				sylcommons++;
			}
			if (p>sylmaxfreq) {
				sylmaxfreq=p;
				mf=p;
			}
			
			double x2=p/(n+0.0);
			sylh1+=Math.log(x2)*log2;

		}
		
		
		//System.out.println(meanSylDiss+" "+n+" "+sylaveraget+" "+sylsingletons+" "+sylrares+" "+sylintermediates+" "+sylcommons+" "+sylmaxfreq+" "+dsim);

		
		//System.out.println("Done with freq measures...");
		sylsingletons=logit(sylsingletons, n);
		sylcommons=logit(sylcommons, n);
		sylintermediates=logit(sylintermediates,n);
		sylrares= logit(sylrares, n);
		sylmaxfreq=logit(sylmaxfreq, n);
		sylaveraget=logit(sylaveraget, n);
		
		//System.out.println(meanSylDiss+" "+n+" "+sylaveraget+" "+sylsingletons+" "+sylrares+" "+sylintermediates+" "+sylcommons+" "+sylmaxfreq+" "+dsim);
                
                double[] out={sylaveraget, sylsingletons, sylrares, sylintermediates, sylcommons, sylmaxfreq, sylh1, sylshzero, sylshsing, sylshmult, sylshall, sylshmean, sylRatio};
                
                return out;
        }
	
//	public void calculateSyllableLD(int[][] x, int[]y, double n) {
		
		int a=x.length;
		
		double pa, pb, exp, act, cord, cor;
		double overall=0;
		double count=0;
		for (int i=0; i<a; i++) {
			//System.out.println(i+" "+a+" "+x[i].length);
			for (int j=0; j<i; j++) {
				boolean found=false;
				for (int k=0; k<x[i].length; k++) {
					if (x[i][k]==j) {
						found=true;
						k=x[i].length;
					}
				}
				if (!found) {
					pa=x[i].length/n;
					pb=x[j].length/n;
					
					exp=pa*pb;
					
					act=0;
					for (int r=0; r<x[i].length; r++) {
						for (int s=0; s<x[j].length; s++) {
							if (y[x[i][r]]==y[x[j][s]]) {
								act++;
								s=x[j].length;
							}
						}
					}
					act/=a+0.0;
					cord=Math.sqrt(pa*(1-pa)+pb*(1-pb));
					cor=(exp-act)/cord;
					overall+=cor;
					count++;
				}
			}
		}
		//sylLD=overall/count;
	}
//	
}
*/
