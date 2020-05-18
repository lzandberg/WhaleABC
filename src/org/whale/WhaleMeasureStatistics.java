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

import java.util.Arrays;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.function.Logistic;
import org.apache.commons.math3.fitting.SimpleCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

public class WhaleMeasureStatistics {
        
    
    

	WhalePopulation population;
        WhalePopulation emppop;
	WhaleParameters param;
        WhaleIndividual[] inds;
        Whale wh;
        WhaleIndividual individual;
        EmpData ed;
	
        int memorylength=100;
	
      
        int maxlearn=10;
        int[] ids;
        int[] syllab;
        int[] indbuffer;
        int[] membuffer;
        int[][][] output1;
        int[][][] output2;
        int[][][] output3;
        int[] subpopsize;
        int[][][] sharedsongs;
        double[][][] output4;
        double[][][][] output6;
        double[][] output7;
        double[] outputavg;
	double overallsongsharing=0;
        	double[][] diss, dissSyl;
	double dsim=0.04;
	
	boolean typeEmpirical=true;
	
	public WhaleMeasureStatistics(WhalePopulation population, WhaleParameters param) {
            this.population=population;
            this.param=param; 
            this.memorylength=param.loglength;
  //          this.individual=individual;
            typeEmpirical=true;
            dsim=param.typeThresh;
            this.subpopsize=population.subpopsize;

            
            //diss=population.calculateEmpDissimilarityMatrix(0);
            //ids=population.calculateEmpIDs();
            //dissSyl=population.calculateEmpSyllDissimilarityMatrix(0);
            //syllab=makeEmpSylLabel();
            //calculateStats();
            //calculateThresholdSpectrum(population.emppop);
            calculateSongSharing();
            calculateSongSharingAvg();
            calculateSharingAcrossPops();
            calculateSummaryAcrossPopulations();
            calculateSongDiversity();
            
            System.out.print(param.mutationVar+" "+param.novbias+" "+param.ntutors+" "+param.problearn1+" ");
            for (int i=0; i<outputavg.length; i++){
                System.out.print(outputavg[i]+" ");
            }
            for (int i=0; i<output7.length; i++){
                for (int j=0; j<output7[i].length; j++){
                    System.out.print(output7[i][j]+" ");
                }
            }
            System.out.println(overallsongsharing);
            output1=null;
            output2=null;
            output3=null;
            output4=null; 
            output6=null;
            output7=null;
                    
            
        }
        
        public WhaleMeasureStatistics(WhalePopulation population, EmpData ed, WhaleParameters param) {
		this.population=population;
		this.param=param;
		dsim=param.typeThresh;
                memorylength=param.loglength;
		ids=population.calculateEmpIDs();
		dissSyl=ed.empSylDiss;
		syllab=ed.syllab;
		diss=ed.empDiss;
                //calculateStats();

                
                
	}
        
       
        
/*
	public void calculateThresholdSpectrum(WhaleIndividual[] pop){
           output1= new int[pop.length][param.memorylength][];
           output2= new int[pop.length][param.memorylength][];
           indbuffer= new int[param.memorylength*pop.length];
           membuffer= new int[param.memorysize*pop.length];
           int ml=param.memorylength;
           for (int i=0; i<pop.length; i++){ //For each individual i in population length
               //if(i%100==0){
               //System.out.println("Individual = " + i);
              // }
               float[] mema=pop[i].getSongLog();
               for (int a=0; a<ml; a++){ //memlength //for each song in total memory
                  
                   int k=0;                    
                  int indexa=pop[i].getMemoryIndex(a); //get the index of song a    
                 // System.out.print(i+" "+pop.length+" ");
                  //for (int j=0; j<pop[i].ns; j++){System.out.print(mema[indexa*pop[i].ns+j]+" ");}
                  //System.out.println();
                  for (int j=0; j<pop.length; j++){ //for each other individual in the population
                      double pdiff=Math.abs(population.subpop[i]-population.subpop[j]);
                      if ((pdiff<2)||(pdiff>8)){
                      float[] memb=pop[j].getSongLog();
                      if(i!=j){ //as long as it is not individual i itself
                          for (int b=0; b<ml; b++){ 
                            int indexb=pop[j].getMemoryIndex(b);
                            //System.out.println(indexa+" "+indexb);
                            if(pop[i].matchSongs(mema, memb, indexa, indexb)){
                              indbuffer[k]=j;
                              membuffer[k]=b;
                              k++;   
                            }                     
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
        */
        
        public void calculateSharingAcrossPops(){
            int years=(memorylength/param.epochsperyear);
            output6=new double[population.subpopsize.length][population.subpopsize.length][years][years];
            for (int i=0; i<population.emppop.length; i++){
               int a=population.emppop[i].subpop; //get the population number
               float[] mema=population.emppop[i].getSongLog();
               for (int j=0; j<years; j++){ //for each year
                   int jj=j*param.epochsperyear+3;  //pick a time shortly after breeding starts/               
                   for (int k=0; k<population.emppop.length; k++){
                       int b=population.emppop[k].subpop;
                       float[] memb=population.emppop[k].getSongLog();
                       for (int l=0; l<years; l++){
                           int ll=l*param.epochsperyear+3;
                           if(population.emppop[i].matchSongs(mema, memb, jj, ll)){
                               output6[a][b][j][l]++;
                           }
                       }
                   }
                   
               }
            }
            /*
            //for (int i=0; i<population.subpopsize.length; i++){
                for (int j=0; j<population.subpopsize.length; j++){
                    for (int a=0; a<param.epochsperyear; a++){
                        for (int b=0; b<param.epochsperyear; b++){
                            System.out.println(" "+j+" "+a+" "+b+" "+output6[0][j][a][b]);
                        }
                    }
                }
            //}
            
               */
        }
        
        public void calculateSummaryAcrossPopulations(){
            output7=new double[3][2];
            for (int i=0; i<population.subpopsize.length; i++){
                int j1=i-1;
                if (j1<0){j1=population.subpopsize.length-1;}
                int j2=i+1;
                if (j2==population.subpopsize.length){j2=0;}
                for (int j=1; j<param.epochsperyear; j++){
                    output7[0][0]+=output6[i][j1][j][j-1];
                    output7[0][1]+=output6[i][j1][j][j];
                    output7[1][0]+=output6[i][i][j][j-1];
                    output7[1][1]+=output6[i][i][j][j];
                    output7[2][0]+=output6[i][j2][j][j-1];
                    output7[2][1]+=output6[i][j2][j][j];
                   // System.out.println(output6[i][j2][j][j]+" "+output6[i][j1][j][j]+" "+j1+" "+j2);
                }
            }   
        }
        
    public void calculateSongSharing(){
        output3= new int[population.emppop.length][memorylength][memorylength];
        for (int i=0; i<population.emppop.length; i++){ //for each individual!! 
          int a=population.emppop[i].subpop; //get the population number
          float[] mema=population.emppop[i].getSongLog();
          for (int j=0; j<population.emppop.length; j++){
              int b=population.emppop[j].subpop;
              if (a==b){
                  float[] memb=population.emppop[j].getSongLog();
                  for (int g=0; g<memorylength; g++){
                      for (int h=0; h<memorylength; h++){
                        if(population.emppop[i].matchSongs(mema, memb, g, h)){
                            output3[a][g][h]++;
                        }
                      }
                  }
              }
          }
        }
          
          
        output4 = new double[subpopsize.length][memorylength][memorylength];
        for (int i=0; i<subpopsize.length; i++){            //for every population
          double x=population.sampleperpop*(population.sampleperpop-0.0); //x=samplesize pop i * samplesize pop i
          
          for (int j=0; j<memorylength; j++){             //for every time point
            for (int k=0; k<memorylength; k++){           //for every other time point 
              //System.out.println("pop = " + i + " t = " + j + "sharing = " + output3[i][j][k]);
                output4[i][j][k]=output3[i][j][k]/x; 
              //System.out.println("pop = "+i+" timepoint = "+j+ " shared song time= " + k + " songsharing " +  output4[i][j][k]); 
              
              
            }
          }
        }
        
        
        }     
    
    public void calculateSongDiversity(){
        int n=population.emppop.length;
        double count=0;
        for (int i=0; i<n; i++){
            float[] mema=population.emppop[i].getSongLog();
            int iter=population.emppop[i].iter;
            
            for (int j=0; j<n; j++){
                float[] memb=population.emppop[j].getSongLog();
                if (population.emppop[i].matchSongs(mema, memb, iter, iter)){
                    count++;
                }
            }
            
        }
        count/=n*n*1.0;
        overallsongsharing=count;
    }
        
        public void calculateSongSharingAvg(){
            outputavg= new double[10];
            double[]count=new double[10];
            
            for(int j=0; j<population.subpopsize.length;j++){    //j=pop
                for(int k=0; k<memorylength;k++){ //timepoint a=k
                    int i=k%10;
                    outputavg[i]+=output4[j][k][k];
                    count[i]++;
                    //System.out.println(j+" "+k+" "+outputavg[i]+" "+count[i]+" "+i+" "+output4[j][k][k]);
                }
            }
            
            for (int i=0; i<maxlearn; i++){
                outputavg[i]/=count[i];
                //System.out.println(i+" "+outputavg[i]/count[i]);
            }
            
            /*

            for(int i=0; i<maxlearn; i++){ //timepoint i in yearly cycle of length maxlearn
                double sum=0;
                int x=0;    
                for(int j=0; j<population.subpopsize.length;j++){    //j=pop
                    for(int k=0; k<memorylength;k++){ //timepoint a=k
                        if(k%10==i){
                            for(int l=0; l<memorylength;l++){ //timepoint b=l
                                if(k==l){
                                sum+=output4[j][k][l];
                                x++;
                                } 
                            }
                        } 
                    }
                }
                outputavg[i]=sum/x;
                System.out.println("Time = " + i + " & sharing avg = " + outputavg[i]);
            }
            */
            
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
        