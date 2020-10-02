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
        
    
        int[][] samplesizes={{2,1,2,1,6,6,6,5,2,1,2},
            {5,6,6,6,5,10,4,6,3,6,6},
            {7,5,6,6,7,7,1,1,9,6,3},
            {1,0,0,3,0,4,2,1,1,5,1},
            {1,6,0,4,0,0,3,6,6,6,6}
        };
        
        int popoffset=6;

	WhalePopulation population;
        WhalePopulation emppop;
	WhaleParameters param;
        WhaleIndividual[] inds;
        Whale wh;
        WhaleIndividual individual;
        EmpData ed;
        int sampleperpop=100;
        int memorylength=100;
        int numSylls;
        int numdims;
        int memorySize;

	
      
        int maxlearn=10;
        int[] ids;
        int[] syllab;
        int[] indbuffer;
        int[] membuffer;
        int[][][] output1;
        int[][][] output2;
        int[][][] output3; //songsharing
        int[] subpopsize;
        int[][][] sharedsongs;
        double[][][] output4;
        double[][][][] output6; //SharingAcrossPops
        double[][] output7;
        double[] outputavg;
	double overallsongsharing=0;
        	double[][] diss, dissSyl;
	double dsim=0.04;
        int ns;
        int[][][] output8;
        float[] buffer;
        int[][][] songfreq;
        int[][][] themefreq; //frequency of each theme type
        int[][] themetot; //total number of themes
        int[][][] themepsong; //themes per song 
        double[][][] themestats; //summary stats for themes
        
        int[][] singlt;
	
	boolean typeEmpirical=true;
	
	public WhaleMeasureStatistics(WhalePopulation population, WhaleParameters param) {
            this.population=population;
            this.param=param; 
            this.memorylength=param.loglength;
            this.numSylls=param.sylsPerSong;
            this.numdims=param.numdims;
  //          this.individual=individual;
            typeEmpirical=true;
            dsim=param.typeThresh;
            this.subpopsize=population.subpopsize;
            memorySize=numSylls*numdims*memorylength;
            ns=numSylls*numdims;
            
            //System.out.println("WhaleMeasureStatistics");
            
            
            //diss=population.calculateEmpDissimilarityMatrix(0);
            //ids=population.calculateEmpIDs();
            //dissSyl=population.calculateEmpSyllDissimilarityMatrix(0);
            //syllab=makeEmpSylLabel();
            //calculateStats();
            //calculateThresholdSpectrum(population.emppop);
            //System.out.println("Starting...");
            WhaleIndividual[][][] sample=sampleIndividuals();
            double[] songfreqstats=calculateSongFreqStats(sample);
            
            double[] themefreqstats=calculateThemeFreqStats(sample);
            
            double[] songlengthstats=calculateSongLengthStats(sample);
            double[][] songsharingstats=calculateSongSharingStats(sample);
            
            System.out.println(param.mutationVar+" "+param.novbias+" "+param.problearn1+" "+param.ntutors+" "+param.dropparam+" "+param.probadd+" "+
                    songfreqstats[0]+" "+songfreqstats[1]+" "+themefreqstats[0]+" "+themefreqstats[1]+" "+songlengthstats[0]+" "+songlengthstats[1]+" "+
                    songlengthstats[2]+" "+songlengthstats[3]+" "+songsharingstats[0][0]+" "+songsharingstats[1][0]+" "+songsharingstats[2][0]+" "+
                    songsharingstats[0][1]+" "+songsharingstats[1][1]+" "+songsharingstats[2][1]);
            
            
            
            /*
            calculateSongSharing();
            System.out.println("Calculate song sharing");
            calculateSongSharingAvg();
            System.out.println("Calculate Song sharing avg");
            calculateSharingAcrossPops();
            System.out.println("Calc sharing across pops");
            calculateSummaryAcrossPopulations();
            System.out.println("Calculate summ across pops");
            calculateSongDiversity();
            System.out.println("Calculate song diversity");
            calculateSongTypes();
            System.out.println("Calculate song types");
            calculateThemeSharing();
            System.out.println("Calculate Theme sharing");
            calculateThemeTypes();
            System.out.println("Calculate Theme types");
            calculateThemeStats();
            System.out.println("Calculate Theme stats");
            
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
            */
            output1=null;
            output2=null;
            output3=null;
            output4=null; 
            output6=null;
            output7=null;
            output8=null;
            buffer=null;
            songfreq=null;
            themefreq=null; //frequency of each theme type
            themetot=null; //total number of themes
            themepsong=null; //themes per song 
            themestats=null; //summary stats for themes
                    
            
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
        
        public WhaleIndividual[][][] sampleIndividuals(){
            population.breeding=false;
            //System.out.println(samplesizes.length+" "+samplesizes[0].length);
            WhaleIndividual[][][] out=new WhaleIndividual[samplesizes.length][samplesizes[0].length][];
            
            for (int i=0; i<samplesizes.length; i++){
                for (int j=0; j<samplesizes[i].length; j++){
                    if (samplesizes[i][j]>0){
                        out[i][j]=population.getWhales(i+popoffset, samplesizes[i][j]);
                    }
                    else{
                        out[i][j]=new WhaleIndividual[0];
                    }
                }
            }
            return out; 
        }
        
        
        
        public double[] calculateSongFreqStats(WhaleIndividual[][][] sample){

            int n=sample.length*sample[0].length;
            double[] singletons=new double[n];
            double[] maxfreq=new double[n];
            double[] numinds=new double[n];
            int k=0;
            for (int i=0; i<sample.length; i++){
            
                for (int j=0; j<sample[i].length; j++){
                    WhaleIndividual[] yp=sample[i][j];
                    numinds[k]=yp.length;
                    int yearsago=10-j;
                    //System.out.println(yearsago);
                    for (int a=0; a<yp.length; a++){
                        float[] mema=yp[a].getSongLog();
                        int x=yp[a].iterlog-yearsago;
                        if (x<0){x+=yp[a].loglength;}
                        int count=0;
                        for (int b=0; b<yp.length; b++){
                            if (a!=b){
                                float[] memb=yp[b].getSongLog();
                                int y=yp[b].iterlog-yearsago;
                                if (y<0){y+=yp[b].loglength;}
                                //if ((x<0)||(y<0)){System.out.println(x+" "+y+" "+i+" "+j+" "+mema.length+" "+memb.length+" "+yearsago);}
                                if(yp[a].matchSongs(mema, memb, x, y)){
                                    count++;
                                }
                            }  
                        }
                        if (count==0){
                            singletons[k]++;
                        }
                        if (count>maxfreq[k]){
                            maxfreq[k]=count;
                        }
                    }
                    k++;
                }
            }
            double[] out=new double[2];
            
            double j=0;
            for (int i=0; i<n; i++){
                //System.out.println(singletons[i]+" "+maxfreq[i]+" "+numinds[i]);
                if (numinds[i]>0){
                    out[0]+=singletons[i];
                    out[1]+=maxfreq[i]/numinds[i];
                    j++;
                }
            }
            out[0]/=j;
            out[1]/=j;
            //System.out.println("SONGFREQS: "+out[0]+" "+out[1]);
            return out;
        }
        
        
        public double[] calculateThemeFreqStats(WhaleIndividual[][][] sample){

            int n=sample.length*sample[0].length;
            double[] singletons=new double[n];
            double[] maxfreq=new double[n];
            double[] numinds=new double[n];
            int k=0;
            for (int i=0; i<sample.length; i++){
            
                for (int j=0; j<sample[i].length; j++){
                    WhaleIndividual[] yp=sample[i][j];
                    int yearsago=10-j;
                    for (int a=0; a<yp.length; a++){
                        float[] mema=yp[a].getSongLog();
                        int x=yp[a].iterlog-yearsago;
                        if (x<0){x+=yp[a].loglength;} 
                        int p=yp[a].calculateSongLength(mema, x);
                        numinds[k]+=p;
                        int[] count=new int[p];
                        
                        for (int b=0; b<yp.length; b++){
                            if (a!=b){
                                float[] memb=yp[b].getSongLog();
                                int y=yp[b].iterlog-yearsago;
                                if (y<0){y+=yp[b].loglength;}
                                int[] count2=yp[a].countSharedThemes(mema, memb, x, y);
                                for (int s=0; s<p; s++){
                                    count[s]+=count2[s];
                                }
                            }  
                        }
                        for (int b=0; b<p; b++){
                            if (count[b]==0){
                                singletons[k]++;
                            }
                            if (count[b]>maxfreq[k]){
                                maxfreq[k]=count[b];
                            }
                        }
                    }
                    k++;
                }
            }
            double[] out=new double[2];
            
            double j=0;
            for (int i=0; i<n; i++){
                //System.out.println(singletons[i]+" "+maxfreq[i]+" "+numinds[i]);
                if (numinds[i]>0){
                    out[0]+=singletons[i];
                    out[1]+=maxfreq[i]/numinds[i];
                    j++;
                }
            }
            out[0]/=j;
            out[1]/=j;
            //System.out.println("THEMEFREQS: "+out[0]+" "+out[1]);
            return out;
        }
        
        public double[] calculateSongLengthStats(WhaleIndividual[][][] sample){
            int n=sample.length;
            int m=sample[0].length;
            int[][][] songlengths=new int[n][m][];
            for (int i=0; i<n; i++){
                for (int j=0; j<m; j++){
                    int yearsago=10-j;
                    int l=sample[i][j].length;
                    songlengths[i][j]=new int[l];
                    for (int k=0; k<l; k++){
                        int x=sample[i][j][k].iterlog-yearsago;
                        if (x<0){x+=sample[i][j][k].loglength;} 
                        songlengths[i][j][k]=sample[i][j][k].calculateSongLength(sample[i][j][k].getSongLog(), x);
                    }
                }
            }
            
            double[][]sampmeans=new double[n][m];
            double[][] sampcounts=new double[n][m];
            int mins=0;
            int maxs=0;
            double gm=0;
            double count=0;
            for (int i=0; i<n; i++){
                for (int j=0; j<m; j++){
                    sampcounts[i][j]=sample[i][j].length;
                    for (int k=0; k<sample[i][j].length; k++){
                        sampmeans[i][j]+=songlengths[i][j][k];
                        if (songlengths[i][j][k]==param.minThemes){mins++;}
                        if (songlengths[i][j][k]==param.maxThemes){maxs++;}
                    }
                    if (sample[i][j].length>0){
                        sampmeans[i][j]/=sample[i][j].length+0.0;
                    }
                    gm+=sampmeans[i][j];
                    if (sampcounts[i][j]>0){count++;}
                }
            }
            
            gm/=count;
            
            double[][] sampvar=new double[n][m];
            for (int i=0; i<n; i++){
                for (int j=0; j<m; j++){
                    if (sampcounts[i][j]>1){
                        for (int k=0; k<sample[i][j].length; k++){
                            double s=songlengths[i][j][k];
                            sampvar[i][j]+=s*s;
                        }
                        sampvar[i][j]*=1/(sampcounts[i][j]-1.0);
                    }
                }
            }
            
            double pv=0;
            double tc=0;
            for (int i=0; i<n; i++){
                for (int j=0; j<m; j++){
                    if (sampcounts[i][j]>1){
                        pv+=sampvar[i][j]*(sampcounts[i][j]-1.0);
                        tc+=sampcounts[i][j]-1.0;
                    }
                }
            }
            
            pv=Math.sqrt(pv/tc);
            
            
            double[]out={gm, pv, mins, maxs}; 
            //System.out.println("SONGLENGTH STATS: "+out[0]+" "+out[1]+" "+out[2]+" "+out[3]);                    
            
            return out;
        }
        
        public double[][] calculateSongSharingStats(WhaleIndividual[][][] sample){
            int n=sample.length;
            int m=sample[0].length;
            
            double[][] out=new double[3][2];
            double[][] count=new double[3][2];
            for (int i=0; i<n; i++){
                for (int j=0; j<m; j++){
                    WhaleIndividual[] q=sample[i][j];
                    if (q.length>0){
                        for (int a=-1; a<=1; a++){
                            int aa=i+a;
                            if ((aa>=0)&&(aa<n)){
                                for (int b=-1; b<=0; b++){
                                    int bb=j+b;
                                    if (bb>=0){
                                        WhaleIndividual[] r=sample[aa][bb];
                                        if (r.length>0){
                                            boolean self=false;
                                            if ((a==0)&&(b==0)){
                                                self=true;
                                            }
                                            if ((!self)||(q.length>1)){
                                                double sh=calcSharing(q,r, j, bb, self);
                                                out[a+1][b+1]+=sh;
                                                count[a+1][b+1]++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }  
                }
            }
            
            for (int i=0; i<3; i++){
                for (int j=0; j<2; j++){
                    if (count[i][j]>0){
                        out[i][j]/=count[i][j];
                    }
                }
            }
            return out; 
        }
        
        public double calcSharing(WhaleIndividual[]aa, WhaleIndividual[]bb, int c, int d, boolean self){
            
           WhaleIndividual[] a=aa;
           WhaleIndividual[] b=bb;
            
            int n=a.length;
            int m=b.length;
            
            
            /*
            if (n<m){
                a=bb;
                b=aa;
                n=a.length;
                m=b.length;
            }
            */
            
            int yearsago1=10-c;
            int yearsago2=10-d;
            double count=0;
            for (int i=0; i<n; i++){
                float[] mema=a[i].getSongLog();
                int x=a[i].iterlog-yearsago1;
                if (x<0){x+=a[i].loglength;}
                
                for (int j=0; j<m; j++){
                    if ((!self)||(i!=j)){
                        float[] memb=b[j].getSongLog();
                        int y=b[j].iterlog-yearsago2;
                        if (y<0){y+=b[j].loglength;} 
                        if (a[i].matchSongs(mema, memb, x,y)){
                            count++;
                            j=m;
                        }
                    }
                }
            }
            
            count=count/(n+0.0);//Jaccard Index of sharing!
            return count; 
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
       

    
        
        
        

    public void calculateSongTypes(){
            buffer= new float[sampleperpop*memorySize];
            songfreq=new int[subpopsize.length][memorylength][sampleperpop*memorylength];
            int bufferindex;
            
        for(int pop=0; pop<subpopsize.length;pop++){
            for(int t=0; t<memorylength;t++){
                buffer= new float[sampleperpop*memorySize];
                bufferindex=0;
                //for (int i=0; i<population.emppop.length; i++){ //for each individual!! (ID-A)
                for (int i=0; i<sampleperpop; i++){ //for each sampled individual...
                  int a=population.emppop[i].subpop; //get the population number of ID-A
                    if(a==pop){
                      float[] mema=population.emppop[i].getSongLog(); //Get its songlog
                       //if(bufferindex==0){
                       //     System.arraycopy(mema, t*ns, buffer, bufferindex*ns, ns); 
                       //     bufferindex++; 
                       // }
                       // else{
                            boolean matched=false;
                            for(int j=0; j<bufferindex; j++){
                                if(population.emppop[i].matchSongs(mema, buffer, t, j)){
                                    matched=true;
                                    songfreq[pop][t][j]++;
                                }
                            }
                            if (!matched){
                                    
                               // }
                                //else{
                                    System.arraycopy(mema, t*ns, buffer, bufferindex*ns, ns); 
                                    bufferindex++;
                                //}
                            //}
                        }      
                    }
                }
            }
        }
    }
    
    
    
        public void calculateThemeTypes(){
            buffer= new float[sampleperpop*memorySize];
            themefreq=new int[subpopsize.length][memorylength][sampleperpop*memorylength*numSylls];
            themetot=new int[subpopsize.length][memorylength];
            themepsong= new int [subpopsize.length][memorylength][sampleperpop];
            //themefreq=new int[subpopsize.length][memorylength][sampleperpop*memorylength*numSylls];
            int bufferindex;
            
        for(int pop=0; pop<subpopsize.length;pop++){
            int x=0;
            for(int t=0; t<memorylength;t++){
                buffer= null;
                bufferindex=0;
                for (int i=0; i<population.emppop.length; i++){ //for each individual!! 
                  int a=population.emppop[i].subpop; //get the subpopulation 
                    if(a==pop){
                      float[] mema=population.emppop[i].getSongLog(); //Get its songlog
                        for(int k=0; k<numSylls; k++){ //go through the themes in its song
                            if(mema[t*ns+k*numdims]!=-1000){
                                themepsong[pop][t][x]++;
                                if(bufferindex>0){
                                    for(int j=0; j<bufferindex; j++){
                                        themetot[pop][t]++; //total number of themes in the (sub)pop at time t
                                        if(population.emppop[i].matchThemes(mema, buffer, t*numSylls+k, j)){
                                        themefreq[pop][t][j]++;
                                        }
                                        else{
                                        System.arraycopy(mema, t*ns+k*numdims, buffer, bufferindex*numdims, numdims); 
                                        bufferindex++;
                                        }
                                    
                                } 
                            }
                            else if(bufferindex==0){
                                
                                    System.arraycopy(mema, t*numdims, buffer, bufferindex*numdims, numdims); 
                                    bufferindex++; 
                                }
                            }
                        }
                        x++;
                    }
                }
            }
        }
    }

    
    //each song can match multiple others!
        
        
    public void calculateSongSharing(){
        output3= new int[population.emppop.length][memorylength][memorylength];
        for (int i=0; i<population.emppop.length; i++){ //for each individual!! (ID-A)
          int a=population.emppop[i].subpop; //get the population number of ID-A
          float[] mema=population.emppop[i].getSongLog(); //Get its songlog
          for (int j=0; j<population.emppop.length; j++){ //go through the all others
              int b=population.emppop[j].subpop; //get its population (ID-B)
              if (a==b){ //if that population is the same
                  float[] memb=population.emppop[j].getSongLog(); //get songLog for ID-b
                  for (int g=0; g<memorylength; g++){ //go through mema: songIDa
                      for (int h=0; h<memorylength; h++){ //go through memb: songIDb
                        if(population.emppop[i].matchSongs(mema, memb, g, h)){ //match songs
                            output3[a][g][h]++; //fill in output: [popa][timepointA(=position in songlog)][timepointB]
                        }
                      }
                  }
              }
          }
        }
          
        //calculate average songsharing per pop  
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
        for (int i=0; i<n; i++){ //for each IDa  in emppop
            float[] mema=population.emppop[i].getSongLog(); //get songlog
            int iter=population.emppop[i].iter; //get iter for each id
            
            for (int j=0; j<n; j++){ //for each IDb
                float[] memb=population.emppop[j].getSongLog(); //
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
        
        public void calculateThemeSharing(){
          for (int i=0; i<population.emppop.length; i++){ //for each individual!! (ID-A)
            int a=population.emppop[i].subpop; //get the population number of ID-A
            float[] mema=population.emppop[i].getSongLog(); //Get its songlog: mema
                for (int j=0; j<population.emppop.length; j++){ //go through all others
                    int b=population.emppop[j].subpop; //get its population (ID-B)
                    if (a==b){ //if that population is the same
                        
                        float[] memb=population.emppop[j].getSongLog(); //get songLog for ID-b: memb
                            
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
        }
                
           public void calculateThemeStats(){
         for(int pop=0; pop<subpopsize.length;pop++){
            for(int t=0; t<memorylength;t++){           
                for(int i=0; i<sampleperpop*memorylength*numSylls; i++){
                    if(themefreq[pop][t][i]==0){
                        break;
                    }
                    else{
                    themestats[pop][t][5]++;
                    }
                }
            }
         }
                
          for(int pop=0; pop<subpopsize.length;pop++){
            for(int t=0; t<memorylength;t++){    
                double c = themestats[pop][t][5]*0.02;   
                for(int i=0; i<sampleperpop*memorylength*numSylls; i++){    
                    if(themefreq[pop][t][i]==0){
                        break;
                    }
                    else{
                        if(themefreq[pop][t][i]==1){
                            themestats[pop][t][1]++; //number of singletons
                        }
                        if(themefreq[pop][t][i]<4){
                            themestats[pop][t][2]++; //number of rare themes
                        }
                        if(themefreq[pop][t][i]>1&&themefreq[pop][t][i]<c){ 
                            themestats[pop][t][3]++;    //number of intermediate themes
                        }
                        if(themefreq[pop][t][i]>c){ 
                            themestats[pop][t][4]++;    //numbers of common themes
                        }
                        if(themefreq[pop][t][i]>themestats[pop][t][6]){
                            themestats[pop][t][6]=themefreq[pop][t][i]; //number of IDs singing the most common theme
                        }  
                        themestats[pop][t][7]=+(themefreq[pop][t][i]/themetot[pop][t])*(Math.log(themefreq[pop][t][i]/themetot[pop][t]) / Math.log(2)); //H index
                        themestats[pop][t][8]=+Math.log(2*(themefreq[pop][t][i]/themetot[pop][t]));
                    }
                   
                }
                // convert to proportion of all themes in subpop pop at time t;
                themestats[pop][t][1]/=themestats[pop][t][5];
                themestats[pop][t][2]/=themestats[pop][t][5];
                themestats[pop][t][3]/=themestats[pop][t][5];
                themestats[pop][t][4]/=themestats[pop][t][5];
                themestats[pop][t][8]=1+themetot[pop][t]*(Math.pow(themestats[pop][t][8],-1)); //Alpha P
                
                c=0;
               
                for(int i=0; i<sampleperpop; i++){
                    themestats[pop][t][9]+=themepsong[pop][t][i];
                    if(themefreq[pop][t][i]==2){
                            themestats[pop][t][10]++; //songs with 2 themes (minimal)
                    }
                    if(themefreq[pop][t][i]==numSylls){
                            themestats[pop][t][11]++; //songs with max themes
                    }
                    if(themefreq[pop][t][i]==Math.round(0.5*numSylls)+1){
                            themestats[pop][t][12]++; //songs with intermediate themes
                    }
                    
                    
                }
                themestats[pop][t][9]/=sampleperpop; //average themes per song
                        
            
            
            
            
            }     
          }
          
          
    }
  
        
        

/* 
           public void calculateSongSimMatrix(){
        output8= new double[subpopsize.length][memorylength][sampleperpop][sampleperpop];
        for (int i=0; i<population.emppop.length; i++){ //for each individual!! (ID-A)
          int a=population.emppop[i].subpop; //get the population number of ID-A
          float[] mema=population.emppop[i].getSongLog(); //Get its songlog
          for (int j=0; j<population.emppop.length; j++){ //go through the all others
              int b=population.emppop[j].subpop; //get its population (ID-B)
              if (a==b){ //if that population is the same
                  float[] memb=population.emppop[j].getSongLog(); //get songLog for ID-b
                  for (int g=0; g<memorylength; g++){ //go through mema: songIDa
                      for (int h=0; h<memorylength; h++){ //go through memb: songIDb
                            output8[a][g][i][j]=population.emppop[i].compareSongsPow(mema, memb, g, h); //compare songs
                            
                        }
                      }
                  }
              }
          }
        }    
        
        
        
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
        
        
        
    public void calculateSingletons(){
        int y;
        double x=population.sampleperpop*(population.sampleperpop-0.0);
        singlt= new int[subpopsize.length][memorylength];
        //rare= new int[subpopsize.length][memorylength];
        
        for (int i=0; i<subpopsize.length; i++){
            for(int j=0; j<memorylength; j++){
                for(int k=0; k<sampleperpop; k++){
                    y=0;
                    for(int l=0; l<sampleperpop; l++){
                        //if(output8[i][j][k][l]<matchthresh){
                            y++;
                        }
                        if(y==1){
                            singlt[i][j]++;
                 }
                }
            }
        }
            
     }
        
            
    public void calculateThemepSongStats(){
        for(int pop=0; pop<subpopsize.length;pop++){
            for(int t=0; t<memorylength;t++){ 
                for(int i=0; i<sampleperpop; i++){
                    themestats[pop][t][9]+=themepsong[pop][t][i];
                    if(themefreq[pop][t][i]==2){
                            themestats[pop][t][10]++; //songs with 2 themes (minimal)
                    }
                    if(themefreq[pop][t][i]==numSylls){
                            themestats[pop][t][11]++; //songs with max themes
                    }
                    if(themefreq[pop][t][i]==Math.round(0.5*numSylls)+1){
                            themestats[pop][t][12]++; //songs with intermediate themes
                    }
                    
                    
                }
                themestats[pop][t][9]/=sampleperpop; //average themes per song
            }
        }
    }  
*/        
        
        

}
                        
        