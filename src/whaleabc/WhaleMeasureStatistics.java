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



import java.util.Arrays;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.function.Logistic;
import org.apache.commons.math3.fitting.SimpleCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

public class WhaleMeasureStatistics {
        
    
        int[][] samplesizes={{2,1,2,1,6,6,6,5,2,1,2},         //RIGHT WAY ROUND!!!
            {5,6,6,6,5,10,4,6,3,6,6},
            {7,5,6,6,7,7,1,1,9,6,3},
            {1,0,0,3,0,4,2,1,1,5,1},
            {1,6,0,4,0,0,3,6,6,6,6}
        };
        
        
        int[][] samplesizesnh={{10,10,10,10,10,10,10,10,10,10},         //RIGHT WAY ROUND!!!
            {10,10,10,10,10,10,10,10,10,10},
            {10,10,10,10,10,10,10,10,10,10},
            {10,10,10,10,10,10,10,10,10,10},
            {10,10,10,10,10,10,10,10,10,10},
            {10,10,10,10,10,10,10,10,10,10}
        };
        
        
        
        
        int popoffset=5;
        int popoffsetnh=0;
        /*
        int[][] samplesizes={{1,6,0,4,0,0,3,6,6,6,6},           //REVERSED!!!
            {1,0,0,3,0,4,2,1,1,5,1},
            {7,5,6,6,7,7,1,1,9,6,3},
            {5,6,6,6,5,10,4,6,3,6,6},
            {2,1,2,1,6,6,6,5,2,1,2}
        };
        int popoffset=1;
    */
	WhalePopulation population;
	WhaleParameters param;
        
        double[] out;
        double[] out1, out2, out3, out4;
       
        
        public WhaleMeasureStatistics(WhalePopulation population, WhaleParameters param, boolean w){
            this.population=population;
            this.param=param;
            if (param.hemisphere==1){
                popoffset=popoffsetnh;
                samplesizes=samplesizesnh;
            }
            
            WhaleIndividual[][][] sample=sampleIndividuals();
            int ns=param.sylsPerSong;
            /*
            //for (int i=0; i<sample.length; i++){
            for (int i=0; i<1; i++){
                System.out.println(population.currentpopsizes[i]+" "+population.tutorpopsizes[i]);
                for (int j=0; j<sample[i].length; j++){
                    for (int k=0; k<sample[i][j].length; k++){
                        System.out.print(i+" "+j+" "+k+" "+sample[i][j][k].age+" "+sample[i][j][k].territory+" "+sample[i][j][k].learningEpochs+" ");
                        float[] x=sample[i][j][k].getSongLog();
                        int y=(sample[i][j][k].iterlog-j-1)*ns;
                        if (y<0){y+=sample[i][j][k].loglength*ns;}
                        
                        //System.out.println("ch: "+k+" "+y+" "+ns+" "+x.length+" "+sample[i][j][k].iterlog);
                        
                        for (int a=0; a<ns; a++){
                            System.out.print(x[a+y]+" ");
                        }
                        System.out.println();
                    }
                }
            }
            */
            
            double[][] themesharingstats=calculateThemeSharingStats(sample);
            double[] ts={themesharingstats[0][0]-themesharingstats[2][0], 
                themesharingstats[1][1]-themesharingstats[1][0], Math.max(themesharingstats[0][0], themesharingstats[2][0]), themesharingstats[1][1]}; 
            out4=ts;
            //System.out.println(themesharingstats[0][0]+" "+themesharingstats[2][0]+" "+themesharingstats[1][1]);
            /*
            sample=sampleIndividualsAll();
            themesharingstats=calculateThemeSharingStatsByPop(sample);
            out1=new double[sample.length];
            out2=new double[sample.length];
            out3=new double[sample.length];
            for (int i=0; i<sample.length; i++){
                out1[i]=themesharingstats[i][0];
                out2[i]=themesharingstats[i][2];
                out3[i]=themesharingstats[i][4];
            }
            */
            //themesharingstats=calculateThemeSharingStats(sample);
            //calculateThemeAgeAndSpread(sample);
            //System.out.println(themesharingstats[0][0]+" "+themesharingstats[2][0]);
            //out=derivedthemesharingstats;
        }

	public WhaleMeasureStatistics(WhalePopulation population, WhaleParameters param) {
            this.population=population;
            this.param=param; 
            if (param.hemisphere==1){
                popoffset=popoffsetnh;
                samplesizes=samplesizesnh;
            }

            WhaleIndividual[][][] sample=sampleIndividuals();
            double[] songfreqstats=calculateSongFreqStats(sample);
            
            double[] themefreqstats=calculateThemeFreqStats(sample);
            
            double[] songlengthstats=calculateSongLengthStats(sample);
            double[][] songsharingstats=calculateSongSharingStats(sample);
            double[][] themesharingstats=calculateThemeSharingStats(sample);
            double[] derivedthemesharingstats={Math.abs(themesharingstats[0][0]-themesharingstats[2][0]), 
                themesharingstats[1][1]-themesharingstats[1][0], Math.max(themesharingstats[0][0], themesharingstats[2][0])};
        
            //System.out.println(param.mutationVar+" "+param.novbias+" "+param.problearn1+" "+param.ntutors+" "+param.dropparam+" "+param.probadd+" "+
            //        songfreqstats[0]+" "+songfreqstats[1]+" "+themefreqstats[0]+" "+themefreqstats[1]+" "+songlengthstats[0]+" "+songlengthstats[1]+" "+
            //        songlengthstats[2]+" "+songlengthstats[3]+" "+songsharingstats[0][0]+" "+songsharingstats[1][0]+" "+songsharingstats[2][0]+" "+
            //        songsharingstats[0][1]+" "+songsharingstats[1][1]+" "+songsharingstats[2][1]+" "+themesharingstats[0][0]+" "+themesharingstats[1][0]+" "+themesharingstats[2][0]+" "+
            //        themesharingstats[0][1]+" "+themesharingstats[1][1]+" "+themesharingstats[2][1]);

            //System.out.println(param.mutationVar+" "+param.novbias+" "+param.problearn1+" "+param.ntutors+" "+param.dropparam+" "+param.probadd+" "+
             //       themefreqstats[0]+" "+" "+songlengthstats[0]+" "+" "+themesharingstats[0][0]+" "+themesharingstats[1][0]+" "+themesharingstats[2][0]+" "+
            //        themesharingstats[0][1]+" "+themesharingstats[1][1]+" "+themesharingstats[2][1]);
       
            //for (int i=0; i<population.currentpopsizes.length; i++){
            //    System.out.print(population.currentpopsizes[i]+" ");
            //}
            //System.out.println();
            out=new double[23];
            System.arraycopy(songfreqstats, 0, out, 0, 2);
            System.arraycopy(themefreqstats, 0, out, 2, 2);
            System.arraycopy(songlengthstats, 0, out, 4, 4);
            System.arraycopy(songsharingstats[0], 0, out, 8, 2);
            System.arraycopy(songsharingstats[1], 0, out, 10, 2);
            System.arraycopy(songsharingstats[2], 0, out, 12, 2);
            System.arraycopy(themesharingstats[0], 0, out, 14, 2);
            System.arraycopy(themesharingstats[1], 0, out, 16, 2);
            System.arraycopy(themesharingstats[2], 0, out, 18, 2);
            System.arraycopy(derivedthemesharingstats, 0, out, 20, 3);
            
            //required stats for inference: 2, 4, 17, 20, 21, 22
            //System.out.println(themefreqstats[0]);
            //System.out.println("Sharing: "+derivedthemesharingstats[0]+" "+derivedthemesharingstats[1]+" "+derivedthemesharingstats[2]+" "+themesharingstats[1][1]);
            if (derivedthemesharingstats[0]>0.3){
                sample=sampleIndividualsAll(11,11);
                calculateThemeAgeAndSpread(sample);
            }
        }
        
        public WhaleIndividual[][][] sampleIndividuals(){
            population.breeding=false;
            //System.out.println(samplesizes.length+" "+samplesizes[0].length);
            WhaleIndividual[][][] out=new WhaleIndividual[samplesizes.length][samplesizes[0].length][];
            
            for (int i=0; i<samplesizes.length; i++){
                for (int j=0; j<samplesizes[i].length; j++){
                    if (samplesizes[i][j]>0){
                        //System.out.println(i+" "+popoffset+" "+samplesizes[i][j]);
                        out[i][j]=population.getWhales(i+popoffset, samplesizes[i][j]);
                        //System.out.println(i+" "+j+" "+out[i][j].length);
                    }
                    else{
                        out[i][j]=new WhaleIndividual[0];
                    }
                }
            }
            return out; 
        }
        
        public WhaleIndividual[][][] sampleIndividualsAll(int a, int b){
            population.breeding=false;
            
            WhaleIndividual[][][] out=new WhaleIndividual[population.subpopsize.length][a][];
            //System.out.println(out.length+" "+out[0].length);
            for (int i=0; i<population.subpopsize.length; i++){
                for (int j=0; j<a; j++){
                    out[i][j]=population.getWhales(i, b);
                }
            }
            return out; 
        }
        
        public void traceThemes(int id){
            WhaleIndividual[][][] sample=sampleIndividualsAll(1, 5);
            for (int i=0; i<sample.length; i++){
                for (int k=0; k<sample[i][0].length; k++){
                    WhaleIndividual w=sample[i][0][k];
                    float[] x=w.getSongLog();
                    int b=w.iterlog;
                    for (int a=0; a<w.loglength; a++){
                        b--;
                        if (b<0){b=w.loglength-1;}
                        for (int c=0; c<w.ns; c++){
                            int d=b*w.ns+c;
                            if (x[d]!=-1000){
                                System.out.println(i+" "+a+" "+k+" "+x[d]);
                            }
                        }
                    }
                }
            }
        }

        
        public double[] calculateSongFreqStats(WhaleIndividual[][][] sample){

            int n=sample.length*sample[0].length;
            double[] singletons=new double[n];
            double[] maxfreq=new double[n];
            double[] numinds=new double[n];
            int k=0;
            for (int i=0; i<sample.length; i++){
            
                for (int j=0; j<sample[i].length; j++){
                    if (sample[i][j].length>1){
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
            }
            double[] out=new double[2];
            
            double j=0;
            for (int i=0; i<k; i++){
                //System.out.println(singletons[i]+" "+maxfreq[i]+" "+numinds[i]);
                if (numinds[i]>0){
                    out[0]+=singletons[i];
                    out[1]+=(1+maxfreq[i])/numinds[i];
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
                    if (sample[i][j].length>1){
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
                            //if (a!=b){
                                float[] memb=yp[b].getSongLog();
                                int y=yp[b].iterlog-yearsago;
                                if (y<0){y+=yp[b].loglength;}
                                int[] count2=yp[a].countSharedThemes(mema, memb, x, y);
                                for (int s=0; s<p; s++){
                                    count[s]+=count2[s];
                                }
                            //}  
                        }
                        for (int b=0; b<p; b++){
                            if (count[b]==1){
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
                            double s=songlengths[i][j][k]-sampmeans[i][j];
                            sampvar[i][j]+=s*s;
                        }
                        //sampvar[i][j]*=1/(sampcounts[i][j]-1.0);
                    }
                }
            }
            
            double pv=0;
            double tc=0;
            for (int i=0; i<n; i++){
                for (int j=0; j<m; j++){
                    if (sampcounts[i][j]>1){
                        //pv+=sampvar[i][j]*(sampcounts[i][j]-1.0);
                        pv+=sampvar[i][j];
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
            
            count=count/(n+0.0);
            return count; 
        }
        
        public double[][] calculateThemeSharingStats(WhaleIndividual[][][] sample){
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
                                                double sh=calcThemeSharing(q,r, j, bb, self);
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
        
        
        
        
        public void calculateThemeAgeAndSpread(WhaleIndividual[][][] sample){
            int n=sample.length;
            int m=sample[0].length;
            
            int np=11;
            
            double[][] shared=new double[11][np];
            double[] count=new double[np];
            
            for (int i=0; i<n; i++){
                for (int j=0; j<sample[i][0].length; j++){
                    WhaleIndividual a=sample[i][0][j];
                    float[]y=a.repertoire;
                    for (int k=0; k<a.ns; k++){
                        if (y[k]!=-1000){
                            for (int ii=0; ii<n; ii++){
                                int dist=i-ii+popoffset;
                                if (dist<0){dist+=np;}
                                if (dist>=np){dist-=np;}
                                for (int kk=0; kk<sample[ii][0].length; kk++){
                                    WhaleIndividual b=sample[ii][0][kk];
                                    float[] x=b.getSongLog();
                                    int maxAge=0;
                                    for (int h=0; h<x.length; h++){
                                        if (y[k]==x[h]){
                                            int p=b.getAge(h);
                                            if (p>maxAge){maxAge=p;}
                                            shared[p][dist]++;
                                        }
                                    }
                                    count[dist]++;
                                }
                            }
                        }   
                    }
                }
            }
            System.out.println();
            for (int i=0; i<shared.length; i++){
                for (int j=0; j<shared[0].length; j++){
                    System.out.println(i+" "+j+" "+(shared[i][j]/count[j]));
                }
                
            }
            
   
        }
        
        
        
        public double[][] calculateThemeSharingStatsByPop(WhaleIndividual[][][] sample){
            int n=sample.length;
            
            double[][] out=new double[n][6];
            double[][] count=new double[n][6];
            for (int i=0; i<n; i++){
                for (int j=1; j<2; j++){
                    WhaleIndividual[] q=sample[i][j];
                    if (q.length>0){
                        for (int a=-1; a<=1; a++){
                            int aa=i+a;
                            if (aa<0){aa+=n;}
                            if (aa==n){aa-=n;}
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
                                            double sh=calcThemeSharing(q,r, j+9, bb+9, self);
                                            out[i][(b+1)*3+a+1]+=sh;
                                            count[i][(b+1)*3+a+1]++;   
                                        }
                                    }
                                }
                            }
                        }
                    }  
                }
            }
            
            for (int i=0; i<n; i++){
                for (int j=0; j<6; j++){
                    if (count[i][j]>0){
                        out[i][j]/=count[i][j];
                    }
                }
            }
            return out; 
        }
        
        public double calcThemeSharing(WhaleIndividual[]aa, WhaleIndividual[]bb, int c, int d, boolean self){
            
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
            //int yearsago1=c;
            //int yearsago2=d;
            
            //int ns=param.sylsPerSong;
            double count=0;
            for (int i=0; i<n; i++){
                float[] mema=a[i].getSongLog();
                int x=(a[i].iterlog-yearsago1);
                if (x<0){x+=a[i].loglength;}
                double c1=0;
                double c2=0;
                for (int j=0; j<m; j++){
                    if ((!self)||(i!=j)){
                        float[] memb=b[j].getSongLog();
                        int y=(b[j].iterlog-yearsago2);
                        if (y<0){y+=b[j].loglength;} 
                        int[] count2=a[i].countSharedThemes(mema, memb, x, y);
                        double s=0;
                        for (int k=0; k<count2.length; k++){
                            s+=count2[k];
                        }
                        c1+=s/(count2.length+0.0);
                        c2++;
                    }
                }
                count+=c1/c2;
            }
            
            count=count/(n+0.0);//Jaccard Index of sharing!
            return count; 
        }
        
        
        
}
                        
        
