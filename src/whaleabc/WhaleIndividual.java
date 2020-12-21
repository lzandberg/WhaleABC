/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package whaleabc;

/**
 *
 * @author robertlachlan lies zandberg
 */


import java.util.Arrays;
import java.util.Random;

//This is the main class for objects representing individuals - including all song learning behaviour

public class WhaleIndividual  {
        int minThemes=2;
        int maxThemes=10; //parameters?
        int age=0;
	float[] repertoire, newRepertoire;
        int subpop=0;
        public float[] songmemory,songlog;
        int memoryActual;
	int territory;
	WhalePopulation population;
	WhaleParameters param;

	double matchThresh;
	double mutationVariance=0.1;
        
        int numSylls;
        int numdims;
        int memorylength, loglength; //memorylength = 50 songs, loglength = 100 songs
        int ns=0;
        int memorysize=0; //numsongs*numsyllables*numdimensions
        int logsize=0;
        int iter=0;
        int iterlog=0;
        
        double probadd=0.01; //parameters? probability of adding a completely new theme
        double novbias = 1;
        public double dropLookUp[];
	
	float[] songBuffer;
        double[] tutsongsim;
        float[] songsim; 
	double[] cumFreq;
	double[] powerLookUp;
        int repSizes[];
        int maxSongLength;
	
	int songtypeCount=0;
        int songthemeCount=0;
        int songtypefreq=0;
        float min=0;
        float max=10;
        int learningEpochs=0;
	        
        public WhaleIndividual(){
            
        }
        
	public WhaleIndividual(int territory, WhaleParameters param){
		//System.out.println("WhaleIndividual");	
		this.territory=territory;
		//this.random=random;
		
		this.param=param;
			
		this.numSylls=param.sylsPerSong;
                this.numdims=param.numdims;
		this.mutationVariance=param.mutationVar;
		this.matchThresh=param.typeThresh*param.typeThresh;
                
                this.memorylength=param.memorylength;
                this.loglength=param.loglength;
                
                this.novbias=param.novbias*0.5; //0.5 because in compareSongsPow, we calculate the squared dissimilarity.
                this.minThemes=param.minThemes;
                this.maxThemes=param.maxThemes;
                
                this.dropLookUp=param.powlookUp;
                this.probadd=param.probadd;
	
                ns=numSylls*numdims;
                memorysize=numSylls*numdims*memorylength;
                logsize=numSylls*numdims*loglength;
		
		repertoire=new float[ns];
		newRepertoire=new float[ns];	
		songBuffer=param.songBuffer;
                tutsongsim=param.tutsongsim;
		cumFreq=param.cumFreq;
		
		
	}
        
        
        
        public void initiate(){
            initiateRepertoire(population.subpop[territory]);
            updateRepertoire();
            initiateMemory();
        }
	
	public void setPopulation(WhalePopulation population) {
		this.population=population;
	}
        
	
	//at beginning of simulation run, repertoires are initiated with randomly selected values.
	public void initiateRepertoire(int subpop){
            int k=0;
            for (int j=0; j<ns; j++) { 
                newRepertoire[k]=subpop; //newRepertoire[k]=ind.subpop number: each song is ns float numbers each float is the same: the subpopulation number?
                k++;
            }
	}
        
        public void initiateMemory(){
            songmemory = new float[memorysize];  //make float[] for memory, with length memorysize (#songs * #dims * #maxthemes)
            for(int i=0; i<memorysize; i++){ //fill up memory with that one song
                songmemory[i] = newRepertoire[i%ns];
            }
            songlog=new float[logsize];
            for(int i=0; i<logsize; i++){ // also fill up songlog with that one song
                songlog[i] = newRepertoire[i%ns];
            }

	}
        
       
	
	//at the end of each year, each replaced individual's song values are updated to their newly calculated ones
	public void updateRepertoire(){
            System.arraycopy(newRepertoire, 0, repertoire, 0, ns);
	}
        
        
        
        public void learnSongs(boolean logcounter){
            WhaleIndividual[] tutors=population.getTutors(territory); //in this case territory is individual ID
            constructMemory(tutors);
            if (songtypeCount>0){
                makeCumDistr();
                pickSongs();
                dropTheme();
                addTheme();
            }
            mutate();
     
            System.arraycopy(newRepertoire, 0, songmemory, iter*ns, ns); 
            iter++;
            if (iter>=memorylength){
                iter=0;
            }
               
            if (logcounter){
                System.arraycopy(newRepertoire, 0, songlog, iterlog*ns, ns);
                iterlog++;
                if (iterlog>=loglength){
                    iterlog=0;
                }
            } 
            age++;
            memoryActual++;
            if (memoryActual>memorylength){memoryActual=memorylength;}
	}


		
    public boolean matchSongs(float[] x, float[] y, int a, int b) {
                double d=0;
                double f;
                double temp=0;
                int m=0; //number of themes in song A
                
                int aa=a*ns;
                int bb=b*ns;
                int a2, b2, a3;
		for (int i=0; i<numSylls; i++) { //Go through the themes of song a
                    
                    a2=aa+i*numdims; //Start of theme[i] of song a
                    if(x[a2]!=-1000){ //if theme in song a is not -1000:
                        f=100000;
                        m++;
                        for (int j=0; j<numSylls; j++) { //Go through themes of song b and find the most similar theme
                            b2=bb+j*numdims; //Start of theme[j] song b
                            if (y[b2]!=1000){
                                a3=a2;
                                temp=0;
                                for(int k=0; k<numdims; k++){ // Go through the dimensions of theme[i] of song a and theme[j] of song b
                                    temp+=(x[a3]-y[b2])*(x[a3]-y[b2]); 
                                    a3++;
                                    b2++;
                                }
                                if (temp<f){
                                    f=temp;                    
                                }
                            }
                        }
                        d+=f;
                        if (d>matchThresh*m){return false;}
                    }
		}
                
        return true;
    }
    
    
    public int getAge(int p){
        int x=(int)Math.floorDiv(p, ns);
        int y=iterlog-x-1;
        if (y<0){y+=loglength;}
        return y;
    }
    
    
        public int[] countSharedThemes(float[] x, float[] y, int a, int b) {
                double f;
                double temp=0;
                
                int aa=a*ns;
                int bb=b*ns;
                int a2, b2, a3;
                
                int p=calculateSongLength(x, a);
                
                int[] count=new int[p];
                p=0;
		for (int i=0; i<numSylls; i++) { //Go through the themes of song a
                    a2=aa+i*numdims; //Start of theme[i] of song a
                    if(x[a2]!=-1000){ //if theme in song a is not -1000:
                        f=100000;
                        for (int j=0; j<numSylls; j++) { //Go through themes of song b and find the most similar theme
                            b2=bb+j*numdims; //Start of theme[j] song b
                            if (y[b2]!=1000){
                                a3=a2;
                                temp=0;
                                for(int k=0; k<numdims; k++){ // Go through the dimensions of theme[i] of song a and theme[j] of song b
                                    temp+=(x[a3]-y[b2])*(x[a3]-y[b2]); 
                                    a3++;
                                    b2++;
                                }
                                if (temp<f){
                                    f=temp;                    
                                }
                            }
                        }
                        if (f<matchThresh){count[p]=1;}
                        p++;
                    }
		}
                
            return count;
        }
        
        
        public int calculateSongLength(float[]x, int a){
            int aa=a*ns;
            int a2;
            int count=0;
            for (int i=0; i<numSylls; i++) { //Go through the themes of song a
                a2=aa+i*numdims; //Start of theme[i] of song a
                if(x[a2]!=-1000){
                    count++;
                }
            }
            return count;
        }

        public boolean matchThemes(float[] x, float[] y, int a, int b) {
                //System.out.println("a = " + a + " & b = " + b + " & ns = " + ns);
		
                double d=0;
		int aa=a*numdims;
                int bb=b*numdims;
		for (int i=0; i<numdims; i++) {
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
                    aa++;
                    bb++;
                    if (d>(matchThresh)){
                        return false;
                    }
		}
		return true;
	}

        public double compareSongs(float[] x, float[] y, int a, int b) { //might need some optimalisation??
                double d=0;
                double f;
                double temp=0;
                int m=0;
                //int aa=a*ns;
                //int bb=b*ns;
                int a2, b2;
                int aa=a*ns;
                int bb=b*ns;
		for (int i=0; i<numSylls; i++) { //Go through the themes of song a  
                    if(x[aa]!=-1000){ //if theme in song a is not -1000:
                        m++;
                        f=100000;
                        b2=bb;
                        for (int j=0; j<numSylls; j++) { //Go through themes of song b and find the most similar theme 
                            if (y[b2]!=-1000){
                                temp=0;
                                a2=aa;
                                for(int k=0; k<numdims; k++){ // Go through the dimensions of theme[i] of song a and theme[j] of song b
                                    temp+=(x[a2]-y[b2])*(x[a2]-y[b2]); 
                                    a2++;
                                    b2++;
                                }
                                if (temp<f){
                                    f=temp;                    
                                }
                            }
                            //b2+=numdims; //Start of theme[j] song b
                        }
                        d+=f;                                                            
                    }
                    aa+=numdims; //Start of theme[i] of song a
		}
                //double p=pow(d/(0.0+m),novbias);                                            
                //System.out.println(d+" "+m);
                
                
                return d/(0.0+m);
                
                //double p=lookUpPower(d/(0.0+m));
                	
                
	}

	
        /*
        public int getMemoryIndex(int a){
            //a=a;
            //System.out.println("iter = " + iter + " & a = " + a);
            index=iter-a;
            if(index<0){
                index=(memorylength)+index;
            }
            //System.out.println("index = " + index);
            return index;
        }
         */           
	
	public float[] getSongMemory(){
            return songmemory;
        }
        

        public float[] getSongLog(){
            return songlog;
        }
        
        
	public void constructMemory(WhaleIndividual[] tutors) { // constructs memory of tutors songs of which one song will eventually be picked
            songtypeCount=0;
            double mindist,x;
  
            int tl=tutors.length; // number of tutors
            
            for (int i=0; i<tl; i++) { //for each tutor               
                mindist=Double.MAX_VALUE;
                for (int k=0; k<memoryActual; k++) {  
                    x=compareSongs(tutors[i].repertoire, songmemory, 0, k); 
                    if (x<mindist){
                        mindist=x;
                    } 
                }
                if (mindist>matchThresh){
                    tutsongsim[songtypeCount]=Math.pow(mindist, novbias);//Change ns
                    System.arraycopy(tutors[i].repertoire, 0, songBuffer, songtypeCount*ns, ns); //add it to the songbuffer
                    songtypeCount++;
                }      
            } 
        }
        
        public void makeCumDistr(){
            cumFreq[0]=tutsongsim[0];
            for (int i=1; i<songtypeCount; i++) {
		cumFreq[i]=tutsongsim[i]+cumFreq[i-1];
            }
        }

	public void pickSongs(){
            int c2=songtypeCount-1;
            double v=param.nextDouble()*cumFreq[c2];
            for (int j=0; j<songtypeCount; j++){
		if (v<cumFreq[j]){
                    System.arraycopy(songBuffer, j*ns, newRepertoire, 0, ns); //copy to repertoire
                }
            }
	}

        public void dropTheme(){
            int m=0;
      
            for (int i=0; i < numSylls; i++) { //Counting the number of themes
                if(newRepertoire[i*numdims]!=-1000){ 
                    m++;
                }
            }
            
            
            for (int i=0; i < numSylls; i++) { //for each theme in a newRepertoire  
                if(m>minThemes){
                    if(newRepertoire[i*numdims]!=-1000){
                        int k=0;
                        int p=numSylls*songtypeCount;
                        //System.out.println(numSylls+" "+songtypeCount);
                        for (int j=0; j<p; j++) {  //go through each theme in the songbuffer (not all tutors repertoires...)
                            if(matchThemes(newRepertoire, songBuffer, i, j)){ 
                                k++;
                            }
                        }
                        //System.out.println(k);
                        double dropprob=dropLookUp[k]; //probability of learning a theme (not dropping the theme)- 
                        // dependent on the frequency of each theme (k) and the number of themes in the tutor song //??
                        if (param.nextDouble()>dropprob){ //
                            //System.out.println("dropTheme: k= " + k + " dropprob= " + dropprob + " chance= " + chance)  ;
                            for(int l=0; l<numdims;l++){
                                newRepertoire[(i*numdims)+l]=-1000; //drop that theme (all dimensions revert to -1000)
                            }
                            m--;
                        }
                    }
                    
                }
            }
        }

                
        public void addTheme(){
            int m=0;
      
            for (int i=0; i < numSylls; i++) { //Counting the number of themes
                if(newRepertoire[i*numdims]!=-1000){ 
                    m++;
                }
            }
            if (m<maxThemes){ 
                for(int i=0; i<numSylls; i++){ //for each theme in new repertoire
                    
                    if(newRepertoire[i*numdims]==-1000){
                        if(param.nextDouble()<probadd){
                            //System.out.println("addTheme");
                            for(int j=0; j<numdims;j++){
                                newRepertoire[i*numdims+j]=param.nextFloat(min, max); 
                            //System.out.println("New theme is [" + j + "] " + newRepertoire[i*numdims+j]);
                            }
                        }
                    }
                
                }
            }
        }

	public void mutate(){
            //float m=0;
            
            //int k=0;
            //for (int i=0; i<newRepSize; i++){ //For each song in newRepertoire
            


                
		for (int j=0; j<newRepertoire.length; j++) { //For each syllable-dimension
                    if(newRepertoire[j]!=-1000){
                    //System.out.println("MutationVar = " + mutationVariance+" "+newRepertoire[k]);
                    newRepertoire[j]+=(float)(param.nextGaussian()*mutationVariance); //Change float k in the repertoire a little bit
                    //k++;
                    }
                }

                
        }
            
   

}