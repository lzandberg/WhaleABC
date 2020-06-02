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
import java.util.Random;


/**
 *
 * @author lieszandberg
 */






//This is the main class for objects representing individuals - including all song learning behaviour

public class WhaleIndividual  {
    
        int age=0;
	int modelType=0;
	int subpop=0;
	float[] repertoire, newRepertoire, tempRep;

	int repSize, newRepSize;
        public int index;
        public float[] songmemory, songlog;
	
	//as with repertoire, there is a current repertoire size, and a value for the newly constructed one
	//after learning. These get merged at the end of each year.
	
	int territory;
	//the id for this individual within the population
	
	//the population object (to which this individual belongs)
	WhalePopulation population;
	WhaleParameters param;
        Whale whale;
	//Below are a list of characteristics that are taken from the defaults object.
	//see defaults for more info
	//Defaults defaults;
	
        boolean isDeadTerritory=false;
	double mortalityRate=0.4;

	double matchThresh;
        double matchThreshTheme; 
	double mutationVariance=0.1;
	double recombinationRate=0.001;
        double missingtheme=0.01; //??

	double confBias=1.1;
        
        int numSylls;
        int numdims;
        int memorylength, loglength; //memorylength = 50 songs, loglength = 100 songs
        int ns=0;
        int memorysize=0; //numsongs*numsyllables*numdimensions
        int logsize=0;
        int iter=0;
        int iterlog=0;
        int maxtheme=10; //parameters?
        double probadd=0.1; //parameters? probability of adding a completely new theme

        double novbias = 1;
         
         
        
	
	//int maxRep=20;	//this is used to limit array sizes
	int maxRep=1;
        int mr;	//maxRep-1 ; to save some simple repetition later
	
	int repType=-1;	
	//for the individuals that correspond to the empirical sample, this parameter is set to their repertoire size
	//and that is what their repertoire size is set to in the simulation too.
	//for the individuals that were not included in the empirical sample, this parameter is set to -1, and then 
	//the program sets their repertoire size sampling from the repertoire size distribution
	
	boolean isDead=true;
	//trigger to let program know that an individual needs to be replaced that year.

	int gens=0;
	
	float[] songBuffer;
        double[] tutsongsim;
	int[] songFreq;
        int[] songLength;
        float[] songsim, newsong; 
	double[] cumFreq, freq;
	double[] powerLookUp;
	int repSizes[];
        int maxSongLength;
	
	int songtypeCount=0;
        int songthemeCount=0;
        int songtypefreq=0;
        float min=0;
        float max=10;
	
	long t1,t2,t3,t4,t5, t6;
	


        
	public WhaleIndividual(int territory, WhaleParameters param, int repType){
		//System.out.println("WhaleIndividual");	
		this.territory=territory;
		//this.random=random;
		
		this.param=param;
		this.repType=repType;
			
		this.modelType=param.modelType;
		this.mortalityRate=param.mortalityRate;
		this.numSylls=param.sylsPerSong;
                this.numdims=param.numdims;
		this.mutationVariance=param.mutationVar;
		this.matchThresh=param.typeThresh*param.typeThresh*param.sylsPerSong;
                this.memorylength=param.memorylength;
                this.loglength=param.loglength;
                this.novbias=param.novbias*0.5; //0.5 because in compareSongsPow, we calculate the squared dissimilarity.

	
                ns=numSylls*numdims;
                memorysize=numSylls*numdims*memorylength;
                logsize=numSylls*numdims*loglength;
                
                //System.out.println(memorysize+" "+memorylength);
                //System.out.println(mutationVariance+" "+recombinationRate);
                
		int maxRep=param.maxRep;
		
		repertoire=new float[maxRep*ns];
		newRepertoire=new float[maxRep*ns];	
		tempRep=new float[ns];
		songBuffer=param.songBuffer;
		songFreq=param.songFreq;
                tutsongsim=param.tutsongsim;
		cumFreq=param.cumFreq;
                freq=param.freq;
		powerLookUp=param.powerLookUp;
		repSizes=param.repSizes;
		mr=repSizes.length;
		
		setRepertoireSize();
		
	}
        
        public void initiate(){
            initiateRepertoire(population.subpop[territory]);
            updateRepertoire();
            initiateMemory();
            initiateSongLength();
        }
	
	public void setPopulation(WhalePopulation population) {
		this.population=population;
	}
        
	
	
	//decides whether individual dies that year or not.
	public void mortality(){
            if(!isDeadTerritory){
		if (param.nextDouble()<mortalityRate){
			isDead=true;
			gens++;
		}
		else{
			isDead=false;
		}
            }
            else{
                isDead=false;
            }
	}
	
	
	
	//at beginning of simulation run, repertoires are initiated with randomly selected values.
	public void initiateRepertoire(int subpop){
            int k=0;
            //System.out.println("newRepSize= " + newRepSize );
            for (int j=0; j<ns; j++) { 
                //newRepertoire[k]=param.nextFloat()+subpop;
                newRepertoire[k]=subpop; //newRepertoire[k]=ind.subpop number: each song is ns float numbers each float is the same: the subpopulation number?
                k++;
            }
	}
        
        
        
        	public void initiateMemory(){
                    //int k=0;
                    songmemory = new float[memorysize];  //make float[] for memory, with length memorysize (#songs * #dims * #maxthemes)
                    for(int i=0; i<memorysize; i++){ //fill up memory with that one song
                        songmemory[i] = newRepertoire[i%ns];
                    }
                    songlog=new float[logsize];
                    for(int i=0; i<logsize; i++){ // also fill up songlog with that one song
                        songlog[i] = newRepertoire[i%ns];
                    }

	}
        
        public void initiateSongLength(){
            songLength = new int[memorylength];
            for(int i=0; i<memorylength; i++){
                songLength[i]=newRepertoire.length;
            }
        }

       /* 	public void initiateSongLength(){
                    //int k=0;
                    songmemory = new int[memorysize]; 
                    for(int i=0; i<memorysize; i++){ //fill up memory with that one song
                        songmemory[i] = newRepertoire[i%ns];
                    }


	}
                */
	
	//at the end of each year, each replaced individual's song values are updated to their newly calculated ones
	public void updateRepertoire(){
            //System.out.println("updateRepertoire");
		if (isDead){
                    
                    System.arraycopy(newRepertoire, 0, repertoire, 0, newRepSize*ns);
                    
                    //System.out.println("rep updated");
		/*	for (int i=0; i<newRepSize; i++){
				//repertoire[i]=newRepertoire[i];
				System.arraycopy(newRepertoire[i], 0, repertoire[i], 0, ns);
				//System.out.println(i+" "+repertoire[i][0]);
			}*/
                    repSize=newRepSize;
		}
	}
        
        
        
        public void learnSongs(){
                   
                
                buildRepertoire3();
                //if (songtypeCount>0){
                mutate(); 
                //}
                
              
                
                System.arraycopy(newRepertoire, 0, songmemory, iter*ns, ns); 
                songLength[iter]=newRepertoire.length;       
                iter++;

                if (iter>=memorylength){
                    iter=0;
                    
                }
                
                System.arraycopy(newRepertoire, 0, songlog, iterlog*ns, ns);
                iterlog++;
                if (iterlog>=loglength){
                    iterlog=0;
                }
                age++;
               
                
	}
	
	//simple function to determine repertoire size for new male after replacement.
	//this depends on whether you are from the original empirically sampled set (repType>=0) or not.
	
	public void setRepertoireSize(){
            //System.out.println("setRepertoireSize");
		if (repType<0){
			newRepSize=repSizes[param.nextInt(mr)];
		}
		else{
                    newRepSize=repType;
		}
		
		if (newRepSize>maxRep) {
			System.out.println(repType+" "+newRepSize);
		}
		
	}

		
	public boolean matchSongs(float[] x, float[] y, int a, int b) {
                //System.out.println("a = " + a + " & b = " + b + " & ns = " + ns);
		double d=0;
		int aa=a*ns;
                int bb=b*ns;
		for (int i=0; i<ns; i++) {
                    if(aa!=-1000){
                        if(bb==-1000){
                            d+=missingtheme;
                        }
                        else{
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
                        }
                    }
                    aa++;
                    bb++;
                    if (d>matchThresh){
                        return false;
                    }
		}
		
		//if (d>matchThresh) {	
		//	return false;
		//}
		//System.out.println(matchThresh+" "+d+" ");
		return true;
		
		
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
                    if (d>(matchThreshTheme/numSylls)){
                        return false;
                    }
		}
		
		//if (d>matchThresh) {	
		//	return false;
		//}
		//System.out.println(matchThresh+" "+d+" ");
		return true;
		
		
	}
        


        

             
        public double compareSyllables(float[] x, int s1, int a, float[] y, int s2, int b){
            
            int aa=a*numdims+s1*ns;
            int bb=b*numdims+s2*ns;
            double d=0;
            double buf=0;
            for (int i=0; i<numdims; i++){
                buf=x[aa+i]-y[bb+i];
                d+=buf*buf;
            }
            
            return(Math.sqrt(d));
        }
	
	public double compareSongs(float[] x, float[] y, int a, int b) {
		double d=0;
                
                int aa=a*ns;
                int bb=b*ns;
                
		for (int i=0; i<ns; i++) {
                    if(aa!=-1000){
                        if(bb==-1000){
                            d+=missingtheme;
                        }
                        else{
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
                        }
                    }
                    aa++;
                    bb++;
			//System.out.println(d+" "+x[i]+" "+y[i]);
		}
                
		return (Math.sqrt(d/(0.0+numSylls)));	
	}
	
	public double[] merge(double[] a, double[]b) {
		double[] c=new double[ns];
		for (int i=0; i<ns; i++) {
			c[i]=0.5*(a[i]+b[i]);
		}
		return c;
	}
        
        public double compareSongsPow(float[] x, float[] y, int a, int b) {
		double d=0;
                
                int aa=a*ns;
                int bb=b*ns;
                
		for (int i=0; i<ns; i++) {
                    if(aa!=-1000){
                        if(bb==-1000){
                            d+=missingtheme;
                        }
                        else{
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
                        }
                    }
                    aa++;
                    bb++;
			//System.out.println(d+" "+x[i]+" "+y[i]);
		}
                //System.out.println("dissim = " + Math.sqrt(d/(0.0+numSylls)) + "novelty bias = " + Math.pow(d/(0.0+numSylls),novbias));	
 
		return (pow(d/(0.0+numSylls),novbias));	
	}
	
	
	public boolean matchSongs(float[] a, float[] b, double x) {
		double d=0;
		double e=0;
		double c=1/x;
		for (int i=0; i<ns; i++) {
			e=a[i]-b[i]*c;
			d+=e*e;			
		}
		
		if (d>matchThresh) {	
			return false;
		}
		//System.out.println(matchThresh+" "+d+" ");
		return true;
	}
	
	public void addSongs(double[] a, double[] b) {
		for (int i=0; i<ns; i++) {
			a[i]+=b[i];
		}
	}
	
        
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
                    
	
	public float[] getSongMemory(){
            return songmemory;
        }
        

        public float[] getSongLog(){
            return songlog;
        }
        
   
       public void constructMemory4(WhaleIndividual[] tutors){
       	
            //System.out.println("constructMemory3"); 
            songtypeCount=0;
            //boolean found;
            double mindist;
            
            
            int r; //repertoire size of each given tutor
            int tl=tutors.length; // number of tutors
            //tutsongsim = new double[tl] ;
                       //list with for each tutor the minimal song similarity comparing its (currently the tutors only) song 
                       //with all songs in the individuals song memory
                       // **Still need to figure out a way how to do this when the song rep size tutor>1**
            for (int i=0; i<tl; i++) { //for each tutor               
                r=tutors[i].repSize; //determine repertoire size (r) for the tutor        
                double x;
           
                for (int j = 0; j < r; j++) { //for each song in a tutors repertoire 
                    mindist=Double.MAX_VALUE;
                    
                    for (int k=0; k<memorylength; k++) {  
                         x=compareSongsPow(tutors[i].repertoire, songmemory, j, k); 
                        if (x<mindist){
                            mindist=x;
                        } 
                    }
                    if (mindist>matchThresh){
                        tutsongsim[i]=mindist;   //add it to the tutor song similarity array
                        //Change ns
                        System.arraycopy(tutors[i].repertoire, 0, songBuffer, songtypeCount*ns, ns); //add it to the songbuffer
                  //System.out.println("copied tutorsong to songBuffer");
                        songtypeCount++;
                    }
                }      
            }
            if (songtypeCount>0){
                makeCumDistr();
            }
            
            // New method to count songtypes

            
            
       }
        
	public void constructMemory3(WhaleIndividual[] tutors) { // constructs memory of tutors songs of which one song will eventually be picked
            //System.out.println("constructMemory3"); 
            songtypeCount=0;
            //boolean found;
            double mindist;
            
            
            int r; //repertoire size of each given tutor
            int tl=tutors.length; // number of tutors
            
            //tutsongsim = new double[tl] ;
                       //list with for each tutor the minimal song similarity comparing its (currently the tutors only) song 
                       //with all songs in the individuals song memory
                       // **Still need to figure out a way how to do this when the song rep size tutor>1**
            for (int i=0; i<tl; i++) { //for each tutor               
                r=tutors[i].repSize; //determine repertoire size (r) for the tutor        
                double x;
           
                for (int j = 0; j < r; j++) { //for each song in a tutors repertoire 
                    mindist=Double.MAX_VALUE;
                    
                    for (int k=0; k<memorylength; k++) {  
                         x=compareSongsPow(tutors[i].repertoire, songmemory, j, k); 
                        if (x<mindist){
                            mindist=x;
                        } 
                    }
                    if (mindist>matchThresh){
                        tutsongsim[i]=mindist;   //add it to the tutor song similarity array
                        //Change ns
                        System.arraycopy(tutors[i].repertoire, 0, songBuffer, songtypeCount*ns, ns); //add it to the songbuffer
                  //System.out.println("copied tutorsong to songBuffer");
                        songtypeCount++;
                    }
                }      
            }
            if (songtypeCount>0){
                makeCumDistr();
            }
            
            // New method to count songtypes

            
        }
        
        
        public void makeFLookUps(){
            for (int i=0; i<songtypeCount; i++){
                freq[i]=powerLookUp[songFreq[i]];
            }
            
            cumFreq[0]=freq[0];
            for (int i=1; i<songtypeCount; i++) {
		cumFreq[i]=freq[i]+cumFreq[i-1];
            }
        }
        
        public void makeCumDistr(){
           //System.out.println("makeCumDistr");
            cumFreq[0]=tutsongsim[0];
            
            for (int i=1; i<songtypeCount; i++) {
		cumFreq[i]=tutsongsim[i]+cumFreq[i-1];
            }
        }
        
	public void pickSongs(){
		int i,j;
		double v, x;
		int t=0;
		int c2=songtypeCount-1;
                if (newRepSize<songtypeCount){
                    newRepSize=songtypeCount;
                }
		for (i=0; i<newRepSize; i++){
			v=param.nextDouble()*cumFreq[c2];
			for (j=0; j<songtypeCount; j++){
				if (v<cumFreq[j]){
					//x=1/(songFreq[j]+0.0);
                                        System.arraycopy(songBuffer, j*ns, newRepertoire, i*ns, ns);
					//for (int k=0; k<ns; k++) {
					//	newRepertoire[i][k]=songBuffer[j][k];
                                            
					//}
					//newRepertoire[i]=songBuffer[j];
					t=j;
					j=songtypeCount;
				}
			}
			for (j=t; j<songtypeCount; j++){
				cumFreq[j]-=freq[t];
			}
			
		}
	}
	
	public void pickSongs2(){
            //System.out.println("PickSongs2");
		int i,j;
		double v;
		int t=0;
		int c2=songtypeCount-1;
                
                if (newRepSize>songtypeCount){
                    //System.out.println(songtypeCount+" "+newRepSize+" "+param.mutationVar+" "+param.confBias);
                    newRepSize=songtypeCount;
                   
                }
                //if (songtypeCount<3){
                  //  System.out.println(songtypeCount+" "+newRepSize+" "+param.mutationVar+" "+param.confBias);
                //}
		for (i=0; i<newRepSize; i++){
			v=param.nextDouble()*cumFreq[c2];
			for (j=0; j<songtypeCount; j++){
				if (v<cumFreq[j]){

                                        // Check that oldsong is replaced by different newSong
                                        float[] newSong = new float[ns];
                                        System.arraycopy(songBuffer, j*ns, newSong, 0, ns);
                                        //System.out.println("newsongOrig = " + Arrays.toString(newSong));
                                        System.arraycopy(songBuffer, j*ns, newRepertoire, i*ns, ns); //copy to repertoire
                                        //System.arraycopy(songBuffer, j*ns, songmemory, iter, ns); //copy to songmemory
                                        //System.out.println("newmemory = " + Arrays.toString(songmemory));
                                        //how to make it overwrite the oldest memory?
                                        
					//newRepertoire[i]=songBuffer[j];
                                        //System.out.println(i+" "+j);
					t=j;
					j=songtypeCount;
				}
			}
			for (j=t; j<songtypeCount; j++){
				cumFreq[j]-=freq[t];
			}
			
		}
                dropTheme();
                addTheme();
	}
            
        
   
        public void buildRepertoire3() {
                //System.out.println("buildRepertoire3");
		//t2=System.nanoTime();
		WhaleIndividual[] tutors=population.getTutors(territory); //in this case territory is individual ID
		//t3=System.nanoTime();
		constructMemory3(tutors);
                if (songtypeCount>0){
                //checkMemory();
		//t4=System.nanoTime();
                    pickSongs2();
                }
		//t5=System.nanoTime();
	}
        
        
        
        public void checkMemory(){
            
            
            for (int i=0; i<songtypeCount; i++){
                for (int j=0; j<i; j++){
                    if(matchSongs(songBuffer, songBuffer, i, j)){
                        System.out.println("ALERT!");
                    }
                    
                }
                
                
            }
            
            
        }
        
        public float[] dropTheme(){
            int m=0;
      
            for (int i=0; i < numSylls; i++) { //for each theme in a newRepertoire  
                            if(newRepertoire[i*numdims]!=-1000){ 
                            m++;
                            }
                        }
            
            
            for (int i=0; i < numSylls; i++) { //for each theme in a newRepertoire  
                        int k=0;

                        for (int j=0; j<(songBuffer.length/numdims); j++) {  //go through each theme in the songbuffer (not all tutors repertoires...)
                            
                            if(matchThemes(newRepertoire, songBuffer, i, j)){ 
                            k++;
                            }
                        }
                    double dropprob=(1/(1+ Math.pow(0.5,(k-1+(5-(0.5*m)))))); //probability of learning a theme - 
                    // dependent on the frequency of each theme (k) and the number of themes in the tutor song //??
                    double chance = param.nextDouble();
                    if (chance>dropprob){ //
                        //System.out.println("dropTheme: k= " + k + " dropprob= " + dropprob + " chance= " + chance)  ;
                        for(int l=0; l<numdims;l++){
                            newRepertoire[(i*numdims)+l]=-1000; //drop that theme (all dimensions revert to -1000)
                        }
                    }
                    
                }
            
                return newRepertoire;
        }
            
        
        
        public float[] dropTheme2(float[] x){
            float oldsong[]=x;
            
            int loc=0;

            loc=param.nextInt((oldsong.length/numdims)); 
            loc=loc*numdims;    //random location in the oldsong array (including 0 and songlength+1)

            for (int i=0; i<oldsong.length; i++) { //remove one theme from the oldsong[] at location loc
                if (i<loc){ 
                newsong[i] = oldsong[i]; 
                }
                else if (i==loc) {
                    i+=numdims;
                }
                else{
                    newsong[i-numdims]=oldsong[i];
                } 

            }
            return newsong;
        }
        
        
                
        public float[] addTheme(){
               
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
                return newRepertoire;
        }
           
        
        
        public float[] addTheme2(float[] x){
            float oldsong[]=x;
            newsong= new float[oldsong.length+numdims];
            float[] newtheme = new float[numdims];
            int loc=0;
            
            loc=param.nextInt((oldsong.length/numdims)+1); 
            loc=loc*numdims;    //random location in the oldsong array (including 0 and songlength+1)
            
            for(int i=0; i<newtheme.length; i++){ //make random new theme
                newtheme[i]=param.nextFloat(); //check if we want this 0-10

            }

            for (int i=0; i<oldsong.length; i++) { //insert the newtheme[] into the oldsong[] at location loc
                if (i<loc){ 
                newsong[i] = oldsong[i]; 
                }
                else if (i==loc) {
                    for(int j=0; j<newtheme.length; j++){
                        newsong[i+j]=newtheme[j];
                    }
                    newsong[i+newtheme.length]=oldsong[i];
                            }
                else{
                    newsong[i+newtheme.length]=oldsong[i];
                } 

            }
            return newsong;
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
            
            
        
        
        
        public void recombine(){
            int r=0;
            int a=0;
            int maxtries=5;
            int j=0;
            for (int i=0; i<newRepSize; i++){ //for each song in the repertoire
		if (param.nextDouble()<recombinationRate) {	 //if probability			
                    boolean found=true;
                    int c=0;
                    while(found){
                        int q=param.nextInt(songtypeCount);         //
                        r=(param.nextInt(numSylls-1)+1)*numdims;    //get a random integer for a syllable and the corresponding dimension
                        int b=r+ns*q;                               //in a certain song in the songBuffer
                        a=0;
                        for (j=r; j<ns; j++) {                      //for sylldim j in a 
                            tempRep[a]=songBuffer[b]+(float)(param.nextGaussian()*mutationVariance);
                            a++;
                            b++;
                        }
                        found=false;
                        //System.out.println("mutate");
                        for (j=0; j<newRepSize; j++){
                            if (matchSongs(tempRep, newRepertoire, 0, j)){
                                found=true;
                            }
                        }
                        c++;
                        if ((found)&&(c==maxtries)){found=false;}
                    }
                    a=r+ns*i;
                    for (j=0; j<ns; j++){
                        newRepertoire[a]=tempRep[j];
                        a++;
                    }
                }
            }
            //System.out.println("newmemory = " + Arrays.toString(songmemory));
            System.arraycopy(newRepertoire, 0, songmemory, iter*ns, ns); //copy mutated song to songmemory
            //System.out.println("iter = " + iter);
	}
        
        
        //Source: https://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c/
        public static double pow(final double a, final double b) {
            final long tmp = Double.doubleToLongBits(a);
            final long tmp2 = (long)(b * (tmp - 4606921280493453312L)) + 4606921280493453312L;
            return Double.longBitsToDouble(tmp2);
        }

	

}