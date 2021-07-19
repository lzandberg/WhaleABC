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


//This is the main class for objects representing individuals - including all song learning behaviour

public class SimpleIndividual  extends WhaleIndividual{
        /*
        int minThemes=2;
        int maxThemes=10; //parameters?
        int age=0;
	int subpop=0;
	float[] repertoire, newRepertoire;

        public int index;
        public float[] songmemory,songlog;
	
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

        double matchThresh=0;
	
        
        int numSylls;
        int memorylength, loglength; //memorylength = 50 songs, loglength = 100 songs
        
        int ns=0;
        int memorysize=0; //numsongs*numsyllables*numdimensions
        int logsize=0;
        int iter=0;
        int iterlog=0;
        
        double probadd=0.01; //parameters? probability of adding a completely new theme         
	
	int maxRep=1;
        int mr;	//maxRep-1 ; to save some simple repetition later
		
	float[] songBuffer;

        int maxSongLength;
	
	int songtypeCount=0;
        int songthemeCount=0;
        int songtypefreq=0;
        float min=0;
        float max=10;
	
        double[] powlookUp;
*/
    
        double mutationRate=0.01;
        int memoryActual;
        
        
	public SimpleIndividual(int territory, WhaleParameters param){
            super();
		//System.out.println("WhaleIndividual");	
		this.territory=territory;
		
		this.param=param;
			
		this.numSylls=param.sylsPerSong;
		this.mutationRate=param.mutationVar;
                this.memorylength=param.memorylength;
                this.loglength=param.loglength;

                this.minThemes=param.minThemes;
                this.maxThemes=param.maxThemes;
                this.matchThresh=param.typeThresh;
                ns=numSylls;
                memorysize=numSylls*memorylength;
                logsize=numSylls*loglength;
                int memoryActual=0;
                
                this.powerLookUp=param.powerLookUp;
                this.tutsongsim=param.tutsongsim;
                this.cumFreq=param.cumFreq;
		
		repertoire=new float[ns];
		newRepertoire=new float[ns];	
		songBuffer=param.songBuffer;
                
                this.probadd=param.probadd;
		this.dropLookUp=param.powlookUp;
	}
        
        public void initiate(){
            //System.out.println("Songs initiating!");
            initiateRepertoire(population.subpop[territory]);
            updateRepertoire();
            initiateMemory();
        }
	
	public void setPopulation(WhalePopulation population) {
		this.population=population;
	}

	//at beginning of simulation run, repertoires are initiated with randomly selected values.
	public void initiateRepertoire(int subpop){
            
            for (int j=0; j<7; j++) { 
                newRepertoire[j]=subpop*0.01f+j*0.001f; //newRepertoire[k]=ind.subpop number: each song is ns float numbers each float is the same: the subpopulation number?
                //System.out.print(newRepertoire[j]+" ");
            }
            //System.out.println();
            for (int j=7; j<ns; j++){
                newRepertoire[j]=-1000;
            }
	}
        
        
        public void initiateMemory(){
            songmemory = new float[memorysize];  //make float[] for memory, with length memorysize (#songs * #dims * #maxthemes)
            for(int i=0; i<memorysize; i++){ //fill up memory with that one song
               songmemory[i] = newRepertoire[i%ns];
            }
            //memoryActual=memorylength;
            songlog=new float[logsize];
            for(int i=0; i<logsize; i++){ // also fill up songlog with that one song
                songlog[i] = newRepertoire[i%ns];
            }
	}

	//at the end of each year, each replaced individual's song values are updated to their newly calculated ones
	public void updateRepertoire(){
            System.arraycopy(newRepertoire, 0, repertoire, 0, ns);
            //System.out.println("repertoire updated"+repertoire.length); 
	}
        
        public void learnSongs(boolean logcounter){
            
            buildRepertoire();
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
                int m=0; //number of themes in song A
                
                int aa=a*ns;
                int bb=b*ns;
                int a2, b2;
		for (int i=0; i<numSylls; i++) { //Go through the themes of song a
                    
                    a2=aa+i; //Start of theme[i] of song a
                    if(x[a2]!=-1000){ //if theme in song a is not -1000:

                        m++;
                        for (int j=0; j<numSylls; j++) { //Go through themes of song b and find the most similar theme
                            b2=bb+j; //Start of theme[j] song b
                            if (y[b2]!= -1000){
                                if (x[a2]==y[b2]){
                                    d++;
                                    j=numSylls;
                                }
                            }
                        }
                        if (m-d>matchThresh*numSylls){return false;}
                    }
		}
                if (m-d>matchThresh*m){return false;}
                
                
        //System.out.println(m+" "+d+" "+numSylls+" "+matchThresh);
        return true;
    }
    
    public int countUnsharedThemes(float[] x, float[] y, int a, int b) {
                int c=0;
                int aa=a*ns;
                int bb=b*ns;
                int a2, b2;
                boolean found;
		for (int i=0; i<numSylls; i++) { //Go through the themes of song a
                    
                    a2=aa+i; //Start of theme[i] of song a
                    if(x[a2]!=-1000){ //if theme in song a is not -1000:
                        found=false;
                        for (int j=0; j<numSylls; j++) { //Go through themes of song b and find the most similar theme
                            b2=bb+j; //Start of theme[j] song b
                            if (y[b2]!= -1000){
                                if (x[a2]==y[b2]){
                                    found=true;
                                    j=numSylls;
                                }
                            }
                        }
                        if (!found){c++;}
                    }
                }
        return c;
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
            a2=aa+i; //Start of theme[i] of song a
            if(x[a2]!=-1000){ //if theme in song a is not -1000:
                for (int j=0; j<numSylls; j++) { //Go through themes of song b and find the most similar theme
                    b2=bb+j; //Start of theme[j] song b
                    if (y[b2]!=1000){
                        if (x[a2]==y[b2]){
                            count[p]=1;
                            j=numSylls;
                        }
                    }
		}
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
            a2=aa+i; //Start of theme[i] of song a
            if(x[a2]!=-1000){
                count++;
            }
        }
        return count;
    }
        

    public boolean matchThemes(float[] x, float[] y, int a, int b) {
        double d=0;
	int aa=a;
        int bb=b;
        
        if (x[aa]==y[bb]){return true;}
        return false;
    }

    public float[] getSongMemory(){
        return songmemory;
    }
        

    public float[] getSongLog(){
        return songlog;
    }
        
   
    public void constructMemory(WhaleIndividual[] tutors) { // constructs memory of tutors songs of which one song will eventually be picked
        songtypeCount=0;
        int tl=tutors.length; // number of tutors
        for (int i=0; i<tl; i++) { //for each tutor 
            boolean found=false;
            for (int k=0; k<memoryActual; k++) {  
                if(matchSongs(tutors[i].repertoire, songmemory, 0, k)){
                    found=true;
                } 
            }    
            if (!found){
                System.arraycopy(tutors[i].repertoire, 0, songBuffer, songtypeCount*ns, ns); //add it to the songbuffer
                songtypeCount++;
            }
        }
        //if (songtypeCount>0){System.out.println(songtypeCount+" "+tutors.length);}
    }
    
      public void constructMemory2(WhaleIndividual[] tutors) { // constructs memory of tutors songs of which one song will eventually be picked
        songtypeCount=0;
        int tl=tutors.length; // number of tutors
        for (int i=0; i<tl; i++) { //for each tutor 
            int minunshared=100;
            for (int k=0; k<memoryActual; k++) {  
                int p=countUnsharedThemes(tutors[i].repertoire, songmemory, 0, k);
                if(p<minunshared){minunshared=p;}
            }
            if (minunshared>0){
                System.arraycopy(tutors[i].repertoire, 0, songBuffer, songtypeCount*ns, ns); //add it to the songbuffer
                tutsongsim[songtypeCount]=powerLookUp[minunshared];
                songtypeCount++;
            }
        }
        //if (songtypeCount>0){System.out.println(songtypeCount+" "+tutors.length);}
    }
    
    public void makeCumDistr(){
        cumFreq[0]=tutsongsim[0];
        for (int i=1; i<songtypeCount; i++) {
            cumFreq[i]=tutsongsim[i]+cumFreq[i-1];
        }
    }

    public void pickSongs2(){
        if (songtypeCount>0){
            int c2=songtypeCount-1;
            double v=param.nextDouble()*cumFreq[c2];
            for (int j=0; j<songtypeCount; j++){
                if (v<cumFreq[j]){
                    System.arraycopy(songBuffer, j*ns, newRepertoire, 0, ns); //copy to repertoire
                }
            }
        }
    }
    
    
        
    public void pickSongs(){
        int a=param.nextInt(songtypeCount);
        System.arraycopy(songBuffer, a*ns, newRepertoire, 0, ns);
    }

            
    public void buildRepertoire() {
        WhaleIndividual[] tutors=population.getTutors(territory); //in this case territory is individual ID
	//constructMemory2(tutors);
        //pickSongs2();
        //makeCumDistr();
        constructMemory(tutors);
        if (songtypeCount>0){
            learningEpochs++;
            pickSongs();
            dropTheme(tutors);
            addTheme();
        }
    }
        
    public void dropTheme(WhaleIndividual[] tutors){
        int m=0;
        for (int i=0; i < numSylls; i++) { //Counting the number of themes
            if(newRepertoire[i]!=-1000){ 
                m++;
            }
        }
        for (int i=0; i < numSylls; i++) { //for each theme in a newRepertoire  
            if(m>minThemes){
                if(newRepertoire[i]!=-1000){
                    int k=0;
                    int p=numSylls*songtypeCount;
                    /*
                    for (int j=0; j<p; j++) {  //go through each theme in the songbuffer (not all tutors repertoires...)
                        if(matchThemes(newRepertoire, songBuffer, i, j)){ 
                            k++;
                        }
                    }               //this compares song to the number of different songs you've memorized which isn't correct...
                    */
                    for (int j=0; j<tutors.length; j++) {  //go through each theme in the tutors' repertoires...
                        for (int a=0; a<numSylls; a++){
                            if(matchThemes(newRepertoire, tutors[j].repertoire, i, a)){ 
                                k++;
                            }
                        }
                    }
                    double dropprob=dropLookUp[k]; //probability of learning a theme (not dropping the theme)- 
                    if (param.nextDouble()>dropprob){ //
                        param.deletionCounter++;
                        newRepertoire[i]=-1000; //drop that theme (all dimensions revert to -1000)
                        m--;
                    }
                }
                    
            }
        }
    }
        
    public void addTheme(){
        for(int i=0; i<numSylls; i++){ //for each theme in new repertoire
            if(newRepertoire[i]==-1000){
                if(param.nextDouble()<probadd){
                    param.insertionCounter++;
                    newRepertoire[i]=param.nextFloat(); 
                    //System.out.println("Add "+newRepertoire[i]);
                }
            }
        }
    }
           
    public void mutate(){
        //
        for (int j=0; j<newRepertoire.length; j++) { //For each syllable-dimension
            if((newRepertoire[j]!=-1000)&&(param.nextDouble()<mutationRate)){
                param.mutationCounter++;
                newRepertoire[j]=param.nextFloat(); //Change float k in the repertoire a little bit
                //System.out.println("mutate " +mutationRate+" "+newRepertoire[j]);
            }
        }
    }
            
}
