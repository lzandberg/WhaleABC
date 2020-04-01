/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.Whale;

import java.util.Arrays;

/**
 *
 * @author lieszandberg
 */






//This is the main class for objects representing individuals - including all song learning behaviour

public class WhaleIndividual extends org.ChaffinchABC.Individual {

	int modelType=0;
	
	float[] repertoire, newRepertoire, tempRep;

	int repSize, newRepSize;
	
	//as with repertoire, there is a current repertoire size, and a value for the newly constructed one
	//after learning. These get merged at the end of each year.
	
	int territory;
	//the id for this individual within the population
	
	//the population object (to which this individual belongs)
	WhalePopulation population;
	WhaleParameters param;
	//Below are a list of characteristics that are taken from the defaults object.
	//see defaults for more info
	//Defaults defaults;
	
        boolean isDeadTerritory=false;
        
	double mortalityRate=0.4;

	double matchThresh=0.1;
	double mutationVariance=0.1;
	double recombinationRate=0.001;

	double confBias=1.1;
        
        	int numSylls=6;
        int numdims=2;
        int memorylength=100;
        int ns=0;
        int memorysize=0;
        double nb;
        float songmemory[];            
        double[] tutsongsim; 
        
	
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
	int[] songFreq;
        float[] songsim; 
	double[] cumFreq, freq;
	double[] powerLookUp;
	int repSizes[];
	
	int songtypeCount=0;
	
	long t1,t2,t3,t4,t5, t6;
	


        
	public WhaleIndividual(int territory, WhaleParameters param, int repType){
			
		this.territory=territory;
		//this.random=random;
		
		this.param=param;
		
		this.repType=repType;
			
		this.modelType=param.modelType;
		this.mortalityRate=param.mortalityRate;
		this.numSylls=param.sylsPerSong;
                this.numdims=param.numdims;
		this.mutationVariance=param.mutationVar;
		this.recombinationRate=param.recomRate;
		this.matchThresh=param.typeThresh*param.typeThresh*param.sylsPerSong;
		this.confBias=param.confBias;
	
                ns=numSylls*numdims;
                memorysize=numSylls*numdims*memorylength;
                
                //System.out.println(mutationVariance+" "+recombinationRate);
                
		int maxRep=param.maxRep;
		
		repertoire=new float[maxRep*ns];
		newRepertoire=new float[maxRep*ns];	
		tempRep=new float[ns];
		songBuffer=param.songBuffer;
		songFreq=param.songFreq;
		cumFreq=param.cumFreq;
                freq=param.freq;
		powerLookUp=param.powerLookUp;
		repSizes=param.repSizes;
		mr=repSizes.length;
		
		setRepertoireSize();
		initiateRepertoire();
		updateRepertoire();
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
	
	
	/*
	//at beginning of simulation run, repertoires are initiated with randomly selected values.
	public void initiateRepertoire(){
            int k=0;
		for (int i=0; i<newRepSize; i++){
			for (int j=0; j<ns; j++) {
				//newRepertoire[i][j]=param.nextFloat()*0.2f+1000*j;
                                newRepertoire[k]=param.nextFloat()*20f;
                                k++;
				//newRepertoire[i][j]=i*0.1+nextInt(3)*0.01;
			}
		}
	}
        */
        
        
        	public void initiateMemory(){
            int k=0;
            songmemory = new float[memorysize]; 
            
		//make a new song? Or move this somewhere else to make population-wide startsong?
			for (int j=0; j<ns; j++) {
                                newRepertoire[j]=param.nextFloat()*20f;
			}
                        
                        for(int i=0; i<memorylength; i++){
                        songmemory[i] = newRepertoire[i%ns];
                        }

	}

	
	
	//at the end of each year, each replaced individual's song values are updated to their newly calculated ones
	public void updateRepertoire(){
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

	//in this function, individuals learn songs. First they are assigned a repertoire size. Then that
	//repertoire is filled according to the learningMethod.
		
	
	/*
	public void learnSongs(){
		if (isDead){
			//t1=System.nanoTime();
			setRepertoireSize();
			if (modelType==0) {
				buildRepertoireSimple();
			}
			//else if (modelType==1){
			//	buildRepertoireConformist();
			//}
			else if (modelType==2){
				buildRepertoireConformist2();
			}
			mutate();
			//t6=System.nanoTime();
			//System.out.println((t2-t1)+" "+(t3-t2)+" "+(t4-t3)+" "+(t5-t4)+" "+(t6-t5)+" "+(t6-t1));
			//double[] out= {(t2-t1), (t3-t2), (t4-t3), (t5-t4), (t6-t5)};
			//return out;
		}
		//return null;
	}
        */
        
        	public void learnSongs(){
		if (isDead){
			//t1=System.nanoTime();
			setRepertoireSize();
			if (modelType==0) {
				buildRepertoireSimple();
			}
			//else if (modelType==1){
			//	buildRepertoireConformist();
			//}
			else if (modelType==2){
				buildRepertoireConformist2();
			}
                        	else if (modelType==3){
				buildRepertoire3();
			}
			mutate();
		}
	}
	
	//simple function to determine repertoire size for new male after replacement.
	//this depends on whether you are from the original empirically sampled set (repType>=0) or not.
	
	public void setRepertoireSize(){
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
		double d=0;
		int aa=a*ns;
                int bb=b*ns;
		for (int i=0; i<ns; i++) {
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
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
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
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
                    d+=(x[aa]-y[bb])*(x[aa]-y[bb]);
                    aa++;
                    bb++;
			//System.out.println(d+" "+x[i]+" "+y[i]);
		}
		return (Math.pow(d/(0.0+numSylls),nb));	
	}
	
        

	
	public void buildRepertoireSimple(){
		
		WhaleIndividual[] tutors=population.getTutors(territory);
		//System.out.println(tutors.length);
		int n=tutors.length; //n=number of tutors
		int checker=maxRep*10; //?
		for (int i=0; i<newRepSize; i++){ //for each available spot in the repertoire (newRepertoireSize)
			boolean found=true; 
			int c2=0;
			while (found){
				found=false;
				int t=param.nextInt(n); //pick random tutor t
				int u=param.nextInt(tutors[t].repSize); //pick random song u from tutor t
				//float[] h=tutors[t].repertoire[u];
				
				for (int j=0; j<i; j++){				
					//if (h==newRepertoire[j]) { 
					//	found=true;	
					//	j=i;
					//}
					if (matchSongs(tutors[t].repertoire, newRepertoire, u, j)){
						found=true;
						//newRepertoire[j]=merge(h, newRepertoire[j]);
						j=i;
					}
				}
				if (!found){
					//newRepertoire[i]=h;
					System.arraycopy(tutors[t].repertoire, u*ns, newRepertoire, i*ns, ns);
				}
				c2++;
				if ((found)&&(c2==checker)){
					found=false;
					newRepSize=i;
				}
			}
		}
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
	
	/*
	public void constructMemory(Individual[] tutors) {
		
		songtypeCount=0;
		float[] x;
		boolean found;
		int i, j, k, a;
		double e,f;
		for (i=0; i<tutors.length; i++) {
			for (j=0; j<tutors[i].repSize; j++) {
				x=tutors[i].repertoire[j];
				found=false;
				for (k=0; k<songtypeCount; k++) {
					f=0;
					//g=songFreq[k];
					for (a=0; a<ns; a++) {
						e=x[a]-songBuffer[k][a];
						f+=e*e;
					}
					
					if (f<matchThresh) {
					//if (matchSongs(x, songBuffer[k], songFreq[k])) {
						for (a=0; a<ns; a++) {
							songBuffer[k][a]=(songBuffer[k][a]*songFreq[k]+x[a])/(songFreq[k]+1.0f);
						}
						//addSongs(songBuffer[k], x);
						songFreq[k]++;
						found=true;
						k=songtypeCount;
					}
				}
				if (!found) {
					System.arraycopy(x, 0, songBuffer[songtypeCount], 0, ns);
					songFreq[songtypeCount]=1;
					songtypeCount++;
				}
			}
		}	
		for (i=0; i<songtypeCount; i++) {
			cumFreq[i]=powerLookUp[songFreq[i]];
		}
		for (i=1; i<songtypeCount; i++) {
			cumFreq[i]+=cumFreq[i-1];
		}
		
		
	}
	*/
	
	
	public void constructMemory2(WhaleIndividual[] tutors) {
		
            songtypeCount=0;
            boolean found;
            //float[] sb;
            int r;
            int tl=tutors.length;
            for (int i=0; i<tl; i++) {
                r=tutors[i].repSize;
                for (int j = 0; j < r; j++) {
                    //sb = tutors[i].repertoire[j];
                    found=false;
                    for (int k=0; k<songtypeCount; k++) {
                        if (matchSongs(tutors[i].repertoire, songBuffer, j, k)) {
                            songFreq[k]++;
                            found=true;
                            k=songtypeCount;
                        }
                    }
                    if (!found) {
                        System.arraycopy(tutors[i].repertoire, j*ns, songBuffer, songtypeCount*ns, ns);
                        songFreq[songtypeCount]=1;
                        songtypeCount++;
                        
                        //songBuffer[songtypeCount] = tutor.repertoire[j];
                        
                    }
                }
            }
            
            makeFLookUps();
            

	}
        
   
	public void constructMemory3(WhaleIndividual[] tutors) {
		
            songtypeCount=0;
            boolean found;
            double mindist;
            
            
            int r; //repertoire size of each given tutor
            int tl=tutors.length; // number of tutors
                       //list with for each tutor the minimal song similarity comparing its (currently the tutors only) song 
                       //with all songs in the individuals song memory
                       // **Still need to figure out a way how to do this when the song rep size tutor>1**
            for (int i=0; i<tl; i++) { //for each tutor
                r=tutors[i].repSize; //determine repertoire size (r) for the tutor
                mindist=Double.MAX_VALUE;
                double x;

                
                double[] sc = new double[memorylength]; // make an array for the comparison of the tutors song with all of songMemory
                for (int j = 0; j < r; j++) { //for each song in a tutors repertoire 
                    
                    for (int k=0; k<memorylength; k++) {  
                                       x=compareSongsPow(tutors[i].repertoire, songmemory, j, k); 
                
                if (x<mindist){
                    mindist=x;
                }
            }
                  
                 // if (minValue>matchThresh){
                  tutsongsim[i]=mindist;   //add it to the tutor song similarity array
                  System.arraycopy(tutors[i].repertoire, i*ns, songBuffer, songtypeCount*ns, ns); //add it to the songbuffer
                  songtypeCount++;
                  
                 makeCumDistr();
          //  }

                }      
            }
            
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
					
					//for (int k=0; k<numSylls; k++) {
					//	newRepertoire[i][k]=songBuffer[j][k];
					//}
                                        
                                        System.arraycopy(songBuffer, j*ns, newRepertoire, i*ns, ns);
                                        
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
	}
	
        	

	
	//public void buildRepertoireConformist() {
	//	WhaleIndividual[] tutors=population.getTutors(territory);
	//	constructMemory(tutors);
	//	pickSongs();
	//}
	
        /*
	public void buildRepertoireConformist2() {
		//t2=System.nanoTime();
		WhaleIndividual[] tutors=population.getTutors(territory);
		//t3=System.nanoTime();
		constructMemory2(tutors);
                
                //checkMemory();
		//t4=System.nanoTime();
		pickSongs2();
		//t5=System.nanoTime();
	}
*/
        
        	public void buildRepertoire3() {
		//t2=System.nanoTime();
		WhaleIndividual[] tutors=population.getTutors(territory); //in this case territory is individual ID
		//t3=System.nanoTime();
		constructMemory3(tutors);
                
                //checkMemory();
		//t4=System.nanoTime();
		pickSongs2();
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
        
	
	public void mutate(){
            int i=0;
            int j=0; 
            float m=0;
            
            int k=0;
            for (i=0; i<newRepSize; i++){
			
		for (j=0; j<ns; j++) {
				
                    newRepertoire[k]+=(float)(param.nextGaussian()*mutationVariance);
                    k++;
		}
            }
            int r=0;
            int a=0;
            int maxtries=5;
            for (i=0; i<newRepSize; i++){
		if (param.nextDouble()<recombinationRate) {				
                    boolean found=true;
                    int c=0;
                    while(found){
                        int q=param.nextInt(songtypeCount);
                        r=(param.nextInt(numSylls-1)+1)*numdims;
                        int b=r+ns*q;
                        a=0;
                        for (j=r; j<ns; j++) {
                            tempRep[a]=songBuffer[b]+(float)(param.nextGaussian()*mutationVariance);
                            a++;
                            b++;
                        }
                        found=false;
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
	}
	
	

}