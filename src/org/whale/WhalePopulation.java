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
import java.util.LinkedList;

import org.apache.commons.math3.distribution.GammaDistribution;

public class WhalePopulation {

    
	WhaleIndividual[] pop, emppop, tutors;
	WhaleIndividual[][] neighbours;
	int ninds=0;
	int nemp=0;
	double[] prob;
	int[] dx, dy, rlu1, rlu2;
	double ba, bb=0;

	//double nthresh=1.8;
        double nthresh=4.1;
	int px, py=0;
        int[][] lookUps;
        int[][] emplocs;
        int[] subpopsize = new int[] {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
        int[] minpops=new int[12];//
        int[] currentpopsizes=new int[12];//
        double growthrate=0.05;
        int minAge=10;
        int[] subpop;
        //int[] tutors;
        int[] subpopstarts;
	WhaleParameters param;
        int[] tutor;
        int ntutors=5;
        double problearn1=0.7; //probability of learning outside population in feeding grounds
        double problearn2=0.0; //probability of learning outside population in breeding grounds
        public boolean breeding; //season (for learning outside or within pop)
        int counter;
        int learnoutside=1; //# potential learning instances outside population
        int sampleperpop=100;

        
        //int[][] rec=new int[100][100];
	
	public WhalePopulation(WhaleIndividual[] pop, WhaleParameters param, int[][] emplocs){
            	//System.out.println("WhalePopulation");

		this.pop=pop;
                this.emppop=emppop;
		this.param=param;
		this.subpopsize=param.popsizes;
                this.minpops=param.minpops;
                this.growthrate=param.popgrowth;
                this.problearn1=param.problearn1;
                this.problearn2=param.problearn2;
                this.ntutors=param.ntutors;
                this.minAge=param.memorylength;
		//this.nthresh=param.neighthresh;
		//this.emplocs=emplocs;
		int popsize=0;
                
		ninds=pop.length;
                
                makeSubPops();
                initiateAges();
//		calculateSpatialArrangement();
//		preparePopulation();
//                makeDeadTerritories();
                
		//makeNeighbours(nthresh);
//		setEmpiricalData();
	}
        
        public void learnSongs(){
            if (!breeding){
                updatePopSize();
            }
            for (int i=0; i<subpopsize.length; i++){        //For every subpop
                for (int j=0; j<currentpopsizes[i]; j++){   //For every individual in the subpop
                    //System.out.println(j);
                    pop[j+subpopstarts[i]].learnSongs(); //Each individual in the currentpopsize goes through learnSongs();
                }
            }
            for (int i=0; i<subpopsize.length; i++){
                for (int j=0; j<currentpopsizes[i]; j++){
                    pop[j+subpopstarts[i]].updateRepertoire();//Each individual in the current popsize does updateRepertoire();
                }
            }
        }
        
        
        public void updatePopSize(){
            for (int i=0; i<subpopsize.length; i++){
                currentpopsizes[i]+=(int)Math.round(currentpopsizes[i]*growthrate*(1-(currentpopsizes[i]/(subpopsize[i]+0.0))));
                //System.out.println(i+" "+currentpopsizes[i]);
            }   
        }
        
        public void initiateAges(){
            for (int i=0; i<subpopsize.length; i++){
                for (int j=0; j<currentpopsizes[i]; j++){
                    pop[j+subpopstarts[i]].age=minAge;
                }
            }
        }
        
        public void makeSubPops(){
         
            currentpopsizes=new int[subpopsize.length];
            
            int ninds = 0; //reset number of inds
            
            for (int i=0; i<subpopsize.length; i++){  //for each subpop
                ninds += subpopsize[i]; //add the number of inds in that subpop to the total ninds
                currentpopsizes[i]=minpops[i]; //currentpopsizes for that subpop=minimalpopsize for that subpop
                //System.out.println(currentpopsizes[i]+" "+subpopsize[i]);
            }
        
         
        subpop=new int[ninds]; //subpop is new array with for each individual the subpop it belongs to
        int k=0;
         for (int i=0; i<subpopsize.length; i++){ //for every subpop
             for(int j=0; j<subpopsize[i]; j++){ //go through that population and specify the subpop
                 subpop[k]=i; //subpopulation in array subpop is specified
                 pop[k].subpop=i; //in WhalePopulation array the subpop is specified
                 k++;
             }
        }
        
        subpopstarts= new int[subpopsize.length];
        int start=0;
        subpopstarts[0]=start; 
        for(int i=1; i<subpopstarts.length; i++){
            int n=subpopsize[i-1];
            start += n;        
            subpopstarts[i]=start;      
        }
}
        
        public void makeEmpPop(){
        int n2=subpopsize.length*sampleperpop; //total number of IDs 10 populations*100
        emppop=new WhaleIndividual[n2]; //New individual array with length n2
        int k=0;
        int a=0; 
        for (int i=0; i<subpopsize.length; i++){ //for every subpop
            for (int j=0; j<sampleperpop; j++){ //for every of the 100 samples of that subpop
                emppop[k]=pop[a+j];
                k++;
        }
       a+=subpopsize[i];
   }
}
        




        
        public void makeDeadTerritories(){
            int[][][]r=param.removes;
            int count2=0;
            for (int i=0; i<r.length; i++){
               //System.out.println(px+" "+py+" "+r[i][0][0]+" "+r[i][0][1]+" "+r[i][1][0]+" "+r[i][1][1]);
                
                for (int j=r[i][0][0]; j<r[i][0][1]; j++){
                    for (int k=r[i][1][0]; k<r[i][1][1]; k++){
                        pop[lookUps[j][k]].isDeadTerritory=true;
                        count2++;
                    }
                }
            }
            /*
            int count=0;
            for (int i=0; i<pop.length; i++){
                if(!pop[i].isDeadTerritory){count++;}
            }
            
            System.out.println("DEAD REPORT: "+pop.length+" "+count2+" "+count);
            */
        }
	
	public void calculateSpatialArrangement() {
            
            
            
		
		GammaDistribution gamma=new GammaDistribution(ba, bb);
		
		prob=new double[ninds];
		dx=new int[ninds];
		dy=new int[ninds];
		
		//If a bird is hatched in territory X, what is the probability that it ends in up territory Y?
		int k=0;
		double sum=0;
		for (int i=0; i<px; i++) {
			for (int j=0; j<py; j++) {
				
				double d=Math.sqrt(i*i+j*j);
				//prob[k]=beta.density(d)/(d*d);
				if (d<1) {d=1;}
				prob[k]=gamma.density(d)/d;		//divided by d because that will be proportional to the circumference at d.
				//System.out.println("SPCHCK: "+i+" "+j+" "+prob[k]);
				//if (k<1000) {System.out.println("CHKD: "+i+" "+j+" "+prob[k]);}
				// In other words, the probability of a bird moving distance d is from beta; prob that it goes to particular territory at distance d is
				//proportional to beta divided by circumference...
				sum+=prob[k];
				dx[k]=i;
				dy[k]=j;
				k++;
			}
		}		
		for (int i=0; i<prob.length; i++) {
			prob[i]/=sum;		// normalizes for quantization effects of territories.
		}
		for (int i=1; i<prob.length; i++) {
			prob[i]+=prob[i-1];		// turns into cumulative.
		}
	
	}
	
	public void preparePopulation() {
		
		lookUps=new int[px][py];
		rlu1=new int[ninds];
		rlu2=new int[ninds];
		int k=0;
		for (int i=0; i<px; i++) {
			for (int j=0; j<py; j++) {
				lookUps[i][j]=k;
				rlu1[k]=i;
				rlu2[k]=j;
				k++;
			}
		}
	}
        
	
	public void setEmpiricalData() {
		nemp=emplocs.length;
		emppop=new WhaleIndividual[emplocs.length];
		int n=0;
		for (int i=0; i<emplocs.length; i++) {
			//System.out.println(i+" "+emplocs[i][0]+ " "+ emplocs[i][1]);
			int p=lookUps[emplocs[i][0]][emplocs[i][1]];
			pop[p].repType=param.repSizes[i];
			pop[p].newRepSize=param.repSizes[i];
			pop[p].initiateRepertoire(i);//specify for which inidvidual (?)
			pop[p].updateRepertoire();
			emppop[i]=pop[p];
                        if (emppop[i].isDeadTerritory){System.out.println("ALERT!!!!: dead emp: "+emplocs[i][0]+" "+emplocs[i][1]);}
			n+=emppop[i].repSize;
		}
		//System.out.println("SONGS: "+n+" "+emplocs.length+" "+param.repSizes.length);
		
	}
	
	public void makeNeighbours(double thresh) {
		neighbours=new WhaleIndividual[ninds][];
		WhaleIndividual[] buf=new WhaleIndividual[ninds];
		double t=thresh*thresh;
		int x1, y1, x2, y2;
		double d;
		for (int i=0; i<ninds; i++) {
                    buf=new WhaleIndividual[ninds];
			x1=rlu1[i];
			y1=rlu2[i];
			int count=0;
			for (int j=0; j<ninds; j++) {
				x2=rlu1[j];
				y2=rlu2[j];
				d=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
				if ((!pop[j].isDeadTerritory)&&(d<=t)) {
					buf[count]=pop[j];
					count++;
				}
			}
			neighbours[i]=new WhaleIndividual[count];
                        for (int j=0; j<count; j++){
                            neighbours[i][j]=buf[j];
                        }
                        //System.out.println(nthresh+" "+count);
			//System.arraycopy(buf, 0, neighbours[i], 0, count);	
		}
	}
        
        public void setSeason(int counter){
            breeding=true;
            if(counter<learnoutside){
                breeding=false;
            }

        }
	
	public WhaleIndividual [] getTutors(int p) {

            int subpopp=subpop[p]; //which subpopulation is ID p from
            int poptut=subpopp;    //tutor population by default is own subpopulation 
            double problearn;
            problearn=problearn1;
           
            if(breeding){
                problearn=problearn2;
            }
            //System.out.println("Breeding = " + breeding + " & problearn= " + problearn);

            
            tutors= new WhaleIndividual [ntutors];
            int[] tut;
            tut= new int [ntutors];
            
        for(int i=0; i<ntutors;i++){
            
            if(param.nextDouble()<problearn){  //only learning from neighbouring population when...
            //System.out.println("learnFromOtherPop");
             if (param.nextInt(2)==1){ //learn from population on the left or learn from the right?
                   poptut=subpopp+1;  
                   if (subpopp==subpopsize.length-1){ 
                      poptut=0;
                   }
               }
               else{
                   if(subpopp==0){
                       poptut=subpopsize.length-1;
                   }
                   else {
                       poptut=subpopp-1;
                   }
               }
            }
          
            
            boolean found=true;
            int a=0;
            while (found){
            //a=param.nextInt(subpopsize[poptut])+subpopstarts[poptut];
            a=param.nextInt(currentpopsizes[poptut])+subpopstarts[poptut];
            found=false;
            if (pop[a].age<minAge){
                found=true;
            }
            else{
                for (int j=0; j<i; j++){
                    if (a==tut[j]){
			found=true;
                    }
                }
                }
            }
            tut[i]=a;
        }
        
        for (int i=0; i<ntutors; i++){
            tutors[i]=pop[tut[i]];
            //if (subpop[tut[i]]!=subpop[p]){
            //    System.out.println("Learn from other pop: "+breeding+" "+tut[i]+" "+p);
            //}
            //3
            //System.out.println(p+" "+tut[i]+" "+subpop[p]+" "+subpop[tut[i]]+" "+breeding+" "+problearn1+" "+problearn2);
 
        }
                             // System.out.println("return Tutors");
                        return tutors;
                    
                    
        }  
            
	
        public int getEmpPop(int a){
                    int subpopsample=-1;
                    int maxsubemppop=0;
                    for(int j=0; j<subpop.length; j++){ //go through the subpopulation to find which subpopulation holds the individual with index a
                        maxsubemppop+=sampleperpop;
                        if(a<maxsubemppop){
                            subpopsample=j;
                            break;
                        }
                    }
            //System.out.println("individual = " + a + " & subpop = " + subpopsample);
            return subpopsample; //returns the number of the subpopulation

        }
	
	public int[] calculateIDs() {
		int n=0;
		for (int i=0; i<pop.length; i++) {
			if (pop[i].repType>=0) {
				n+=pop[i].repSize;
			}
		}
		
		int[] out=new int[n];
		
		n=0;
		for (int i=0; i<pop.length; i++) {
			if (pop[i].repType>=0) {
				for (int j=0; j<pop[i].repSize; j++) {
					out[n]=i;
					n++;
				}
			}
		}
		return out;
	}
	
	public int[] calculateEmpIDs() {
		int n=0;
		for (int i=0; i<emppop.length; i++) {
			n+=emppop[i].repSize;
		}
		
		int[] out=new int[n];
		
		n=0;
		for (int i=0; i<emppop.length; i++) {
			//System.out.println(i+" "+emppop[i].repSize);
			for (int j=0; j<emppop[i].repSize; j++) {
				out[n]=i;
				
				n++;
			}
		}
		
		//System.out.println(n+" "+emppop.length);
		
		return out;
	}
	
	public double[][] calculateDissimilarityMatrix(){
		int n=0;
		for (int i=0; i<pop.length; i++) {
			if (pop[i].repType>=0) {
				n+=pop[i].repSize;
			}
		}
		
		double[][] out=new double[n][n];
		int x=0;
		int y=0;
		double s=0;
		for (int i=0; i<pop.length; i++) {
			if (pop[i].repType>=0) {
				for (int j=0; j<pop[i].repSize; j++) {
					y=0;
					for (int a=0; a<pop.length; a++) {
						if (pop[a].repType>=0) {
							for (int b=0; b<pop[a].repSize; b++) {
								out[x][y]=pop[i].compareSongs(pop[i].repertoire, pop[a].repertoire, j, b);
								s+=out[x][y];
								y++;
							}
						}
					}
					x++;
				}
			}
		}
		
		//System.out.println("Mean Dissimilarity: "+(s/(n*n*1.0)));
		
		return out;
		
	}
        
        public void reportEmpPop(){
            for (int i=0; i<emppop.length; i++){
                for (int j=0; j<emppop[i].repSize; j++){
                    for (int k=0; k<emppop[i].ns; k++){
                        System.out.print(emppop[i].repertoire[j*emppop[i].ns+k]+" ");
                    }
                    System.out.println();
                }
            }
        }
	
	public double[][] calculateEmpDissimilarityMatrix(double adj){
		int n=0;
		for (int i=0; i<nemp; i++) {
			n+=emppop[i].repSize;
		}
		//System.out.println("SONGS: "+n);
		
		double[][] out=new double[n][n];
		int x=0;
		int y=0;
		double s=0;
		double ps=0;
		for (int i=0; i<emppop.length; i++) {
			for (int j=0; j<emppop[i].repSize; j++) {
				y=0;
				for (int a=0; a<emppop.length; a++) {
					for (int b=0; b<emppop[a].repSize; b++) {
						out[x][y]=emppop[i].compareSongs(emppop[i].repertoire, emppop[a].repertoire, j, b)+adj;
						s+=out[x][y];
						if (out[x][y]<param.typeThresh) {ps++;}
						
						y++;
					}
				}
				x++;
			}
		}
		
		//System.out.println("Mean Dissimilarity: "+(s/(n*n*1.0)));
		//System.out.println("Prop shared: "+(ps/(n*n*1.0))+" "+param.typeThresh);
	
		return out;
		
	}
	
	public double[][] calculateEmpSyllDissimilarityMatrix(double adj){
		int n=0;
		for (int i=0; i<nemp; i++) {
			n+=emppop[i].repSize;
		}
		//System.out.println("SONGS: "+n);
		n*=param.sylsPerSong;
		double[][] out=new double[n][n];
		int x=0;
		int y=0;
		double s=0;
		double ps=0;
		double t=0;
		//double d1, d2;
		
		for (int i=0; i<emppop.length; i++) {
			for (int j=0; j<emppop[i].repSize; j++) {
				for (int k=0; k<param.sylsPerSong; k++) {
					//d1=emppop[i].repertoire[j][k];
					y=0;
					for (int a=0; a<emppop.length; a++) {
						for (int b=0; b<emppop[a].repSize; b++) {
							for (int c=0; c<param.sylsPerSong; c++) {
								if (c==k) {
									//d2=emppop[a].repertoire[b][c]-d1;
									//out[x][y]=Math.sqrt(d2*d2);
                                                                        out[x][y]=emppop[i].compareSyllables(emppop[i].repertoire, j, k, emppop[a].repertoire, b, c)+adj;
									s+=out[x][y];
									t++;
									if (out[x][y]<param.typeThresh) {ps++;}
								}
								else {
									out[x][y]=1;
								}
								y++;
							}
						}
					}
					x++;
				}
			}
		}
		
		//System.out.println("Mean Dissimilarity: "+(s/(t)));
		//System.out.println("Prop shared: "+(ps/(n*n*1.0))+" "+param.typeThresh);
	
		return out;
		
	}
	
	public double getDistance(int a, int b) {
		int ax=rlu1[a];
		int bx=rlu1[b];
		int ay=rlu2[a];
		int by=rlu2[b];
		
		double d=Math.sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by));
		return d;
	}
	
	public double getEmpDistance(int a, int b) {
		int ax=emplocs[a][0];
		int bx=emplocs[b][0];
		int ay=emplocs[a][1];
		int by=emplocs[b][1];
		
		double d=Math.sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by));
		return d;
	}
	/*
        public void printrec(){
            for (int i=0; i<100; i++){
                for (int j=0; j<100; j++){
                    System.out.print(rec[i][j]+" ");
                }
                System.out.println();
            }
        }
	*/
	
}