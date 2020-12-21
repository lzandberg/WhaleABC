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

import java.util.Arrays;
import java.util.LinkedList;

import org.apache.commons.math3.distribution.GammaDistribution;

public class WhalePopulation {

    /*
        int[][] NHInteractionTable={{1,1,1,0,0,0},
            {1,1,1,1,0,0},
            {1,1,1,1,0,0},
            {0,1,1,1,0,0},
            {0,0,0,0,1,1},
            {0,0,0,0,1,1}
        };
        */
        
        int[][] NHInteractionTable={{1,2},
            {0,2,3},
            {0,1,3},
            {2,3},
            {5},
            {4}
        };
        
        /*int[][] SHInteractionTable={{1,1,0,0,0,0,0,0,0,0,0,1},
            {1,1,1,0,0,0,0,0,0,0,0,0},
            {0,1,1,1,0,0,0,0,0,0,0,0},
            {0,0,1,1,1,0,0,0,0,0,0,0},
            {0,0,0,1,1,1,0,0,0,0,0,0},
            {0,0,0,0,1,1,1,0,0,0,0,0},
            {0,0,0,0,0,1,1,1,0,0,0,0},
            {0,0,0,0,0,0,1,1,1,0,0,0},
            {0,0,0,0,0,0,0,1,1,1,0,0},
            {0,0,0,0,0,0,0,0,1,1,1,0},
            {0,0,0,0,0,0,0,0,0,1,1,1},
            {1,0,0,0,0,0,0,0,0,0,1,1}  
        };*/
        
         int[][] SHInteractionTable={{1,11},
            {0,2},
            {1,3},
            {2,4},
            {3,5},
            {4,6},
            {5,7},
            {6,8},
            {7,9},
            {8,10},
            {9,11},
            {10,0}  
        };
        
        int[][] interactionTable;
    
	WhaleIndividual[] pop, emppop, tutors;
	WhaleIndividual[][] neighbours;
	int ninds=0;
	int nemp=0;
	double[] prob;
	int[] dx, dy, rlu1, rlu2;
	double ba, bb=0;
        
        boolean agecounter=false;
        int agecounterval=9;

	//double nthresh=1.8;
        double nthresh=4.1;
	int px, py=0;
        int[][] lookUps;
        int[][] emplocs;
        int[] subpopsize = new int[] {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
        int[] minpops=new int[12];//
        int[] kpops=new int[12];
        int[] currentpopsizes=new int[12];//
        int[] tutorpopsizes=new int[12];
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
        double[] growthrates=new double[12];
        int hemisphere=0;
        //int[][] rec=new int[100][100];
	
	public WhalePopulation(WhaleIndividual[] pop, WhaleParameters param, int[][] emplocs){
            	//System.out.println("WhalePopulation");

		this.pop=pop;
                this.emppop=emppop;
		this.param=param;
		this.subpopsize=param.popsizes;
                this.minpops=param.minpops;
                this.kpops=param.kpops;
                this.growthrate=param.popgrowth;
                this.problearn1=param.problearn1;
                this.problearn2=param.problearn2;
                this.ntutors=param.ntutors;
                //this.minAge=param.memorylength;
                this.minAge=param.loglength*10;
                this.hemisphere=param.hemisphere;
                
                if (hemisphere==0){
                    interactionTable=SHInteractionTable;
                }
                else{
                    interactionTable=NHInteractionTable;
                }
                
		//this.nthresh=param.neighthresh;
		//this.emplocs=emplocs;
		int popsize=0;
                
		ninds=pop.length;
                calculateLogisticParameters();
                makeSubPops();
                initiateAges();
//		calculateSpatialArrangement();
//		preparePopulation();
//                makeDeadTerritories();
                
		//makeNeighbours(nthresh);
//		setEmpiricalData();
	}
        
        
        public void calculateLogisticParameters(){
            int n=kpops.length;
            double nyears=50;
            growthrates=new double[n];
            for (int i=0; i<n; i++){
                
                double c1=minpops[i]/(kpops[i]-minpops[i]+0.0);
                double c2=(kpops[i]-subpopsize[i])/(subpopsize[i]+0.0);
                
                growthrates[i]=-1*Math.log(c1*c2)/nyears;
            }
        }
        
        
        public void learnSongs(){
            if (!breeding){
                //System.out.println("updating pop size");
                updatePopSize();
            }
            for (int i=0; i<subpopsize.length; i++){        //For every subpop
                for (int j=0; j<currentpopsizes[i]; j++){   //For every individual in the subpop
                    //System.out.println("a "+i+" "+j+" "+subpopsize.length+" "+currentpopsizes[i]);
                    pop[j+subpopstarts[i]].learnSongs(agecounter); //Each individual in the currentpopsize goes through learnSongs();
                }
            }
            for (int i=0; i<subpopsize.length; i++){
                for (int j=0; j<currentpopsizes[i]; j++){
                    pop[j+subpopstarts[i]].updateRepertoire();//Each individual in the current popsize does updateRepertoire();
                }
            }
            //System.out.println("Done");
        }
        
        
        public void updatePopSize(){
            for (int i=0; i<subpopsize.length; i++){
                //currentpopsizes[i]+=(int)Math.round(currentpopsizes[i]*growthrate*(1-(currentpopsizes[i]/(subpopsize[i]+0.0))));
                //currentpopsizes[i]+=(int)Math.round(currentpopsizes[i]*growthrates[i]*(1-(currentpopsizes[i]/(subpopsize[i]+0.0))));
                currentpopsizes[i]+=(int)Math.ceil(currentpopsizes[i]*growthrates[i]*(1-(currentpopsizes[i]/(kpops[i]+0.0))));
                //System.out.println(i+" "+currentpopsizes[i]);
            }   
            for (int i=0; i<subpopsize.length; i++){
                int count=0;
                for (int j=0; j<currentpopsizes[i]; j++){
                    //System.out.println(pop[j+subpopstarts[i]].age);
                    if(pop[j+subpopstarts[i]].age>=minAge){count++;}
                }
                tutorpopsizes[i]=count;
                //System.out.println("s: "+i+" "+count+" "+currentpopsizes[i]);
            }
        }
        
        public void initiateAges(){
            for (int i=0; i<subpopsize.length; i++){
                for (int j=0; j<currentpopsizes[i]; j++){
                    pop[j+subpopstarts[i]].age=minAge;
                    pop[j+subpopstarts[i]].memoryActual=pop[j+subpopstarts[i]].memorylength;
                }
            }
        }
        
        public void makeSubPops(){
         
            currentpopsizes=new int[subpopsize.length];
            
            int ninds = 0; //reset number of inds
            
            for (int i=0; i<subpopsize.length; i++){  //for each subpop
                //ninds += subpopsize[i]; //add the number of inds in that subpop to the total ninds
                ninds += kpops[i]; //add the number of inds in that subpop to the total ninds
                currentpopsizes[i]=minpops[i]; //currentpopsizes for that subpop=minimalpopsize for that subpop
                tutorpopsizes[i]=minpops[i];
                //System.out.println(currentpopsizes[i]+" "+subpopsize[i]+" "+tutorpopsizes[i]);
            }
        
         
        subpop=new int[ninds]; //subpop is new array with for each individual the subpop it belongs to
        int k=0;
         for (int i=0; i<subpopsize.length; i++){ //for every subpop
            //for(int j=0; j<subpopsize[i]; j++){ //go through that population and specify the subpop
            for(int j=0; j<kpops[i]; j++){ //go through that population and specify the subpop
                 subpop[k]=i; //subpopulation in array subpop is specified
                 pop[k].subpop=i; //in WhalePopulation array the subpop is specified
                 k++;
             }
        }
        
        subpopstarts= new int[subpopsize.length];
        int start=0;
        subpopstarts[0]=start; 
        for(int i=1; i<subpopstarts.length; i++){
            //int n=subpopsize[i-1];
            int n=kpops[i-1];
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
        
	
        public void setSeason(int counter){
            breeding=true;
            if(counter<learnoutside){
                breeding=false;
            }
            
            agecounter=false;
            if (counter==agecounterval){
                agecounter=true;
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

            //System.out.println(ntutors+" "+subpopsize.length+" "+tutorpopsizes[poptut]+" "+currentpopsizes[poptut]);
            
            int[] tut;
            int ntut=ntutors;
            
                        
            if (ntutors>currentpopsizes[poptut]){
                ntut=currentpopsizes[poptut];
                //tutors= new WhaleIndividual [ntut];
                //for (int i=0; i<ntut; i++){
                //    tutors[i]=pop[subpopstarts[poptut]+i];
                //}
                //return tutors;
            }  
            //tutors= new WhaleIndividual [ntut];
            tut=new int[ntut];
            
        //boolean flagged=false;
        for(int i=0; i<ntut;i++){
            boolean found=true;
            int count=0;
            while (found){
                
                //System.out.println(i+" "+ntut);
                if(param.nextDouble()<problearn){  //only learning from neighbouring population when...
                    //System.out.println("learnFromOtherPop "+subpopp);
                    int oc=param.nextInt(interactionTable[subpopp].length);
                    poptut=interactionTable[subpopp][oc];
                    
                    /*
                    if (param.nextInt(2)==1){ //learn from population on the left or learn from the right?
                        
                        int npoptut=subpopp+1;  
                        if (subpopp==subpopsize.length-1){ 
                            npoptut=0;
                        }
                        //if (param.nextDouble()<currentpopsizes[npoptut]/(currentpopsizes[npoptut]+currentpopsizes[subpopp]+0.0)){
                            poptut=npoptut;
                        //}
                    }
                    else{
                        int npoptut=subpopp-1;
                        if(subpopp==0){
                            npoptut=subpopsize.length-1;
                        }
                        
                        //if (param.nextDouble()<currentpopsizes[npoptut]/(currentpopsizes[npoptut]+currentpopsizes[subpopp]+0.0)){
                            poptut=npoptut;
                        //}
                    }
                    */
                }
                //System.out.println(poptut+" "+subpopp+" "+currentpopsizes[poptut]+" "+subpopstarts[poptut]);
                //a=param.nextInt(subpopsize[poptut])+subpopstarts[poptut];
                int a=param.nextInt(currentpopsizes[poptut])+subpopstarts[poptut];
                found=false;
                //if (pop[a].age<minAge){
                if (pop[a].age==0){
                    found=true;
                    //System.out.println("here");
                }
                else{
                    for (int j=0; j<i; j++){
                        if (a==tut[j]){
                            found=true;
                            //System.out.println("or here");
                        }
                    }
                }
                tut[i]=a;
                count++;
                if (count>50){
                    //System.out.println(count+" "+ntut+" "+currentpopsizes[poptut]);
                    ntut--;
                    found=false;
                    //flagged=true;
                    //System.out.println(ntut);
                }
            }
            
        }
        tutors=new WhaleIndividual[ntut];
        for (int i=0; i<ntut; i++){
            //if (flagged){System.out.println(tut[i]);}
            tutors[i]=pop[tut[i]];
        }
        //if (flagged){ System.out.println("return Tutors"+tutors.length;);}
                        return tutors;
                    
                    
        }  
        
        public WhaleIndividual[] getWhales(int subpopp, int ntut){
            int poptut=subpopp;    //tutor population by default is own subpopulation 

            //System.out.println(subpopp+" "+ntut+" "+tutorpopsizes[poptut]+" "+currentpopsizes[poptut]);
            
            
            int[] tut;            
                        
            if (ntut>tutorpopsizes[poptut]){
                ntut=tutorpopsizes[poptut]; 
            }  
            
            
            tutors= new WhaleIndividual [ntut];
            tut=new int[ntut];
            
            
            for(int i=0; i<ntut;i++){
                //System.out.println(subpopp+" "+tutorpopsizes[poptut]+" "+i+" "+ntut);
                boolean found=true;
                int count=0;
                while (found){
                    count++;
                    //if (count>50){System.out.println(count+" "+ntut+" "+tutorpopsizes[poptut]);}
                    int a=param.nextInt(tutorpopsizes[poptut])+subpopstarts[poptut];
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
                    tut[i]=a;
                }
            
            }
        
            for (int i=0; i<ntut; i++){
                tutors[i]=pop[tut[i]];
            }
            return tutors;  
        }
        
        
        
        public WhaleIndividual [] getWhalesX(int subpopp, int nind) {

            int poptut=subpopp;    //tutor population by default is own subpopulation 
            
            int[] tut;
            int ntut=nind;
            if (nind>tutorpopsizes[poptut]){
                ntut=tutorpopsizes[poptut];
                tutors= new WhaleIndividual [ntut];
                for (int i=0; i<ntut; i++){
                    tutors[i]=pop[subpopstarts[poptut]+i];
                }
                return tutors;
            }  
            tut=new int[ntut];
            tutors= new WhaleIndividual [ntut];
            
            for(int i=0; i<ntut;i++){
                boolean found=true;
                int count=0;
                while (found){
                    count++;
                    if (count>50){System.out.println(count+" "+ntut+" "+tutorpopsizes[poptut]);}
                    int a=param.nextInt(currentpopsizes[poptut])+subpopstarts[poptut];
                    found=false;
                    if (pop[a].age<minAge){
                        found=true;
                        //System.out.println("here");
                    }
                    else{
                        for (int j=0; j<i; j++){
                            if (a==tut[j]){
                                found=true;
                            }
                        }
                    }
                    tut[i]=a;
                }
            
            }
        
            for (int i=0; i<ntut; i++){
                tutors[i]=pop[tut[i]];
            }
            
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
	
}