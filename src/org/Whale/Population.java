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
import java.util.LinkedList;

import org.apache.commons.math3.distribution.GammaDistribution;

public class Population extends org.ChaffinchABC.Population{

	Individual[] pop, emppop;
	Individual[][] neighbours;
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
        int[] subpopsize = new int[] {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
        int[] subpop;
        int[] tutors;
        int[] subpopstarts;
	Parameters param;
        int[] tutor;
        int ntutors=5;
        double problearn=0.7;
        
        
        //int[][] rec=new int[100][100];
	
	public Population(Individual[] pop, Parameters param, int[][] emplocs){
		this.pop=pop;
		px=param.nx;
		py=param.ny;
		this.param=param;
		this.ba=param.betaa;
		this.bb=param.betab/param.betaa;
		this.nthresh=param.neighthresh;
		this.emplocs=emplocs;
		
		ninds=pop.length;
                
                makeSubPops();
//		calculateSpatialArrangement();
//		preparePopulation();
//                makeDeadTerritories();
                
		makeNeighbours(nthresh);
//		setEmpiricalData();
	}
        
        public void makeSubPops(){
         int ninds = 0;
         
         for (int i=0; i<subpopsize.length; i++){  
            ninds += subpopsize[i];
    }
        
         
        subpop=new int[ninds];
        int k=0;
         for (int i=0; i<subpopsize.length; i++){ 
             for(int j=0; j<subpopsize[i]; j++){
                 subpop[k]=i;
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
		emppop=new Individual[emplocs.length];
		int n=0;
		for (int i=0; i<emplocs.length; i++) {
			//System.out.println(i+" "+emplocs[i][0]+ " "+ emplocs[i][1]);
			int p=lookUps[emplocs[i][0]][emplocs[i][1]];
			pop[p].repType=param.repSizes[i];
			pop[p].newRepSize=param.repSizes[i];
			pop[p].initiateRepertoire();
			pop[p].updateRepertoire();
			emppop[i]=pop[p];
                        if (emppop[i].isDeadTerritory){System.out.println("ALERT!!!!: dead emp: "+emplocs[i][0]+" "+emplocs[i][1]);}
			n+=emppop[i].repSize;
		}
		//System.out.println("SONGS: "+n+" "+emplocs.length+" "+param.repSizes.length);
		
	}
	
	public void makeNeighbours(double thresh) {
		neighbours=new Individual[ninds][];
		Individual[] buf=new Individual[ninds];
		double t=thresh*thresh;
		int x1, y1, x2, y2;
		double d;
		for (int i=0; i<ninds; i++) {
                    buf=new Individual[ninds];
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
			neighbours[i]=new Individual[count];
                        for (int j=0; j<count; j++){
                            neighbours[i][j]=buf[j];
                        }
                        //System.out.println(nthresh+" "+count);
			//System.arraycopy(buf, 0, neighbours[i], 0, count);	
		}
	}
	
	public org.Whale.Individual[] getTutors(int p) {
            int subpopp=subpop[p];
            int poptut=subpopp;

            double randomprob=param.nextDouble();
            
            if(randomprob>problearn){
             if (param.nextInt(2)==1){
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
            

            tutors=new int[ntutors];
            
        for(int i=0; i<ntutors;i++){
            boolean found=true;
            int a=0;
            while (found){
            a=param.nextInt(subpopsize[poptut])+subpopstarts[poptut];
            found=false;
            for (int j=0; j<i; j++){
		if (a==tutors[j]){
			found=true;
		}
             }
            }
            tutors[i]=a;
        }
            
            
/*
                int a=param.nextInt(subpopsize[poptut])+subpopstarts[poptut];
                for(int j=0;j<i;j++){
                    if(tutors[j]==a){
                    a=param.nextInt(subpopsize[poptut])+subpopstarts[poptut];
                    j=-1;
                    }
                }
                tutors[i]=a;
               }
 */           
                                   
            return tutors;
                    
                    
        }           

            
            
            /*
            boolean fail=true;
		
		int tx=0;
		int ty=0;
		int ox=rlu1[p];
		int oy=rlu2[p];
		
		while (fail) {
			fail=false;
			double x=param.nextDouble();
			int y=Arrays.binarySearch(prob, x);
			if (y<0) {y=(y+1)*-1;}
		
			int sx=dx[y];
			int sy=dy[y];
			
			
			if (param.nextLong()<0) {
				tx=ox-sx;
				if (tx<0) {
					tx=ox+sx;
					if (tx>=px) {
						fail=true;
					}
				}
			}
			else {
				tx=ox+sx;
				if (tx>=px) {
					tx=ox-sx;
					if (tx<0) {
						fail=true;
					}
				}
			}
			if (param.nextLong()<0) {
				ty=oy-sy;
				if (ty<0) {
					ty=oy+sy;
					if (ty>=py) {
						fail=true;
					}
				}
			}
			else {
				ty=oy+sy;
				if (ty>=py) {
					ty=oy-sy;
					if (ty<0) {
						fail=true;
					}
				}
			}
                        if (!fail){
                            if (pop[lookUps[tx][ty]].isDeadTerritory){fail=true;}
                        }
		}
		int q=lookUps[tx][ty];
                
                /*
                int rx=ox-tx+50;
                if (rx<0){rx=0;}
                if (rx>99){rx=99;}
                int ry=oy-ty+50;
                if (ry<0){ry=0;}
                if (ry>99){ry=99;}
                
                rec[rx][ry]++;
                
		//System.out.println("SELCHECK: "+p+" "+q+" "+ox+" "+tx+" "+oy+" "+ty);
		
                //System.out.println(getDistance(p,q)+" "+ba+" "+bb);
                
		Individual[] x=neighbours[q];
                
		int m=x.length;
		
                Individual[] y=new Individual[m];
                for (int j=0; j<m; j++) {
			int a=param.nextInt(m-j)+j;
			Individual b=x[a];
			x[a]=x[j];
			x[j]=b;
                        //System.out.println(p+" "+j+" "+b.territory);
                        //System.out.println(p+" "+a+" "+x[a].territory);
		}
		
		return x;
	*/
	
	
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
