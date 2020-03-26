package org.ChaffinchABC;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.LinkedList;



//This is a class to read the empirical data from a csv file and process it. It's a bit ugly in places.

public class EmpData {
	
    double adj=0.03;
    double adjsyl=0.03;
	
	double[][] repSizes;	//a double[][] that includes the frequencies of different repertoire sizes (with index 0 equivalent to rep size 1) for the different populations
	String[] populations;	//array of names of populations
	
	int[][][] songFreqs;
	
	int currentPopulation=0;	//switch indicating population-specific simulations

	//Parameters param;
	
	int[][] locs;
	int[] syllab;
	String[][] ids;
	int[][] idssyll;
	double[][] empDiss, empSylDiss;
	int[] repSizeD;
	int maxRepSize;
        int[][] popSizes;
	String fileLocation="/data/home/btw774/";
	
        
        public EmpData(){}
        
	public EmpData(String fileLocation){
		this.fileLocation=fileLocation;
		//this.param=param;
		
		locs=parseGeogDataFile(fileLocation+"chaffgeog.csv");
		ids=parseLabelFile(fileLocation+"chafflabels.csv");
		idssyll=parseLabelFileSyl(ids, fileLocation+"chafflabelssyl.csv");
                popSizes=parsePopulationSize(fileLocation+"chaffpopsize.csv");
		//locs=parseGeogDataFile(fileLocation+"tengeog.csv");
		//ids=parseLabelFile(fileLocation+"tenlabels.csv");
		//locs=parseGeogDataFile(fileLocation+"bluegeog.csv");
		//ids=parseLabelFile(fileLocation+"bluelabels.csv");
		//idssyll=parseLabelFileSyl(ids, fileLocation+"bluelabelssyl.csv");
		syllab=makeSylLabel(idssyll);
		repSizeD=new int[ids.length];
		maxRepSize=0;
		//int repcount=0;
		for (int i=0; i<ids.length; i++) {
			repSizeD[i]=ids[i].length;
			//repcount+=repSizeD[i];
			if (repSizeD[i]>maxRepSize) {maxRepSize=repSizeD[i];}
			//System.out.println(param.repSizes[i]+" "+maxRepSize);
		}
		//System.out.println(repcount+" "+maxRepSize);
		maxRepSize++;
		empDiss=parseDissFile(ids, fileLocation+"chaffdiss.csv");
		empSylDiss=parseDissFileSyl(idssyll, fileLocation+"chaffdisssyl.csv");
		//empDiss=parseDissFile(ids, fileLocation+"tendiss.csv");
		//empDiss=parseDissFile(ids, fileLocation+"bluediss.csv");
		//empSylDiss=parseDissFileSyl(idssyll, fileLocation+"bluedisssyl.csv");
		
		//System.out.println(empDiss.length);

	}
	
	public int[][] parseGeogDataFile(String fileLocation){
		File file=new File(fileLocation);
		int[][] out=null;
		try{
			String cvsSplitBy = ",";
			BufferedReader reader=new BufferedReader(new FileReader(file));
			String line=null;
			LinkedList<int[]> locs=new LinkedList<int[]>();
			

			while((line=reader.readLine())!=null){
				
				String[] s=line.split(cvsSplitBy);
				
				int[] t=new int[s.length];
				for (int i=0; i<s.length; i++){
					t[i]=Integer.parseInt(s[i]);
				}
				locs.add(t);
			}
			
			reader.close();
			
			out=new int[locs.size()][2];
			for (int i=0; i<out.length; i++) {
				int[] a=locs.get(i);
				//System.out.println(i+" "+a[0]+" "+a[1]);
				out[i][0]=a[0];
				out[i][1]=a[1];
			}
			
			
		}
		
		catch(Exception e){
			e.printStackTrace();
		}
	
		return out;
	}
	
        public int[][] parsePopulationSize(String fileLocation){
            File file=new File(fileLocation);
            int[][] out=null;
            try{
		String cvsSplitBy = ",";
		BufferedReader reader=new BufferedReader(new FileReader(file));
		String line=null;
                LinkedList<int[]> d=new LinkedList<int[]>();
		while((line=reader.readLine())!=null){	
                    String[] s=line.split(cvsSplitBy);
                    if (s.length==2) {
			int[] x=new int[2];
                        x[0]=Integer.parseInt(s[0]);
                        x[1]=Integer.parseInt(s[1]);
                        d.add(x);
                        System.out.println("POPULATION SIZE: "+x[0]+" "+x[1]);
                    }
		}
                out=new int[d.size()][];
                for (int i=0; i<out.length; i++){
                    out[i]=d.get(i);
                }
                
		reader.close();			
            }	
            catch(Exception e){
                e.printStackTrace();
            }
	
            return out;
	}
        
        
	public String[][] parseLabelFile(String fileLocation){
		File file=new File(fileLocation);
		String[][] out=null;
		try{
			String cvsSplitBy = ",";
			BufferedReader reader=new BufferedReader(new FileReader(file));
			String line=null;
			LinkedList<String[]> locs=new LinkedList<String[]>();
			

			while((line=reader.readLine())!=null){	
				String[] s=line.split(cvsSplitBy);
				//System.out.println(s.length);
				if (s.length==2) {
					locs.add(s);	
				}
			}
			
			reader.close();
			
			int n=0;
			String last=" ";
			for (String[] s : locs) {
				if (!s[0].equals(last)) {
					n++;
				}
				last=s[0];
			}
			
			out=new String[n][];
			n=0;
			int m=0;
			last=locs.get(0)[0];
			for (String[] s : locs) {
				if (!s[0].equals(last)) {
					out[n]=new String[m];
					m=1;
					n++;
				}
				else {
					m++;
				}
				last=s[0];
			}
			out[n]=new String[m];
			
			n=0;
			m=-1;
			last=locs.get(0)[0];
			for (String[] s : locs) {
				if (!s[0].equals(last)) {
					m=0;
					n++;
				}
				else {
					m++;
				}
				
				out[n][m]=s[1];
				last=s[0];
			}		
		}
		
		catch(Exception e){
			e.printStackTrace();
		}
		
		for (int i=0; i<out.length; i++) {
			System.out.print(i+" ");
			for (int j=0; j<out[i].length; j++) {
				System.out.print(out[i][j]+" ");
			}
			System.out.println();
		}
		
                
            System.out.println("NumInds "+out.length);
		return out;
	}
	
	public int[][] parseLabelFileSyl(String[][] songNames, String fileLocation){
		File file=new File(fileLocation);
		int[][] out=new int[songNames.length][];
		for (int i=0; i<songNames.length; i++) {
			out[i]=new int[songNames[i].length];
		}
		try{
			String cvsSplitBy = ",";
			BufferedReader reader=new BufferedReader(new FileReader(file));
			String line=null;			
			
			while((line=reader.readLine())!=null){	
				String[] s=line.split(cvsSplitBy);
				
				if (s.length>=2) {
					boolean found=false;
					for (int i=0; i<songNames.length; i++) {
						for (int j=0; j<songNames[i].length; j++) {
							if (songNames[i][j].equals(s[0])) {
								out[i][j]++;	
                                                                found=true;
							}
						}
					}
                                        if (!found){System.out.println("MISSING SONG:"+s[0]+" "+s[1]);}
				}
			}
			
			reader.close();
			
		}
		
		catch(Exception e){
			e.printStackTrace();
		}
		/*
		for (int i=0; i<out.length; i++) {
			System.out.print(i+" ");
			for (int j=0; j<out[i].length; j++) {
				System.out.print(out[i][j]+" ");
			}
			System.out.println();
		}
		*/
		return out;
	}
	
	public int[] makeSylLabel(int[][] labels) {
		int n=0;
		for (int i=0; i<labels.length; i++) {
			for (int j=0; j<labels[i].length; j++) {
				n+=labels[i][j];
			}
		}
		
		int[] out=new int[n];
		n=0;
		int m=0;
		for (int i=0; i<labels.length; i++) {
			for (int j=0; j<labels[i].length; j++) {
				for (int k=0; k<labels[i][j]; k++) {
					out[n]=m;
					n++;
				}
				m++;
			}
		}
		return out;
	}
	
	public double[][] parseDissFile(String[][] labels, String fileLocation){
		
		int n=0;
		for (int i=0; i<labels.length; i++) {
			n+=labels[i].length;
		}
		String[] x=new String[n];
		n=0;
		for (int i=0; i<labels.length; i++) {
			for (int j=0; j<labels[i].length; j++) {
				x[n]=labels[i][j];
				n++;
			}
		}
		//System.out.println(n);
		double[][] out=new double[n][n];
		//for (int i=0; i<n; i++) {
		//	out[i]=new double[i+1];
		//}
		
		
		
		File file=new File(fileLocation);
		try{
			String cvsSplitBy = ",";
			BufferedReader reader=new BufferedReader(new FileReader(file));
			String line=null;			

			while((line=reader.readLine())!=null){	
				String[] s=line.split(cvsSplitBy);
				
				int a=0;
				int b=0;
				for (int i=0; i<n; i++) {
					if (x[i].equals(s[2])) {a=i;}
					if (x[i].equals(s[3])) {b=i;}
				}
				//System.out.println(a+" "+b+" "+s[0]+" "+s[1]+" "+s[2]+" "+s[3]+" "+s[4]);
				out[a][b]=Double.parseDouble(s[4])-adj;
                                if (out[a][b]<0){out[a][b]=0;}
				out[b][a]=out[a][b];
			}
			
			reader.close();
			
		}
		
		catch(Exception e){
			e.printStackTrace();
		}
		/*
		for (int i=0; i<out.length; i++) {
			System.out.print(i+" ");
			for (int j=0; j<out[i].length; j++) {
				System.out.print(out[i][j]+" ");
			}
			System.out.println();
		}
		*/
		return out;
	}
	
	public double[][] parseDissFileSyl(int[][] labels, String fileLocation){
		
		int n=0;
		for (int i=0; i<labels.length; i++) {
			for (int j=0; j<labels[i].length; j++) {
				n+=labels[i][j];
			}
		}
		double[][] out=new double[n][n];
		System.out.println(fileLocation+" "+n);
		int a=0;
		int b=0;
		
		File file=new File(fileLocation);
		try{
			String cvsSplitBy = ",";
			BufferedReader reader=new BufferedReader(new FileReader(file));
			String line=reader.readLine();			

			while((line=reader.readLine())!=null){	
				String[] s=line.split(cvsSplitBy);
				int a1=Integer.parseInt(s[4]);
                                int a2=Integer.parseInt(s[5]);
                                double b1=1;
                                if (a1==a2){b1=Double.parseDouble(s[6]);}
                                out[a][b]=b1-adjsyl;
                                out[b][a]=out[a][b];

				b++;
				if (b>a) {
					b=0;
					a++;
				}
				
			}
			
			reader.close();
			
		}
		
		catch(Exception e){
			e.printStackTrace();
		}
		/*
		for (int i=0; i<out.length; i++) {
			System.out.print(i+" ");
			for (int j=0; j<out[i].length; j++) {
				System.out.print(out[i][j]+" ");
			}
			System.out.println();
		}
		*/
		return out;
	}


}
