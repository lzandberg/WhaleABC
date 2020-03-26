/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SimulationRunner;

import org.ChaffinchABC.Priors;
import org.ChaffinchABC.ChaffinchABC;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedList;
import org.ABCRunner.DocumentSave;
import org.ChaffinchABC.EmpData;

/**
 *
 * @author rflachlan
 */
public class SimulationRunner {
        
    int n, id, round, numpops;
    long seed;
    String[] fileloc;
    
    public double[][] results;
    public double[][] stats;
    public double[][] params;
    public double[] score;
    
    double[][] pvals1, cov;
    double[] weights;
    
    Priors p;
    ChaffinchABC[] cemp;
    int[][][] locs;
    int[][] repSizes;
    int[][][] popSizes;
    String currentpop="";
    
    public SimulationRunner(ChaffinchABC[] cemp, int numpops, int n, int id, int round, long seed, String[] fileloc, double[] weights, double[][] pvals1, double[][]cov){
    
	this.n=n;
        this.numpops=numpops;
	this.id=id;
	this.round=round;
	this.results=new double[n][];
        this.seed=seed;
        this.fileloc=fileloc;
        p=new Priors(seed);
        
        this.cemp=cemp;
        //double[] x=p.sampleFromPriors();
        //cemp=new ChaffinchABC(fileloc, x, p.variables, System.currentTimeMillis());
	locs=new int[numpops][][];
        repSizes=new int[numpops][];
        popSizes=new int[numpops][][];
        for (int i=0; i<numpops; i++){
            locs[i]=cemp[i].locs;
            repSizes[i]=cemp[i].repSizes;
            popSizes[i]=cemp[i].popSizes;
        }
        
        if (weights!=null){
            this.weights=new double[weights.length];
            System.arraycopy(weights, 0, this.weights, 0, weights.length);
        }
        
        if (pvals1!=null){
            this.pvals1=new double[pvals1.length][];
            for (int i=0; i<pvals1.length; i++){
                this.pvals1[i]=new double[pvals1[1].length];
                System.arraycopy(pvals1[i], 0, this.pvals1[i], 0, pvals1[i].length);
            }
        }
        
        if (cov!=null){
            this.cov=new double[cov.length][];
            for (int i=0; i<cov.length; i++){
                this.cov[i]=new double[cov[1].length];
                System.arraycopy(cov[i], 0, this.cov[i], 0, cov[i].length);
            }
        }
        
    }
    
    public void readRound(){
        System.out.println("Reading round...");
        File file=new File(fileloc+"/cvmat.csv");
        cov=null;
	try{
            String cvsSplitBy = ",";
            BufferedReader reader=new BufferedReader(new FileReader(file));
            String line=null;
			
            LinkedList<double[]> locs=new LinkedList<double[]>();
            while((line=reader.readLine())!=null){
				
                String[] s=line.split(cvsSplitBy);
				
		double[] t=new double[s.length-1];
		for (int i=0; i<s.length-1; i++){
                    t[i]=Double.parseDouble(s[i]);
		}
		locs.add(t);
            }
			
            reader.close();
			
            cov=new double[locs.size()][];
            for (int i=0; i<cov.length; i++) {
                cov[i]=locs.get(i);	
            }
        }
		
	catch(Exception e){
            e.printStackTrace();
	}
        
        file=new File(fileloc+"/pvals.csv");
        pvals1=null;
        weights=null;
	try{
            String cvsSplitBy = ",";
            BufferedReader reader=new BufferedReader(new FileReader(file));
            String line=null;
			
            LinkedList<double[]> locs=new LinkedList<double[]>();
            while((line=reader.readLine())!=null){
				
                String[] s=line.split(cvsSplitBy);
				
		double[] t=new double[s.length-1];
		for (int i=0; i<s.length-1; i++){
                    t[i]=Double.parseDouble(s[i]);
		}
		locs.add(t);
            }
			
            reader.close();
			
            pvals1=new double[locs.size()][];
            weights=new double[locs.size()];
            for (int i=0; i<weights.length; i++) {
                double[] x=locs.get(i);
                pvals1[i]=new double[x.length-1];
                System.arraycopy(x, 0, pvals1[i], 0, x.length-1);
                weights[i]=x[x.length-1];
            }
        }
		
	catch(Exception e){
            e.printStackTrace();
	}
        
        
       // for (int i=0; i<cov.length; i++){
        //    System.out.println("COV: "+i+" "+cov[i][0]+" "+cov[i].length);
        //}
        //for (int i=0; i<weights.length; i++){
        //    System.out.println("WEI: "+i+" "+weights[i]);
        //}
        
    }
		
    public  void runSimulation(){
                    
	stats=new double[n][];
	params=new double[n][];
	score=new double[n];
	ChaffinchABC cabc;
	
			
        //if (round>0){
          //  readRound();
        //}
        
	for (int i=0; i<n; i++) {
            //System.out.println("Start: "+id+" "+round+" "+i);
				
            double[] x;
            if (round==0) {
		x=p.sampleFromPriors();		
            }
            else {
		double[] y=pvals1[pickParams()];
		x=p.drawFromProposal(cov, y);
            }
				
            params[i]=x;
            int xx=cemp[0].stats.length;
            stats[i]=new double[xx*numpops];
            score[i]=0;
            for (int j=0; j<numpops; j++){
                currentpop=fileloc[j];
                cabc=new ChaffinchABC(fileloc[j], p.nextLong(), locs[j], repSizes[j], popSizes[j], x, p.variables);		
                System.arraycopy(cabc.stats, 0, stats[i], j*xx, xx);
                score[i]+=calculateDifference(cabc.stats, cemp[j].stats);
            }		
            score[i]/=numpops+0.0;
            System.out.println("End: "+id+" "+round+" "+i+" "+score[i]+" "+x[0]+" "+x[1]+" "+x[2]+" "+x[3]+" "+x[4]);
            cabc=null;
            System.gc();
        }
       // writeRound();
     
    }
                
    public int pickParams() {
        int n=weights.length;
        double x=p.nextDouble()*weights[n-1];
        int loc=0;
        for (int i=0; i<n; i++){
            if (weights[i]>x){
		loc=i;
		i=n;
            }
        }	
        return loc;
		
    }
                
    public double calculateDifference(double[] x1, double[] x2) {
		
        double[] compsx=calculatePLSComponents(x1);
        double[] compsy=calculatePLSComponents(x2);
        for (int i=0; i<x1.length; i++){
            if (Double.isNaN(x1[i])){System.out.println("NaN error: xs: "+i);}
            if (Double.isNaN(x2[i])){System.out.println("NaN error: ys: "+i);}
            if (Double.isInfinite(x1[i])){System.out.println("Inf error: xs: "+i+" "+currentpop);}
            if (Double.isInfinite(x2[i])){System.out.println("Inf error: ys: "+i+" "+currentpop);}
        }
		
        double d[]=new double[compsx.length];
        double r=0;
        for (int i=0; i<compsx.length; i++){
			//d[i]=Math.abs(compsx[i]-compsy[i]);
            d[i]=compsx[i]-compsy[i];
			//System.out.println(i+" "+compsx[i]+" "+compsy[i]);
            if (Double.isNaN(d[i])){
		System.out.println("NaN error: "+i+" "+currentpop);
		d[i]=1000000000;
            }
            r+=d[i]*d[i];
        }
		//System.out.println(r);
        return Math.sqrt(r);	
    }
	
    public double[] calculatePLSComponents(double[]x) {
        
        
        double[] meansEur={-5.01123446549048,
0.620488261164743,
-2.25868867109894,
-4.3816361484787,
-4.5932997785268,
-4.10344067272516,
-3170.21252339761,
2.03905044281808,
0.180481813430034,
0.056756324584284,
-4.62349164008799,
-2.33239023280552,
-1.67143651376686,
1.27603578070916,
-2.70298492815375,
-4.19024254542585,
-4.44128548194377,
-4.79520589484783,
-1221.16135747739,
2.09773941767954,
0.10162726139512,
0.063559315072075,
-2.31841609951652,
-2.32654902050727,
0.366139672279454,
0.396921158635796,
0.279012362172349,
0.175466128432252,
0.109095400873605,
0.021137602947223,
0.071915163372294,
0.161055638171774,
0.197686116273492,
0.21718806938804,
0.30014098765264,
0.110153678309971,
-4.76123046896868,
-1.94140866459996,
-1.73159592337588,
-1.01366814391866,
-4.40214742359431,
-3.37331014229546,
-20635.5160277545,
0.619223732256057,
0.140277231982254,
0.240499035761728,
0.056806197932046,
0.855831114036807,
0.690979011722273
    };
			
	double[] meansten={-5.12400003386634,
0.812020351210565,
-1.42262513694081,
-4.00452951649113,
-5.31686952029794,
-4.14553452257471,
-2461.46324759867,
2.11000660053601,
0.136360622467409,
0.074458110055308,
-4.20279003627876,
-1.65006771396764,
-2.32869477956419,
0.966749713357919,
-2.0213148865739,
-3.56154093510974,
-3.97934624216939,
-4.39003787998439,
-1271.53884212326,
2.1561409195868,
0.126929460425612,
0.100614741724562,
-2.14940690616875,
-1.61614153627425,
0.160488967496936,
0.443689377599039,
0.30402873944658,
0.157083636016491,
0.075300211947095,
0.004690661718906,
0.060725459777789,
0.199846330332282,
0.242018597542398,
0.240264158192006,
0.248973013365408,
0.057492370613989,
-5.03183200472383,
-1.06525267862521,
-0.922008117164706,
-1.7394545226785,
-4.63105754913971,
-3.53008613685035,
-16573.1188774749,
0.666132403183601,
0.145307662294364,
0.188559934522023,
0.025137574203404,
0.649203454325181,
0.804648160582998,
};
        
        double[] means={-3.648851222745,
0.56207909481315,
-1.84646697405124,
-4.1259101185633,
-3.03112891148099,
-2.81362748929108,
-619.849553066165,
1.98028440494438,
0.156001179978834,
0.038759479245486,
-4.47195609306871,
-3.38034783819188,
-0.51990091771468,
0.79751768641915,
-2.2321526962907,
-3.30160453014194,
-3.14718929514628,
-3.34060497692709,
-198.833900962382,
1.99625704004,
0.095676760926729,
0.026470936980676,
-2.28952128000962,
-3.50425210292109,
0.623793259166247,
0.36807055733703,
0.281629660157898,
0.181872609803064,
0.106543035446961,
0.024722846932467,
0.065835043593322,
0.198978996175324,
0.226796760831086,
0.202944701882568,
0.242153937288592,
0.101254636080494,
-4.25385435808852,
-0.657566417081893,
-1.40521858645804,
-2.96746141702557,
-3.17112073549762,
-3.1106430233401,
-3021.57252922586,
0.722839296073216,
0.10317832501665,
0.173982378910143,
0.124669237647957,
0.662851155188058,
-1.19244627576924,
    };
        
        double[] sdEur={1.45022326103438,
2.81682872010525,
1.77194566165274,
2.19483075201863,
3.11553347315334,
1.38120510151514,
813.69212067217,
0.478840966928311,
0.17942250257961,
0.140122378862232,
2.90143838730985,
2.08845230824882,
0.443267666657538,
1.76899101397925,
1.55955277102211,
1.90262659045755,
2.26625297944275,
0.699994077676105,
240.542342898902,
0.600849252160861,
0.090829363604302,
0.112597000907675,
2.49789048014121,
3.03096710079162,
0.624028455006547,
0.265463130518279,
0.184613073383505,
0.16529963228303,
0.177559000212162,
0.057516859743295,
0.055278800156175,
0.13657251229218,
0.13804113884476,
0.127613995889303,
0.19032513064508,
0.170746411779942,
1.68338763417259,
2.2120957518041,
1.66082666882401,
2.60257486703841,
3.67829400184601,
1.46627228627455,
5040.00796685785,
0.388323620935189,
0.154298742280684,
0.288610303366468,
0.115024530338603,
0.893681227068193,
4.15427552466775
    };


        double[] sdten={0.715233616783147,
1.95899161327385,
1.41294001456335,
2.11258290367711,
1.40828236101682,
0.811334392227518,
302.069758999391,
0.374744935217507,
0.088811987422266,
0.12234337389111,
1.78057185206335,
1.68008309810733,
0.297945699977022,
1.4001232594843,
1.206691891812,
1.59455173554663,
2.01394569476323,
0.586949158265868,
161.446918858657,
0.419660893312224,
0.080348456007583,
0.098176210382283,
1.73559440748499,
2.08381990194041,
0.223312260113263,
0.228692177479827,
0.152856239575118,
0.135470125899722,
0.115011956888659,
0.01868845654684,
0.050875053992851,
0.133385413207911,
0.138600321311987,
0.131657890120888,
0.170962796443147,
0.123379600993575,
1.52931067004568,
1.30776369794955,
0.714472458033387,
1.94887690524881,
3.41227269465387,
1.41096346508033,
2857.55041619034,
0.378507179319225,
0.160242947186846,
0.258804438220394,
0.04075269116346,
0.72469990209225,
3.42971862243186
};
        
        
        double[] sd={1.33856304445208,
2.12884964908644,
1.4496687388267,
0.966148777238149,
2.5720501540236,
1.41916391012216,
172.595561203644,
0.534754908413121,
0.160489698856514,
0.127024419384134,
2.47459052537992,
2.3710470079981,
0.89512084546277,
1.69437109398182,
1.12406879836061,
1.31718731662705,
1.97398202673517,
0.806522972382953,
73.4415041606147,
0.630546139346848,
0.083910609287055,
0.113710932097979,
2.26753902461797,
2.88120241259806,
0.913599155968361,
0.254285887861418,
0.195906942655129,
0.172934586538636,
0.158033258996358,
0.074438606814368,
0.053073266695958,
0.176029540087753,
0.17578809152087,
0.159215389525663,
0.190204842947808,
0.16602635917295,
1.32020291130261,
1.61317196296432,
1.37045145965958,
2.14671354958667,
2.96025365453465,
1.23930312913513,
705.903680259474,
0.306520860979434,
0.16415776516175,
0.244820134676151,
0.217567396781175,
0.854612129834398,
4.08425433298923
    };
			
        double[][] loadingsEur={{-0.177708694035116,0.135043143980217,-0.014829163922633,-0.045630365756478,0.031447658259743,-0.001252346851566,-0.074486028747279,-0.084650914060972},
            {0.262788987728988,0.06644976229574,0.105055824841973,0.056972579236845,0.075666842751083,0.195407520052299,0.260708002440891,0.192075631959133},
            {-0.145073248241628,-0.388987855308478,-0.110142623859636,0.182864871223062,-0.001021575768633,0.036433342031581,-0.058106268190271,-0.126541856741582},
            {-0.249761405514343,-0.21963146918551,-0.118715256459125,0.095788478634302,0.169295583840668,0.029635636477928,-0.11723687479828,-0.124795666800199},
            {-0.115167890857104,0.229788032971897,0.026693664962345,-0.089826686748995,-0.021524153334275,-0.046801084965997,-0.062156404041742,-0.053958030363264},
            {-0.20179734984289,0.087463687235456,-0.005616723713995,0.022484152600132,0.041839228714516,0.015454157969257,-0.065391910977295,-0.076371321800993},
            {-0.173111340411226,0.1416008647945,-0.019356851456101,-0.074362122243914,0.008938650144593,-0.039703530285907,-0.082522785082472,-0.055310692515499},
            {0.136548544363576,-0.075816479246767,0.027472271482808,0.006178273081728,-0.18486451652181,-0.009338012070355,0.0598722809462,-0.134920024318776},
            {-0.124162480140835,0.124569783872691,-0.014163173533087,-0.153151397649396,-0.20056510897829,-0.288366432621786,-0.07692323657199,0.214917965435886},
            {0.125518447927742,-0.010948946917164,0.076058646636984,0.087696626849626,-0.159480013479776,-0.02837212127489,0.070684370765586,0.068676471162524},
            {-0.235653403689134,0.004155460641808,-0.059442317872394,0.003470225842474,-0.00696764715264,-0.010935751839211,-0.062340795434181,-0.037122011661628},
            {-0.082294021571036,0.119178103849667,0.086543035193519,-0.003820937748503,-0.17805113769043,0.050860068759232,0.162477300151541,-0.245517689272361},
            {-0.182326283390689,0.130098293929415,0.018433173412701,-0.013266090074159,0.019452168842352,0.018154281975701,0.130650976548179,0.113711558816968},
            {0.232308559681699,-0.023926230370846,0.000144833940219,-0.015853625876053,0.048404312793288,0.079983798497062,-0.002895013245245,-0.047168370669678},
            {-0.256965825592154,-0.169699305227333,0.011564884228748,0.114510464072536,-0.154563743585158,-0.057750972678883,-0.094534718604133,-0.035312752813463},
            {-0.239293061706232,-0.035816782886913,-0.006457116884198,0.088635918192088,-0.006948064336597,-0.129365833851763,-0.105851763970699,-0.099324902193897},
            {-0.194609343003836,0.09384970840506,0.008813676456996,0.017258821306469,0.037528101516706,-0.008896726503977,0.045104178572926,-0.011392449144234},
            {-0.228140372182563,0.016131210489905,0.057942053959117,0.10929162068196,-0.028757960193883,0.082943815370511,0.104798613912841,0.054757224223203},
            {-0.185650761973312,0.123176323336095,0.011321787243369,-0.02413117714149,0.008639503250856,-0.004520447720623,0.114106287936015,0.112027039706437},
            {-0.055199356755607,-0.151776159089233,0.043604770286374,-0.001931896954907,-0.309465501399977,0.244629333737483,0.10506241051456,0.033233766970767},
            {-0.161123010803586,-0.158622801980008,0.080923536510744,0.254513860790246,-0.029225479995191,0.011507963921498,0.126640924536774,0.236097891435171},
            {-0.040956870894175,-0.260759676939584,0.054939723075993,0.254566595475972,-0.034359690075322,-0.040361012362427,-0.041914484638777,0.075083883730915},
            {-0.241156043672186,-0.004476957653316,0.01511916580505,0.05905457252868,-0.052582813479738,-0.011420819477726,0.04661382661798,0.063275935828957},
            {-0.157066779629747,-0.029099152687009,0.078282239793002,0.025404800360506,-0.206721193086341,0.257209873388094,0.172640857731083,0.011574974550967},
            {-0.143364389152545,0.191430852490749,0.022106750523618,-0.059556189386752,9.32176391712876E-05,0.051506195308934,0.189214795629324,0.212844521912264},
            {-0.019161767065651,-0.10796437644887,0.072237133704279,-0.011658069825651,-0.326752929835455,0.103694660574659,-0.025909118201798,-0.038343905087668},
            {-0.019069720096728,-0.048295497371556,0.057842841803863,0.051411603024778,-0.118440404376361,0.06878090379446,-0.089388996198522,0.067004818645374},
            {-0.029035689856421,0.013414437431402,-0.051373498208272,0.030216533031735,0.144108595188792,-0.191239876710567,0.27556913779987,-0.288028979232349},
            {0.032342601809572,0.098400147073267,-0.124322358359127,-0.055133124209705,0.280679327522307,-0.170332711055493,-0.128622233758979,0.008627629559586},
            {-0.039255431288671,0.12298770515141,-0.095234275106617,-0.117422790820178,0.051912487696246,0.10469798075006,-0.154039125486193,0.127418986687067},
            {-0.014047367972775,0.203433603442228,0.523299158453023,0.193395964887235,0.140389011664121,0.029639658785607,-0.197845359905801,-0.322952652735514},
            {0.000147981316638,0.064132248010152,0.201142416592048,-0.103520806666602,-0.220807946610423,0.152236273375498,-0.162849986944895,-0.064496939930732},
            {0.021352751258096,0.058027549359592,0.204429525009925,-0.0936311723777,-0.247112314438774,0.008118321999127,-0.274596505136373,0.103736368194326},
            {-0.018198318601356,-0.004403374413734,0.108313036065263,0.111707872077669,-0.148590711687356,-0.423212305495223,0.268801000462068,-0.222580371652846},
            {0.033049484494831,-0.042829881310946,-0.174644942762658,0.098439685416916,0.089952057918536,-0.188696174322908,0.326214342524283,-0.029988134280124},
            {-0.084201129755078,-0.094675122854168,-0.239625122897583,-0.057347775932782,0.246788570866858,0.32974421242293,-0.254108615757003,0.00691651385509},
            {0.032116629872059,0.185388071799933,-0.159026490448425,0.260888306886599,-0.070783791008575,0.011327320649335,-0.056157433194074,0.076190671347662},
            {0.143199072523585,-0.127635340946332,0.117251001382556,-0.090639162745952,0.070712424533847,0.094007406939028,0.149495626588078,0.058721265986256},
            {0.049736518454777,-0.287813214215554,-0.014992844484936,0.120219572271104,0.080050887055911,0.091326666067218,-0.015149609391283,-0.048461298389553},
            {-0.114573074656267,0.154193776225673,0.175826958765384,0.362058525747316,0.193507735066319,0.010710661485098,-0.003134563940788,0.192219259776301},
            {0.082944544436336,0.104155684943538,-0.249288361443096,0.129569758845194,-0.184653925155898,-0.018554156320356,-0.06773917463319,0.048780391304934},
            {0.068514567868857,0.167504823304842,-0.085042952031524,0.399307883945382,0.029787704916163,0.12921009300292,0.056144871366476,0.220955071662361},
            {-0.034177587930368,0.204358176717521,-0.154524141549259,0.171190461023482,-0.084133355031459,0.023104662499369,-0.007409771161171,0.0411613254604},
            {-0.104533087029805,-0.104258280232026,0.252505758041688,-0.146399081549549,0.189626256265273,0.031680523507786,0.112105977263436,0.139565018580437},
            {0.140214168058834,0.06877514810434,-0.078589011418042,0.263877710011884,-0.129784639885312,-0.211738754911092,-0.225289200994209,0.151634240575059},
            {0.065686487489914,0.103509936031546,-0.297729372959147,0.055903141675748,-0.18575479506564,0.070575193385764,-0.030392396093379,-0.268851683648604},
            {-0.065297597929542,0.212914861434966,0.003899440813155,0.086539651666944,0.09993401782827,0.385336635882524,0.18571764254019,-0.264643954669812},
            {0.041835803037498,0.16491745118681,-0.23840097130653,0.106167975631537,-0.121424407217751,0.172153591160413,0.05392826684486,-0.211042259304935},
            {0.216290450831419,0.0268237415177,0.225634275242838,0.303453941101196,0.082849508144976,-0.031568031837992,-0.213066261173697,-0.149002665439142}};
    
    
    
        double[][] loadingsten={{-0.183782362062397,0.15361726576304,-0.063334970720672,-0.040146346956658,-0.016632738487994,0.00454284237725,0.012323635066188},
            {0.24634019837153,0.032256404647491,0.027272257638151,-0.028961133380933,0.09492960776028,0.095019275684254,0.086640075389439},
            {-0.228168284809058,-0.268498710009128,0.047137056654792,0.178690466073811,0.000985114733009,0.022137556862042,-0.0040085659879},
            {-0.22409018050005,-0.003042614588222,-0.091472628416238,0.031070311711212,0.111921440660652,0.038603752779815,-0.112596445636087},
            {-0.05960817113279,0.329592545963153,-0.014257574079545,-0.097932347387124,-0.086047205433842,0.045114237028269,0.138898782728471},
            {-0.199889621004176,0.112787657230528,-0.003305981809251,0.032774849750874,-0.027957407101142,0.022511046326509,0.018374293446186},
            {-0.175475688228897,0.163957317446917,-0.068845890615115,-0.062153017469602,-0.036122640483173,-0.020235094546918,0.026379578123156},
            {0.101359479760105,-0.100312020312014,0.084915888402651,-0.018634864047209,-0.234231309583948,-0.042351215972969,-0.142760475430993},
            {-0.055923778978547,0.094993399896085,0.087744273108761,-0.043196217787374,-0.387116630888251,-0.292891008412985,0.206835270271475},
            {0.098101333690991,-0.114295237457966,0.160266212722287,0.100293682216327,-0.177571621141939,-0.088856408905627,0.184389458590808},
            {-0.236645741934172,-0.00171571995489,-0.01311324757161,0.059576168373177,-0.052943988550443,-0.005290455687513,0.033472360874969},
            {-0.099647945104639,0.03053371127447,0.100026094366716,0.039579743274542,-0.175664369054206,0.199709530918234,-0.000311621997514},
            {-0.196486652326058,0.118208069963899,-0.052412661093236,-0.011675959848935,-0.0292776873907,-0.06498367912049,-0.110765927771065},
            {0.242381174951015,0.015653882116785,0.014800615089718,-0.055123790236571,0.064447914419873,0.077773944534337,0.10485080435643},
            {-0.224391243274866,-0.161757352681048,0.077163270396288,0.142002619672938,-0.118266815886616,0.049267342970276,0.032519701318288},
            {-0.236578658041011,-0.083789226518113,0.006400376971799,0.13671573945067,0.026307244824523,-0.044954276203383,-0.052224748518022},
            {-0.211045507480315,0.055789561618755,-0.059674280725296,0.035159589620044,0.033219969988879,-0.05562337050627,-0.167383101963586},
            {-0.213992530140863,0.047572990210761,0.039046012578911,0.095600328586302,-0.04804496961275,0.056689367888109,-0.007719112005185},
            {-0.197566522300924,0.116243594470773,-0.054350030915792,-0.016073598250541,-0.024166945732967,-0.053046168017894,-0.063776769433201},
            {0.013688096756108,-0.115873601192689,0.081864158661474,-0.024169064295483,-0.289125944149424,0.18280593715995,-0.055712095699485},
            {-0.074848043117956,-0.098126393830483,0.172424647017842,0.145075973766405,-0.197272339176829,-0.040627486769651,0.307096272407667},
            {0.012360354739438,-0.158765159903823,0.192968846183186,0.170774011464242,-0.122311415547479,0.031546368645375,0.310492125133286},
            {-0.238237732135796,-0.015312838663312,0.009814488460246,0.08854527243211,-0.061989435210442,0.009155050453487,-0.011994221379663},
            {-0.125248288021972,-0.016710712534397,0.093912875107359,0.041447099587707,-0.200512680144228,0.265895950400184,0.01271568962551},
            {-0.154662557234438,0.202491329190405,-0.056215355485775,-0.062287718347765,-0.06247798647108,-0.034194454586162,-0.00373801452509},
            {-0.00103593020937,-0.129796803194193,0.018167157535052,-0.033308886843385,-0.291575086828856,0.082197576218254,-0.25188637700776},
            {-0.000407998600837,-0.009100015458565,0.071740611884796,0.04113232630423,-0.056286941366123,-0.009277393867563,-0.173687035096881},
            {-0.034681587187162,0.062937013264732,0.006459707192011,0.051026933858516,0.159365556752329,-0.13078081166902,0.287225491972247},
            {-0.022174646877701,0.107008342671839,-0.074065981818198,-0.036481924586342,0.189581350542989,-0.0749133370583,0.133690542977909},
            {-0.026344666075313,0.174086070180301,-0.104429381504834,-0.20582704286252,-0.203609224057182,-0.135135591739889,0.074310945351494},
            {0.035095072525579,0.31367742514936,0.561150122557538,0.058902772368942,0.102796312109504,0.16009939735182,0.109593011941761},
            {0.064885441406185,0.05772686166911,0.083156828084764,-0.140244191805694,-0.178045746988959,0.300152563645859,-0.14729922924535},
            {0.082508387140864,0.033053934873747,0.031280539773741,-0.165055868742779,-0.208062451463129,0.131540347065371,-0.186449400131038},
            {-0.01083930083593,0.004361693685818,0.079620961032097,0.085319769365838,-0.151967476651572,-0.416513577113793,0.130591589880716},
            {-0.077793795743224,-0.087022517876952,-0.04811424220989,0.218775298146621,0.168116921031766,-0.226945760218853,-0.030253310109943},
            {-0.092347105122809,-0.02608764692787,-0.109743570506221,-0.072111678489899,0.138433190344658,0.216893433084064,0.098148561591818},
            {0.099462507507051,0.117999706404683,-0.146155242715382,0.283519931393636,-0.058389107953814,0.04546390084051,-0.010064865587462},
            {0.077924724633062,-0.170722364671267,0.113242266642474,-0.174212328611935,0.112108604516319,0.038288376568302,0.236844986628074},
            {-0.136797145240543,-0.366555899523748,0.056956715472915,0.144383946143864,0.201801814793518,0.247533799481453,0.028910278281119},
            {-0.062104500757232,0.334611115486623,0.219432699951327,0.30036408498653,0.130283178052087,0.14966049847749,-0.088020900767684},
            {0.11385525308362,-0.018476194078346,-0.20211922903624,0.205840726238599,-0.125958321107532,-0.104540588121414,-0.145911355615818},
            {0.115850669633567,0.124102578899373,-0.054305462163013,0.360195955965427,0.024966418635163,0.113625241888575,0.02294459565783},
            {0.059929987399336,0.16969796441669,-0.172305661166196,0.232567474924276,-0.088791174480338,0.098200893299385,0.075984222432162},
            {-0.12723470445796,-0.003395209351557,0.223220142680081,-0.19258164421607,0.107996233911011,-0.059250475345226,-0.060531860002825},
            {0.145583527186833,0.066098595329403,-0.028959465114839,0.280828585232197,-0.060896302977118,-0.107635784404478,-0.101851646522229},
            {0.095943160049197,-0.035960444378902,-0.308533644489189,0.107775334073703,-0.120241936447848,0.15329955654855,0.151592267400569},
            {-0.052013957225348,0.203336660940078,-0.114195263438709,0.007236597834457,-0.023246304366643,0.235110434003837,0.194504310621183},
            {0.09705454000899,0.021188765069992,-0.279205368028283,0.146284643338188,-0.111247782906267,0.157013930276187,0.162891670669553},
            {0.219307864248627,0.103561672345657,0.254690544353281,0.242896145046886,0.04263592055797,-0.1633493130975,-0.344057538802728}};			
		
        double[][] loadings={{-0.166620369116212,0.128895420220523,-0.061784768685505,0.004040075192428,-0.061887669484508,0.015705753303001,-0.135541405676834,-0.038726211867736},
                {0.279746632208979,0.056998625948644,0.056131260302331,0.087160782613944,0.071566137454795,0.177434536625853,0.147686332627191,0.127146767879547},
                {-0.187452112728977,-0.302110308576685,-0.015607668034225,-0.041773687601175,0.272624661168793,-0.084539965943007,0.036716948046652,-0.051420338739709},
                {-0.195383647681802,-0.136850333304669,0.031271257161807,0.080994601273359,0.01181033591935,-0.111216693366873,0.119448082514422,-0.255431078857862},
                {-0.165484222867283,0.128501201512541,-0.026060737656012,-0.016973254057765,-0.187373618636919,-0.089499022956233,-0.133490110590397,-0.167766478484296},
                {-0.182420642844356,0.112345365229347,-0.020324022953551,0.007768014258091,-0.019687674833043,-0.043995306849365,-0.204193996949199,-0.095170264205983},
                {-0.172558425914371,0.123121864019993,-0.057298332606485,-0.014835044976662,-0.117450341063667,-0.016401546221827,-0.116480868349203,-0.019645720025326},
                {0.05160358274625,-0.113561571159411,0.018615822717013,-0.227105463922672,0.120331017140105,-0.008254499098966,0.013855595159914,0.183836816858142},
                {-0.145744742290852,0.083612850720384,-0.003749517544625,-0.052966659491467,-0.100689911384356,-0.119963833660815,-0.205318309883738,0.277619476580677},
                {0.053440888574283,-0.071623501088434,-0.020322411985499,0.036549130804956,0.372213169849781,0.060463315165044,-0.270772056583283,0.14137693598703},
                {-0.234221674235222,0.032283035567858,-0.062032495611903,-0.013316283918674,0.036682609434008,-0.007701180245352,-0.122327146716598,-0.029515475169268},
                {-0.012542389526164,0.003411963459651,0.165996973804779,-0.207670592332764,0.021776984344568,-0.03796533641182,0.051703967674287,0.278892617948932},
                {-0.170656339522651,0.114865017429063,0.03926486832919,0.044299022192297,-0.044178070995551,0.134540084625251,0.133447029114732,0.041786235974551},
                {0.256644937807026,0.003332559958981,-0.023504835426641,0.066554022461523,0.053918721294277,0.050819184582924,0.004870351689455,0.002871926295122},
                {-0.247494458434039,-0.165725855124926,0.067317428630211,-0.163850051683005,0.126011250013683,-0.072909633518743,-0.053741378317896,0.067540411695941},
                {-0.247095892890199,-0.043300536806101,0.003564663399179,-0.017879680074261,-0.035014592342203,-0.149065019433126,-0.017020087239969,0.008879737908201},
                {-0.20254481602317,0.081712575289887,0.013387457896961,0.024650995850414,-0.096269857221509,-0.026497612633377,0.001561132979679,-0.12094692862129},
                {-0.227085782837202,0.04695801732269,0.080316294371452,-0.02961854514931,0.006786399363423,0.0372659726381,0.015991453843088,-0.025614672820814},
                {-0.205287447066714,0.085489348177031,0.020245725968447,0.012081660348663,-0.086102168512284,0.037343856381916,0.072562417866455,0.014619721126797},
                {-0.059900056683411,-0.124084460003568,0.074385738862973,-0.368836948574589,0.071300843217549,0.170266069414763,0.035736930618874,0.083895664871741},
                {0.184750967912749,-0.052721764399532,0.078287701672727,0.033896503760636,0.212214886896933,-0.045252074542245,-0.140845962514859,0.139215581077227},
                {0.012544435879144,-0.176720117776923,0.082477656674902,0.007030229133674,0.253199630615684,-0.138709604551654,-0.150388830001867,0.103967674888587},
                {-0.242971688552699,0.024122767503028,0.037205571648995,-0.017161891848385,0.019052427571079,0.047238348634714,0.048725120164735,0.01976218554194},
                {-0.100348430344508,-0.059023440502887,0.156134437805396,-0.22887904748851,0.104256492196955,0.122065169218288,-0.011884301273742,0.053474388103788},
                {-0.155654036187184,0.147324725338388,0.004652130780744,0.036196061995071,-0.049820578142114,0.130488799988095,0.062102834241589,0.055798874351039},
                {-0.017518454502675,-0.08504999694307,0.098347961523264,-0.285096885276162,0.037289063303173,0.075611177388305,-0.067520696405225,0.047667138723792},
                {-0.02336878175342,-0.045804927960362,0.075972103359209,-0.188723706165377,-0.017537669986237,0.10405552850842,0.010692039387871,-0.155377006187783},
                {-0.063636111645484,0.010960007551115,0.012132463365252,0.093645031251555,-0.089319130738608,-0.212886320174271,0.21891908907476,0.168000469128409},
                {-0.002276027984419,0.047537285400558,-0.171749620710046,0.224475401732943,0.035654774496469,-0.161017210644236,-0.075983273599726,0.103485964876557},
                {-0.014260883698654,0.080990564377653,-0.219823425641009,0.086469504920265,0.166625088045113,0.103720561302019,-0.251961570729343,0.08787641146787},
                {-0.014489491965012,0.19049682574609,0.582761242258047,0.305845798632518,0.021456703373894,-0.081134503648739,-0.353112331898722,0.157512792292645},
                {0.032192851310357,0.033251914532303,0.152417325454769,-0.220218683407384,-0.108528229246674,0.196586846744548,-0.116330361826877,-0.057088171378766},
                {0.043660492408815,0.035154414476142,0.184165650877017,-0.229316401248577,-0.127030721882958,0.159322717841983,-0.101987362447223,-0.154305862149647},
                {-0.068475975196447,-0.025154172132004,0.082580489714317,-0.022985345352009,-0.098362995448359,-0.333110878124149,0.189703185942304,0.329441445681429},
                {-0.061377554642163,-0.077840070992461,-0.151790218392694,0.119625233022235,0.147914884645498,-0.213722813826072,0.153890914854548,0.102065585866957},
                {-0.041314836093956,-0.02767903212385,-0.268737409966972,0.1292589576985,0.213459451927808,0.119733423826442,-0.162947375729128,-0.057029554044069},
                {-0.000661199985369,0.252843349250827,-0.068602072255113,-0.145601779659645,0.146600954411495,-0.102741914407044,0.096370348310868,-0.025314891834753},
                {0.159009211982334,-0.109731411323987,-0.055651890336802,0.047079772993288,-0.125169833786819,0.057527474596309,-0.198705854620988,-0.020334553102954},
                {-0.028669945967961,-0.25544162867312,0.022770850166911,0.011496159737917,0.245124653847855,-0.09441403609055,-0.004432820848447,-0.305558027986685},
                {-0.1105081210731,0.044555659335849,0.340627015777117,0.203630829188983,0.220046708445777,0.011077133780623,0.310175684526789,-0.354254509290075},
                {0.031808025330311,0.28105565647291,-0.052282176528261,-0.148061147916376,0.078044772774438,-0.084111122127253,0.102483229206012,0.071937747723693},
                {0.058397063516856,0.283679406919544,0.069594264018011,-0.046792443671583,0.3208613332752,0.011070032655021,0.271304645092224,0.059679012245067},
                {-0.064524954557761,0.238665759571723,-0.030605800970305,-0.081637456571908,0.105146238076795,0.013091397268965,0.15067586563143,0.128950395979741},
                {-0.047497510928023,-0.25288021023573,0.122135700978381,0.19055630038043,-0.129577092001791,0.193515697733349,0.118193310195063,0.174164867655428},
                {0.160113485667188,0.131379875593219,-0.134868915234525,-0.326214883348362,0.006156976920607,-0.446946686987764,-0.076489559011845,-0.183444925842742},
                {-0.047891869885833,0.228518920943788,-0.062484078467528,-0.019846304925577,0.158105322029002,0.057402022994037,-0.096692864615613,-0.095054502605769},
                {-0.09234985503902,0.187072803257845,-0.043611779225768,0.086220529990819,0.195736898727236,0.297097532822489,0.054864737209237,0.065269446092855},
                {-0.044733885693107,0.225308901834641,-0.107241797932954,-0.044972998920748,0.173695291214132,0.100843364565521,-0.022574890052296,-0.045867537657268},
                {0.243222369873139,0.150593021440182,0.349006690421896,0.000563630363581,0.124093824611887,-0.270980299515074,-0.141339473666601,-0.173336834206536}};
        
        
	double[] out=new double[loadings[0].length];
		
		//System.out.println(loadings.length+" "+loadings[0].length);
		//System.out.println(x.length);
		
		//for (int j=0; j<loadings.length; j++){
		//	System.out.println(j+" "+x[j]+" "+means[j]+" "+sd[j]+" "+(x[j]-means[j])/sd[j]);
		//}
		
		
		
	for (int i=0; i<loadings[0].length; i++){
            for (int j=0; j<loadings.length; j++){
		out[i]+=((x[j]-means[j])/sd[j])*loadings[j][i];
                
                if(Double.isNaN(x[j])){System.out.println(currentpop+" "+j);}
				//System.out.println(i+" "+j+" "+out[i]+" "+x[j]+" "+means[j]+" "+sd[j]+" "+loadings[j][i]);
            }
        }
                
                /*
                for (int i=0; i<x.length; i++){
                    System.out.println("IN: "+x[i]);
                }
                for (int i=0; i<out.length; i++){
                    System.out.println("OUT: "+out[i]);
                }
                */
	return out;
    }
    
    
    public int passThresh(double threshold){
        int count=0;
        for (int i=0; i<score.length; i++){
            if (score[i]<threshold){count++;}
        }
        
        return count;
    }
    
    
     public void writeRound(){
            	DocumentSave ds=new DocumentSave(fileloc+"/outputstats"+id+".csv", ", ");

                for (int i=0; i<stats.length; i++){
                    for (int j=0; j<stats[i].length; j++){
                        ds.writeDouble(stats[i][j]);
                    }
                    ds.writeDouble(score[i]);
                    ds.writeLine();      
                }
		
		ds.finishWriting();
                
                ds=new DocumentSave(fileloc+"/outputparams"+id+".csv", ", ");
		
                for (int i=0; i<params.length; i++){
                    for (int j=0; j<params[i].length; j++){
                        ds.writeDouble(params[i][j]);
                    }
                    ds.writeLine();   
                }	
		ds.finishWriting();       
        }
    
    
    public static void main (String args[]) {
        
        int a=Integer.parseInt(args[0]);
        int b=Integer.parseInt(args[1]);
        int c=Integer.parseInt(args[2]);
        long d=Long.parseLong(args[3]);
        String[] e=new String[1];
        e[0]=args[4];
       int f=1;
        Priors p=new Priors(System.currentTimeMillis());
        double[] x=p.sampleFromPriors();
        ChaffinchABC[] cemp=new ChaffinchABC[1];
        cemp[0]=new ChaffinchABC(e[0], x, p.variables, System.currentTimeMillis());
        
        
	SimulationRunner sr=new SimulationRunner(cemp, 1, a,b,c,d,e, null, null, null);
	sr.runSimulation();
        
        
    }
		
	
}
        
        
        
    

