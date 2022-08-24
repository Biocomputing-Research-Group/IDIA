/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.PeptidePeakClusterDetection;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.ScanCollection;
import MSUmpire.BaseDataStructure.ScanData;
import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.BaseDataStructure.XYZData;
import MSUmpire.DIA.CorrCalcCluster2CurveUnit;
import MSUmpire.LCMSPeakStructure.LCMSPeakBase;
import MSUmpire.MathPackage.GaussianFilter;
import MSUmpire.PeakDataStructure.IonChargeHashSet;
import MSUmpire.PeakDataStructure.PeakCluster;
import MSUmpire.PeakDataStructure.PeakCurve;
import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;

import java.io.*;
import java.lang.reflect.Array;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;

import net.sf.javaml.core.kdtree.KDTree;
import net.sf.javaml.core.kdtree.KeySizeException;
import umontreal.iro.lecuyer.util.Num;
import weka.core.parser.java_cup.internal_error;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.log4j.Logger;
//import org.jfree.chart.ChartFactory;
//import org.jfree.chart.ChartUtilities;
//import org.jfree.chart.JFreeChart;
//import org.jfree.chart.plot.PlotOrientation;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import com.compomics.util.experiment.biology.ions.ElementaryIon;

//import org.eclipse.collections.impl.list.mutable.primitive.BooleanArrayList;
import GridMass.module.GridMassModel;
import GridMass.datamodule.IsotopicTraces;

/**
 * Peak detection processing parent class
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PDHandlerBase {
    
//    protected HashSet<String> IncludedHashMap;
    protected HashMap<XYData,Boolean> InclusionFound=new HashMap();
    protected XYPointCollection InclusionRT=new XYPointCollection();
    protected KDTree InclusionRange=new KDTree(2);    
    protected int NoCPUs = 4;
    public float minSNR;
    public TreeMap<Float, XYData>[] IsotopePatternMap;
    public TreeMap<Float, XYData>[] IsotopePatternFragMap;
    protected LCMSPeakBase LCMSPeakBase;
    protected InstrumentParameter parameter;
    protected boolean ReleaseScans = true;
//    protected boolean TargetedOnly = false;
    protected float PPM;
    public int MSlevel=1;

    public PDHandlerBase() {
    }

//    public void SetTargetedDetectionOnly(){
//        TargetedOnly=true;
//    }
    
    protected String InclusionCheckInfo(){
        int count=0;
        
        for(XYData point : InclusionFound.keySet()){
            boolean value = InclusionFound.get(point);
            if(value){
                count++;
            }
            else{
                //Logger.getRootLogger().warn("Missing signals: mz="+point.getX()+", RT="+point.getY());
            }
        }
        return count+"/"+InclusionFound.size();                
    }
    
    //Add mz and RT coordinate to inclusion list
     public void AddToInclusionList(float mz, float rt){
         XYData point=new XYData(mz, rt);
        InclusionFound.put(point, false);
        InclusionRT.AddPoint(rt,mz);        
        try {
            InclusionRange.insert(new double[]{rt,mz}, point);
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
          
    public void ClearAllPeaks() {
        LCMSPeakBase.BaseClearAllPeaks();
    }

    public void ClearRawPeaks() {
        for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
            peakCurve.CalculateMzVar();
            peakCurve.StartRT();
            peakCurve.EndRT();
            peakCurve.ReleaseRawPeak();
        }
    }

    //Detect all m/z trace peak curves given a list of ScanCollection
    public void FindAllMzTracePeakCurvesForScanCollections(ArrayList<ScanCollection> scanCollections, GridMassModel gridmass,float strt, float endrt,String msLevel) throws IOException {
        //Read peptide isotope pattern table
        ReadPepIsoMS1PatternMap();
        LCMSPeakBase.UnSortedPeakCurves = new ArrayList<>();
        
        for (ScanCollection scanCollection : scanCollections) {
            //Detect mz trace peak curves for each ScanCollection
            FindAllMzTracePeakCurves(scanCollection, gridmass, strt, endrt,"ms1");
        }
        gridmass = null;
        Logger.getRootLogger().info("Inclusion mz values found: "+InclusionCheckInfo());        
    }
    
    static private int int_id(final int[] idx_arr,final int scannum_idx, final int peakindex){
        final int ret = idx_arr[scannum_idx]+peakindex;
        assert ret<idx_arr[scannum_idx+1]:"indexing error";
        return ret;
    }
    static int step = -1;
    //Detect all m/z trace / peak curves
    protected void FindAllMzTracePeakCurves(ScanCollection scanCollection, GridMassModel gridmass, float scanstrt, float scanendrt, String msLevel) throws IOException {
        Logger.getRootLogger().info("Processing all scans to detect possible m/z peak curves and");
        Logger.getRootLogger().info("Smoothing detected signals......");
        float preRT = 0f;
        
        //Loop for each scan in the ScanCollection
        final ArrayList<ForkJoinTask<ArrayList<PeakCurve>>> ftemp = new ArrayList<>();
        final ForkJoinPool fjp=new ForkJoinPool(NoCPUs);
        final int idx_end = scanCollection.GetScanNoArray(MSlevel).size();


        final int[] ia=new int[idx_end+1];
        ia[0]=0;
        for (int idx = 0; idx < idx_end; idx++) {
            final int scanNO = scanCollection.GetScanNoArray(MSlevel).get(idx);
            final ScanData sd=scanCollection.GetScan(scanNO);
            ia[idx+1]=sd.Data.size()+ia[idx];
        }

        final boolean[] included = new boolean[ia[ia.length-1]];
        if (step == -1)
            step = fjp.getParallelism() * 32;
        long peakCurvesCount = 0;
        for (int idx = 0; idx < idx_end; idx++) {
            int scanNO = scanCollection.GetScanNoArray(MSlevel).get(idx);
            ScanData scanData = scanCollection.GetScan(scanNO);
            
//            //If we are doing targeted peak detection and the RT of current scan is not in the range of targeted list, jump to the next scan 
//            if(TargetedOnly && !FoundInInclusionRTList(scanData.RetentionTime)){
//                continue;
//            }
            if (idx == 0) {
                preRT = scanData.RetentionTime - 0.01f;
            }
            for (int i = 0; i < scanData.PointCount(); i++) {
                XYData peak = scanData.Data.get(i);
                //If we are doing targeted peak detection and the RT and m/z of current peak is not in the range of targeted list, jump to the next peak 
//                if (TargetedOnly && !FoundInInclusionMZList(scanData.RetentionTime,peak.getX())) {
//                    continue;
//                }
                
                if(peak.getX()<parameter.MinMZ){
                    continue;
                }
                
                //Check if the current peak has been included in previously developed peak curves
//                if (!IncludedHashMap.contains(scanNO + "_" + peak.getX())) {//The peak hasn't been included
                final int id_scanNO_peak=int_id(ia, idx, i);
                if (!included[id_scanNO_peak]) {//The peak hasn't been included
                    //The current peak will be the starting peak of a new peak curve
                    //Add it to the hash table

//                    IncludedHashMap.add(scanNO + "_" + peak.getX());
                    included[id_scanNO_peak] = true;

                    float startmz = peak.getX();
                    float startint = peak.getY();
                    
                   //Find the maximum peak within PPM window as the starting peak
                    for (int j = i + 1; j < scanData.PointCount(); j++) {
                        XYData currentpeak = scanData.Data.get(j);
                        final int id_scanNO_currentpeak = int_id(ia,idx,j);
                        if (!included[id_scanNO_currentpeak]) {
//                        if (!IncludedHashMap.contains(scanNO + "_" + currentpeak.getX())) {
                            if (InstrumentParameter.CalcPPM(currentpeak.getX(), startmz) <= PPM) {
                                included[id_scanNO_currentpeak]=true;
//                                IncludedHashMap.add(scanNO + "_" + currentpeak.getX());

                                if (currentpeak.getY() >= startint) {
                                    startmz = currentpeak.getX();
                                    startint = currentpeak.getY();
                                }
                            } else {
                                break;
                            }
                        }
                    }

                    //Initialize a new peak curve
                    PeakCurve Peakcurve = new PeakCurve(parameter);
                    //Add a background peak
                    Peakcurve.AddPeak(preRT, startmz, scanData.background);
                    //Add the starting peak
                    Peakcurve.AddPeak(scanData.RetentionTime, startmz, startint);
                    Peakcurve.StartScan = scanNO;

                    int missedScan = 0;
                    float endrt=scanData.RetentionTime;
                    int endScan=scanData.ScanNum;
                    float bk=0f;

                     //Starting from the next scan, find the following peaks given the starting peak
                    for (int idx2 = idx + 1; idx2 < scanCollection.GetScanNoArray(MSlevel).size() && (missedScan < parameter.NoMissedScan /*|| (TargetedOnly && Peakcurve.RTWidth()<parameter.MaxCurveRTRange)*/); idx2++) {
                        int scanNO2 = scanCollection.GetScanNoArray(MSlevel).get(idx2);
                        ScanData scanData2 = scanCollection.GetScan(scanNO2);

                        endrt=scanData2.RetentionTime;
                        endScan=scanData2.ScanNum;
                        bk=scanData2.background;
                        float currentmz = 0f;
                        float currentint = 0f;

                        //If the scan is empty
                        if (scanData2.PointCount() == 0) {                            
                            if (parameter.FillGapByBK) {
                                Peakcurve.AddPeak(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background);
                            }
                            missedScan++;
                            continue;
                        }

                        //Find the m/z index 
                        int mzidx = scanData2.GetLowerIndexOfX(Peakcurve.TargetMz);
                        for (int pkidx = mzidx; pkidx < scanData2.Data.size(); pkidx++) {
                            XYData currentpeak = scanData2.Data.get(pkidx);
                            if (currentpeak.getX() < parameter.MinMZ) {
                                continue;
                            }
                            //Check if the peak has been included or not
                            final int int_id_scanNO2_currentpeak=int_id(ia,idx2,pkidx);
//                            if (!included.get(int_id_scanNO2_currentpeak)) {
                            if (!included[int_id_scanNO2_currentpeak]) {
                                if (InstrumentParameter.CalcPPM(currentpeak.getX(), Peakcurve.TargetMz) > PPM) {
                                    if (currentpeak.getX() > Peakcurve.TargetMz) {
                                        break;
                                    }
                                } else {
                                    //////////The peak is in the ppm window, select the highest peak
                                    included[int_id_scanNO2_currentpeak] = true;
//                                    IncludedHashMap.add(scanNO2 + "_" + currentpeak.getX());
                                    if (currentint < currentpeak.getY()) {
                                        currentmz = currentpeak.getX();
                                        currentint = currentpeak.getY();
                                    }
                                }
                            }
                        }
                        
                        //No peak in the PPM window has been found
                        if (currentmz == 0f) {
                            if (parameter.FillGapByBK) {
                                Peakcurve.AddPeak(scanData2.RetentionTime, Peakcurve.TargetMz, scanData2.background);
                            }
                            missedScan++;
                        } else {
                            missedScan = 0;
                            Peakcurve.AddPeak(scanData2.RetentionTime, currentmz, currentint);
                        }
                    }
                    Peakcurve.AddPeak(endrt, Peakcurve.TargetMz, bk);
                    Peakcurve.EndScan=endScan;
                    								

                    //First check if the peak curve is in targeted list
                    if (FoundInInclusionList(Peakcurve.TargetMz, Peakcurve.StartRT(), Peakcurve.EndRT()) ) {
//                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
                        ++peakCurvesCount;
                        ftemp.add(fjp.submit(new PeakCurveSmoothingUnit(Peakcurve, parameter,msLevel)));
                    //Then check if the peak curve passes the criteria
                    } else if (Peakcurve.GetRawSNR() > LCMSPeakBase.SNR && Peakcurve.GetPeakList().size() >= parameter.MinPeakPerPeakCurve + 2 )  {
//                        LCMSPeakBase.UnSortedPeakCurves.add(Peakcurve);
                        ++peakCurvesCount;
                        ftemp.add(fjp.submit(new PeakCurveSmoothingUnit(Peakcurve, parameter,msLevel)));
                    } else {
                        Peakcurve = null;
                    }
                }
            }                               
            preRT = scanData.RetentionTime;
            /** the if statement below does PeakCurveSmoothing() and ClearRawPeaks()
             */
            final boolean last_iter = idx + 1 == idx_end;
            if (ftemp.size() == step || last_iter) {
                final List<ForkJoinTask<ArrayList<PeakCurve>>> ftemp_sublist_view = 
                        last_iter?
                        ftemp
                        : ftemp.subList(0, step/2);
                for(final Future<ArrayList<PeakCurve>> f : ftemp_sublist_view){
                    try{LCMSPeakBase.UnSortedPeakCurves.addAll(f.get());}
                    catch(InterruptedException|ExecutionException e){throw new RuntimeException(e);}
                }
                ftemp_sublist_view.clear();
                if(!last_iter && fjp.getActiveThreadCount()<fjp.getParallelism()){
                    step *= 2;
                }
            }
        }
        
        assert ftemp.isEmpty();
        
//        Logger.getRootLogger().info("diaumpire"+LCMSPeakBase.UnSortedPeakCurves.size());
        
        
        ArrayList<IsotopicTraces> isotopicTraces = gridmass.startrun((double) scanstrt,(double) scanendrt, MSlevel);
//        Logger.getRootLogger().info("gridmass"+isotopicTraces.size());

        
        //Assign peak curve index
        
        ArrayList<PeakCurve> tmppeakcurves = new ArrayList();
        TreeMap<Float, Integer> peaksTargetMZ = new TreeMap();
        int ii = 0;
        for(PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
        	peaksTargetMZ.put(peakCurve.TargetMz, ii);
        	ii++;
        }

        int count_overlapping =0;
        for(IsotopicTraces isotp : isotopicTraces) {
        	SortedMap<Float, Integer> rangeTargetMZ = peaksTargetMZ.subMap((float) (isotp.getTargetMZ()-0.02), (float) (isotp.getTargetMZ()+0.02));
//        	SortedMap<Float, Integer> rangeTargetMZ = peaksTargetMZ.subMap((float) (isotp.getTargetMZ()-0.02), (float) (isotp.getTargetMZ()+0.02));
        	PeakCurve tmpcurve = transformdataformat(isotp);
        	if (tmpcurve==null) {
        		continue;
        	}
        	boolean newpeakcurve = true;
        	for(float i : rangeTargetMZ.keySet()) {
        		float diffstartRT = Math.abs(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).StartRT() - (float) (isotp.getRawDataPointsRTRange().lowerEndpoint()*60/60));
        		float diffendRT = Math.abs(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).EndRT() - (float) (isotp.getRawDataPointsRTRange().upperEndpoint()*60/60));
        		float diffApex = Math.abs(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).ApexRT - (float) isotp.getRT());
        		float tmpstart = Math.max(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).StartRT() , (float) (isotp.getRawDataPointsRTRange().lowerEndpoint()*60/60));
        		float tmpend = Math.min(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).EndRT(),(float) (isotp.getRawDataPointsRTRange().upperEndpoint()*60/60));
        		float tmpstartl = Math.min(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).StartRT() , (float) (isotp.getRawDataPointsRTRange().lowerEndpoint()*60/60));
        		float tmpendl = Math.max(LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)).EndRT(),(float) (isotp.getRawDataPointsRTRange().upperEndpoint()*60/60));
        		if( diffstartRT<=0.5 && diffendRT <= 0.5 && diffApex<=0.2 && tmpend-tmpstart>0 && (tmpend - tmpstart)/(tmpendl - tmpstartl)>0.3) {
        			float corr =-1f;
                    try {
                        corr = PeakCurveCorrCalc.CalPeakCorr(tmpcurve, LCMSPeakBase.UnSortedPeakCurves.get(rangeTargetMZ.get(i)), parameter.NoPeakPerMin);
                    } catch (IOException ex) {
                    	Logger.getRootLogger().info("error");
                    }
                    if(corr>parameter.similarityTwoMethods) {
                    	count_overlapping += 1;
                    	newpeakcurve = false;        			
                    	break;
                    }
        		}        		
        	}
        	if(newpeakcurve && tmpcurve!=null) {
        		tmppeakcurves.add(tmpcurve);
        		
        	}
        	
        }
//        Logger.getRootLogger().info("number of overlapping of diaumpire and gridmass:"+String.valueOf(count_overlapping));
        gridmass = null;
        
        
        int i = 1;
        for (PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
            peakCurve.Index = i++;
        }
        
        for(PeakCurve peakCurve : tmppeakcurves) {
        	peakCurve.Index = i++;
        	LCMSPeakBase.UnSortedPeakCurves.add(peakCurve);
        }
        tmppeakcurves = null;
/*--------------------------------------------------------------*/
//  	   String csvFile = "/home/UNT/jl0734/Desktop/DIA/myeclipseproject/diaumpire.csv";
//  	   try {
//  		   BufferedWriter bw2 = Files.newBufferedWriter(Paths.get(csvFile), StandardCharsets.ISO_8859_1);
//  		   bw2.append("mz,rt,height,startrt,endrt\n");
//  		   for(PeakCurve peakCurve : LCMSPeakBase.UnSortedPeakCurves) {
//  			   bw2.append(String.valueOf(peakCurve.TargetMz)+","+String.valueOf(peakCurve.ApexRT)+","+String.valueOf(peakCurve.ApexInt)+","+String.valueOf(peakCurve.StartRT())+","+String.valueOf(peakCurve.EndRT())+"\n");
//  		   }
//  		   bw2.close();
//  	   }
//  	   catch (Exception e) {
// 		e.printStackTrace();
//  	   }
//  	   
//  	   String csvFile1 = "/home/UNT/jl0734/Desktop/DIA/myeclipseproject/diaumpiresame.csv";
//  	   try {
//  		   BufferedWriter bw21 = Files.newBufferedWriter(Paths.get(csvFile1), StandardCharsets.ISO_8859_1);
//  		   bw21.append("mz,startrt,endrt,mz,startrt,endrt\n");
//  		   for(float[] tt: samepeak) {
//  			   bw21.append(String.valueOf(tt[0])+","+String.valueOf(tt[1])+","+String.valueOf(tt[2])+","+String.valueOf(tt[3])+","+String.valueOf(tt[4])+","+String.valueOf(tt[5])+"\n");
//  		   }
//  		   bw21.close();
//  	   }
//  	   catch (Exception e) {
// 		e.printStackTrace();
//  	   }
//        Logger.getRootLogger().info("combined"+LCMSPeakBase.UnSortedPeakCurves.size());
        System.gc();
//        Logger.getRootLogger().info(LCMSPeakBase.UnSortedPeakCurves.size() + " Peak curves found (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
    }
    
    public PeakCurve transformdataformat(IsotopicTraces isotopicTrace) {
         
         PeakCurve resultcurve = new PeakCurve(parameter);
         resultcurve.RegionRidge = null;     
         Set<Float> checkset = new HashSet<>();
         
         for (Integer index:isotopicTrace.getDataPoint().keySet()) {
        	 resultcurve.AddPeak((float) isotopicTrace.getDataPoint().get(index).getRT(), (float) isotopicTrace.getDataPoint().get(index).getMZ(), (float) isotopicTrace.getDataPoint().get(index).getIntensity());
         }

        	 resultcurve.DoBspline();
        	 resultcurve.GetSmoothedList().Data.Finalize();
        	 resultcurve.CalculateMzVar();
        	 resultcurve.StartRT();
        	 resultcurve.EndRT();
        	 
             GaussianFilter gau = new GaussianFilter(resultcurve.GetPeakCollection(), resultcurve.StartRT(), resultcurve.EndRT(),"gridmass");
             gau.gausfilter();
             resultcurve.SmoothData = gau.getsmoothing();
//             resultcurve.newdata = gau.getsmoothing();
             gau=null;
        	 
        	 resultcurve.ReleaseRawPeak();
        	 resultcurve.types="gridmass";
//         }
//         else {
//        	 resultcurve=null;
//         }
//         checkset=null;
         return resultcurve;
    }
        
    
    private boolean FoundInInclusionMZList(float rt, float mz) {
        if(InclusionRT.PointCount()==0){
            return false;
        }
        float lowrt = rt - parameter.MaxCurveRTRange;
        float highrt = rt + parameter.MaxCurveRTRange;
        float lowmz = InstrumentParameter.GetMzByPPM(mz, 1, PPM);
        float highmz = InstrumentParameter.GetMzByPPM(mz, 1, -PPM);
       
        Object[] found=null;
        try {
            found = InclusionRange.range(new double[]{lowrt,lowmz}, new double[]{highrt,highmz});
        } catch (KeySizeException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        if(found!=null && found.length>0){
            return true;
        }
        return false;        
    }
    
    private boolean FoundInInclusionList(float mz, float startrt, float endrt){       
        if(InclusionRT.PointCount()==0){
            return false;
        }
        float lowmz = InstrumentParameter.GetMzByPPM(mz, 1, PPM);
        float highmz = InstrumentParameter.GetMzByPPM(mz, 1, -PPM);
        float lowrt=startrt-parameter.RTtol ;
        float highrt=endrt+parameter.RTtol;
                
        Object[] found=null;
        try {
            found = InclusionRange.range(new double[]{lowrt,lowmz}, new double[]{highrt,highmz});
        } catch (KeySizeException ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
        if(found!=null && found.length>0){
            for(Object point: found){
                InclusionFound.put((XYData) point,true);
            }
            return true;
        }
        return false;     
    }


    //Load pre-built peptide isotope pattern table
    protected void ReadPepIsoMS1PatternMap() throws FileNotFoundException, IOException {

        InputStream is = this.getClass().getClassLoader().getResourceAsStream("resource/IsotopicPatternRange.csv");
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        IsotopePatternMap = new TreeMap[Math.max(2,LCMSPeakBase.MaxNoPeakCluster - 1)];
        for (int i = 0; i < IsotopePatternMap.length; i++) {
            IsotopePatternMap[i] = new TreeMap<>();
        }
        String line = "";
        while ((line = reader.readLine()) != null) {
            float MW = Float.parseFloat(line.split(",")[0]);

            for (int i = 0; i < IsotopePatternMap.length; i++) {
                float Mean = Float.parseFloat(line.split(",")[1 + (i * 2)]);
                float SD = Float.parseFloat(line.split(",")[2 + (i * 2)]);

                if (!Float.isNaN(Mean)) {
                    IsotopePatternMap[i].put(MW, new XYData(Mean + 3.3f * SD, Mean - 3.3f * SD));
                }
            }
        }
        reader.close();
    }

    static private int step_pccc=-1;
    //Group peak curves based on peak profile correlation of isotope peaks
    protected void PeakCurveCorrClustering(XYData mzRange) throws IOException{
        Logger.getRootLogger().info("Grouping isotopic peak curves........");

        LCMSPeakBase.PeakClusters = new ArrayList<>();
        
        //Thread pool
        final ForkJoinPool fjp = new ForkJoinPool(NoCPUs);
//        ArrayList<PeakCurveClusteringCorrKDtree> ResultList = new ArrayList<>();
        final ArrayList<ForkJoinTask<ArrayList<PeakCluster>>> ftemp = new ArrayList<>();
        final int end_idx=LCMSPeakBase.UnSortedPeakCurves.size();
        final ArrayList<PeakCluster> resultClusters=new ArrayList<>();
        //For each peak curve
//        for (PeakCurve Peakcurve : LCMSPeakBase.UnSortedPeakCurves) {
        for (int i=0;i<end_idx;++i) {
            final PeakCurve Peakcurve = LCMSPeakBase.UnSortedPeakCurves.get(i);
            if (Peakcurve.TargetMz >= mzRange.getX() && Peakcurve.TargetMz <= mzRange.getY()) {
                //Create a thread unit for doing isotope clustering given a peak curve as the monoisotope peak
                PeakCurveClusteringCorrKDtree unit = new PeakCurveClusteringCorrKDtree(Peakcurve, LCMSPeakBase.GetPeakCurveSearchTree(), parameter, IsotopePatternMap, LCMSPeakBase.StartCharge, LCMSPeakBase.EndCharge, LCMSPeakBase.MaxNoPeakCluster, LCMSPeakBase.MinNoPeakCluster);
//                ResultList.add(unit);
                ftemp.add(fjp.submit(unit));
            }
            if(step_pccc==-1)
                step_pccc=fjp.getParallelism()*32;
            final boolean last_iter = i + 1 == end_idx;
            if(ftemp.size()==step_pccc || last_iter){
                final List<ForkJoinTask<ArrayList<PeakCluster>>> ftemp_sublist_view=
                        last_iter?
                        ftemp:
                        ftemp.subList(0, step_pccc/2);
                for (final ForkJoinTask<ArrayList<PeakCluster>> fut : ftemp_sublist_view)
                    try {resultClusters.addAll(fut.get());}
                    catch (InterruptedException | ExecutionException ex) {throw new RuntimeException(ex);}
                ftemp_sublist_view.clear();
                if(!last_iter && fjp.getActiveThreadCount()<fjp.getParallelism()){
//                    System.out.println("PeakCurveSmoothingUnit: fjp.getActiveThreadCount()\t"+fjp.getActiveThreadCount()+"\t"+step_pccc);
                    step_pccc *= 2;
                }
            }
        }

        assert ftemp.isEmpty():"temp storage for futures should be empty by end of loop";
        fjp.shutdown();

        try {
            fjp.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Logger.getRootLogger().info("interrupted..");
        }

        for (final PeakCluster peakCluster : resultClusters) {
            //Check if the monoistope peak of cluster has been grouped in other isotope cluster, if yes, remove the peak cluster
            if (!parameter.RemoveGroupedPeaks ||
//                    !peakCluster.MonoIsotopePeak.ChargeGrouped.contains(peakCluster.Charge)
                    !IonChargeHashSet.contains(peakCluster.MonoIsotopePeak.ChargeGrouped, peakCluster.Charge)
                    ) {
                peakCluster.Index = LCMSPeakBase.PeakClusters.size() + 1;
                peakCluster.GetConflictCorr();
                LCMSPeakBase.PeakClusters.add(peakCluster);
            }
        }
        
 /*------------------------------------------------------------------------------------------------------------*/      
//        int ii=0;
//        for(PeakCluster tmp: resultClusters) {
//        	if(ii%500==0) {
//        		XYSeries series3 = new XYSeries("xySeries3");
//                for(int aaa=0; aaa<tmp.MonoIsotopePeak.GetPeakList().size();aaa++) {
//                	series3.add(tmp.MonoIsotopePeak.GetPeakList().getXat(aaa),tmp.MonoIsotopePeak.GetPeakList().getZat(aaa));
//                }
//        		XYSeriesCollection dataset3 = new XYSeriesCollection();
//            	dataset3.addSeries(series3);
//            	JFreeChart chart3 = ChartFactory.createScatterPlot(
//            			"mz/rt",
//            			"rt",
//            			"int",
//            			dataset3, // data
//            			PlotOrientation.VERTICAL,
//            			false, // include legend
//            			false, // tooltips
//            			false // urls
//           		);
//             
//
//            	OutputStream out3 = null;
//
//                try {
//                	out3 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/"+ii+"_org.png");
//                	ChartUtilities.writeChartAsPNG(out3, chart3, 520, 520);
//                	
//            	} catch (Exception e) {
//            	    e.printStackTrace();
//            	}
//                finally {
//        	    	try {
//        	    		out3.flush();
//        	    		out3.close();
//        	    	} catch (Exception e) {
//        				e.printStackTrace();
//        			}
//                }
//        		
//                XYSeries series1 = new XYSeries("xySeries1");
//                for(int aaa=0; aaa<tmp.MonoIsotopePeak.newdata.Data.size();aaa++) {
//                	series1.add(tmp.MonoIsotopePeak.newdata.Data.get(aaa).getX(),tmp.MonoIsotopePeak.newdata.Data.get(aaa).getY());
//                }
//        		XYSeriesCollection dataset1 = new XYSeriesCollection();
//            	dataset1.addSeries(series1);
//            	JFreeChart chart1 = ChartFactory.createScatterPlot(
//            			"mz/rt",
//            			"rt",
//            			"int",
//            			dataset1, // data
//            			PlotOrientation.VERTICAL,
//            			false, // include legend
//            			false, // tooltips
//            			false // urls
//           		);
//             
//
//            	OutputStream out1 = null;
//
//                try {
//                	out1 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/"+ii+"_new.png");
//                	ChartUtilities.writeChartAsPNG(out1, chart1, 520, 520);
//                	
//            	} catch (Exception e) {
//            	    e.printStackTrace();
//            	}
//                finally {
//        	    	try {
//        	    		out1.flush();
//        	    		out1.close();
//        	    	} catch (Exception e) {
//        				e.printStackTrace();
//        			}
//                }
//                
//                XYSeries series2 = new XYSeries("xySeries2");
//                for(int aaa=0; aaa<tmp.MonoIsotopePeak.GetSmoothedList().Data.size();aaa++) {
//                	series2.add(tmp.MonoIsotopePeak.GetSmoothedList().Data.get(aaa).getX(),tmp.MonoIsotopePeak.GetSmoothedList().Data.get(aaa).getY());
//                }
//        	    	XYSeriesCollection dataset2 = new XYSeriesCollection();
//        	    	dataset2.addSeries(series2);
//        	    	JFreeChart chart2 = ChartFactory.createScatterPlot(
//        	    			"mz/rt",
//        	    			"rt",
//        	    			"int",
//        	    			dataset2, // data
//        	    			PlotOrientation.VERTICAL,
//        	    			false, // include legend
//        	    			false, // tooltips
//        	    			false // urls
//        	   		);
//        	     
//
//        	    	OutputStream out2 = null;
//
//        	        try {
//        	        	out2 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/"+ii+"smoothing.png");
//        	        	ChartUtilities.writeChartAsPNG(out2, chart2, 520, 520);
//        	        	
//        	    	} catch (Exception e) {
//        	    	    e.printStackTrace();
//        	    	}
//        	        finally {
//        		    	try {
//        		    		out2.flush();
//        		    		out2.close();
//        		    	} catch (Exception e) {
//        					e.printStackTrace();
//        				}
//        	        }
//        		
//        	}
//        	ii++;
//        }

        System.gc();
        Logger.getRootLogger().info("No of ion clusters:" + LCMSPeakBase.PeakClusters.size() + " (Memory usage:" + Math.round((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1048576) + "MB)");
    }
}
