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

import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.MathPackage.PearsonCorr;
import MSUmpire.PeakDataStructure.PeakCurve;
import cern.colt.list.DoubleArrayList;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.eclipse.collections.impl.list.mutable.primitive.FloatArrayList;

/**
 * Calculate peak profile peak correlation given two peak curves
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurveCorrCalc {

    public static float CalPeakCorr_Overlap(PeakCurve peakA, PeakCurve peakB, int Astart, int Aend, int Bstart, int Bend, int NoPeakPerMin) throws IOException {
        return CalPeakCorr_Overlap(peakA, peakB, Astart, Aend, Bstart, Bend, false, NoPeakPerMin);
    }

    public static float CalPeakCorr(PeakCurve peakA, PeakCurve peakB, int NoPointPerMin) throws IOException {        
        PearsonCorr corr = new PearsonCorr();
        float startRT = Math.max(peakA.StartRT(), peakB.StartRT());
        float endRT = Math.min(peakA.EndRT(), peakB.EndRT());
        XYPointCollection PeakACollection = peakA.GetSmoothPeakCollection(startRT, endRT);
        XYPointCollection PeakBCollection = peakB.GetSmoothPeakCollection(startRT, endRT);
        float corre = 0f;
        
        //double corre2 = 0f;
        if (PeakACollection.Data.size() > 0 && PeakBCollection.Data.size() > 0) {
            corre = corr.CalcCorr(PeakACollection, PeakBCollection, NoPointPerMin);   
        }
        return corre;
    }
    
    public static float CalPeakCorrOrg(PeakCurve peakA, PeakCurve peakB, int NoPointPerMin) throws IOException{
//    	Regression regression = new Regression();
    	PearsonsCorrelation regression = new PearsonsCorrelation();
    	float corre = -1f;
    	float startRT = Math.max(peakA.GetPeakList().getXat(0), peakB.GetPeakList().getXat(0));
    	float endRT = Math.min(peakA.GetPeakList().getXat(peakA.GetPeakList().size()-1), peakB.GetPeakList().getXat(peakB.GetPeakList().size()-1));
    	DoubleArrayList listA = new DoubleArrayList();
    	DoubleArrayList listB = new DoubleArrayList();
    	
    	int startA = 0;
    	int startB =0;
    	for(int i =0; i<peakA.GetPeakList().size(); i++) {
    		if(peakA.GetPeakList().getXat(i)>startRT ) {
    			break;
    		}
    		startA=i;
    	}
    	for(int i=0; i<peakB.GetPeakList().size(); i++) {
    		if(peakB.GetPeakList().getXat(i)>startRT) {
    			break;
    		}
    		startB=i;
    	}
    	
    	while(true) {
    		if (startA>=peakA.GetPeakList().size() || startB>=peakB.GetPeakList().size() ) {
    			break;
    		}
    		if(peakB.GetPeakList().getXat(startB)>endRT && peakA.GetPeakList().getXat(startA)>endRT) {
    			break;
    		}
    		if(peakB.GetPeakList().getXat(startB)>endRT && peakA.GetPeakList().getXat(startA)>endRT) {
    			break;
    		}
    		if(Math.abs(peakA.GetPeakList().getXat(startA)-peakB.GetPeakList().getXat(startB))<0.02) {
    			listA.add(Double.parseDouble(String.valueOf(peakA.GetPeakList().getZat(startA))));
    			listB.add(Double.parseDouble(String.valueOf(peakB.GetPeakList().getZat(startB))));
    			startA++;
    			startB++;
    		}
    		else {
    			if(Math.abs(peakA.GetPeakList().getXat(startA)-startRT)>Math.abs(peakB.GetPeakList().getXat(startB)-startRT)) {
    				int tmpIndex=0;
    				if(startA-1>0) {
    					tmpIndex=startA-1;
    				}
    				listA.add((Double.parseDouble(String.valueOf(peakA.GetPeakList().getZat(startA)))+Double.parseDouble(String.valueOf(peakA.GetPeakList().getZat(tmpIndex))))/2);
    				listB.add(Double.parseDouble(String.valueOf(peakB.GetPeakList().getZat(startB))));
    				startB++;
    			}
    			else {
    				int tmpIndex=0;
    				if(startB-1>0) {
    					tmpIndex=startB-1;
    				}
    				listA.add(Double.parseDouble(String.valueOf(peakA.GetPeakList().getZat(startA))));
    				listB.add((Double.parseDouble(String.valueOf(peakB.GetPeakList().getZat(startB)))+Double.parseDouble(String.valueOf(peakB.GetPeakList().getZat(tmpIndex))))/2);
    				startA++;
    			}
    		}
    	}
    	
    	listA.trimToSize();
    	listB.trimToSize();
        double[] arrayA = new double[listA.size()];
        double[] arrayB = new double[listB.size()];
        for(int i=0; i<listA.size(); i++) {
        	arrayA[i]=listA.get(i);
        }
        for(int i=0; i<listB.size();i++) {
        	arrayB[i]=listB.get(i);
        }
        if(arrayA.length != arrayB.length) {
        	System.out.println("not equal");
        	System.exit(1);
        }
        if(arrayA.length>4) {
        	corre = (float) regression.correlation(arrayA, arrayB);
        }
    	return corre;
    }
    
    public static float CalPeakCorr_Overlap(PeakCurve peakA, PeakCurve peakB, int Astart, int Aend, int Bstart, int Bend, boolean output, int NoPeakPerMin) throws IOException {
        PearsonCorr corr = new PearsonCorr();
        float startRT = Math.max(peakA.GetPeakRegionList().get(Astart).getX(), peakB.GetPeakRegionList().get(Bstart).getX());
        float endRT = Math.min(peakA.GetPeakRegionList().get(Aend).getZ(), peakB.GetPeakRegionList().get(Bend).getZ());
        XYPointCollection PeakACollection = peakA.GetSmoothPeakCollection(startRT, endRT);
        XYPointCollection PeakBCollection = peakB.GetSmoothPeakCollection(startRT, endRT);
        float corre = 0f;
        if (PeakACollection.Data.size() > 0 && PeakBCollection.Data.size() > 0) {
            corre = corr.CalcCorr(PeakACollection, PeakBCollection, NoPeakPerMin);
            if (output) {
                FileWriter writer = new FileWriter("PeakA.csv");
                for (int i = 0; i < PeakACollection.PointCount(); i++) {
                    writer.write(PeakACollection.Data.get(i).getX() + "," + PeakACollection.Data.get(i).getY() + "\n");
                }
                writer.close();
                writer = new FileWriter("PeakB.csv");
                for (int i = 0; i < PeakBCollection.PointCount(); i++) {
                    writer.write(PeakBCollection.Data.get(i).getX() + "," + PeakBCollection.Data.get(i).getY() + "\n");
                }
                writer.close();
            }
        }
        return corre;
    }
}
