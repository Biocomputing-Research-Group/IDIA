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
package MSUmpire.MathPackage;

import MSUmpire.BaseDataStructure.XYPointCollection;

import java.io.FileOutputStream;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PearsonCorr {

    
//    public double CalcCorrV2(XYPointCollection CollectionA, XYPointCollection CollectionB, int NoPointPerInterval) {
//        SpearmansCorrelation pearsonsCorrelation = new SpearmansCorrelation();
//        
//        int num = Math.max(CollectionA.PointCount(), CollectionB.PointCount()) / 2;
//        float timeinterval = 2f / (float) NoPointPerInterval;
//        if (num < 6) {
//            return 0f;
//        }
//
//        double[] arrayA = new double[num];
//        double[] arrayB = new double[num];
//
//        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());
//
//        int i = 0;
//        float low = start;
//        float up = start + timeinterval;
//
//        for (int j = 0; j < CollectionA.PointCount(); j++) {
//            while (CollectionA.Data.get(j).getX() > up) {
//                i++;
//                low = up;
//                up = low + timeinterval;
//            }
//            if (i >= num) {
//                break;
//            }
//            if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
//                if (CollectionA.Data.get(j).getY() > arrayA[i]) {
//                    arrayA[i] = CollectionA.Data.get(j).getY();
//                }
//            }
//        }
//        i = 0;
//        low = start;
//        up = start + timeinterval;
//        for (int j = 0; j < CollectionB.PointCount(); j++) {
//            while (CollectionB.Data.get(j).getX() > up) {
//                i++;
//                low = up;
//                up = low + timeinterval;
//            }
//            if (i >= num) {
//                break;
//            }
//            if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
//                if (CollectionB.Data.get(j).getY() > arrayB[i]) {
//                    arrayB[i] = CollectionB.Data.get(j).getY();
//                }
//            }
//        }
//
//        if(arrayA[0]==0f){
//            arrayA[0]=arrayA[1];
//        }
//        if(arrayB[0]==0f){
//            arrayB[0]=arrayB[1];
//        }
//        for (int idx = 1; idx < num - 1; idx++) {
//            if (arrayA[idx] == 0f) {
//                arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
//            }
//            if (arrayB[idx] == 0f) {
//                arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
//            }
//        }
//        
//        if(arrayA[num - 1]==0f){
//            arrayA[num - 1]=arrayA[num - 2];
//        }
//        if(arrayB[num - 1]==0f){
//            arrayB[num - 1]=arrayB[num - 2];
//        }
//        double R2 =pearsonsCorrelation.correlation(arrayA, arrayB); 
//        return R2;
//    }
    
    public float CalcCorr(XYPointCollection CollectionA, XYPointCollection CollectionB, int NoPointPerInterval) {
        Regression regression = new Regression();

        int num = Math.max(CollectionA.PointCount(), CollectionB.PointCount()) / 2;
        float timeinterval = 2f / (float) NoPointPerInterval;
        if (num < 6) {
            return 0f;
        }

        float[] arrayA = new float[num];
        float[] arrayB = new float[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        int i = 0;
        float low = start;
        float up = start + timeinterval;

        for (int j = 0; j < CollectionA.PointCount(); j++) {
            while (CollectionA.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                if (CollectionA.Data.get(j).getY() > arrayA[i]) {
                    arrayA[i] = CollectionA.Data.get(j).getY();
                }
            }
        }
        i = 0;
        low = start;
        up = start + timeinterval;
        for (int j = 0; j < CollectionB.PointCount(); j++) {
            while (CollectionB.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                if (CollectionB.Data.get(j).getY() > arrayB[i]) {
                    arrayB[i] = CollectionB.Data.get(j).getY();
                }
            }
        }

        for (int idx = 1; idx < num - 1; idx++) {
            if (arrayA[idx] == 0f) {
                arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
            }
            if (arrayB[idx] == 0f) {
                arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
            }
        }

        XYPointCollection pointset = new XYPointCollection();
        for (int idx = 0; idx < num; idx++) {
            if (arrayA[idx] > 0 && arrayB[idx] > 0) {
                pointset.AddPoint(arrayA[idx], arrayB[idx]);
            }
        }

        float R2 = 0f;
        if (pointset.PointCount() > 5) {
            regression.SetData(pointset);
            if (regression.equation.Mvalue > 0) {
                R2 = regression.GetR2();
            }
        }
        return R2;
    }
    public float CalcCorrOut(XYPointCollection CollectionA, XYPointCollection CollectionB, int NoPointPerInterval, String filepath) {
        Regression regression = new Regression();
//    	PearsonsCorrelation regression = new PearsonsCorrelation();
//    	OverlappingCorrelation regression = new OverlappingCorrelation();
        int num = Math.max(CollectionA.PointCount(), CollectionB.PointCount()) / 2;
        float timeinterval = 2f / (float) NoPointPerInterval;
        if (num < 6) {
            return 0f;
        }

        float[] arrayA = new float[num];
        float[] arrayB = new float[num];
//        double[] arrayA = new double[num];
//        double[] arrayB = new double[num];

        float start = Math.max(CollectionA.Data.get(0).getX(), CollectionB.Data.get(0).getX());

        int i = 0;
        float low = start;
        float up = start + timeinterval;

        for (int j = 0; j < CollectionA.PointCount(); j++) {
            while (CollectionA.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionA.Data.get(j).getX() >= low && CollectionA.Data.get(j).getX() < up) {
                if (CollectionA.Data.get(j).getY() > arrayA[i]) {
                    arrayA[i] = CollectionA.Data.get(j).getY();
//                	arrayA[i] = Double.parseDouble(String.valueOf(CollectionA.Data.get(j).getY()));
                }
            }
        }
        i = 0;
        low = start;
        up = start + timeinterval;
        for (int j = 0; j < CollectionB.PointCount(); j++) {
            while (CollectionB.Data.get(j).getX() > up) {
                i++;
                low = up;
                up = low + timeinterval;
            }
            if (i >= num) {
                break;
            }
            if (CollectionB.Data.get(j).getX() >= low && CollectionB.Data.get(j).getX() < up) {
                if (CollectionB.Data.get(j).getY() > arrayB[i]) {
                    arrayB[i] = CollectionB.Data.get(j).getY();
//                	arrayB[i] = Double.parseDouble(String.valueOf(CollectionB.Data.get(j).getY()));
                }
            }
        }

        for (int idx = 1; idx < num - 1; idx++) {
            if (arrayA[idx] == 0f) {
                arrayA[idx] = (arrayA[idx - 1] + arrayA[idx + 1]) / 2;
            }
            if (arrayB[idx] == 0f) {
                arrayB[idx] = (arrayB[idx - 1] + arrayB[idx + 1]) / 2;
            }
        }

        XYPointCollection pointset = new XYPointCollection();
        for (int idx = 0; idx < num; idx++) {
            if (arrayA[idx] > 0 && arrayB[idx] > 0) {
                pointset.AddPoint(arrayA[idx], arrayB[idx]);
            }
        }
        

        float R2 = 0f;
        if (pointset.PointCount() > 5) {
            regression.SetData(pointset);
            if (regression.equation.Mvalue > 0) {
                R2 = regression.GetR2();
//            }
//        if(arrayA.length>5) {
//        	regression.setData(pointset);
//        	R2=(float) regression.getOverlapping();        	
        	System.out.println(arrayA.length);
        	System.out.println(arrayB.length);
            
            
            XYSeries series1 = new XYSeries("xySeries1");
            for (int iii=0; iii<arrayA.length; iii++) {
            	series1.add((float) iii,(float) arrayA[iii]);
            }
            XYSeriesCollection dataset1 = new XYSeriesCollection();
        	dataset1.addSeries(series1);
        	JFreeChart chart1 = ChartFactory.createScatterPlot(
        			"ms1",
        			"rt",
        			"int",
        			dataset1, // data
        			PlotOrientation.VERTICAL,
        			false, // include legend
        			false, // tooltips
        			false // urls
       		);
         

            try {
            	FileOutputStream out1 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/"+filepath+"ms1.png");
            	ChartUtilities.writeChartAsPNG(out1, chart1, 520, 520);
        		out1.flush();
        		out1.close();
        	} catch (Exception e) {
        	    e.printStackTrace();
        	}
            XYSeries series2 = new XYSeries("xySeries2");
            for (int iii=0; iii<arrayB.length; iii++) {
            	series2.add((float) iii,(float) arrayB[iii]);
            }
            XYSeriesCollection dataset2 = new XYSeriesCollection();
        	dataset2.addSeries(series2);
        	JFreeChart chart2 = ChartFactory.createScatterPlot(
        			"ms2",
        			"rt",
        			"int",
        			dataset2, // data
        			PlotOrientation.VERTICAL,
        			false, // include legend
        			false, // tooltips
        			false // urls
       		);
            try {
            	FileOutputStream out2 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/"+filepath+"ms2.png");
            	ChartUtilities.writeChartAsPNG(out2, chart2, 520, 520);
            	out2.flush();
        		out2.close();
            	
        	} catch (Exception e) {
        	    e.printStackTrace();
        	}

        	System.out.println(R2);
//        	System.exit(1);
        }
        }
        return R2;
    }
}
