package MSUmpire.MathPackage;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.TreeMap;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.eclipse.collections.impl.list.mutable.primitive.FloatArrayList;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import MSUmpire.BaseDataStructure.XYData;
import umontreal.iro.lecuyer.functionfit.SmoothingCubicSpline;
import weka.core.parser.java_cup.internal_error;

public class test{
	
	public static void main(String[] args) {
//		int[] aa=new int[2];
//		aa[0]=0;
//		aa[1]=1;
//		TreeMap<Integer, Integer> bbMap = new TreeMap<>();
//		for(int i =0; i<aa.length;i++) {
//			bbMap.put(aa[i], aa[i]);
//		}
//		System.out.println(bbMap.size());
//		FloatArrayList xyzdata=new FloatArrayList();
//		System.out.println(xyzdata.size());
		
		
		Random uniform = new Random();
		NormalDistribution normal = new NormalDistribution();
		
		
		double[] unlist = new double[20];
		double[] norlist = new double[20];
		
		for(int i=0; i<20; i++) {
//			unlist[i]=uniform.nextDouble();
			unlist[i]=uniform.nextGaussian();
//			norlist[i]=normal.density(unlist[i])+uniform.nextDouble()*0.1;
		}
		Arrays.sort(unlist);
		for(int i=0; i<20; i++) {
//			unlist[i]=uniform.nextDouble();
//			unlist[i]=uniform.nextGaussian();
			norlist[i]=normal.density(unlist[i])+uniform.nextDouble()*0.1;
		}
		try {
//			loessinterpolator loess = new loessinterpolator(0.5, 5);
//			double[] result = loess.smooth(unlist, norlist);
//			polynomialsplinefunction f = loess.interpolate(unlist, norlist);
			SmoothingCubicSpline fit=new SmoothingCubicSpline(unlist, norlist, 0.7);
//			XYSeries series2 = new XYSeries("xySeries2");
//	        for(int i=0; i<20; i++) {
//	        	series2.add(unlist[i],result[i]);
//	        }
//	        
//	        XYSeriesCollection dataset2 = new XYSeriesCollection();
//	    	dataset2.addSeries(series2);
//	    	JFreeChart chart2 = ChartFactory.createScatterPlot(
//	    			"mz/rt",
//	    			"rt",
//	    			"int",
//	    			dataset2, // data
//	    			PlotOrientation.VERTICAL,
//	    			false, // include legend
//	    			false, // tooltips
//	    			false // urls
//	   		);
//	     
//
//	    	OutputStream out2 = null;
//
//	        try {
//	        	out2 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/test1.png");
//	        	ChartUtilities.writeChartAsPNG(out2, chart2, 520, 520);
//	        	
//	    	} catch (Exception e) {
//	    	    e.printStackTrace();
//	    	}
//	        finally {
//		    	try {
//		    		out2.flush();
//		    		out2.close();
//		    	} catch (Exception e) {
//					e.printStackTrace();
//				}
//			}
	        
	        XYSeries series3 = new XYSeries("xySeries3");
	        double range = (unlist[19]-unlist[0])/40;
	        for(int i=0; i<40; i++) {
	        	series3.add(unlist[0]+i*range,fit.evaluate(unlist[0]+i*range));
	        }

	        XYSeriesCollection dataset3 = new XYSeriesCollection();
	    	dataset3.addSeries(series3);
	    	JFreeChart chart3 = ChartFactory.createScatterPlot(
	    			"mz/rt",
	    			"rt",
	    			"int",
	    			dataset3, // data
	    			PlotOrientation.VERTICAL,
	    			false, // include legend
	    			false, // tooltips
	    			false // urls
	   		);
	     

	    	OutputStream out3 = null;

	        try {
	        	out3 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/test2.png");
	        	ChartUtilities.writeChartAsPNG(out3, chart3, 520, 520);
	        	
	    	} catch (Exception e) {
	    	    e.printStackTrace();
	    	}
	        finally {
		    	try {
		    		out3.flush();
		    		out3.close();
		    	} catch (Exception e) {
					e.printStackTrace();
				}
			}
	        
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		 XYSeries series1 = new XYSeries("xySeries1");
		 for(int i=0; i<20; i++) {
			 series1.add(unlist[i],norlist[i]);
		 }
	        
	        
	        XYSeriesCollection dataset1 = new XYSeriesCollection();
	    	dataset1.addSeries(series1);
	    	JFreeChart chart1 = ChartFactory.createScatterPlot(
	    			"mz/rt",
	    			"rt",
	    			"int",
	    			dataset1, // data
	    			PlotOrientation.VERTICAL,
	    			false, // include legend
	    			false, // tooltips
	    			false // urls
	   		);
	     

	    	OutputStream out1 = null;

	        try {
	        	out1 = new FileOutputStream("/home/UNT/jl0734/Desktop/DIA/myeclipseproject/result/test.png");
	        	ChartUtilities.writeChartAsPNG(out1, chart1, 520, 520);
	        	
	    	} catch (Exception e) {
	    	    e.printStackTrace();
	    	}
	        finally {
		    	try {
		    		out1.flush();
		    		out1.close();
		    	} catch (Exception e) {
					e.printStackTrace();
				}
	        }
	        
	 
//	FloatArrayList tmpArrayList = new FloatArrayList();
//	tmpArrayList.add(1.0f);
//	tmpArrayList.add(10.0f);
//	tmpArrayList.add(5.0f);
//	for(int aa=0; aa<tmpArrayList.size(); aa++) {
//		System.out.println(tmpArrayList.get(aa));
//	}
//	tmpArrayList.trimToSize();
//	tmpArrayList.sortThis().reverseThis();
//	for(int aa=0; aa<tmpArrayList.size(); aa++) {
//		System.out.println(tmpArrayList.get(aa));
//	}
//		
	}
	
}