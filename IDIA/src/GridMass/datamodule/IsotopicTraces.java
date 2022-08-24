package GridMass.datamodule;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import javax.annotation.Nonnull;
import com.google.common.collect.Range;
import com.google.common.primitives.Ints;

import weka.core.parser.java_cup.internal_error;

public class IsotopicTraces{

  private Hashtable<Integer, DataPoint> dataPointsMap;

  private double mz, rt, height, area;
//  private Double fwhm = null, tf = null, af = null;

  // scan of max intensity
  private int representativeScan = -1;

  private Range<Double> rawDataPointsIntensityRange, rawDataPointsMZRange, rawDataPointsRTRange;
  private double mzSum = 0;
  private int mzN = 0;

//  private final int scanNumbers[];
  private double targetmz;

  public IsotopicTraces( ) {
//    this.scanNumbers = scanNumbers;

    dataPointsMap = new Hashtable<Integer, DataPoint>();
  }

  public void addMzPeak(int scanNumber, DataPoint mzValue) {
    dataPointsMap.put(scanNumber, mzValue);
    mzSum += mzValue.getMZ();
    mzN++;
    mz = mzSum / mzN;

  }

  public int getRepresentativeScan() {
	  return representativeScan;
  }
  public DataPoint getDataPoint(int scanNumber) {
    return dataPointsMap.get(scanNumber);
  }

//  public int getDataPointSize() {
//	  return dataPointsMap.size();
//  }
  public Hashtable<Integer, DataPoint> getDataPoint(){
	  return dataPointsMap;
  }
  public double getMZ() {
    return mz;
  }

  public double getArea() {
    return area;
  }

  public double getHeight() {
    return height;
  }

  public double getRT() {
    return rt;
  }
  
  public double getStartRT() {
	  return this.rawDataPointsRTRange.lowerEndpoint();
  }
  
  public double getEndRT() {
	  return this.rawDataPointsRTRange.upperEndpoint();
  }
  
  public double getMaxMZ() {
	  return this.rawDataPointsMZRange.upperEndpoint();
  }
  
  public double getMinMZ() {
	  return this.rawDataPointsMZRange.lowerEndpoint();
  }
  
  public @Nonnull Range<Double> getRawDataPointsIntensityRange() {
    return rawDataPointsIntensityRange;
  }

  public @Nonnull Range<Double> getRawDataPointsMZRange() {
    return rawDataPointsMZRange;
  }

  public @Nonnull Range<Double> getRawDataPointsRTRange() {
    return rawDataPointsRTRange;
  }

  
  public void finishChromatogram(HashMap scantort) {

    int allScanNumbers[] = Ints.toArray(dataPointsMap.keySet());
    Arrays.sort(allScanNumbers);

    // Calculate median m/z
    double allMzValues[] = new double[allScanNumbers.length];
    for (int i = 0; i < allScanNumbers.length; i++) {
      allMzValues[i] = dataPointsMap.get(allScanNumbers[i]).getMZ();
    }
    mz = calcQuantile(allMzValues, 0.5f);

    // Update raw data point ranges, height, rt and representative scan
    height = Double.MIN_VALUE;
    for (int i = 0; i < allScanNumbers.length; i++) {

      DataPoint mzPeak = dataPointsMap.get(allScanNumbers[i]);

      // Replace the MzPeak instance with an instance of SimpleDataPoint,
      // to reduce the memory usage. After we finish this Chromatogram, we
      // don't need the additional data provided by the MzPeak

      dataPointsMap.put(allScanNumbers[i], mzPeak);

      if (i == 0) {
        rawDataPointsIntensityRange = Range.singleton(mzPeak.getIntensity());
        rawDataPointsMZRange = Range.singleton(mzPeak.getMZ());
      } else {
        rawDataPointsIntensityRange =
            rawDataPointsIntensityRange.span(Range.singleton(mzPeak.getIntensity()));
        rawDataPointsMZRange = rawDataPointsMZRange.span(Range.singleton(mzPeak.getMZ()));
      }

      if (height < mzPeak.getIntensity()) {
        height = mzPeak.getIntensity();
//        rt = dataFile.getScan(allScanNumbers[i]).getRetentionTime();
        rt = (double) scantort.get(allScanNumbers[i]);
        representativeScan = allScanNumbers[i];
      }
    }

    // Update area
    area = 0;
    for (int i = 1; i < allScanNumbers.length; i++) {
      // For area calculation, we use retention time in seconds
      double previousRT = (double) scantort.get(allScanNumbers[i-1]);
      double currentRT = (double) scantort.get(allScanNumbers[i]);
//      double previousRT = dataFile.getScan(allScanNumbers[i - 1]).getRetentionTime() * 60d;
//      double currentRT = dataFile.getScan(allScanNumbers[i]).getRetentionTime() * 60d;
      double previousHeight = dataPointsMap.get(allScanNumbers[i - 1]).getIntensity();
      double currentHeight = dataPointsMap.get(allScanNumbers[i]).getIntensity();
      area += (currentRT - previousRT) * (currentHeight + previousHeight) / 2;
    }



    rawDataPointsRTRange = null;

    for (int scanNum : allScanNumbers) {
//      double scanRt = dataFile.getScan(scanNum).getRetentionTime();
    	double scanRt = (double) scantort.get(scanNum);
      DataPoint dp = getDataPoint(scanNum);

      if ((dp == null) || (dp.getIntensity() == 0.0))
        continue;

      if (rawDataPointsRTRange == null)
        rawDataPointsRTRange = Range.singleton(scanRt);
      else
        rawDataPointsRTRange = rawDataPointsRTRange.span(Range.singleton(scanRt));
    }

  }
  
  public static double calcQuantile(double[] values, double q) {

	    if (values.length == 0)
	      return 0;

	    if (values.length == 1)
	      return values[0];

	    if (q > 1)
	      q = 1;

	    if (q < 0)
	      q = 0;

	    double[] vals = values.clone();

	    Arrays.sort(vals);

	    int ind1 = (int) Math.floor((vals.length - 1) * q);
	    int ind2 = (int) Math.ceil((vals.length - 1) * q);

	    return (vals[ind1] + vals[ind2]) / 2;

	  }
  
  public void calTargeMZ() {
	  double sumMZ = 0;
	  double sumIntensity = 0;
	  for(Integer keys: this.dataPointsMap.keySet()) {
		  sumMZ += dataPointsMap.get(keys).getIntensity() * dataPointsMap.get(keys).getIntensity() * dataPointsMap.get(keys).getMZ();
		  sumIntensity += dataPointsMap.get(keys).getIntensity() * dataPointsMap.get(keys).getIntensity();
	  }
	  if (sumIntensity == 0) {
		  this.targetmz = 0;
	  }
	  else {
		  this.targetmz = sumMZ/sumIntensity;
	  }
  }
  
  public double getTargetMZ() {
	  return this.targetmz;
  }
  
}
