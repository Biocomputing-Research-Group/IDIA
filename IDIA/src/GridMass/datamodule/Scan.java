package GridMass.datamodule;

import com.google.common.collect.Range;

import weka.core.parser.java_cup.internal_error;

import javax.annotation.Nonnull;
import java.util.Vector;


public class Scan{

  private int scanNumber = -1;
  private int msLevel = -1;
  private DataPoint dataPoints[] = null;
  private double retentionTime = -1;
  private Range<Double> mzRange = null;
  private DataPoint basePeak = null;    // max intensity datapoint
  private double totalIonCurrent = -1;  // sum all intensity
  private boolean centroided = true;
  private String scanType = null;
  private int peakCount = -1;
  private double startMZ = -1;
  private double endMZ = -1;
//  private int precusorScanNum = -1;
  private double precusorMZ = -1;
  private String activationMethod = null;
//  private int windowWideness = 25;
  
//  public Scan(Scan sc) {
//    this(sc.getScanNumber(), sc.getMSLevel(), sc.getRetentionTime(),
//            sc.getDataPoints());
//  }
  
  public Scan() {
	  this.scanNumber = -1;
	  this.msLevel = -1;
	  this.dataPoints = null;
	  this.retentionTime = -1;
	  this.mzRange = null;
	  this.basePeak = null;    // max intensity datapoint
	  this.totalIonCurrent = -1;  // sum all intensity
	  this.centroided = true;
	  this.scanType = null;
	  this.peakCount = -1;
	  this.startMZ = -1;
	  this.endMZ = -1;
//	  this.precusorScanNum = -1;
	  this.activationMethod = null;
//	  this.windowWideness = 25;
  }


  public @Nonnull DataPoint[] getDataPoints() {
    return dataPoints;
  }

  public @Nonnull DataPoint[] getDataPointsByMass(@Nonnull Range<Double> mzRange) {

    int startIndex, endIndex;
    for (startIndex = 0; startIndex < dataPoints.length; startIndex++) {
      if (dataPoints[startIndex].getMZ() >= mzRange.lowerEndpoint())
        break;
    }

    for (endIndex = startIndex; endIndex < dataPoints.length; endIndex++) {
      if (dataPoints[endIndex].getMZ() > mzRange.upperEndpoint())
        break;
    }

    DataPoint pointsWithinRange[] = new DataPoint[endIndex - startIndex];
    System.arraycopy(dataPoints, startIndex, pointsWithinRange, 0, endIndex - startIndex);

    return pointsWithinRange;
  }

  public @Nonnull DataPoint[] getDataPointsOverIntensity(double intensity) {
    int index;
    Vector<DataPoint> points = new Vector<DataPoint>();

    for (index = 0; index < dataPoints.length; index++) {
      if (dataPoints[index].getIntensity() >= intensity)
        points.add(dataPoints[index]);
    }

    DataPoint pointsOverIntensity[] = points.toArray(new DataPoint[0]);

    return pointsOverIntensity;
  }

  public void setDataPoints(DataPoint[] dataPoints) {

    this.dataPoints = dataPoints;
    mzRange = Range.singleton(0.0);
    basePeak = null;
    totalIonCurrent = 0;

    if (dataPoints.length > 0) {

      basePeak = dataPoints[0];
      mzRange = Range.singleton(dataPoints[0].getMZ());

      for (DataPoint dp : dataPoints) {

        if (dp.getIntensity() > basePeak.getIntensity())
          basePeak = dp;

        mzRange = mzRange.span(Range.singleton(dp.getMZ()));
        totalIonCurrent += dp.getIntensity();

      }
    }
  }

  public int getNumberOfDataPoints() {
    return dataPoints.length;
  }

  public int getScanNumber() {
    return scanNumber;
  }

  public int getMSLevel() {
    return msLevel;
  }

  public double getRetentionTime() {
    return retentionTime;
  }

  public @Nonnull Range<Double> getDataPointMZRange() {
    return mzRange;
  }

  public DataPoint getHighestDataPoint() {
    return basePeak;
  }

  public double getTIC() {
    return totalIonCurrent;
  }
  
  public String getscanType() {
	  return this.scanType;
  }
  
  public int getpeakCount() {
	  return this.peakCount;
  }
  
  public double getPrecusorMZ() {
	  return this.precusorMZ;
  }
  
  public void setScanNumber(int scanNumber) {
	    this.scanNumber = scanNumber;
  }
  
  public void setMSLevel(int msLevel) {
	    this.msLevel = msLevel;
  }
  
  public void setRetentionTime(double retentionTime) {
	  this.retentionTime = retentionTime;
  }
  
  public void setMzRange(Range mzRnage) {
	  this.mzRange = mzRnage;
  }
  
  public void setbasePeak(DataPoint basePeak) {
	  this.basePeak = basePeak;
  }
  
  public void settotalIonCurrent(double totalIonCurrent) {
	  this.totalIonCurrent = totalIonCurrent;
  }
  
  public void setcentroided(boolean centroided) {
	  this.centroided = centroided;
  }
  
  public void setscanType(String scanType) {
	  this.scanType = scanType;
  }
  
  public void setpeakCount(int peakCount) {
	  this.peakCount = peakCount;
  }
  public void setStartMZ(double startMZ) {
	  this.startMZ = startMZ;
  }
  
  public void setEndMZ(double endMZ) {
	  this.endMZ = endMZ;
  }
  
  public void setPrecusorMZ(double precursorMZ) {
	  this.precusorMZ = precursorMZ;
  }
//  public void setprecusorScanNum(int precusorScanNum) {
//	  this.precusorScanNum = precusorScanNum;
//  }
  
  public void setactivationMethod(String activationMethod) {
	  this.activationMethod = activationMethod;
  }
  
//  public void setwindowWideness(int windowWideness) {
//	  this.windowWideness = windowWideness;
//  }
 
}
