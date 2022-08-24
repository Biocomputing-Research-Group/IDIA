package GridMass.datamodule;

public class DataPoint{

  private double mz, intensity,rt;


  public DataPoint(DataPoint dp) {
    this.mz = dp.getMZ();
    this.intensity = dp.getIntensity();
    this.rt = dp.getRT();
  }

  public DataPoint(double mz, double intensity, double rt) {
    this.mz = mz;
    this.intensity = intensity;
    this.rt = rt;
  }

  public double getIntensity() {
    return intensity;
  }

  public double getMZ() {
    return mz;
  }
  
  public double getRT() {
	  return rt;
  }
  
  public boolean equals(Object obj) {
    if (!(obj instanceof DataPoint))
      return false;
    DataPoint dp = (DataPoint) obj;
    return (this.mz == dp.getMZ()) && (this.intensity == dp.getIntensity()) && (this.rt == dp.getRT());
  }

  public int hashCode() {
    return (int) (this.mz + this.intensity+this.rt);
  }

}
