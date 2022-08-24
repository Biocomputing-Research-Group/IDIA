package GridMass.datamodule;


public class IndexedDataPoint {
  DataPoint datapoint;
  int index;

  public IndexedDataPoint(int index, DataPoint dp) {
    this.datapoint = dp;
    this.index = index;
  }
  public DataPoint getDatapoint() {
	  return datapoint;
  }
  public int getIndex() {
	  return index;
  }
}
