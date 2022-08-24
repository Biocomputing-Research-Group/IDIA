package MSUmpire.MathPackage;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
//import MSUmpire.MathPackage.Regression.Equation;


public class conCC {
//	public Equation equation;
	protected XYPointCollection pointset;
//	private float SigXY = 0;
	private float SigX = 0;
	private float SigY = 0;
	private float SX2 = 0;
	private float SY2 = 0;
	private float SXY=0;
	private float numCount;
	
	public void SetData(XYPointCollection pointset) {
		this.pointset = pointset;
		this.numCount = pointset.PointCount();
		calAverage();
		calOthers();
	}
	
	public void calAverage() {
		for(int i=0; i<pointset.PointCount();i++) {
			XYData point = pointset.Data.get(i);
//			SigXY += point.getX()*point.getY();
			SigX += point.getX();
			SigY += point.getY();
		}
		SigX = SigX/numCount;
		SigY = SigY/numCount;
	}
	
	public void calOthers() {
		for(int i=0; i<pointset.PointCount(); i++) {
			XYData point = pointset.Data.get(i);
			SX2 += (point.getX()-SigX)*(point.getX()-SigX);
			SY2 += (point.getY()-SigY)*(point.getY()-SigY);
			SXY += (point.getX()-SigX)*(point.getY()-SigY);
		}
		SX2 = SX2/numCount;
		SY2 = SY2/numCount;
		SXY = SXY/numCount;
	}
	public float calCondCC() {
		float cc=-1;
		cc=2*SXY/(SX2+SY2+(SigX-SigY)*(SigX-SigY));
		return cc;
	}

}