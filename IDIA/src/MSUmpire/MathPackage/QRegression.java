package MSUmpire.MathPackage;

import MSUmpire.BaseDataStructure.XYData;
import MSUmpire.BaseDataStructure.XYPointCollection;
//import MSUmpire.MathPackage.Regression.Equation;

public class QRegression {
//	public Equation equation;
	protected XYPointCollection pointset;
//	private float SigXY = 0;
	private float SigX = 0;
	private float SigY = 0;
	private float SigX2 = 0;
//	private float SigY2 = 0;
//	private float sigX3 = 0;
//	private float sigY3 = 0;
//	private float sigX4 = 0;
//	private float sigY4 = 0;
//	private float sigX2Y = 0;
	private float SXX=0;
	private float SXY=0;
	private float SXX2=0;
	private float SX2X2=0;
	private float SX2Y=0;
	private float SSE=0;
	private float SST=0;
	private int numCount=0;
	private float bValue=0;
	private float cValue=0;
	private float aValue=0;
	
	public void setdata(XYPointCollection pointset) {
		this.pointset = pointset;
		this.numCount = pointset.PointCount();
		calAverage();
		calSomeValues();
	}
	
	public void calAverage() {
		for(int i=0; i<pointset.PointCount();i++) {
			XYData point = pointset.Data.get(i);
//			SigXY += point.getX()*point.getY();
			SigX += point.getX();
			SigY += point.getY();
			SigX2 += point.getX()*point.getX();
		}
	}
	
	public void calSomeValues() {
		for(int i=0; i<pointset.PointCount(); i++) {
			XYData point = pointset.Data.get(i);
			SXX += (point.getX()-SigX/numCount)*(point.getX()-SigX/numCount);
			SXY += (point.getX()-SigX/numCount)*(point.getY()-SigY/numCount);
			SXX2 += (point.getX()-SigX/numCount)*(point.getX()*point.getX()-SigX2/numCount);
			SX2X2 += (point.getX()*point.getX()-SigX2/numCount)*(point.getX()*point.getX()-SigX2/numCount);
			SX2Y += (point.getX()*point.getX()-SigX2/numCount)*(point.getY()-SigY/numCount);
		}
		bValue = (SXY*SX2X2-SX2Y*SXX2)/(SXX*SX2X2-SXX2*SXX2);
		cValue = (SX2Y*SXX-SXY*SXX2)/(SXX*SX2X2-SXX2*SXX2);
		aValue = SigY/numCount-bValue*SigX/numCount-cValue*SigX2/numCount;
	}
	
	public float getY(float x) {
		return aValue+bValue*x+cValue*x*x;
	}
	
	public void calSSE() {
		for(int i=0; i<pointset.PointCount(); i++) {
			SSE += (pointset.Data.get(i).getY()-getY(pointset.Data.get(i).getX()))*(pointset.Data.get(i).getY()-getY(pointset.Data.get(i).getX()));
		}
	}
	
	public void calSST() {
		for(int i=0; i<pointset.PointCount(); i++) {
			SST += (pointset.Data.get(i).getY()-SigY/numCount)*(pointset.Data.get(i).getY()-SigY/numCount);
		}
	}
	
	public float getR2() {
		float r2=0;
		calSSE();
		calSST();
		r2 = 1-SSE/SST;
		return r2;
	}
}