package MSUmpire.MathPackage;

import java.util.Arrays;
import java.util.HashMap;

import org.eclipse.collections.impl.list.mutable.primitive.FloatArrayList;
import MSUmpire.BaseDataStructure.XYPointCollection;
import weka.core.parser.java_cup.internal_error;

public class GaussianFilter{
	private XYPointCollection pointSet;
	private FloatArrayList RTList;
	private FloatArrayList intyList;
	private XYPointCollection resultData ;
	private float curveStart = 0.0f;
	private float curveEnd = 0.0f;
	private float halfwithd=0f;
	private float sigma = 0f;
	private int size=0;
	private String type = "";

	public  GaussianFilter(XYPointCollection pointset,float curveStart, float curveEnd, String type) {
		this.pointSet = pointset;	
		this.curveStart = curveStart;
		this.curveEnd = curveEnd;
		this.resultData = new XYPointCollection();
		this.size=pointSet.Data.size();
		this.type = type;

		
	}
	

	public void gausfilter() {
		XYPointCollection extendPoinstSet = new XYPointCollection();
		FloatArrayList tmpWindowRT ;
		FloatArrayList tmpWindowInty ;
		float maxValue=0.0f;
		float halfVlaue=0.0f;
		float startposition=Float.POSITIVE_INFINITY;
		float endposition=0.0f;
		int windowStart = 0;
		int windowEnd = 0;
		float tmpCenter = 0f;
		float diffRT=0f;
		

		for(int i = 0; i<pointSet.Data.size(); i++) {
			if(maxValue<pointSet.Data.get(i).getY()) {
				maxValue=pointSet.Data.get(i).getY();
			}
		}
//		halfVlaue = 1*maxValue/2;
		halfVlaue = 3*maxValue/4;
//		halfVlaue = 2*maxValue/3;
		for(int i =0; i<pointSet.Data.size();i++) {
			if(pointSet.Data.get(i).getY()>halfVlaue) {
				if(pointSet.Data.get(i).getX()<startposition) {
					startposition=pointSet.Data.get(i).getX();
				}
				if(pointSet.Data.get(i).getX()>endposition) {
					endposition=pointSet.Data.get(i).getX();
				}
			}	
		}
//		float tmphalfwidth = (endposition - startposition);
//		this.sigma = tmphalfwidth/2.57f;
//		this.halfwithd = 3*sigma;
		
		float tmpwidth = 1*(endposition - startposition);
		this.sigma = tmpwidth/2.355f;
		this.halfwithd = 6*sigma;
		
		
//		this.halfwithd = (endposition - startposition)/2f;
//		this.sigma = 2*this.halfwithd/2.57f;
			
		tmpWindowRT = new FloatArrayList();
		tmpWindowInty = new FloatArrayList();
		for (int i = 0; i < pointSet.Data.size(); i++) {
			if(pointSet.Data.get(i).getX()<=(pointSet.Data.get(0).getX()+halfwithd)) {
				windowStart = i;
				tmpWindowRT.add(pointSet.Data.get(i).getX());
				tmpWindowInty.add(pointSet.Data.get(i).getY());
			}
			else {
				break;
			}
		}
		float oldValue=0f;
		for (int i = tmpWindowRT.size()-1; i >= 0; i--) {
			if(i==tmpWindowRT.size()-1) {
				diffRT = pointSet.Data.get(1).getX()-pointSet.Data.get(0).getX();				
				oldValue = tmpWindowRT.get(i);
				tmpWindowRT.set(i, pointSet.Data.get(0).getX()-diffRT);
			}
			else {
				diffRT = oldValue - tmpWindowRT.get(i);			
				oldValue = tmpWindowRT.get(i);
				tmpWindowRT.set(i, tmpWindowRT.get(i+1)-diffRT);
			}
			
		}
		for (int i = 0; i <tmpWindowRT.size(); i++) {
			extendPoinstSet.AddPoint(tmpWindowRT.get(i),tmpWindowInty.get(i));
//			extendPoinstSet.AddPoint(tmpWindowRT.get(i),0);
		}
		tmpWindowRT = null;
		tmpWindowInty = null;


		for(int i=pointSet.Data.size()-1; i>=0; i--) {
			if(pointSet.Data.get(i).getX()>=(pointSet.Data.get(pointSet.Data.size()-1).getX()-halfwithd)) {

				windowEnd = i;
			}
			else {
				break;
			}
		}
		diffRT=pointSet.Data.get(size-1).getX()-pointSet.Data.get(size-2).getX();
		for(int i=windowEnd, j =0; i<size; i++,j++) {
			pointSet.AddPoint(pointSet.Data.get(size-1+j).getX()+diffRT,pointSet.Data.get(i).getY());
//			pointSet.AddPoint(pointSet.Data.get(size-1+j).getX()+diffRT,0);
			diffRT =  pointSet.Data.get(i+1).getX()- pointSet.Data.get(i).getX();
		}
		for (int i = 0; i < pointSet.Data.size(); i++) {
			extendPoinstSet.AddPoint(pointSet.Data.get(i).getX(),pointSet.Data.get(i).getY());
		}
		extendPoinstSet.Data.Finalize();
		

		for (int i = windowStart+1; i <size+windowStart; i++) {
			tmpWindowRT = new FloatArrayList();
			tmpWindowInty = new FloatArrayList();
			int tmpstartindex = extendPoinstSet.GetClosetIndexOfX(extendPoinstSet.Data.get(i).getX()-halfwithd);
			int tmpendindex = extendPoinstSet.GetClosetIndexOfX(extendPoinstSet.Data.get(i).getX()+halfwithd);
			
			if((i-tmpstartindex)>(tmpendindex-i)) {
				for(int j =2*i-tmpendindex; j<=tmpendindex; j++ ) {
					tmpWindowRT.add(extendPoinstSet.Data.get(j).getX());
					tmpWindowInty.add(extendPoinstSet.Data.get(j).getY());
				}
			}
			else {
				for(int j = tmpstartindex; j <= 2*i-tmpstartindex; j++) {
					tmpWindowRT.add(extendPoinstSet.Data.get(j).getX());
					tmpWindowInty.add(extendPoinstSet.Data.get(j).getY());
				}
			}
			FloatArrayList tmpKernal=gauskernal(tmpWindowRT, extendPoinstSet.Data.get(i).getX());
			resultData.AddPoint(extendPoinstSet.Data.get(i).getX(),converlution(tmpWindowInty, tmpKernal));
			tmpKernal = null;
			tmpWindowRT = null;
			tmpWindowInty = null;							
		}

			
	}
	
	public XYPointCollection getsmoothing() {
		return resultData;
	}
	
	public float converlution(FloatArrayList orgData, FloatArrayList kernal) {
		float result = 0f;
		for (int i = 0; i < orgData.size(); i++) {
			result += orgData.get(i)*kernal.get(orgData.size()-1-i);			
		}
		return result;
	}
	
	public FloatArrayList gauskernal(FloatArrayList rt, float centerTime){
		FloatArrayList resultKernal = new FloatArrayList();
		FloatArrayList kernal = new FloatArrayList();
		float kernalSum = 0f;
		float tmpValue =0f;
		for (int i = 0; i < rt.size(); i++) {
			tmpValue = (float)Math.exp(-((rt.get(i)-centerTime)*(rt.get(i)-centerTime))/(2f*sigma*sigma)) / ((float)(Math.sqrt(2*Math.PI)*sigma));
			resultKernal.add(tmpValue);		
			kernalSum += tmpValue;
		}
		for (int i = 0; i <resultKernal.size(); i++) {
			kernal.add(resultKernal.get(i)/kernalSum);
			
		}
		
		return kernal;
	}
}