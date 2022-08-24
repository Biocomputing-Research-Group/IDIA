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

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.PeakDataStructure.PeakCurve;
import jehep.actions.newAction;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import MSUmpire.MathPackage.GaussianFilter;

/**
 * Peak shape smoothing process thread unit
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class PeakCurveSmoothingUnit implements Callable<ArrayList<PeakCurve>> {

    PeakCurve curve;
    boolean export;
    String msLevelString;

    InstrumentParameter parameter;

    public PeakCurveSmoothingUnit(PeakCurve curve, InstrumentParameter para, String msLevelString) {
        this.curve = curve;
        this.parameter = para;
        this.msLevelString=msLevelString;
    }

    @Override
    public ArrayList<PeakCurve> call() {
        final ArrayList<PeakCurve> ResultCurves;
        //If we want to split multimodal peak curve by CWT
        if (parameter.DetectByCWT) {
        	if(msLevelString.equals("ms1")) {
        		curve.DoBspline();
//        		curve.newDoBspline();
        	}
        	else {
//        		curve.newDoBspline();
        		curve.DoBspline();
        	}

            curve.DetectPeakRegion();
            ResultCurves = curve.SeparatePeakByRegion(parameter.SNThreshold);
            curve=null;
        }
        else{
            curve.DoBspline();
            ResultCurves=new ArrayList<>();
            ResultCurves.add(curve);
        }

        for (final PeakCurve peakCurve : ResultCurves) {
            peakCurve.CalculateMzVar();
            peakCurve.StartRT();
            peakCurve.EndRT();
            
            
            GaussianFilter gau = new GaussianFilter(peakCurve.GetPeakCollection(), peakCurve.StartRT(), peakCurve.EndRT(),"diaumpire");
            gau.gausfilter();
            peakCurve.SmoothData = gau.getsmoothing();
//            peakCurve.newdata = gau.getsmoothing();
            gau=null;
            
            
            peakCurve.ReleaseRawPeak();
            peakCurve.types="diaumpire";
        }

        return ResultCurves;
    }
}
