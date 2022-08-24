

package GridMass.module;

import java.text.Format;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import org.apache.log4j.Logger;

import GridMass.datamodule.DataPoint;
import GridMass.datamodule.Scan;
import bsh.This;
import gnu.trove.impl.hash.TIntShortHash;
import GridMass.datamodule.IndexedDataPoint;
import GridMass.datamodule.IsotopicTraces;
import weka.core.parser.java_cup.internal_error;
import MSUmpire.BaseDataStructure.InstrumentParameter;

public class GridMass{

  private HashMap<Integer, DataPoint[]> dpCache = null;


  private HashMap scantort;
  private int totalScans;
  private Scan[] scans;
  private int[] scanNumbers;
  Datum[] roi[];
  double retentiontime[];

  // User parameters
//  private String suffix;
  private double mzTol = 0.01;
  private double intensitySimilarity = 0.5;
  private double minimumTimeSpan = 0.5, maximumTimeSpan = 5;
  private double smoothTimeSpan = 0.3, smoothTimeMZ = 0.05, smoothMZ = 0.05;
  private double additionTimeMaxPeaksPerScan;
  private double minimumHeight = 10;
  private double rtPerScan;
  private int tolScans;
  private int maxTolScans;
//  private int debug = 0;

  private double minMasa = 0;
  private double maxMasa = 0;

  private ArrayList<IsotopicTraces> resultIsotopicTraces = new ArrayList();
  private int newPeakID=0;

  private String ignoreTimes = "";

  public GridMass(HashMap scantort, Scan[] scans, InstrumentParameter parameter) {
	  this.scans = scans;
	  this.scantort = scantort;
	  int i = 0;
	  scanNumbers = new int[scantort.size()];
	  for(Object keys: scantort.keySet()) {
		  scanNumbers[i] = (int) keys;
		  i++;
	  }
	  this.totalScans = scans.length;
	  this.mzTol = parameter.mzTol;
	  this.minimumTimeSpan = parameter.minimumTimeSpan;
	  this.maximumTimeSpan = parameter.maximumTimeSpan;
	  this.smoothTimeSpan = parameter.smoothTimeSpan;
	  this.smoothTimeMZ = parameter.smoothTimeMZ;
	  this.smoothMZ = parameter.smoothMZ;
	  this.minimumHeight = parameter.minimumHeight;
  }


  public void run() {

    // Check if we have any scans
    if (totalScans == 0) {
    	Logger.getRootLogger().error("No scans match the selected criteria");
    	System.exit(1);
    }

    // Check if the scans are properly ordered by RT
    double prevRT = Double.NEGATIVE_INFINITY;
    for (Scan s : scans) {
      if (s.getRetentionTime() < prevRT) {
        final String msg = "Retention time of scan #" + s.getScanNumber()
            + " is smaller then the retention time of the previous scan."
            + " Please make sure you only use scans with increasing retention times."
            + " You can restrict the scan numbers in the parameters, or you can use the Crop filter module";
        Logger.getRootLogger().error(msg);
        System.exit(1);
      }
      prevRT = s.getRetentionTime();
    }

    int j;
    // minimumTimeSpan
    Scan scan = scans[0];
    double minRT = scan.getRetentionTime();
    double maxRT = scan.getRetentionTime();
    retentiontime = new double[totalScans];
    int i;
    for (i = 0; i < totalScans; i++) {
      scan = scans[i];
      double irt = scan.getRetentionTime();
      if (irt < minRT)
        minRT = irt;
      if (irt > maxRT)
        maxRT = irt;
      retentiontime[i] = irt;
    }
    rtPerScan = (maxRT - minRT) / i;
    // "tolerable" units in scans
    tolScans = Math.max(2, (int) ((minimumTimeSpan / rtPerScan)));
    maxTolScans = Math.max(2, (int) ((maximumTimeSpan / rtPerScan)));

    // Algorithm to find masses:
    // (1) copy masses:intensity > threshold
    // (2) sort intensities descend
    // (3) Find "spot" for each intensity
    // (3.1) if they have not spot ID assigned
    // (3.1.1) Extend mass in mass and time while > 70% pixels > threshold
    // (3.1.2) If extension > mintime ==> mark all pixels with the spot ID
    // (3.1.3) if extension < mintime ==> mark all pixels with spot ID = -1
    // (4) Group spots within a time-tolerance and mass-tolerance

    roi = new Datum[totalScans][];
    ArrayList<Datum> roiAL = new ArrayList<Datum>();
    long passed = 0, nopassed = 0;
    minMasa = Double.MAX_VALUE;
    maxMasa = 0;
    int maxJ = 0;
    boolean[] scanOk = new boolean[totalScans];
    Arrays.fill(scanOk, true);

    IndexedDataPoint[][] data = smoothDataPoints(smoothTimeSpan, smoothTimeMZ, 0, smoothMZ, 0, minimumHeight);

    for (i = 0; i < totalScans; i++) {
      scan = scans[i];
      IndexedDataPoint mzv[] = data[i]; // scan.getDataPoints();
      double prev = (mzv.length > 0 ? mzv[0].getDatapoint().getMZ() : 0);
      double massSum = 0;
      for (j = 0; j < mzv.length; j++) {
        if (mzv[j].getDatapoint().getIntensity() >= minimumHeight)
          massSum += mzv[j].getDatapoint().getMZ() - prev;
        prev = mzv[j].getDatapoint().getMZ();
        if (mzv[j].getDatapoint().getMZ() < minMasa)
          minMasa = mzv[j].getDatapoint().getMZ();
        if (mzv[j].getDatapoint().getMZ() > maxMasa)
          maxMasa = mzv[j].getDatapoint().getMZ();
      }
      double dm = 100.0 / (maxMasa - minMasa);
      if (scanOk[i]) {
        if (!scanOk[i]) {
          // Disable neighbouring scans, how many ?
          for (j = i; j > 0
              && retentiontime[j] + additionTimeMaxPeaksPerScan > retentiontime[i]; j--) {
            scanOk[j] = false;
          }
          for (j = i; j < totalScans
              && retentiontime[j] - additionTimeMaxPeaksPerScan < retentiontime[i]; j++) {
            scanOk[j] = false;
          }
        }       
      }
    }

    String[] it = ignoreTimes.trim().split(", ?");
    for (j = 0; j < it.length; j++) {
      String itj[] = it[j].split("-");
      if (itj.length == 2) {
        Double a = Double.parseDouble(itj[0].trim());
        Double b = Double.parseDouble(itj[1].trim());
        for (i = Math.abs(Arrays.binarySearch(retentiontime, a)); i < totalScans
            && retentiontime[i] <= b; i++) {
          if (retentiontime[i] >= a) {
            scanOk[i] = false;
          }
        }
      }
    }

    passed = 0;
    nopassed = 0;
    for (i = 0; i < totalScans; i++) {
      if (scanOk[i]) {
        scan = scans[i];
        IndexedDataPoint mzv[] = data[i];
        DataPoint mzvOriginal[] = scan.getDataPoints();
        ArrayList<Datum> dal = new ArrayList<Datum>();
        for (j = 0; j < mzv.length; j++) {
          if (mzv[j].getDatapoint().getIntensity() >= minimumHeight) {
            dal.add(new Datum(mzv[j].getDatapoint(), i, mzvOriginal[mzv[j].getIndex()]));
            passed++;
          } else {
            nopassed++;
          }
        }
        if (j > maxJ)
          maxJ = j;
        roi[i] = dal.toArray(new Datum[0]);
        roiAL.addAll(dal);
      }
    }

    // New "probing" algorithm
    // (1) Generate probes all over chromatograms
    // (2) Move each probe to their closest maximum until it cannot find a
    // new maximum
    // (3) assign spot id to each "center" using all points within region

    // (1) Generate probes all over
    double byMZ = Math.max(mzTol * 2, 1e-6);
    int byScan = Math.max(1, tolScans / 4);
    double m;
    int ndata = (int) Math
        .round((((double) totalScans / (double) byScan) + 1) * ((maxMasa - minMasa + byMZ) / byMZ));
    Probe probes[] = new Probe[ndata];
    int idata = 0;
    for (i = 0; i < totalScans; i += byScan) {
      for (m = minMasa - (i % 2) * byMZ / 2; m <= maxMasa; m += byMZ) {
        probes[idata++] = new Probe(m, i);
      }
    }

    double mzR = byMZ / 2;
    int scanR = Math.max(byScan - 1, 2);
    int okProbes = 0;
    for (i = 0; i < idata; i++) {
      moveProbeToCenter(probes[i], scanR, mzR);
      if (probes[i].intensityCenter < minimumHeight) {
        probes[i] = null;
      } else {
        okProbes++;
      }
    }
    if (okProbes > 0) {
      Probe[] pArr = new Probe[okProbes];
      for (okProbes = i = 0; i < idata; i++) {
        if (probes[i] != null) {
          pArr[okProbes++] = probes[i];
        }
      }
      probes = pArr;
      pArr = null;
    }
    // (3) Assign spot id to each "center"
    Arrays.sort(probes);
    SpotByProbes sbp = new SpotByProbes();
    ArrayList<SpotByProbes> spots = new ArrayList<SpotByProbes>();
    double mzA = -1;
    int scanA = -1;
    for (i = 0; i < probes.length; i++) {
      if (probes[i] != null && probes[i].intensityCenter >= minimumHeight) {
        if (probes[i].mzCenter != mzA || probes[i].scanCenter != scanA) {
          if (sbp.size() > 0) {
            spots.add(sbp);
            sbp.assignSpotId();

          }
          sbp = new SpotByProbes();
          mzA = probes[i].mzCenter;
          scanA = probes[i].scanCenter;
        }
        sbp.addProbe(probes[i]);
      }
    }
    if (sbp.size() > 0) {
      spots.add(sbp);
      sbp.assignSpotId();
    }
    i = 0;
    for (SpotByProbes sx : spots) {
      if (sx.size() > 0) {
        assignSpotIdToDatumsFromScans(sx, scanR, mzR);
      }
    }

    // (4) Join Tolerable Centers
    int criticScans = Math.max(1, tolScans / 2);
    int joins = 0;
    for (i = 0; i < spots.size() - 1; i++) {
      SpotByProbes s1 = spots.get(i);
      if (s1.center != null && s1.size() > 0) {
        for (j = i; j > 0 && j < spots.size() && spots.get(j - 1).center != null
            && spots.get(j - 1).center.mzCenter + mzTol > s1.center.mzCenter; j--);
        for (; j < spots.size(); j++) {
          SpotByProbes s2 = spots.get(j);
          if (i != j && s2.center != null) {
            if (s2.center.mzCenter - s1.center.mzCenter > mzTol)
              break;
            int l = Math.min(Math.abs(s1.minScan - s2.minScan), Math.abs(s1.minScan - s2.maxScan));
            int r = Math.min(Math.abs(s1.maxScan - s2.minScan), Math.abs(s1.maxScan - s2.maxScan));
            int d = Math.min(l, r);
            boolean overlap = !(s2.maxScan < s1.minScan || s2.minScan > s1.maxScan);
            if ((d <= criticScans || overlap) && (intensityRatio(s1.center.intensityCenter,
                s2.center.intensityCenter) > intensitySimilarity)) {
              assignSpotIdToDatumsFromSpotId(s1, s2, scanR, mzR);
              s1.addProbesFromSpot(s2, true);
              j = i; // restart
              joins++;
            }
          }
        }
      }
    }

    // (5) Remove "Large" spanned masses
    for (i = 0; i < spots.size() - 1; i++) {
      SpotByProbes s1 = spots.get(i);
      if (s1.center != null && s1.size() > 0) {
        int totalScans = s1.maxScan - s1.minScan + 1;
        int lScan = s1.minScan;
        int rScan = s1.maxScan;
        ArrayList<Integer> toRemove = new ArrayList<Integer>();
        toRemove.add(i);
        for (j = i; j > 0 && j < spots.size() && spots.get(j - 1).center != null
            && spots.get(j - 1).center.mzCenter + mzTol > s1.center.mzCenter; j--);
        for (; j < spots.size(); j++) {
          SpotByProbes s2 = spots.get(j);
          if (i != j && s2.center != null) {
            if (s2.center.mzCenter - s1.center.mzCenter > mzTol)
              break;
            if (intensityRatio(s1.center.intensityCenter,
                s2.center.intensityCenter) > intensitySimilarity) {
              int dl = Math.min(Math.abs(lScan - s2.minScan), Math.abs(lScan - s2.maxScan));
              int dr = Math.min(Math.abs(rScan - s2.minScan), Math.abs(rScan - s2.maxScan));
              int md = Math.min(dl, dr);
              if (md <= maxTolScans || !(s2.maxScan < lScan || s2.minScan > rScan)) {
                // distancia tolerable o intersectan
                totalScans += s2.maxScan - s2.minScan + 1;
                toRemove.add(j);
                lScan = Math.min(lScan, s2.minScan);
                rScan = Math.max(rScan, s2.maxScan);
              }
            }
          }
        }
        if (totalScans * rtPerScan > maximumTimeSpan) {
          for (Integer J : toRemove) {
            spots.get(J).clear();
          }
        }
      }
    }

    // Build peaks from assigned datums
//    i = 0;
    for (SpotByProbes sx : spots) {
      if (sx.size() > 0 && sx.maxScan - sx.minScan + 1 >= tolScans) {
        sx.buildMaxDatumFromScans(roi, minimumHeight);
        if (sx.getMaxDatumScans() >= tolScans && (sx.getContigousMaxDatumScans() >= tolScans
            || sx.getContigousToMaxDatumScansRatio() > 0.5)) {
          IsotopicTraces peak = new IsotopicTraces();
          if (addMaxDatumFromScans(sx, peak) > 0) {
            peak.finishChromatogram(scantort);
            if (peak.getArea() > 1e-6) {
//            	resultIsotopicTraces[newPeakID] = peak;
            	peak.calTargeMZ();
            	this.resultIsotopicTraces.add(peak);
            	newPeakID++;
            }
          }
        }
      }
    }

  }

  public double intensityRatio(double int1, double int2) {
    return Math.min(int1, int2) / Math.max(int1, int2);
  }


  public IndexedDataPoint[][] smoothDataPoints(double timeSpan,
      double timeMZSpan, int scanSpan, double mzTol, int mzPoints, double minimumHeight) {
    int[] scanNumbers = this.scanNumbers;
    int totalScans = scanNumbers.length;
    DataPoint mzValues[][] = null; // [relative scan][j value]
    DataPoint mzValuesJ[] = null;
    int mzValuesScan[] = null;
    int mzValuesMZidx[] = null;
    IndexedDataPoint newMZValues[][] = null;
    IndexedDataPoint tmpDP[] = new IndexedDataPoint[0];
    newMZValues = new IndexedDataPoint[totalScans][];
    int i, j, si, sj, ii, k, ssi, ssj, m;
    double timeSmoothingMZtol = Math.max(timeMZSpan, 1e-6);

    int modts = Math.max(1, totalScans / 10);

    for (i = 0; i < totalScans; i++) {

      // Smoothing in TIME space
      Scan scan = scans[i];
      double rt = retentiontime[i];
      DataPoint[] xDP = null;
      IndexedDataPoint[] iDP = null;
      sj = si = i;
      ssi = ssj = i;
      int t = 0;
      if (timeSpan > 0 || scanSpan > 0) {
        if (scan != null) {
          for (si = i; si > 1; si--) {
            if (retentiontime[si - 1] < rt - timeSpan / 2) {
              break;
            }
          }
          for (sj = i; sj < totalScans - 1; sj++) {
            if (retentiontime[sj + 1] >= rt + timeSpan / 2) {
              break;
            }
          }
          ssi = i - (scanSpan - 1) / 2;
          ssj = i + (scanSpan - 1) / 2;
          if (ssi < 0) {
            ssj += -ssi;
            ssi = 0;
          }
          if (ssj >= totalScans) {
            ssi -= (ssj - totalScans + 1);
            ssj = totalScans - 1;
          }
          if (sj - si + 1 < scanSpan) {
            si = ssi;
            sj = ssj;
          }
        }
        if (scan != null && sj > si) {
          // Allocate
          if (mzValues == null || mzValues.length < sj - si + 1) {
            mzValues = new DataPoint[sj - si + 1][];
            mzValuesScan = new int[sj - si + 1];
            mzValuesMZidx = new int[sj - si + 1];
          }
          // Load Data Points
          for (j = si; j <= sj; j++) {
            int jsi = j - si;
            if (mzValues[jsi] == null || jsi >= mzValuesScan.length - 1
                || mzValuesScan[jsi + 1] != scanNumbers[j]) {
              Scan xscan = scans[j];
              mzValues[jsi] = xscan.getDataPoints();
              mzValuesScan[jsi] = scanNumbers[j];
            } else {
              mzValues[jsi] = mzValues[jsi + 1];
              mzValuesScan[jsi] = mzValuesScan[jsi + 1];
            }
            mzValuesMZidx[jsi] = 0;
          }
          // Estimate Averages
          ii = i - si;
          if (tmpDP.length < mzValues[ii].length)
            tmpDP = new IndexedDataPoint[mzValues[ii].length * 3 / 2];
          for (k = 0; k < mzValues[ii].length; k++) {
            DataPoint dp = mzValues[ii][k];
            double mz = dp.getMZ();
            double intensidad = 0;
            if (dp.getIntensity() > 0) { // only process those > 0
              double a = 0;
              short c = 0;
              int f = 0;
              for (j = 0; j <= sj - si; j++) {
                for (mzValuesJ = mzValues[j]; mzValuesMZidx[j] < mzValuesJ.length - 1
                    && mzValuesJ[mzValuesMZidx[j] + 1].getMZ() < mz
                        - timeSmoothingMZtol; mzValuesMZidx[j]++);

                f = mzValuesMZidx[j];

                for (m = mzValuesMZidx[j] + 1; m < mzValuesJ.length
                    && mzValuesJ[m].getMZ() < mz + timeSmoothingMZtol; m++) {
                  if (Math.abs(mzValuesJ[m].getMZ() - mz) < Math.abs(mzValuesJ[f].getMZ() - mz)) {
                    f = m;
                  } else {
                    // siempre debe ser mas cercano porque
                    // están ordenados por masa, entonces
                    // parar la búsqueda
                    break;
                  }
                }
                if (f > 0 && f < mzValuesJ.length
                    && Math.abs(mzValuesJ[f].getMZ() - mz) <= timeSmoothingMZtol
                    && mzValuesJ[f].getIntensity() > 0) { // >=
                  // minimumHeight
                  // ?
                  // System.out.println("mz="+mz+"; Closer="+mzValuesJ[f].getMZ()+", f="+f+",
                  // Intensity="+mzValuesJ[f].getIntensity());
                  a += mzValuesJ[f].getIntensity();
                  c++;
                }
              }
              intensidad = c > 0 ? a / c : 0;
              if (intensidad >= minimumHeight) {
                tmpDP[t++] = new IndexedDataPoint(k,new DataPoint(mz, intensidad,dp.getRT()));
              }
            }
          }

        }
      } else if (scan != null) {
        xDP = scan.getDataPoints();
        if (tmpDP.length < xDP.length)
          tmpDP = new IndexedDataPoint[xDP.length];
        for (k = 0; k < xDP.length; k++) {
          if (xDP[k].getIntensity() >= minimumHeight) {
            tmpDP[t++] = new IndexedDataPoint(k, xDP[k]);
          }
        }
      }
      iDP = new IndexedDataPoint[t];
      for (k = 0; k < t; k++) {
        iDP[k] = tmpDP[k];
      }
      newMZValues[i] = iDP;
    }

    return newMZValues;
  }


  int addMaxDatumFromScans(SpotByProbes s, IsotopicTraces peak) {

    int i, j;
    int adds = 0;
    for (i = s.minScan; i <= s.maxScan; i++) {
      Datum[] di = roi[i];
      if (di != null && di.length > 0) {
        Datum max = new Datum(new DataPoint(0, -1, 0), 0, new DataPoint(0, -1, 0));
        int idx = findFirstMass(s.minMZ, di);
        for (j = idx; j < di.length && di[j].mz <= s.maxMZ; j++) {
          Datum d = di[j];
          if (d.spotId == s.spotId) {
            if (d.intensity > max.intensity && d.mz >= s.minMZ && d.intensity > minimumHeight) {
              max = d;
            }
          }
        }
        if (max.intensity > 0) {
          adds++;
          peak.addMzPeak(scans[i].getScanNumber(),
              new DataPoint(max.mzOriginal, max.intensityOriginal, max.rt));
        }
      }
    }
    return adds;
  }

  void assignSpotIdToDatumsFromScans(SpotByProbes s, int sRadius, double mzRadius) {

    int i, j;
    for (i = s.minScan; i <= s.maxScan; i++) {
      Datum[] di = roi[i];
      if (di != null && di.length > 0) {
        int idx = findFirstMass(s.minMZ - mzRadius, di);
        for (j = idx; j < di.length && di[j].mz <= s.maxMZ + mzRadius; j++) {
          Datum d = di[j];
          if (d.mz >= s.minMZ - mzRadius) {
            if (d.spotId != 0) {
              // Some spot already assigned this to it. Check
              // exactly who is the winner
              Probe p = new Probe(d.mz, d.scan);
              moveProbeToCenter(p, sRadius, mzRadius);
              if (p.mzCenter == s.center.mzCenter && p.scanCenter == s.center.scanCenter) {
                // This datum is actually MINE (s) !!!, this
                // will happen to datums close to spot borders
                // and that compete with other spot
                // System.out.println("Reassigning spot to Id="+s.spotId+" from
                // Spot:"+d.toString());
                s.setSpotIdToDatum(d);
              }
            } else {
              s.setSpotIdToDatum(d);
            }
          }
        }
      }
    }
  }

  void assignSpotIdToDatumsFromSpotId(SpotByProbes s, SpotByProbes s2, int sRadius,
      double mzRadius) {

    int i, j;
    int oldSpotId = s2.spotId;
    int mxScan = Math.max(s.maxScan, s2.maxScan);
    double minMZ = Math.min(s.minMZ, s2.minMZ);
    double maxMZ = Math.max(s.maxMZ, s2.maxMZ);
    for (i = Math.min(s.minScan, s2.minScan); i <= mxScan; i++) {
      Datum[] di = roi[i];
      if (di != null && di.length > 0) {
        int idx = findFirstMass(minMZ - mzRadius, di);
        for (j = idx; j < di.length && di[j].mz <= maxMZ + mzRadius; j++) {
          Datum d = di[j];
          if (d.spotId == oldSpotId) {
            s.setSpotIdToDatum(d);
          }
        }
      }
    }
  }

  void moveProbeToCenter(Probe p, int sRadius, double mzRadius) {

    int i, j, k;
    double maxMZ, minMZ;
    boolean move = true;
    Datum max = new Datum(new DataPoint(0, -1, 0), 0, new DataPoint(0, -1, 0));
    while (move) {
      k = Math.min(totalScans - 1, p.scanCenter + sRadius);
      for (i = Math.max(p.scanCenter - sRadius, 0); i <= k; i++) {
        Datum[] di = roi[i];
        if (di != null && di.length > 0) {
          minMZ = p.mzCenter - mzRadius;
          int idx = findFirstMass(minMZ, di);
          maxMZ = p.mzCenter + mzRadius;
          for (j = idx; j < di.length && di[j].mz <= maxMZ; j++) {
            Datum d = di[j];
            if (d.intensity > max.intensity && d.mz >= minMZ) {
              max = d;
            }
          }
        }
      }
      if (max.intensity >= 0 && (max.mz != p.mzCenter || max.scan != p.scanCenter)) {
        p.mzCenter = max.mz;
        p.scanCenter = max.scan;
        p.intensityCenter = max.intensity;
        // p.moves++;
      } else {
        move = false;
      }
    }
  }

  static int findFirstMass(double mass, DataPoint mzValues[]) {
    int l = 0;
    int r = mzValues.length - 1;
    int mid = 0;
    while (l < r) {
      mid = (r + l) / 2;
      if (mzValues[mid].getMZ() > mass) {
        r = mid - 1;
      } else if (mzValues[mid].getMZ() < mass) {
        l = mid + 1;
      } else {
        return mid;
      }
    }
    while (l > 0 && mzValues[l].getMZ() > mass)
      l--;
    return l;
  }

  static int findFirstMass(double mass, Datum mzValues[]) {
    return findFirstMass(mass, mzValues, 0, mzValues.length - 1);
  }

  static int findFirstMass(double mass, Datum mzValues[], int l, int r) {
    int mid = 0;
    while (l < r) {
      mid = (r + l) / 2;
      if (mzValues[mid].mz > mass) {
        r = mid - 1;
      } else if (mzValues[mid].mz < mass) {
        l = mid + 1;
      } else {
        return mid;
      }
    }
    while (l > 0 && mzValues[l].mz > mass)
      l--;
    return l;
  }

  public ArrayList<IsotopicTraces> getFinalresult(){
	  return this.resultIsotopicTraces;
  }
  

}
