package GridMass.module;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.io.StringReader;
//import java.lang.reflect.Array;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
//import java.util.List;
//import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;
//import java.util.concurrent.ForkJoinPool;
//import java.util.concurrent.ForkJoinTask;
//import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.codec.binary.Base64;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
//import org.apache.commons.net.telnet.WindowSizeOptionHandler;
//import org.apache.derby.tools.sysinfo;
import org.apache.log4j.Logger;
import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import com.google.common.collect.Range;

import org.w3c.dom.Document;
import org.w3c.dom.Node;


import GridMass.datamodule.DataPoint;
import GridMass.datamodule.IsotopicTraces;
import GridMass.datamodule.Scan;
import GridMass.module.GridMass;
import bsh.This;
import weka.core.parser.java_cup.internal_error;

import MSUmpire.BaseDataStructure.InstrumentParameter;

public final class GridMassModel {
	private int windowSize=25;
//	private HashMap<Range, ArrayList<Integer>> ms2WindowsRange = new HashMap();
	private HashMap<String, ArrayList<Integer>> ms2WindowsRange = new HashMap();
	private String filename =""; //"/home/UNT/jl0734/Desktop/DIA/myeclipseproject/napedro_L120417_001_SW_400AQUA_human_2ul_dilution_10.mzXML";
	public HashMap<Integer,Double> scantort = new HashMap();
	public HashMap<Double, Integer> rttoscan = new HashMap();
	private TreeMap<Integer, Long> ScanIndex = new TreeMap();
	private HashMap<Integer, Integer> msLevels = new HashMap();
//	private IntArrayList ms1ScansNum = new IntArrayList();
	private ArrayList<Integer> ms1ScansNum = new ArrayList();

	
	private InstrumentParameter parameter;
	
	public  GridMassModel(InstrumentParameter parameter, String Filename) {
		this.parameter = parameter;
		this.filename=Filename;
		getOffset();
		elutionImfo();
	}
	
	public ArrayList<IsotopicTraces> startrun(double startMZ, double endMZ,int msLevel) {
		ArrayList<IsotopicTraces> resultTraces;
//		Range tmpwindowRange = Range.closed(startMZ, endMZ);
		String tmpwindowRange = String.valueOf(startMZ).split("[.]")[0]+"-"+String.valueOf(endMZ).split("[.]")[0];
		HashMap< Integer, Double> taskScantort = new HashMap();
		Scan[] msScans = null;
		if (msLevel == 1) {
			try {		
				msScans = readLCMSData(ms1ScansNum);		
			}
			catch (Exception e) {
				new RuntimeException(e);
			}
//			System.out.println("set Parameters!");
			for(Integer tmpkey: scantort.keySet()){
				if (msLevels.get(tmpkey) == 1) {
					taskScantort.put(tmpkey,scantort.get(tmpkey));
				}
			}
			
		}
		else {
			try {			
				msScans = readLCMSData(ms2WindowsRange.get(tmpwindowRange));		
			}
			catch (Exception e) {
				System.out.println(tmpwindowRange);
				System.out.println("---------");
				for(String aaa: ms2WindowsRange.keySet()) {
					System.out.println(aaa);
				}
//				new RuntimeException(e);
			}
//			System.out.println("set Parameters!");
			for(Integer tmpkey: ms2WindowsRange.get(tmpwindowRange)) {
				taskScantort.put(tmpkey, scantort.get(tmpkey));
			}

		}
		GridMass gridmassprocess = new GridMass(taskScantort,msScans,parameter);
//		Logger.getRootLogger().info("create gridmass object!");
		taskScantort = null;
		gridmassprocess.run();
		resultTraces = gridmassprocess.getFinalresult();
		gridmassprocess = null;
//		Logger.getRootLogger().info("finish!");
		return resultTraces;
	}
	
	public Scan readScan(String XMLtext)throws ParserConfigurationException, SAXException, IOException, DataFormatException {
		String readCompressionType = "zlib";
		int readPrecision = 32;
		double readbasePeakMZ = 0;
		double readbasePeakIntensity = 0;
		Scan oneScan = new Scan();
		
        if (XMLtext.replaceFirst("</scan>", "").contains("</scan>")) {
            XMLtext = XMLtext.replaceFirst("</scan>", "");
        }
        if (!XMLtext.contains("</scan>")) {
            XMLtext += "</scan>";
        }
        
        Document doc = null;
        try {
        	final DocumentBuilder docBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            
            docBuilder.reset();
            InputSource input = new InputSource(new StringReader(XMLtext));
            
            doc = docBuilder.parse(input);
        } catch (Exception ex) {
            Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            Logger.getRootLogger().error(XMLtext);
        }
        Node root = doc.getFirstChild();
        for (int i = 0; i < root.getAttributes().getLength(); i++) {
            switch (root.getAttributes().item(i).getNodeName()) {
                case ("num"):
                    oneScan.setScanNumber(Integer.parseInt(root.getAttributes().item(i).getNodeValue()));
                    break;
                case ("centroided"): {
                    if ("1".equals(root.getAttributes().item(i).getNodeValue())) {
                    	oneScan.setcentroided(true);
                    } else {
                    	oneScan.setcentroided(false);
                    }
                    break;
                }
                case ("msLevel"):
                	oneScan.setMSLevel(Integer.parseInt(root.getAttributes().item(i).getNodeValue()));
                    break;
                case ("scanType"):
                	oneScan.setscanType(root.getAttributes().item(i).getNodeValue());
                    break;
                case ("peaksCount"):
                	oneScan.setpeakCount(Integer.parseInt(root.getAttributes().item(i).getNodeValue()));
                    break;
                case ("retentionTime"):
                	oneScan.setRetentionTime(Double.parseDouble(root.getAttributes().item(i).getNodeValue().substring(2, root.getAttributes().item(i).getNodeValue().indexOf("S"))) / 60d);
                    break;
//                case ("lowMz"): {
//                    String value = root.getAttributes().item(i).getNodeValue();
//                    if ("inf".equals(value)) {
//                        value = String.valueOf(Double.MIN_VALUE);
//                    }
//                    oneScan.setStartMZ(Double.parseDouble(value));
//                    break;
//                }
//                case ("highMz"): {
//                    String value = root.getAttributes().item(i).getNodeValue();
//                    if ("inf".equals(value)) {
//                        value = String.valueOf(Double.MAX_VALUE);
//                    }
//                    oneScan.setEndMZ(Double.parseDouble(value));
//                    break;
//                }
//                case ("startMz"):{
//                	String value = root.getAttributes().item(i).getNodeValue();
//                	if ("inf".equals(value)) {
//                        value = String.valueOf(Double.MIN_VALUE);
//                    }
//                	oneScan.setStartMZ(Double.parseDouble(value));
//                    break;
//                }
//                case ("endMz"): {
//                    String value = root.getAttributes().item(i).getNodeValue();
//                    if ("inf".equals(value)) {
//                        value = String.valueOf(Double.MAX_VALUE);
//                    }
//                    oneScan.setEndMZ(Double.parseDouble(value));
//                    break;
//                }
                case ("basePeakMz"): {
                    String value = root.getAttributes().item(i).getNodeValue();
                    if ("inf".equals(value)) {
                        value = String.valueOf(Double.MAX_VALUE);
                    }
                    if (!"null".contentEquals(value)) {
                       readbasePeakMZ = Double.parseDouble(value);}
                    break;
                }
                case ("basePeakIntensity"):
                    final String value = root.getAttributes().item(i).getNodeValue();
                    if (!"null".contentEquals(value)) {
                        readbasePeakIntensity = Double.parseDouble(value);}
                    break;
                case ("totIonCurrent"):
                	oneScan.settotalIonCurrent(Double.parseDouble(root.getAttributes().item(i).getNodeValue()));
                    break;
            }
        }
        for (int i = 0; i < root.getChildNodes().getLength(); i++) {
            Node childNode = root.getChildNodes().item(i);
            switch (childNode.getNodeName()) {
                case ("precursorMz"): {
                    oneScan.setPrecusorMZ(Double.parseDouble(childNode.getTextContent()));
                    for (int j = 0; j < childNode.getAttributes().getLength(); j++) {
                        switch (childNode.getAttributes().item(j).getNodeName()) {
//                            case ("precursorScanNum"):
//                                int readPrecursorScanNum = Integer.parseInt(childNode.getAttributes().item(j).getNodeValue());
//                                break;
//                            case ("precursorIntensity"):
//                                double readPrecursorIntensity = Double.parseDouble(childNode.getAttributes().item(j).getNodeValue());
//                                break;
//                            case ("precursorCharge"):
//                                int readPrecursorCharge = Integer.parseInt(childNode.getAttributes().item(j).getNodeValue());
//                                break;
                            case ("activationMethod"):
                            	oneScan.setactivationMethod(childNode.getAttributes().item(j).getNodeValue());
                                break;
//                            case ("windowWideness"):
//                            	oneScan.setwindowWideness(Float.parseFloat(childNode.getAttributes().item(j).getNodeValue()));
//                                break;
                        }
                    }
                    break;
                }
                case ("peaks"): {
                    for (int j = 0; j < childNode.getAttributes().getLength(); j++) {
                        switch (childNode.getAttributes().item(j).getNodeName()) {
                            case ("compressionType"):
                                readCompressionType = childNode.getAttributes().item(j).getNodeValue();
                                break;
                            case ("precision"):
                                readPrecision = Integer.parseInt(childNode.getAttributes().item(j).getNodeValue());
                                break;
                        }
                    }
//                    ParsePeakString(scan, childNode.getTextContent());
                    oneScan.setDataPoints(decodevalues(childNode.getTextContent(), readPrecision, readCompressionType, oneScan.getpeakCount(), oneScan.getRetentionTime()));
                    break;
                }
            }
            childNode = null;
        }
        if ("calibration".equals(oneScan.getscanType())) {
            oneScan.setMSLevel(-1);
        }
        oneScan.setbasePeak(new DataPoint(readbasePeakMZ, readbasePeakIntensity,oneScan.getRetentionTime()));
        XMLtext = null;
        return oneScan;
	}
	private DataPoint[] decodevalues(String stringpeak, int precision, String compressionType, int peaksCount, double rt) throws IOException, DataFormatException  {
		int offset;
		DataPoint[] allPeaks = new DataPoint[peaksCount];

		stringpeak = stringpeak.replaceAll("\n", "");
        byte[] decoded = Base64.decodeBase64(stringpeak.getBytes());

        if ("zlib".equals(compressionType)) {
            decoded = ZlibUncompressBuffer(decoded);
        }
        switch (precision) {
            case (32): {
                offset = 0;
                for (int i = 0; i < peaksCount; i++) {
                    byte[] mz = new byte[]{decoded[offset], decoded[offset + 1], decoded[offset + 2], decoded[offset + 3]};
                    byte[] intensity = new byte[]{decoded[offset + 4], decoded[offset + 5], decoded[offset + 6], decoded[offset + 7]};
                    ByteBuffer mzBuffer = ByteBuffer.wrap(mz);
                    ByteBuffer intBuffer = ByteBuffer.wrap(intensity);
                    double intensityd = (double) intBuffer.getFloat();
                    double mzd = (double) mzBuffer.getFloat();
                    if (intensityd > 0d) {
                        allPeaks[i] = new DataPoint(mzd, intensityd,rt);
                    }
                    mz = null;
                    intensity = null;
                    mzBuffer.clear();
                    intBuffer.clear();
                    mzBuffer = null;
                    intBuffer = null;
                    offset += 8;
                }
                break;
            }
            case (64): {
                offset = 0;
                for (int i = 0; i < peaksCount; i++) {
                    byte[] mz = new byte[]{decoded[offset], decoded[offset + 1], decoded[offset + 2], decoded[offset + 3], decoded[offset + 4], decoded[offset + 5], decoded[offset + 6], decoded[offset + 7]};
                    byte[] intensity = new byte[]{decoded[offset + 8], decoded[offset + 9], decoded[offset + 10], decoded[offset + 11], decoded[offset + 12], decoded[offset + 13], decoded[offset + 14], decoded[offset + 15]};
                    ByteBuffer mzBuffer = ByteBuffer.wrap(mz);
                    ByteBuffer intBuffer = ByteBuffer.wrap(intensity);
                    double intensityd = intBuffer.getDouble();
                    double mzd =  mzBuffer.getDouble();
                    if (intensityd > 0d) {
                    	allPeaks[i] = new DataPoint(mzd, intensityd,rt);
                    }
                    mz = null;
                    intensity = null;
                    mzBuffer.clear();
                    intBuffer.clear();
                    mzBuffer = null;
                    intBuffer = null;
                    offset += 16;
                }
                break;
            }
        }
        stringpeak = null;
        decoded = null;
        return allPeaks;
	}
	
    private byte[] ZlibUncompressBuffer(byte[] compressed) throws IOException, DataFormatException {

        Inflater decompressor = new Inflater();
        decompressor.setInput(compressed);

        ByteArrayOutputStream bos = null;
        try {

            bos = new ByteArrayOutputStream(compressed.length);

            // Decompress the data
            byte[] buf = new byte[decompressor.getRemaining() * 2];
            while (decompressor.getRemaining() > 0) {
                int count = decompressor.inflate(buf);
                bos.write(buf, 0, count);
            }

        } finally {
            try {
                bos.close();
            } catch (Exception nope) { /* This exception doesn't matter */ }
        }
        decompressor.end();
        compressed = null;
        decompressor = null;
        byte[] result = bos.toByteArray();
        bos = null;
        return result;
    }
	
	public void getOffset() {
		final LongArrayList offsets = new LongArrayList();
        long pos = 0;
        long posend = -1;
		try(InputStream is = Files.newInputStream(Paths.get(this.filename))){	        
	        BufferedReader bis = new BufferedReader(new InputStreamReader(is, StandardCharsets.ISO_8859_1));	        
	        while (true) {
	            final String line = bis.readLine();
	            if (line == null) {
	                break;
	            }
	            final int idx = line.indexOf("<scan");
	            final int idxend = line.lastIndexOf("</scan>");
	            if (idxend != -1) {
	                posend = pos + idxend + "</scan>".length();
	            }
	            if (idx != -1) {
	                offsets.add(pos + idx);
	            }
	            pos += line.length() + 1;
	        }
	        if (posend < offsets.getLast()) {
	            throw new RuntimeException("mzXML parse error");
	        }
	        offsets.add(posend);
	        for (int i=0;i<offsets.size()-1;++i){
                this.ScanIndex.put(i+1,offsets.get(i));
            }
            this.ScanIndex.put(Integer.MAX_VALUE, offsets.getLast());
            
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		if(ScanIndex==null | ScanIndex.isEmpty()){
			System.out.println(Paths.get(this.filename));
            throw new RuntimeException("mzXML file wrong!");
            
        }
	}
	
	public void elutionImfo(){
//		int NoMS1Scans = 0;
		try(RandomAccessFile fileHandler = new RandomAccessFile(this.filename, "r")){
            Iterator<Entry<Integer, Long>> iter = this.ScanIndex.entrySet().iterator();
            Long currentIdx = iter.next().getValue();
            while (iter.hasNext()) {
                long startposition = currentIdx;
                long nexposition = iter.next().getValue();
                currentIdx = nexposition;
                fileHandler.seek(startposition);
                
                byte[] bufr = new byte[(int) (nexposition - startposition)];
                fileHandler.readFully(bufr, 0, (int) (nexposition - startposition));

                String temp = new String(bufr);

                double rt = 0f;
                int scanno = 0;
                int mslevel = 0;
                //float precursorF=0f;
                if (!temp.contains("<scan")) {
                    fileHandler.close();
                    return;
                }

                if (temp.contains("<scan num=") && (temp.contains("retentionTime=\"PT"))) {
                    String substr = temp.substring(temp.indexOf("<scan num=") + 11);
                    scanno = Integer.parseInt(substr.substring(0, substr.indexOf("\"")));
                    
                    rt = Double.parseDouble(temp.substring(temp.indexOf("retentionTime=\"PT") + 17).split("S\"")[0]);
                    rt=rt/60d;
                    mslevel = Integer.parseInt(temp.substring(temp.indexOf("msLevel=") + 9, temp.indexOf("msLevel=") + 10));
                    if (temp.contains("scanType=\"calibration\"")) {
                        mslevel = -1;
                    }
//                    if (mslevel == 1) {
//                        NoMS1Scans++;                        
//                        if (temp.contains("scanType=\"SIM\"") && datatype == SpectralDataType.DataType.WiSIM) {
//                            int startidx = temp.indexOf("lowMz=\"") + 7;
//                            int stopidx = startidx + 1;
//                            for (int i = startidx + 1; i < temp.length(); i++) {
//                                if (temp.charAt(i) == '\"') {
//                                    stopidx = i;
//                                    break;
//                                }
//                            }
//                            float lowmz = Float.parseFloat(temp.substring(startidx, stopidx));
//                            startidx = temp.indexOf("highMz=\"") + 8;
//                            stopidx = startidx + 1;
//                            for (int i = startidx + 1; i < temp.length(); i++) {
//                                if (temp.charAt(i) == '\"') {
//                                    stopidx = i;
//                                    break;
//                                }
//                            }
//                            float highmz = Float.parseFloat(temp.substring(startidx, stopidx));
//                            for (XYData MS1win : dIA_Setting.MS1Windows.keySet()) {
//                                if (MS1win.getX() <= lowmz && MS1win.getY() >= highmz) {
//                                    dIA_Setting.MS1Windows.get(MS1win).add(scanno);
//                                }
//                            }
//                        }                        
//                    }
                    if (mslevel == 2) {
                    	int stopidx = temp.indexOf("</precursorMz>");
                        if (stopidx == -1) {
                        	Logger.getRootLogger().error("Parsing </precursorMz> failed. scan number :" + scanno);                                    
                            System.exit(3);
                        }
                    int startidx = 0;
                    for (int i = stopidx; i > 0; i--) {
                    	if (temp.charAt(i) == '>') {
                    		startidx = i + 1;
                            break;
                        }
                    }
                    double precursormz = Double.parseDouble(temp.substring(startidx, stopidx));
                    //By default, assuming it's 5600 data, 
                    //and assume the precursor m/z is at 0.25 * window size Da to the lower bound of isolation window
                    double Loffset = (this.windowSize + 1) * 0.2d;
                    double Roffset = (this.windowSize + 1) * 0.8d;
                                
                    //If the scan contains "windowWideness", then it is a Thermo data, overwrite the isolation window ranges
                    if (temp.contains("windowWideness=\"")) {
                    	startidx = temp.indexOf("windowWideness=\"") + 16;
                        stopidx = startidx + 1;
                        for (int i = startidx + 1; i < temp.length(); i++) {
                        	if (temp.charAt(i) == '\"') {
                            	stopidx = i;
                            	break;
                            }
                        }
                        double windowwideness = Double.parseDouble(temp.substring(startidx, stopidx));
                        //Again assume the precursor m/z is at the center of isolation window, because it is a Thermo data
                        Loffset = windowwideness / 2f;
                        Roffset = windowwideness / 2f;
                   }
//                   Range tmpwindowRange = Range.closed(precursormz - Loffset, precursormz + Roffset);
//                    System.out.println(String.valueOf(precursormz - Loffset));
                    String tmpwindowRange = String.valueOf(precursormz - Loffset).split("[.]")[0]+"-"+String.valueOf(precursormz + Roffset).split("[.]")[0];
                   if (!this.ms2WindowsRange.containsKey(tmpwindowRange)) {
                	   ArrayList<Integer> tmpScanList = new ArrayList();
                	   this.ms2WindowsRange.put(tmpwindowRange, tmpScanList);
                   }
                   this.ms2WindowsRange.get(tmpwindowRange).add(scanno);
                           
                   }
                } else {
                    Logger.getRootLogger().error("index of mzXML error");
                    System.exit(1);
                }
                this.rttoscan.put(rt, scanno);
                this.scantort.put(scanno,rt);
                this.msLevels.put(scanno,mslevel);
            }
            fileHandler.close();       
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		for(int ScanNum : this.msLevels.keySet()){
            if(this.msLevels.get(ScanNum)==1){
                this.ms1ScansNum.add(ScanNum);
            }
//            if(this.msLevels.get(ScanNum)==2){
//                this.ms2ScansNum.add(ScanNum);
//            }
        }
        
	}
	
	public Scan[] readLCMSData(ArrayList mslevels) {
		Scan[] msScan = new Scan[mslevels.size()];
		
       Iterator<Entry<Integer, Long>> iter = this.ScanIndex.entrySet().iterator();        
       Entry<Integer, Long> ent = iter.next(); // first element of scanIndex
       long currentIdx = ent.getValue();
       int nextScanNo = ent.getKey();
       final RandomAccessFile fileHandler;
       try{
    	   fileHandler = new RandomAccessFile(filename, "r");
       }
       catch(FileNotFoundException e){
    	   throw new RuntimeException(e);
       }
       byte[] buffer = new byte[1<<10];
       BitSet labelMsLevel = new BitSet();
       for(int i = 0; i < mslevels.size(); i++) {
    	   int scannum = (int) mslevels.get(i);
    	   labelMsLevel.set(scannum,true);
       }
       int i = 0;
       while (iter.hasNext()) {
           ent = iter.next();
           long startposition = currentIdx;
           long nexposition = ent.getValue();
           int currentScanNo = nextScanNo;
           nextScanNo = ent.getKey();
           currentIdx = nexposition;

           if (labelMsLevel.get(currentScanNo)) {
        	   try { 
        		   final int bufsize =  (int) (nexposition - startposition);
        		   if(buffer.length<bufsize)
        			   buffer = new byte[Math.max(bufsize,buffer.length<<1)];
        		   fileHandler.seek(startposition);
        		   fileHandler.read(buffer, 0, bufsize);
        		   String xmltext = new String(buffer,0,bufsize,StandardCharsets.ISO_8859_1);
        		   if (ent.getKey() == Integer.MAX_VALUE) {
        			   xmltext = xmltext.replaceAll("</msRun>", "");
        		   }
        		   Scan tmpScan = readScan(xmltext);
        		   msScan[i] = tmpScan;
        		   i++;
        		   tmpScan = null;

        	   }
               catch (Exception e) {
            	   e.printStackTrace();
               }

           }
       }
       return msScan;       
	}
	
	public void outputCSV(ArrayList<IsotopicTraces> resultTraces) {
 	   String csvFile = FilenameUtils.getFullPath(this.filename)+"gridmasstarget.csv";
 	   try {
 		   BufferedWriter bw2 = Files.newBufferedWriter(Paths.get(csvFile), StandardCharsets.ISO_8859_1);
 		   bw2.append("mz,rt,height,areaa,startrt,endrt,minmz,maxmz,startscan,endscan\n");
 		   for(IsotopicTraces oneTrace:resultTraces) {
 			   bw2.append(String.valueOf(oneTrace.getTargetMZ())+","+String.valueOf(oneTrace.getRT())+","+String.valueOf(oneTrace.getHeight())+","+String.valueOf(oneTrace.getArea())+","+String.valueOf(oneTrace.getStartRT())+","+String.valueOf(oneTrace.getEndRT())+","+String.valueOf(oneTrace.getMinMZ())+","+String.valueOf(oneTrace.getMaxMZ())+","+String.valueOf(rttoscan.get(oneTrace.getStartRT()))+","+String.valueOf(rttoscan.get(oneTrace.getEndRT()))+"\n");
 		   }
 		   bw2.close();
 	   }
 	   catch (Exception e) {
		e.printStackTrace();
 	   }
    }
//	public ArrayList<IsotopicTraces> getResult() {
//		return this.resultTraces;
//	}
	
	
}