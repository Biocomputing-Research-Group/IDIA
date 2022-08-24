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
package MSUmpire.DIA;

import MSUmpire.BaseDataStructure.InstrumentParameter;
import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
//import MSUmpire.BaseDataStructure.XYPointCollection;
import MSUmpire.LCMSPeakStructure.LCMSPeakDIAMS2;
import MSUmpire.LCMSPeakStructure.LCMSPeakMS1;
//import MSUmpire.PeakDataStructure.PeakCluster;
//import MSUmpire.PeakDataStructure.PrecursorFragmentPairEdge;
import MSUmpire.SpectrumParser.DIA_Setting;
import MSUmpire.SpectrumParser.SpectrumParserBase;
//import MSUmpire.SpectrumParser.mzXMLParser;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.zip.DataFormatException;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.derby.tools.sysinfo;
import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import GridMass.module.GridMassModel;
/**
 * Main data structure which represents a DIA file
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIAPack {

    private InstrumentParameter parameter;
    public String Filename;
    private int NoCPUs = 4;
    public LCMSPeakMS1 MS1FeatureMap;
    public ArrayList<LCMSPeakDIAMS2> DIAWindows;
    public DIA_Setting dIA_Setting = new DIA_Setting();
    private SpectrumParserBase SpectrumParser;
//    public LCMSID IDsummary;
    public HashMap<Integer, Integer> ScanClusterMap_Q1;
    public HashMap<Integer, Integer> ScanClusterMap_Q2;
    public HashMap<Integer, String> ScanClusterMap_Q3;
//    public ArrayList<String> iProphPepXMLs = new ArrayList<String>();
    public boolean ExportPrecursorPeak = false;
    public boolean ExportFragmentPeak = false;
    public HashMap<Integer, Double> FactorialTable;
//    TargetMatchScoring TScoring;
    public boolean Resume;
    
    public int Q1Scan = 0;
    public int Q2Scan = 0;
    public int Q3Scan = 0;
  
    //Whether to use IDs from targeted re-extraction for quantification analysis
    public boolean UseMappedIon = false;
    //Whether to use probability threshold for targeted re-extraction to filter IDs
    public boolean FilterMappedIonByProb = true;
    //Probability threshold for targeted re-extraction
    public float MappedIonProbThreshold = 0.95f;

    
    private void RemoveIndexFile() {
        new File(FilenameUtils.removeExtension(Filename) + ".ScanPosFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".RTidxFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".ScanRTFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".ScanidxFS").delete();
        new File(FilenameUtils.removeExtension(Filename) + ".DIAWindowsFS").delete();
    }
    
    public void FixScanidx() {
        RemoveIndexFile();
        GetSpectrumParser();
    }

    public DIAPack(String Filename, int NoCPUs) throws FileNotFoundException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, DataFormatException {
        super();
        this.Filename = Filename;
        this.NoCPUs = NoCPUs;
    }

    public String GetBaseName() {
        return FilenameUtils.getBaseName(Filename);
    }

    public void SetNoCPUs(int cpu) {
        this.NoCPUs = cpu;
    }

    public void SetParameter(InstrumentParameter parameter) {
        this.parameter = parameter;
    }

    public void SetDataType(SpectralDataType.DataType datatype) {
        this.dIA_Setting.dataType = datatype;
    }

    public void AddVariableWindow(XYData window) {
        if (dIA_Setting.DIAWindows == null) {
            dIA_Setting.DIAWindows = new TreeMap<>();
        }
        dIA_Setting.DIAWindows.put(window, new ArrayList<Integer>());
    }
    
    public String GetQ1Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q1";
        return FilenameUtils.getBaseName(Filename) + "_Q1";
    }

    public String GetQ2Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q2";
        return FilenameUtils.getBaseName(Filename) + "_Q2";
    }

    public String GetQ3Name() {
        //return FilenameUtils.getBaseName(ScanCollectionName) + ".Q3";
        return FilenameUtils.getBaseName(Filename) + "_Q3";
    }


    public void SetWindowSize(float size) {
        dIA_Setting.F_DIA_WindowSize = size;
    }

    public SpectrumParserBase GetSpectrumParser() {
        if (SpectrumParser == null) {
            try {
                SpectrumParser = SpectrumParserBase.GetInstance(Filename, parameter, dIA_Setting.dataType, dIA_Setting, NoCPUs);
            } catch (Exception ex) {
                Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
                Logger.getRootLogger().error("Read spectral file:" + Filename + " failed.");
                System.exit(2);
            }
            dIA_Setting = SpectrumParser.dIA_Setting;
//            SaveDIASetting();            
        }
        return SpectrumParser;
    }
 
    //Entry of processing of signal extraction and generating pseudo MS/MS spectra
    public void process() throws SQLException, IOException, InterruptedException, ExecutionException, ParserConfigurationException, SAXException, FileNotFoundException, DataFormatException, Exception {
   
    	GridMassModel gridmass = new GridMassModel(parameter,Filename);
    	BuildDIAWindows();
        MS1PeakDetection(gridmass);
        DIAMS2PeakDetection(gridmass);
    }

    //Building empty data structure for MS2 feature maps 
    public void BuildDIAWindows() throws IOException, DataFormatException, IOException, IOException, IOException, InterruptedException {
        if (dIA_Setting.DIAWindows == null || dIA_Setting.DIAWindows.isEmpty()) {
            GetSpectrumParser();
        }
        DIAWindows = new ArrayList<>();
        Object[] WindowRange = dIA_Setting.DIAWindows.keySet().toArray();
        for (int i = 0; i < WindowRange.length; i++) {
            XYData DiaWinMz = (XYData) WindowRange[i];
            XYData LastWinMz = null;
            if (i < WindowRange.length - 1) {
                LastWinMz = (XYData) WindowRange[i + 1];
            }
            
            LCMSPeakDIAMS2 diawindow = new LCMSPeakDIAMS2(Filename, this, parameter, DiaWinMz, LastWinMz, GetSpectrumParser(), NoCPUs);
            diawindow.Resume=Resume;
            //pass the settings through to MS2 feature map
            diawindow.datattype = dIA_Setting.dataType;
            diawindow.ExportPeakCurveTable = ExportFragmentPeak;
            diawindow.ExportPeakClusterTable = ExportPrecursorPeak;
            DIAWindows.add(diawindow);
        }
    }

   

        
    public void ClearStructure() {
        MS1FeatureMap = null;
        DIAWindows = null;
        dIA_Setting = null;
        SpectrumParser = null;
    }


    //Perform MS1 feature detection
    private void MS1PeakDetection(GridMassModel gridmass) throws SQLException, InterruptedException, ExecutionException, IOException, ParserConfigurationException, SAXException, FileNotFoundException, Exception {
        //Remove existing pseudo MS/MS MGF files
        RemoveMGF();
                
        MS1FeatureMap = new LCMSPeakMS1(Filename, NoCPUs);
        MS1FeatureMap.datattype = dIA_Setting.dataType;
        MS1FeatureMap.SetParameter(parameter);
        MS1FeatureMap.Resume=Resume;
        //Assign MS1 feature maps
        MS1FeatureMap.SetMS1Windows(dIA_Setting.MS1Windows);
//        MS1FeatureMap.CreatePeakFolder();
        MS1FeatureMap.ExportPeakCurveTable = false;
        
        MS1FeatureMap.SetSpectrumParser(GetSpectrumParser());
        Logger.getRootLogger().info("Processing MS1 peak detection");        
        MS1FeatureMap.ExportPeakClusterTable = ExportPrecursorPeak;
        //Start MS1 feature detection
        MS1FeatureMap.PeakClusterDetection(gridmass,0.0f,0.0f);
        gridmass = null;

        Logger.getRootLogger().info("==================================================================================");
    }

    //Remove pseudo MS/MS spectra MGF files
    private void RemoveMGF() {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf";
        File file = new File(mgffile);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.delete();
        }
        mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf.temp";
        mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf.temp";
        mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        file = new File(mgffile);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.delete();
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.delete();
        }
    }

    public static enum OutputFile{
        ScanClusterMapping_Q1,
        Mgf_Q1,
        ScanClusterMapping_Q2,
        Mgf_Q2,
        ScanClusterMapping_Q3,
        Mgf_Q3;
    }
    private static final EnumMap<OutputFile, BufferedWriter> file_handlers = new EnumMap<>(OutputFile.class);

    public static BufferedWriter get_file(final OutputFile f,final String filename) {
        final BufferedWriter bw = DIAPack.file_handlers.get(f);
        if (bw == null) {
            final BufferedWriter bw2;
            try {
                bw2 = Files.newBufferedWriter(Paths.get(filename), StandardCharsets.ISO_8859_1);
            } catch (final IOException ex) {
                throw new UncheckedIOException(ex);
            }
            DIAPack.file_handlers.put(f, bw2);
            return bw2;
        }
        return bw;
    }

    //Perform MS2 fragment feature detection
    public void DIAMS2PeakDetection(GridMassModel gridmass) throws SQLException, IOException, InterruptedException, ExecutionException, FileNotFoundException, Exception {
        int count = 1;
        //CreateSWATHTables();
        for (LCMSPeakDIAMS2 DIAwindow : DIAWindows) {
            Logger.getRootLogger().info("Processing DIA MS2 (mz range):" + DIAwindow.DIA_MZ_Range.getX() + "_" + DIAwindow.DIA_MZ_Range.getY() + "( " + (count++) + "/" + GetSpectrumParser().dIA_Setting.DIAWindows.size() + " )");
            DIAwindow.ExportPeakCurveTable = ExportFragmentPeak;
            DIAwindow.ExportPeakClusterTable = ExportPrecursorPeak;
            DIAwindow.PeakDetectionPFGrouping(MS1FeatureMap, gridmass, DIAwindow.DIA_MZ_Range.getX(), DIAwindow.DIA_MZ_Range.getY());
            
            DIAwindow.ClearAllPeaks();
            Logger.getRootLogger().info("==================================================================================");
        }
        gridmass = null;
        for(final BufferedWriter bw : DIAPack.file_handlers.values())
            bw.close();
        RenameMGF("");
//        Renameft2("");
        convertMGF("");
    }

    private void RenameMGF(String tag) {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".mgf.temp";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".mgf.temp";
//        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        File file = new File(mgffile);

        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
//        file = new File(mgffile3);
//        if (file.exists()) {
//            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
//        }
    }
    
    private void Renameft2(String tag) {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + ".ft2.temp";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + ".ft2.temp";
        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + ".mgf.temp";
        File file = new File(mgffile);

        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".ft2.temp", tag + ".ft2")));
        }
        file = new File(mgffile2);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".ft2.temp", tag + ".ft2")));
        }
        file = new File(mgffile3);
        if (file.exists()) {
            file.renameTo(new File(file.getAbsolutePath().replace(".mgf.temp", tag + ".mgf")));
        }
    }

    private void convertMGF(String tag) throws java.security.NoSuchAlgorithmException, IOException {
        String mgffile = FilenameUtils.getFullPath(Filename) + GetQ1Name() + tag + ".mgf";
        String mgffile2 = FilenameUtils.getFullPath(Filename) + GetQ2Name() + tag + ".mgf";
//        String mgffile3 = FilenameUtils.getFullPath(Filename) + GetQ3Name() + tag + ".mgf";
        {
            final Path file = Paths.get(mgffile);
            if (Files.exists(file)) {
                final String fn = file.getFileName().toString();
                try (final BufferedWriter bw = Files.newBufferedWriter(file.getParent().resolve(fn.substring(0, fn.length() - 4) + ".mzML"), StandardCharsets.US_ASCII)) {
                    MSUmpire.SpectrumParser.mzXMLParser.to_mzML(file, bw, NoCPUs);
                }
            }
        }
        {
            final Path file = Paths.get(mgffile2);
            if (Files.exists(file)) {
                final String fn = file.getFileName().toString();
                try (final BufferedWriter bw = Files.newBufferedWriter(file.getParent().resolve(fn.substring(0, fn.length() - 4) + ".mzML"), StandardCharsets.US_ASCII)) {
                    MSUmpire.SpectrumParser.mzXMLParser.to_mzML(file, bw, NoCPUs);
                }
            }
        }
//        {
//            final Path file = Paths.get(mgffile3);
//            if (Files.exists(file)) {
//                final String fn = file.getFileName().toString();
//                try (final BufferedWriter bw = Files.newBufferedWriter(file.getParent().resolve(fn.substring(0, fn.length() - 4) + ".mzML"), StandardCharsets.US_ASCII)) {
//                    MSUmpire.SpectrumParser.mzXMLParser.to_mzML(file, bw, NoCPUs);
//                }
//            }
//        }


    }

//    public void SaveParams() {
//        parameter.WriteParamSerialization(Filename);
//    }

//    public void SaveDIASetting() {
//        dIA_Setting.WriteDIASettingSerialization(Filename);
//    }

    public boolean LoadParams() {
        parameter = InstrumentParameter.ReadParametersSerialization(Filename);
        return parameter != null;
    }

    public boolean LoadDIASetting() {
        dIA_Setting = DIA_Setting.ReadDIASettingSerialization(Filename);
        return dIA_Setting != null;
    }

    public InstrumentParameter GetParameter() {
        return parameter;
    }
}
