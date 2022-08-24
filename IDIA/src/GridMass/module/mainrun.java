package GridMass.module;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.io.StringReader;
import java.lang.reflect.Array;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.codec.binary.Base64;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.net.telnet.WindowSizeOptionHandler;
import org.apache.log4j.Logger;
import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import com.google.common.collect.Range;

import org.w3c.dom.Document;
import org.w3c.dom.Node;


import GridMass.datamodule.DataPoint;
import GridMass.datamodule.Scan;
import GridMass.module.GridMass;

public final class mainrun {
	private int windowSize;
	private HashMap<Range, ArrayList<Integer>> windowsRange;
	private ArrayList<Scan> ms1scans;
//	private Scan[] ms2scans;
	private String filename;
	private HashMap<Integer,Double> scantort = new HashMap();
	private HashMap<Double, Integer> rttoscan = new HashMap();
	private TreeMap<Integer, Long> ScanIndex = new TreeMap();
	private HashMap<Integer, Integer> msLevels = new HashMap();
	private IntArrayList ms1ScansNum = new IntArrayList();
	private IntArrayList ms2ScansNum = new IntArrayList();
	public static void main(String[] args) {
//		GridMassModel newgridmass = new GridMassModel();
//		newgridmass.startrun();
	}
}