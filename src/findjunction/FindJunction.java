/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.span.SimpleMutableSeqSpan;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symloader.BAM.SeqSymmetryIterator;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import com.affymetrix.genometryImpl.util.SynonymLookup;
import com.affymetrix.igb.IGB;
import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Anuj
 */
public class FindJunction {

    /**
     * @param args the command line arguments
     */
    private static final int default_threshold = 5; // maybe make all caps
    float totalLength = 0,currentLength = 0;
    int absPercentage = 0;
    public static boolean DEBUG = true;
    
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        File directory = new File(".");
        String home = System.getProperty("user.home");
        int threshold = default_threshold;
        boolean twoTracks = false;
        boolean uniqueness = false;
        String input = args[args.length-1];
        String last_argument_message = "Last argument should be a comma-separated list of one or more BAM files";
        if(input.lastIndexOf(".") < 0){
            System.err.println(last_argument_message);
            System.exit(1);   
        }
        if(!(input.substring(input.lastIndexOf(".")+1, input.length()).equals("bam"))){
            System.err.println(last_argument_message);
            System.exit(1);
        }
        String output = "";
        String twoBit = "";
        String thresh = getArg("-n", args);
        String unique = getArg("-u", args);
        String debug_option = getArg("-d",args);
        if (debug_option != null) {
            DEBUG = true;
        }
        
        if(thresh != null)
            threshold = Integer.parseInt(thresh);
        if(unique != null)
            uniqueness = true;
        output = getArg("-o", args);
        if((getArg("-s", args) != null) && (getArg("-b" , args) != null)){
            System.err.println("Both -s and -b cannot be given together");
            return;
        }
        if((getArg("-s", args) == null) && (getArg("-b", args) == null)){
            System.err.println("Provide either -s or -b option to decide strands");
            return;
        }
        twoBit = getArg("-b", args);
        if(getArg("-s", args) != null){
            twoTracks = true;
            twoBit = null;
        }        
        fJ.init(input, output, threshold, twoTracks, twoBit, uniqueness);
    }
    
    public static String getArg(String label, String[] args) {
        String to_return = null;
        boolean got_it = false;
        if (label != null && args != null) {
            for (String item : args) {
                if (got_it) {
                    to_return = item;
                    break;
                }
                if (item.equals(label)) {
                    got_it = true;
                }
            }
        }
        if (got_it && to_return == null) {
            to_return = "true";
        }
        return to_return;
    }
    
    //This is the method where the control of the program gets started
    public void init(String input, String output, int threshold, boolean twoTracks, String twoBit, boolean uniqueness) throws URISyntaxException, Exception{
        File inputFile, twoBitFile;
        URI inputURI, twoBitURI = null; 
        if(!(input.startsWith("file:") && !(input.startsWith("http:")) && !(input.startsWith("ftp:")))){
            inputFile = new File(input);
            inputFile = inputFile.getAbsoluteFile();
            inputURI = inputFile.toURI();
        }
        else
            inputURI = new URI(input);
        if(twoBit != null){
            if(!(twoBit.startsWith("file:") && !(twoBit.startsWith("http:")) && !(twoBit.startsWith("ftp:")))){
                twoBitFile = new File(twoBit);
                twoBitFile = twoBitFile.getAbsoluteFile();
                twoBitURI = twoBitFile.toURI();
            }
            else
                twoBitURI = new URI(twoBit);
                
        }
        if(output != null){
            File outputFile = new File(output);
            outputFile = outputFile.getAbsoluteFile();
            output = outputFile.getAbsolutePath();
        }
        convertBAMToBed(inputURI , output, threshold, twoTracks, twoBitURI, uniqueness);        
    }
    
    //Takes BAM file in the given path as an input and filters it with the Simple Filter Class
    public void convertBAMToBed(URI input , String output, int threshold, boolean twoTracks, URI twoBit, boolean uniqueness) throws URISyntaxException, Exception{
        URI uri = input;
        URI twoBitURI;
        String twoBitFileName = null;
        TwoBit twoBitFile = null;
        InputStream isreader = IGB.class.getResourceAsStream("/chromosomes.txt");
        SynonymLookup.getChromosomeLookup().loadSynonyms(isreader, true) ;
        AnnotatedSeqGroup group = new AnnotatedSeqGroup("Find Junctions");
        String fileName = uri.toString().substring(uri.toString().lastIndexOf("/")+1, uri.toString().lastIndexOf("."));
        BAM bam = new BAM(uri,fileName,group);
        if(twoBit != null){
            twoBitURI = twoBit;
            twoBitFileName = twoBitURI.toString().substring(twoBitURI.toString().lastIndexOf("/")+1, twoBitURI.toString().lastIndexOf("."));
            twoBitFile = new TwoBit(twoBitURI, twoBitFileName, group);
        }
        FindJunctionOperator operator = new FindJunctionOperator(threshold, twoTracks, twoBitFile, uniqueness);
        List<BioSeq> list = bam.getChromosomeList();
        OutputStream os;
        DataOutputStream dos;
        if(output != null){    
            os  = new FileOutputStream(output);
            dos = new DataOutputStream(os);
        }
        else{
            dos = new DataOutputStream(System.out);
            os = null;
        }
        for(BioSeq bioseq : list)
            totalLength = totalLength + bioseq.getLength(); 
        for(BioSeq bioSeq : list)
            writeJunctions(bioSeq, bam, uri, twoBitFile, operator, dos);
        if(isreader  != null)
            isreader.close();
        if(dos != null && os != null){
            dos.close();
            os.close();
        }
    }    
    
    public void writeJunctions(BioSeq bioseq, BAM bam, URI bamURI, TwoBit twoBitFile, FindJunctionOperator operator, DataOutputStream dos) throws FileNotFoundException, Exception {
        SeqSpan currentSpan = new SimpleMutableSeqSpan(bioseq.getMin(), bioseq.getMax(), bioseq);
        List<SeqSymmetry> syms = new ArrayList<SeqSymmetry>();
        SeqSymmetryIterator iter = bam.getIterator(bioseq, bioseq.getMin(), bioseq.getMax(), false);
        if(twoBitFile != null)
            operator.setResidueString(twoBitFile.getRegionResidues(currentSpan));
        if (iter != null) {
            System.err.print(bioseq.getID()+": ");
            int currentProgress = (int)(iter.getProgress()*100);
            int prevProgress = currentProgress;
            while (iter.hasNext()) {
                syms.add(iter.next());
                if (syms.size() >= operator.offset) {
                    write(bioseq, syms, operator, dos);
                    currentProgress = (int)(iter.getProgress()*100);
                    for(int i=0; i<currentProgress - prevProgress; i++){
                        System.err.print("|");
                    }
                    prevProgress = currentProgress;
                    syms.clear();
                }
            }
            write(bioseq, syms, operator, dos);
            for(int i=0; i<currentProgress - prevProgress; i++){
                System.err.print("|");
            }
            System.err.println("100%");
            iter.close();
        }
    }
    
    private void write(BioSeq bioseq, List<SeqSymmetry> syms, FindJunctionOperator operator, DataOutputStream dos) throws IOException {
        List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
        SeqSymmetry container = operator.operate(bioseq, syms);
        int children = container.getChildCount();
        for (int k = 0; k < children; k++) {
            junctions.add(container.getChild(k));
        }
        BedParser.writeBedFormat(dos, junctions, bioseq);
        junctions.clear();
    }
}