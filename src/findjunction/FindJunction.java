/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.operator.FindJunctionOperator;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import com.affymetrix.genometryImpl.util.SynonymLookup;
import com.affymetrix.igb.IGB;
import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Anuj
 */
public class FindJunction {

    /**
     * @param args the command line arguments
     */
    private static final int DEFAULT_THRESHOLD = 5; // maybe make all caps
    public static boolean DEBUG = false;
    
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        File directory = new File(".");
        String home = System.getProperty("user.home");
        int threshold = DEFAULT_THRESHOLD;
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
        if(DEBUG)
            System.err.println("Initial Heap Memory: "+Runtime.getRuntime().freeMemory());
        URI inputURI, twoBitURI = null; 
        if(!(input.startsWith("file:") && !(input.startsWith("http:")) && !(input.startsWith("ftp:")))){
            inputURI = relativeToAbsolute(input);
        }
        else
            inputURI = new URI(input);
        if(twoBit != null){
            if(!(twoBit.startsWith("file:") && !(twoBit.startsWith("http:")) && !(twoBit.startsWith("ftp:")))){
                twoBitURI = relativeToAbsolute(twoBit);
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
    
    public URI relativeToAbsolute(String path){
        File tempFile = new File(path);
        tempFile = tempFile.getAbsoluteFile();
        return tempFile.toURI();
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
        FindJunctionOperator operator = new FindJunctionOperator();
        HashMap<String, Object> paraMeters = new HashMap<String, Object>();
        paraMeters.put("threshold", threshold);
        paraMeters.put("twoTracks", twoTracks);
        paraMeters.put("uniqueness", uniqueness);
        operator.setParameters(paraMeters);
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
        WriteJunctionsThread thread = new WriteJunctionsThread(bam, twoBitFile, operator, dos, DEBUG);
        for(BioSeq bioSeq : list){
            thread.run(bioSeq);
        }
        if(isreader  != null)
            isreader.close();
        if(dos != null && os != null){
            dos.close();
            os.close();
        }
    }    
}