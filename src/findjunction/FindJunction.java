/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import findjunction.filters.SimpleFilter;
import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import findjunction.filters.ThresholdFilter;
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
    private static int threshold;
    public FindJunction() {
    }
    
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        String home = System.getProperty("user.home");
        if(args.length == 2){
            if(args[1].equals("-n")){
                fJ.setThreshold(5);
                fJ.init(args[0],home+"/test2.bed");
            }
            else
                fJ.init(args[0],args[1]);
        }
        else if(args.length == 3){
            if(args[1].equals("-n")){
                fJ.setThreshold(Integer.parseInt(args[2]));
                fJ.init(args[0],"/Users/auser/Desktop/test2.bed");
            }
            else{
                fJ.setThreshold(5);
                fJ.init(args[0], args[1]);
            }
        }
        else if(args.length == 4){
            fJ.setThreshold(Integer.parseInt(args[3]));
            fJ.init(args[0],args[1]);
        }
        else if(args.length == 1){
            fJ.setThreshold(5);
            fJ.init(args[0],"/Users/auser/Desktop/test2.bed");
        }
        else{
            fJ.setThreshold(5);
            fJ.init("/Users/auser/Desktop/genoviz/branches/igb_6_4/genometryImpl/test/data/bam/small_hits.bam", "/Users/auser/Desktop/test2.bed");
        }
    }
    
    //This is the method where the control of the program gets started
    public void init(String input, String output) throws URISyntaxException, Exception{
        if(input.startsWith("file:") || input.startsWith("http:") || input.startsWith("ftp:"))
            convertBAMToBed(input , output);
        else
            convertBAMToBed("file:"+input , output);
        
    }
    
    //Takes BAM file in the given path as an input and filters it with the Simple Filter Class
    public void convertBAMToBed(String input , String output) throws URISyntaxException, Exception{
       ThresholdFilter filter = new ThresholdFilter();
       URI uri = new URI(input);
       BAM bam = new BAM(uri,"small_hits",new AnnotatedSeqGroup("small_hits")); 
       OutputStream os  = new FileOutputStream(output);
       DataOutputStream dos = new DataOutputStream(os);
       List<BioSeq> list = bam.getChromosomeList();
       BedParser parser = new BedParser();
       List<SeqSymmetry> result = new ArrayList<SeqSymmetry>();
       for(int i=0;i<list.size();i++){
            System.out.println(list.get(i).getID());
            List<SeqSymmetry> syms = bam.getChromosome(list.get(i));
            for(SeqSymmetry sym : syms)
                if(filter.filterSymmetry(list.get(i), sym))
                    result.add(sym);
             parser.writeBedFormat(dos, result, list.get(i));
       }    
    }
    
    public int getThreshold(){
        return threshold;
    }
    
    public void setThreshold(int threshold){
        this.threshold = threshold;
    }
}
