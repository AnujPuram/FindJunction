/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symmetry.BAMSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import com.affymetrix.genometryImpl.symmetry.SimpleMutableSeqSymmetry;
import com.affymetrix.genometryImpl.util.SeqUtils;
import findjunction.filters.ChildThresholdFilter;
import findjunction.filters.NoIntronFilter;
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
        fJ.setThreshold(5);
        String input = "";
        String output = home+"/test2.bed";
        if(args.length % 2 == 0 && args.length>0){
            for(int i=0;i<args.length;i=i+2){
                if(args[i].equals("-n"))
                    fJ.setThreshold(Integer.parseInt(args[i+1]));
                else if(args[i].equals("-i"))
                    input = args[i+1];
                else if(args[i].equals("-o"))
                    output = args[i+1];
            }
            if(!input.equals(""))
                fJ.init(input,output);
            else
                System.out.println("Please give the input file path");
        }
        else
            System.out.println("Invalid Number of Arguments");
        
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
       NoIntronFilter noIntronFilter = new NoIntronFilter();
       ChildThresholdFilter childFilter = new ChildThresholdFilter();
       CreateBedSymOperator operator = new CreateBedSymOperator();
       URI uri = new URI(input);
       BAM bam = new BAM(uri,"small_hits",new AnnotatedSeqGroup("small_hits")); 
       OutputStream os  = new FileOutputStream(output);
       DataOutputStream dos = new DataOutputStream(os);
       List<BioSeq> list = bam.getChromosomeList();
       BedParser parser = new BedParser();
       for(int i=0;i<list.size();i++){
            System.out.println(list.get(i).getID());
            List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
            List<SeqSymmetry> nonZeroIntronSyms = new ArrayList<SeqSymmetry>();
            List<SeqSymmetry> syms = bam.getChromosome(list.get(i));
            List<SeqSymmetry> result = new ArrayList<SeqSymmetry>();
            if(syms.size()>0){
                for(SeqSymmetry sym : syms){
                    if(noIntronFilter.filterSymmetry(list.get(i), sym))
                        nonZeroIntronSyms.add(sym);
                }
                for(SeqSymmetry sym : nonZeroIntronSyms){
                    if(filter.filterSymmetry(list.get(i), sym))
                        result.add(SeqUtils.getIntronSym(sym, list.get(i)));
                    else{
                        SeqSymmetry intronSym = SeqUtils.getIntronSym(sym, list.get(i));
                        int childCount = sym.getChildCount();
                        List<SeqSymmetry>eligibleIntrons = new ArrayList<SeqSymmetry>();
                        for(int j=0;j<intronSym.getChildCount();j++)
                            eligibleIntrons.add(intronSym.getChild(j));
                        Object array[] = eligibleIntrons.toArray();
                        for(int j=0; j<childCount; j++){
                            if(!(childFilter.filterSymmetry(list.get(i), sym.getChild(j)))){
                                if(j-1 >= 0)
                                    array[j-1] = null;
                                if(j < childCount-1)
                                    array[j] = null;
                            }
                        }
                        ((SimpleMutableSeqSymmetry)intronSym).removeChildren();
                        for(int j=0;j<array.length;j++){
                            if(array[j] != null)
                                ((SimpleMutableSeqSymmetry)intronSym).addChild((SeqSymmetry)array[j]);
                        }
                        eligibleIntrons.clear();
                        if(intronSym.getChildCount() > 0)
                            result.add(intronSym);
                    }
                }
                SeqSymmetry container =  operator.operate(list.get(i), result);
                int children = container.getChildCount();
                for(int j=0;j<children;j++){
                    junctions.add(container.getChild(j));
                }
            }
            parser.writeBedFormat(dos, junctions, list.get(i));
       }    
    }
    
    public int getThreshold(){
        return threshold;
    }
    
    public void setThreshold(int threshold){
        this.threshold = threshold;
    }
}