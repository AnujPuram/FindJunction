/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.span.SimpleSeqSpan;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symmetry.BAMSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import com.affymetrix.genometryImpl.symmetry.SimpleMutableSeqSymmetry;
import com.affymetrix.genometryImpl.util.SeqUtils;
import findjunction.filters.ChildThresholdFilter;
import findjunction.filters.NoIntronFilter;
import findjunction.filters.DuplicateSymFilter;
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
    private static boolean twoTracks;
    private static DuplicateSymFilter duplicateSymFilter = new DuplicateSymFilter();
    public FindJunction() {
    }
    
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        String home = System.getProperty("user.home");
        fJ.setThreshold(5);
        fJ.setTwoTracks(false);
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
                else if(args[i].equals("-t")){
                    if(args[i+1].equals("true"))
                        fJ.setTwoTracks(true);
                    else if(args[i+1].equals("false"))
                        fJ.setTwoTracks(false);
                    else{
                        System.out.println("Invalid input for -t option");
                        break;
                    }
                }
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
       FindJunctionOperator operator = new FindJunctionOperator();
       URI uri = new URI(input);
       BAM bam = new BAM(uri,"small_hits",new AnnotatedSeqGroup("small_hits")); 
       OutputStream os  = new FileOutputStream(output);
       DataOutputStream dos = new DataOutputStream(os);
       List<BioSeq> list = bam.getChromosomeList();
       BedParser parser = new BedParser();
       SeqSpan currentSpan;
       List<SeqSymmetry> uniqueSyms = new ArrayList<SeqSymmetry>();
       List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
       int off = 50000;
       for(int i=0;i<list.size();i++){
           for(int j=list.get(i).getMin(); j < list.get(i).getMax(); j= j+off){
               int start =j;
               int end;
               if((start + off) < list.get(i).getMax())
                   end = start + off;
               else
                   end = list.get(i).getMax();
               currentSpan = new SimpleSeqSpan(start, end, list.get(i));
               List<SeqSymmetry> syms = bam.getRegion(currentSpan);
               if(syms.size()>0){
                   duplicateSymFilter.setParam(start);
                   operator.setFilter(duplicateSymFilter);
                   SeqSymmetry container =  operator.operate(list.get(i), syms);
                   int children = container.getChildCount();
                   for(int k=0;k<children;k++){
                       junctions.add(container.getChild(k));
                   }
               }
               parser.writeBedFormat(dos, junctions, list.get(i));
               junctions.clear();
          }
      }
   }    
    
    public int getThreshold(){
        return threshold;
    }
    
    public void setThreshold(int threshold){
        this.threshold = threshold;
    }
    
    public void setTwoTracks(boolean twoTracks){
        this.twoTracks = twoTracks;
    }
    
    public boolean isTwoTracks(){
        return twoTracks;
    }
}