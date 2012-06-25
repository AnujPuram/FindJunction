/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.AnnotatedSeqGroup;
import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
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
    private static final int offset = 100000; 
    public FindJunction() {
    }
    
    
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        String home = System.getProperty("user.home");
        int threshold = 5;
        boolean twoTracks = false;
        String input = "";
        String output = home+"/test2.bed";
        if(args.length % 2 == 0 && args.length>0){
            for(int i=0;i<args.length;i=i+2){
                if(args[i].equals("-n"))
                    threshold = Integer.parseInt(args[i+1]);
                else if(args[i].equals("-i"))
                    input = args[i+1];
                else if(args[i].equals("-o"))
                    output = args[i+1];
                else if(args[i].equals("-t")){
                    if(args[i+1].equals("true"))
                        twoTracks = true;
                    else if(args[i+1].equals("false"))
                        twoTracks = false;
                    else{
                        System.out.println("Invalid input for -t option");
                        break;
                    }
                }
            }
            if(!input.equals(""))
                fJ.init(input, output, threshold, twoTracks);
            else
                System.out.println("Please give the input file path");
        }
        else
            System.out.println("Invalid Number of Arguments");
        
       }
    
    //This is the method where the control of the program gets started
    public void init(String input, String output, int threshold, boolean twoTracks) throws URISyntaxException, Exception{
        if(input.startsWith("file:") || input.startsWith("http:") || input.startsWith("ftp:"))
            convertBAMToBed(input , output, threshold, twoTracks);
        else
            convertBAMToBed("file:"+input , output, threshold, twoTracks);
        
    }
    
    //Takes BAM file in the given path as an input and filters it with the Simple Filter Class
    public void convertBAMToBed(String input , String output, int threshold, boolean twoTracks) throws URISyntaxException, Exception{
        URI uri = new URI(input);
        FindJunctionOperator operator = new FindJunctionOperator(threshold, twoTracks);
        BAM bam = new BAM(uri,"small_hits",new AnnotatedSeqGroup("small_hits")); 
        List<BioSeq> list = bam.getChromosomeList();
        BedParser parser = new BedParser();
        OutputStream os  = new FileOutputStream(output);
        DataOutputStream dos = new DataOutputStream(os);
        for(BioSeq bioSeq : list)
            writeJunctions(bioSeq, parser, bam, operator, dos);
        dos.close();
        os.close();
    }    
    
    public void writeJunctions(BioSeq bioseq, BedParser parser, BAM bam, FindJunctionOperator operator, DataOutputStream dos) throws FileNotFoundException, Exception{
        SymmetryFilterI duplicateSymFilter = new DuplicateSymFilter();
        SeqSpan currentSpan;
        List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
        for(int j=bioseq.getMin(); j < bioseq.getMax(); j= j+offset){
            int start =j;
            int end;
            if((start + offset) < bioseq.getMax())
                end = start + offset;
            else
                end = bioseq.getMax();
            currentSpan = new SimpleSeqSpan(start, end, bioseq);
            List<SeqSymmetry> syms = bam.getRegion(currentSpan);
            if(syms.size()>0){
                duplicateSymFilter.setParam(start);
                operator.setFilter(duplicateSymFilter);
                SeqSymmetry container =  operator.operate(bioseq, syms);
                int children = container.getChildCount();
                for(int k=0;k<children;k++)
                    junctions.add(container.getChild(k));
                parser.writeBedFormat(dos, junctions, bioseq);
                junctions.clear();
            }            
        }
    }
}