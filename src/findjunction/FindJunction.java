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
import com.affymetrix.genometryImpl.symloader.TwoBit;
import com.affymetrix.genometryImpl.symmetry.BAMSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import com.affymetrix.genometryImpl.symmetry.SimpleMutableSeqSymmetry;
import com.affymetrix.genometryImpl.util.SeqUtils;
import com.affymetrix.genometryImpl.util.SynonymLookup;
import findjunction.filters.ChildThresholdFilter;
import findjunction.filters.NoIntronFilter;
import findjunction.filters.DuplicateSymFilter;
import com.affymetrix.igb.IGB;
import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Anuj
 */
public class FindJunction {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        String home = System.getProperty("user.home");
        int threshold = 5;
        boolean twoTracks = false;
        String input = "";
        String output = "";
        String twoBit = "";
        String thresh = getArg("-n", args);
        if(thresh != null)
            threshold = Integer.parseInt(thresh);
        output = getArg("-o", args);
        twoBit = getArg("-b", args);
        if(getArg("-s", args) != null)
            twoTracks = true;
        for(int i=0;i<args.length;i++){
            if(!args[i].startsWith("-")){
                input = args[i];
                break;
            }
            if(args[i].equals("-n") || args[i].equals("-o") || args[i].equals("-b"))
                i++;
        }
        fJ.init(input, output, threshold, twoTracks, twoBit);
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
    public void init(String input, String output, int threshold, boolean twoTracks, String twoBit) throws URISyntaxException, Exception{
                if(!(input.startsWith("file:") || input.startsWith("http:") || input.startsWith("ftp:")))
            input = "file:"+input;
        if(!(twoBit.startsWith("file:") || twoBit.startsWith("http:") || twoBit.startsWith("ftp:")))
            twoBit = "file:"+twoBit;
        convertBAMToBed(input , output, threshold, twoTracks, twoBit);        
    }
    
    //Takes BAM file in the given path as an input and filters it with the Simple Filter Class
    public void convertBAMToBed(String input , String output, int threshold, boolean twoTracks, String twoBit) throws URISyntaxException, Exception{
        URI uri = new URI(input);
        URI twoBitURI = new URI(twoBit);
        InputStream isreader = IGB.class.getResourceAsStream("/chromosomes.txt");
        SynonymLookup.getChromosomeLookup().loadSynonyms(isreader, true) ;
        AnnotatedSeqGroup group = new AnnotatedSeqGroup("Find Junctions");
        String fileName = uri.toString().substring(uri.toString().lastIndexOf("/")+1, uri.toString().lastIndexOf("."));
        String twoBitFileName = twoBitURI.toString().substring(twoBitURI.toString().lastIndexOf("/")+1, twoBitURI.toString().lastIndexOf("."));
        BAM bam = new BAM(uri,fileName,group);
        TwoBit twoBitFile = new TwoBit(twoBitURI, twoBitFileName, group);
        FindJunctionOperator operator = new FindJunctionOperator(threshold, twoTracks, twoBitFile);
        List<BioSeq> list = bam.getChromosomeList();
        BedParser parser = new BedParser();
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
        for(BioSeq bioSeq : list)
            writeJunctions(bioSeq, parser, bam, twoBitFile, operator, dos);
        if(dos != null && os != null){
            dos.close();
            os.close();
        }
    }    
    
    public void writeJunctions(BioSeq bioseq, BedParser parser, BAM bam, TwoBit twoBitFile, FindJunctionOperator operator, DataOutputStream dos) throws FileNotFoundException, Exception{
        SymmetryFilterI duplicateSymFilter = new DuplicateSymFilter();
        SeqSpan currentSpan;
        List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
        List<BioSeq> seq = twoBitFile.getChromosomeList();
        for(int j=bioseq.getMin(); j < bioseq.getMax(); j= j+operator.offset){
            int start =j;
            int end;
            if((start + operator.offset) < bioseq.getMax())
                end = start + operator.offset;
            else
                end = bioseq.getMax();
            currentSpan = new SimpleSeqSpan(start, end, bioseq);
            List<SeqSymmetry> syms = bam.getRegion(currentSpan);
            //synonyms = lookup.getSynonyms(bioseq.getID());
            if(syms.size()>0){
                duplicateSymFilter.setParam(start);
                operator.setFilter(duplicateSymFilter);
                operator.setResidueString(twoBitFile.getRegionResidues(currentSpan));
                SeqSymmetry container =  operator.operate(bioseq, syms);
                int children = container.getChildCount();
                for(int k=0;k<children;k++){
                    junctions.add(container.getChild(k));
                }
                parser.writeBedFormat(dos, junctions, bioseq);
                junctions.clear();
            }            
        }
    }
}