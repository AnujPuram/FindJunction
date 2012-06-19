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
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author auser
 */
public class FindJunction {

    /**
     * @param args the command line arguments
     */
    public FindJunction(){
    }
    public static void main(String[] args)throws FileNotFoundException,IOException, URISyntaxException, Exception {
        FindJunction fJ = new FindJunction();
        SimpleFilter filter = new SimpleFilter();
        URI uri = new URI("file:"+args[0]);
        BAM bam = new BAM(uri,"small_hits",new AnnotatedSeqGroup("small_hits"));
        OutputStream os  = new FileOutputStream("/Users/auser/Desktop/test2.bed");
        DataOutputStream dos = new DataOutputStream(os);
        BedParser parser = new BedParser();
        List<BioSeq> list = bam.getChromosomeList();
        for(int i=0;i<list.size();i++){
            System.out.println(list.get(i).getID());
            List<SeqSymmetry> syms = bam.getChromosome(list.get(i));
            List<SeqSymmetry> result = new ArrayList<SeqSymmetry>();
            for(SeqSymmetry sym : syms)
                if(filter.filterSymmetry(list.get(i), sym))
                    result.add(sym);
            parser.writeBedFormat(dos, result, list.get(i));
            for(SeqSymmetry sym : syms){
                SeqSpan span = sym.getSpan(list.get(i));
                System.out.println("Child Count: "+sym.getChildCount());
                
                System.out.println("Start: "+span.getStart());
                System.out.println("End: "+span.getEnd());
                for(int j=0;j<sym.getChildCount();j++){
                    SeqSymmetry child = sym.getChild(j);
                    span = child.getSpan(list.get(i));
                    System.out.println("Child "+(j+1)+"Start: "+span.getStart());
                    System.out.println("Child "+(j+1)+" End: "+span.getEnd());
                    System.out.println("Length of Child "+(j+1)+":"+(span.getEnd()-span.getStart()));
                }
               System.out.println();
            }
        }
        // TODO code application logic here
    }
}
