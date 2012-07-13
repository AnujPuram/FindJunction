/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.span.SimpleMutableSeqSpan;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.URI;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author auser
 */
public class WriteJunctionsThread implements Runnable{   
    
    BioSeq bioseq;
    BAM bam;
    URI bamURI;
    TwoBit twoBitFile;
    FindJunctionOperator operator;
    DataOutputStream dos;
    boolean DEBUG;
    public WriteJunctionsThread(BioSeq bioseq, BAM bam, URI bamURI, TwoBit twoBitFile, FindJunctionOperator operator, DataOutputStream dos, boolean DEBUG){
        this.bioseq = bioseq;
        this.bam = bam;
        this.bamURI = bamURI;
        this.twoBitFile = twoBitFile;
        this.operator = operator;
        this.dos = dos;
        this.DEBUG = DEBUG;
    }
    
    @Override
    public void run() {
        SeqSpan currentSpan = new SimpleMutableSeqSpan(bioseq.getMin(), bioseq.getMax(), bioseq);
        List<SeqSymmetry> syms = new ArrayList<SeqSymmetry>();
        List<SeqSymmetry> junctions;
        BAM.SeqSymmetryIterator iter = null;
        try {
            iter = bam.getIterator(bioseq, bioseq.getMin(), bioseq.getMax(), false);
            if(twoBitFile != null)
                operator.setResidueString(twoBitFile.getRegionResidues(currentSpan));
        } catch (Exception ex) { 
            System.err.println("Error1: "+ex.getMessage());
        }
        if (iter != null) {
            System.err.print(bioseq.getID()+": ");
            int currentProgress = (int)(iter.getProgress()*100);
            int prevProgress = currentProgress;
            while (iter.hasNext()) {
                syms.add(iter.next());
                if (syms.size() >= operator.offset) {
                    if(DEBUG){
                        System.err.println("Available Heap Memory: "+ Runtime.getRuntime().freeMemory());
                    }
                    try {
                        junctions = createJunctions(syms);
                        write(junctions);
                    } catch (IOException ex) {
                        System.err.println("Error2 "+ex.getMessage());
                    }
                    currentProgress = (int)(iter.getProgress()*100);
                    for(int i=0; i<currentProgress - prevProgress; i++){
                        System.err.print("|");
                    }
                    prevProgress = currentProgress;
                    syms.clear();
                }
            }
            try {
                junctions = createJunctions(syms);
                write(junctions);
            } catch (IOException ex) {  
                System.err.println("Error3: "+ex.getMessage());
            }
            for(int i=0; i<currentProgress - prevProgress; i++){
                System.err.print("|");
            }
            System.err.println("100%");
            iter.close();
        }
    }
    
    private List<SeqSymmetry> createJunctions(List<SeqSymmetry> syms){
        List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
        SeqSymmetry container = operator.operate(bioseq, syms);
        int children = container.getChildCount();
        for (int k = 0; k < children; k++) {
            junctions.add(container.getChild(k));
        }
        return junctions;
    }
    private void write(List<SeqSymmetry> junctions) throws IOException {
        BedParser.writeBedFormat(dos, junctions, bioseq);
        junctions.clear();
    } 
    
}
