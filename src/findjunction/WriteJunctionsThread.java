/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.operator.FindJunctionOperator;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.span.SimpleMutableSeqSpan;
import com.affymetrix.genometryImpl.symloader.BAM;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import com.affymetrix.genometryImpl.symmetry.JunctionUcscBedSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author auser
 */
public class WriteJunctionsThread{   
    
    public static final int SIZE = 30000;
    BAM bam;
    TwoBit twoBitFile;
    FindJunctionOperator operator;
    DataOutputStream dos;
    boolean DEBUG;
    public WriteJunctionsThread(BAM bam, TwoBit twoBitFile, FindJunctionOperator operator, DataOutputStream dos, boolean DEBUG){
        this.bam = bam;
        this.twoBitFile = twoBitFile;
        this.operator = operator;
        this.dos = dos;
        this.DEBUG = DEBUG;
    }

    public void run(BioSeq bioseq) throws IOException {
        SeqSpan currentSpan = new SimpleMutableSeqSpan(bioseq.getMin(), bioseq.getMax(), bioseq);
        List<SeqSymmetry> syms = new ArrayList<SeqSymmetry>();
        List<SeqSymmetry> junctions;
        HashMap<String, JunctionUcscBedSym> map = new HashMap<String, JunctionUcscBedSym>();
        BAM.SeqSymmetryIterator iter = null;
        try {
            iter = bam.getIterator(bioseq, bioseq.getMin(), bioseq.getMax(), false);
            if(twoBitFile != null)
                BioSeq.addResiduesToComposition(bioseq, twoBitFile.getRegionResidues(currentSpan), currentSpan);
        } catch (Exception ex) { 
            System.err.println("Error1: "+ex.getMessage());
        }
        if (iter != null) {
            while (iter.hasNext()) {
                syms.add(iter.next());
                if (syms.size() >= SIZE) {
                    if(DEBUG){
                        System.err.println("Available Heap Memory: "+ Runtime.getRuntime().freeMemory());
                    }
                    operator.subOperate(bioseq, syms, map);
                    syms.clear();
                }
            }
            operator.subOperate(bioseq, syms, map);
            syms.clear();
            System.err.println(bioseq.getID()+": done");
            write(bioseq, new ArrayList<SeqSymmetry>(map.values()));
            map.clear();
            iter.close();
        }
    }
    
    private void write(BioSeq bioseq, List<SeqSymmetry> junctions) throws IOException {
        BedParser.writeBedFormat(dos, junctions, bioseq);
        junctions.clear();
    } 
    
}
