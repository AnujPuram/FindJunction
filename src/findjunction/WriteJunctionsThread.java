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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author auser
 */
public class WriteJunctionsThread{   
    
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
                operator.setResidueString(twoBitFile.getRegionResidues(currentSpan));
        } catch (Exception ex) { 
            System.err.println("Error1: "+ex.getMessage());
        }
        if (iter != null) {
            while (iter.hasNext()) {
                syms.add(iter.next());
                if (syms.size() >= operator.offset) {
                    if(DEBUG){
                        System.err.println("Available Heap Memory: "+ Runtime.getRuntime().freeMemory());
                    }
                    syms.clear();
                }
            }
            junctions = createJunctions(bioseq, syms, map);
            write(bioseq, junctions);
            iter.close();
        }
    }
    
    private List<SeqSymmetry> createJunctions(BioSeq bioseq, List<SeqSymmetry> syms, HashMap<String, JunctionUcscBedSym> map){
        List<SeqSymmetry> junctions = new ArrayList<SeqSymmetry>();
        operator.subOperate(bioseq, syms, map);
        Collection<JunctionUcscBedSym> junctionValues = map.values(); 
        Object junctionsArray[] = junctionValues.toArray();
        for (int k = 0; k < junctionsArray.length; k++) {
            junctions.add((JunctionUcscBedSym)junctionsArray[k]);
        }
        return junctions;
    }
    private void write(BioSeq bioseq, List<SeqSymmetry> junctions) throws IOException {
        BedParser.writeBedFormat(dos, junctions, bioseq);
        junctions.clear();
    } 
    
}
