/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.operator.FindJunctionOperator;
import com.affymetrix.genometryImpl.parsers.BedParser;
import com.affymetrix.genometryImpl.symloader.BAM;
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
    private final BAM bam;
    private final DataOutputStream dos;
    private final boolean DEBUG;
    private final FindJunctionOperator operator;
    
    public WriteJunctionsThread(BAM bam, int threshold, boolean twoTracks, boolean uniqueness, DataOutputStream dos, boolean DEBUG){
        this.bam = bam;
        this.dos = dos;
        this.DEBUG = DEBUG;
        this.operator = new FindJunctionOperator();
        HashMap<String, Object> paraMeters = new HashMap<String, Object>();
        paraMeters.put(FindJunctionOperator.THRESHOLD, threshold);
        paraMeters.put(FindJunctionOperator.TWOTRACKS, twoTracks);
        paraMeters.put(FindJunctionOperator.UNIQUENESS, uniqueness);
        operator.setParameters(paraMeters);
    }

    public void run(BioSeq bioseq) throws IOException {
        List<SeqSymmetry> syms = new ArrayList<SeqSymmetry>();
        HashMap<String, SeqSymmetry> map = new HashMap<String, SeqSymmetry>();
        BAM.SeqSymmetryIterator iter = null;
        try {
            iter = bam.getIterator(bioseq, bioseq.getMin(), bioseq.getMax(), false);
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
