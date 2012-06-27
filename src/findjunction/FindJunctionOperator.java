/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.filter.SymmetryFilterI;
import com.affymetrix.genometryImpl.operator.Operator;
import com.affymetrix.genometryImpl.parsers.FileTypeCategory;
import com.affymetrix.genometryImpl.symloader.TwoBit;
import com.affymetrix.genometryImpl.symmetry.*;
import com.affymetrix.genometryImpl.util.SeqUtils;
import findjunction.filters.ChildThresholdFilter;
import findjunction.filters.DuplicateSymFilter;
import findjunction.filters.NoIntronFilter;
import findjunction.filters.UniqueLocationFilter;
import java.util.*;


/**
 *
 * @author Anuj
 */
public class FindJunctionOperator implements Operator{
    public static final int offset = 100000;
    private static final int default_threshold = 5;
    private static SymmetryFilterI noIntronFilter = new NoIntronFilter();
    private static SymmetryFilterI childThresholdFilter = new ChildThresholdFilter();
    private static SymmetryFilterI duplicateSymFilter = new DuplicateSymFilter();
    private static SymmetryFilterI uniqueLocationFilter = new UniqueLocationFilter();
    private static int threshold = default_threshold;
    private static boolean twoTracks, uniqueness;
    private static TwoBit twoBit;
    private static String residueString;
    public FindJunctionOperator(int threshold, boolean twoTracks, TwoBit twoBit, boolean uniqueness){
        this.threshold = threshold;
        this.twoTracks = twoTracks;
        this.twoBit = twoBit;
        this.uniqueness = uniqueness;
    }   
    
    @Override
    public String getName() {
        return "findJunction";
    }

    @Override
    public String getDisplay() {
        return "Find Junction";
    }

    public void setFilter(SymmetryFilterI filter){
        duplicateSymFilter = filter;
    }
    
    public void setResidueString(String residueString){
        this.residueString = residueString;
    }
    
    @Override
    public SeqSymmetry operate(BioSeq bioseq, List<SeqSymmetry> list) {
        TypeContainerAnnot container = new TypeContainerAnnot("test", "bed");
        HashMap<String, SpecificUcscBedSym> map = new HashMap<String , SpecificUcscBedSym>();
        
        for(SeqSymmetry sym : list){
            if(noIntronFilter.filterSymmetry(bioseq, sym) && duplicateSymFilter.filterSymmetry(bioseq, sym) && 
                    (!uniqueness) || (uniqueness && uniqueLocationFilter.filterSymmetry(bioseq, sym))){
                updateIntronHashMap(sym , bioseq, map);
            }

        }
        Collection<SpecificUcscBedSym> symmetrySet = map.values();
        Object syms[] = symmetrySet.toArray();
        for(int i=0;i<syms.length;i++){
            container.addChild((SpecificUcscBedSym)syms[i]);
        }
        return container;
    }
        

    
    @Override
    public int getOperandCountMin(FileTypeCategory ftc) {
        return ftc == FileTypeCategory.Alignment ? 1 : 0;
    }

    @Override
    public int getOperandCountMax(FileTypeCategory ftc) {
        return ftc == FileTypeCategory.Alignment ? 1 : 0;
    }

    @Override
    public Map<String, Class<?>> getParameters() {
        return null;
    }

    @Override
    public boolean setParameters(Map<String, Object> map) {
        if(map.size() == 1 && map.get(0) instanceof Integer){
            threshold = (Integer)map.get(0);
            return true;
        }
        return false;
    }

    @Override
    public boolean supportsTwoTrack() {
        return true;
    }

    @Override
    public FileTypeCategory getOutputCategory() {
        return FileTypeCategory.Annotation;
    }
    
    //This method splits the given Sym into introns and filters out the qualified Introns
    private static void updateIntronHashMap(SeqSymmetry sym , BioSeq bioseq, HashMap<String, SpecificUcscBedSym> map){
        List<Integer> childIntronIndices = new ArrayList<Integer>();
        int childCount = sym.getChildCount();
        SeqSymmetry intronChild, intronSym;
        childThresholdFilter.setParam(threshold);
        for(int i=0;i<childCount - 1;i++){
            if(childThresholdFilter.filterSymmetry(bioseq, sym.getChild(i)) && childThresholdFilter.filterSymmetry(bioseq, sym.getChild(i+1))){
                childIntronIndices.add(i);
            }
        }
        if(childIntronIndices.size() > 0){
            intronSym = SeqUtils.getIntronSym(sym, bioseq);
            for(Integer i : childIntronIndices){
                intronChild = intronSym.getChild(i);
                addToMap(intronChild, map, bioseq);
            }
        }
    }
    
    private static void addToMap(SeqSymmetry intronSym , HashMap<String, SpecificUcscBedSym> map, BioSeq bioseq){
        int blockMins[] = new int[2];
        int blockMaxs[] = new int[2];
        SeqSpan span = intronSym.getSpan(bioseq);
        blockMins[0] = span.getMin() - threshold;
        blockMins[1] = span.getMax();
        blockMaxs[0] = span.getMin();
        blockMaxs[1] = span.getMax() + threshold;
        String name;
        boolean currentForward = true;
        SpecificUcscBedSym tempSym;
        int minimum = span.getMin();
        int maximum = span.getMax();
        while(minimum > offset)
            minimum = minimum-offset;
        while(maximum > offset)
            maximum = maximum-offset;
        if(maximum < minimum)
            maximum = maximum+offset;
        String leftResidues = residueString.substring(minimum, minimum+2);
        String rightResidues = residueString.substring(maximum-2,maximum);
        if(twoTracks){
            if(leftResidues.equals("GT") && rightResidues.equals("AG"))
                currentForward = true;
            else if(leftResidues.equals("CA") && rightResidues.equals("TC"))
                currentForward = false;
            else
                currentForward = span.isForward();
            if(currentForward)
                name = "J:"+bioseq.getID()+":"+span.getMin()+"-"+span.getMax()+":+";
            else
                name = "J:"+bioseq.getID()+":"+span.getMin()+"-"+span.getMax()+":-";
        }
        else{
            name = "J:"+bioseq.getID()+":"+span.getMin()+"-"+span.getMax()+":+";
            currentForward = true;
        }
        if(map.containsKey(name)){
            float score = map.get(name).getScore();
            map.get(name).setScore(++score);                    
        }
        else{
            tempSym = new SpecificUcscBedSym("test", bioseq, span.getMin()-threshold,
               span.getMax()+threshold, name, 1, currentForward, 0, 0, blockMins, blockMaxs);
            map.put(name,tempSym);
        }
    }
    
}
class SpecificUcscBedSym extends UcscBedSym{
    
    public SpecificUcscBedSym(String type, BioSeq seq, int txMin, int txMax, String name, float score,boolean forward, 
            int cdsMin, int cdsMax, int[] blockMins, int[] blockMaxs){
        super(type, seq, txMin, txMax, name, score, forward, cdsMin, cdsMax, blockMins, blockMaxs);
    }
    
    public void setScore(float score){
        this.score = score;
    }
}