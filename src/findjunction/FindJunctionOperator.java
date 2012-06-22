/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.operator.Operator;
import com.affymetrix.genometryImpl.parsers.FileTypeCategory;
import com.affymetrix.genometryImpl.symmetry.*;
import com.affymetrix.genometryImpl.util.SeqUtils;
import findjunction.filters.ChildThresholdFilter;
import findjunction.filters.NoIntronFilter;
import java.util.*;


/**
 *
 * @author Anuj
 */
public class FindJunctionOperator implements Operator{
    private final int default_threshold = 5;
    int threshold = default_threshold;
    
    @Override
    public String getName() {
        return "findJunction";
    }

    @Override
    public String getDisplay() {
        return "Find Junction";
    }

    @Override
    public SeqSymmetry operate(BioSeq bioseq, List<SeqSymmetry> list) {
        int count;
        NoIntronFilter noIntronFilter = new NoIntronFilter();
        ChildThresholdFilter childFilter = new ChildThresholdFilter();
        List<SeqSymmetry> filteredIntrons = new ArrayList<SeqSymmetry>(); 
        SeqSymmetry currentIntron;
        for(SeqSymmetry sym : list){
            if(noIntronFilter.filterSymmetry(bioseq, sym)){
                List<SeqSymmetry> splitFilteredIntrons = split(sym , bioseq);
                filteredIntrons.addAll(splitFilteredIntrons);
            }
        }
        TypeContainerAnnot container = createJunctions(filteredIntrons, bioseq);
        return container;
    }
    
    public TypeContainerAnnot createJunctions(List<SeqSymmetry> filteredIntrons, BioSeq bioseq){
        TypeContainerAnnot container = new TypeContainerAnnot("test", "bed");
        HashMap<String, SpecificUcscBedSym> map = new HashMap<String , SpecificUcscBedSym>();
        for(SeqSymmetry intronSym : filteredIntrons){
            addToMap(intronSym, map, bioseq);
        }
        Set<String> keySet = map.keySet();
        Object keys[] = keySet.toArray();
        for(int i=0;i<keys.length;i++){
            container.addChild(map.get((String)keys[i]));
        }
        return container;
    }
    
    public void addToMap(SeqSymmetry intronSym , HashMap<String, SpecificUcscBedSym> map, BioSeq bioseq){
        int blockMins[] = new int[2];
        int blockMaxs[] = new int[2];
        FindJunction junction = new FindJunction();
        SeqSpan span = intronSym.getSpan(bioseq);
        blockMins[0] = span.getMin() - junction.getThreshold();
        blockMins[1] = span.getMax();
        blockMaxs[0] = span.getMin();
        blockMaxs[1] = span.getMax() + junction.getThreshold();
        String name = "J:"+bioseq.getID()+":"+span.getMin()+"-"+span.getMax();
        SpecificUcscBedSym tempSym = new SpecificUcscBedSym("test", bioseq, span.getMin()-junction.getThreshold(),
            span.getMax()+junction.getThreshold(), name, 1, true, 0, 0, blockMins, blockMaxs);
        if(map.containsKey(name)){
            float score = ((SpecificUcscBedSym)(map.get(name))).getScore();
            map.get(name).setScore(++score);                    
        }
        else
            map.put(name,tempSym);
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
    public List<SeqSymmetry> split(SeqSymmetry sym , BioSeq bioseq){
        List<SeqSymmetry> splitFilteredIntrons = new ArrayList<SeqSymmetry>();
        List<Integer> childIntronIndices = new ArrayList<Integer>();
        int childCount = sym.getChildCount();
        SeqSymmetry intronChild;
        ChildThresholdFilter filter = new ChildThresholdFilter();
        for(int i=0;i<childCount - 1;i++){
            if(filter.filterSymmetry(bioseq, sym.getChild(i)) && filter.filterSymmetry(bioseq, sym.getChild(i+1))){
                childIntronIndices.add(i);
            }
        }
        if(childIntronIndices.size() > 0 && childIntronIndices.size() != (sym.getChildCount()-1)){
            for(Integer i : childIntronIndices){
                SeqSymmetry currentIntron = SeqUtils.getIntronSym(sym, bioseq).getChild(i);
                splitFilteredIntrons.add(currentIntron);
            }
        }
        else if(childIntronIndices.size() == (sym.getChildCount() - 1)){
            SeqSymmetry intronSym = SeqUtils.getIntronSym(sym, bioseq);
            for(int i=0;i<intronSym.getChildCount() ;i++)
                splitFilteredIntrons.add(intronSym.getChild(i));
        }
        return splitFilteredIntrons;
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