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
import findjunction.filters.ThresholdFilter;
import java.util.*;


/**
 *
 * @author Anuj
 */
public class FindJunctionOperator implements Operator{

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
        HashMap<String, SpecificUcscBedSym> map = new HashMap<String , SpecificUcscBedSym>();
        int count;
        int blockMins[] = new int[2];
        int blockMaxs[] = new int[2];
        TypeContainerAnnot container = new TypeContainerAnnot("test", "bed");
        FindJunction junction = new FindJunction();
        ThresholdFilter filter = new ThresholdFilter();
        NoIntronFilter noIntronFilter = new NoIntronFilter();
        ChildThresholdFilter childFilter = new ChildThresholdFilter();
        SeqSymmetry intronSym;
        for(SeqSymmetry sym : list){
            if(noIntronFilter.filterSymmetry(bioseq, sym)){
                if(filter.filterSymmetry(bioseq, sym))
                    intronSym = SeqUtils.getIntronSym(sym, bioseq);
                else{
                    intronSym = SeqUtils.getIntronSym(sym, bioseq);
                    int childCount = sym.getChildCount();
                    List<SeqSymmetry>eligibleIntrons = new ArrayList<SeqSymmetry>();
                    for(int j=0;j<intronSym.getChildCount();j++)
                        eligibleIntrons.add(intronSym.getChild(j));
                    Object array[] = eligibleIntrons.toArray();
                    for(int j=0; j<childCount; j++){
                        if(!(childFilter.filterSymmetry(bioseq, sym.getChild(j)))){
                            if(j-1 >= 0)
                                array[j-1] = null;
                            if(j < childCount-1)
                                array[j] = null;
                        }
                    }
                    ((SimpleMutableSeqSymmetry)intronSym).removeChildren();
                    for(int j=0;j<array.length;j++){
                        if(array[j] != null)
                            ((SimpleMutableSeqSymmetry)intronSym).addChild((SeqSymmetry)array[j]);
                    }
                    eligibleIntrons.clear();
                }
                count = intronSym.getChildCount();
                for(int i=0;i<count;i++){
                    SeqSpan span = intronSym.getChild(i).getSpan(bioseq);
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
            }      
        }
        Set<String> keySet = map.keySet();
        Object keys[] = keySet.toArray();
        for(int i=0;i<keys.length;i++){
            container.addChild(map.get((String)keys[i]));
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