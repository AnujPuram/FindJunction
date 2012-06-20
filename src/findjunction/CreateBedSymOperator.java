/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package findjunction;

import com.affymetrix.genometryImpl.BioSeq;
import com.affymetrix.genometryImpl.SeqSpan;
import com.affymetrix.genometryImpl.operator.Operator;
import com.affymetrix.genometryImpl.parsers.FileTypeCategory;
import com.affymetrix.genometryImpl.symmetry.BAMSym;
import com.affymetrix.genometryImpl.symmetry.SeqSymmetry;
import com.affymetrix.genometryImpl.symmetry.TypeContainerAnnot;
import com.affymetrix.genometryImpl.symmetry.UcscBedSym;
import com.affymetrix.genometryImpl.util.SeqUtils;
import java.util.*;

/**
 *
 * @author auser
 */
public class CreateBedSymOperator implements Operator{

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
        TypeContainerAnnot container = new TypeContainerAnnot("test", "bed");
        FindJunction junction = new FindJunction();
        HashMap<String, SeqSymmetry> map = new HashMap<String , SeqSymmetry>();
        SeqSymmetry intronSym;
        int count;
        int blockMins[] = new int[2];
        int blockMaxs[] = new int[2];
        for(SeqSymmetry sym : list){
            intronSym = SeqUtils.getIntronSym(sym, bioseq);
            count = intronSym.getChildCount();
            for(int i=0;i<count;i++){
                SeqSpan span = intronSym.getChild(i).getSpan(bioseq);
                blockMins[0] = span.getMin() - junction.getThreshold();
                blockMins[1] = span.getMax();
                blockMaxs[0] = span.getMin();
                blockMaxs[1] = span.getMax() + junction.getThreshold();
                String name = "J:"+bioseq.getID()+":"+span.getMin()+"-"+span.getMax();
                SeqSymmetry tempSym = new SpecificUcscBedSym("test", bioseq, span.getMin()-junction.getThreshold(),
                        span.getMax()+junction.getThreshold(), name, 1, true, 0, 0, blockMins, blockMaxs);
                if(map.containsKey(name)){
                    float score = ((SpecificUcscBedSym)(map.get(name))).getScore();
                    ((SpecificUcscBedSym)(map.get(name))).setScore(++score);                    
                }
                else
                    map.put(name,tempSym);
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
