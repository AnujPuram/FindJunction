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
import findjunction.filters.NoIntronFilter;
import findjunction.filters.UniqueLocationFilter;
import java.util.*;


/**
 *
 * @author Anuj
 */
public class FindJunctionOperator implements Operator{
    public static final int offset = 200000;
    private static final int default_threshold = 5;
    private SymmetryFilterI noIntronFilter = new NoIntronFilter();
    private SymmetryFilterI childThresholdFilter = new ChildThresholdFilter();
    private SymmetryFilterI uniqueLocationFilter = new UniqueLocationFilter();
    private int threshold = default_threshold;
    private boolean twoTracks, uniqueness;
    private TwoBit twoBit;
    private String residueString;
    private int random = 1;
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
    
    public void setResidueString(String residueString){
        this.residueString = residueString;
    }
    
    @Override
    public SeqSymmetry operate(BioSeq bioseq, List<SeqSymmetry> list) {
        TypeContainerAnnot container = new TypeContainerAnnot("test", "bed");
        HashMap<String, JunctionUcscBedSym> map = new HashMap<String , JunctionUcscBedSym>();
        subOperate(bioseq, list, map);
        Collection<JunctionUcscBedSym> symmetrySet = map.values();
        Object syms[] = symmetrySet.toArray();
        for(int i=0;i<syms.length;i++){
            container.addChild((JunctionUcscBedSym)syms[i]);
        }
        map.clear();
        symmetrySet.clear();
        return container;
    }
    
    public void subOperate(BioSeq bioseq, List<SeqSymmetry> list, HashMap<String, JunctionUcscBedSym> map){
      for(SeqSymmetry sym : list){
            if(noIntronFilter.filterSymmetry(bioseq, sym) && ((!uniqueness) || (uniqueness && uniqueLocationFilter.filterSymmetry(bioseq, sym)))){
                updateIntronHashMap(sym , bioseq, map);
            }

        }
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
    private void updateIntronHashMap(SeqSymmetry sym , BioSeq bioseq, HashMap<String, JunctionUcscBedSym> map){
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
                if(intronChild != null)
                    addToMap(intronChild, map, bioseq);
            }
        }
    }
    
    private void addToMap(SeqSymmetry intronSym , HashMap<String, JunctionUcscBedSym> map, BioSeq bioseq){
        int blockMins[] = new int[2];
        int blockMaxs[] = new int[2];
        boolean canonical = true;
        String rightResidues= "",leftResidues= "";
        SeqSpan span = intronSym.getSpan(bioseq);
        blockMins[0] = span.getMin() - threshold;
        blockMins[1] = span.getMax();
        blockMaxs[0] = span.getMin();
        blockMaxs[1] = span.getMax() + threshold;
        String name;
        boolean currentForward = false;
        JunctionUcscBedSym tempSym;
        int minimum = span.getMin();
        int maximum = span.getMax();
        if(!twoTracks){
            if(minimum >= residueString.length() || maximum >= residueString.length()){
                for(int j=residueString.length();j<maximum;j++)
                    residueString = residueString.concat("-");
            }
            leftResidues = residueString.substring(minimum, minimum+2);
            rightResidues = residueString.substring(maximum-2,maximum);
            boolean c;
            if(leftResidues.equalsIgnoreCase("GT") && rightResidues.equalsIgnoreCase("AG")){
                canonical = true;
                currentForward = true;
            }
            else if(leftResidues.equalsIgnoreCase("CT") && rightResidues.equalsIgnoreCase("AC")){
                canonical = true;
                currentForward = false;
            }
            else if((leftResidues.equalsIgnoreCase("AT") && rightResidues.equalsIgnoreCase("AC")) || 
                    (leftResidues.equalsIgnoreCase("GC") && rightResidues.equalsIgnoreCase("AG"))){
                canonical = false;
                currentForward = true;
            }
            else if((leftResidues.equalsIgnoreCase("GT") && rightResidues.equalsIgnoreCase("AT")) || 
                    (leftResidues.equalsIgnoreCase("CT") && rightResidues.equalsIgnoreCase("GC"))){
                canonical = false;
                currentForward = false;
            }
            else{
                canonical = false;
                currentForward = span.isForward();
            }
        }
        else{
            currentForward = span.isForward();
        }
        name = "J:" + bioseq.getID() + ":" + span.getMin() + "-" + span.getMax() + ":";
        String key = name;
        if(map.containsKey(key)){
            map.get(key).updateScore(currentForward);
        }
        else{
            tempSym = new JunctionUcscBedSym("test", bioseq, span.getMin()-threshold,
               span.getMax()+threshold, name, 1, currentForward, 0, 0, blockMins, blockMaxs, currentForward?1:0, currentForward?0:1, canonical);
            map.put(key,tempSym);
        }
    }
    
}
class JunctionUcscBedSym extends UcscBedSym{
    
    int positiveScore, negativeScore;
    float localScore = 1;
    boolean canonical;
    public JunctionUcscBedSym(String type, BioSeq seq, int txMin, int txMax, String name, float score,boolean forward, 
            int cdsMin, int cdsMax, int[] blockMins, int[] blockMaxs, int positiveScore, int negativeScore, boolean canonical){
        super(type, seq, txMin, txMax, name, score, forward, cdsMin, cdsMax, blockMins, blockMaxs);
        this.positiveScore = positiveScore;
        this.negativeScore = negativeScore;
        this.canonical = canonical;
    }
    
    public void updateScore(boolean isForward){
       localScore++;
        if(!canonical){
            if(isForward)
                this.positiveScore++;
            else
               this.negativeScore++;
        }
    }
    @Override
     public float getScore(){
        return localScore;
     }
    
     @Override
     public String getName(){
         return getID();
     }
     
     @Override
     public String getID(){
         return super.getID() + (isForward()? "+" : "-");
     }
     
     @Override
     public boolean isForward(){
         return canonical ? super.isForward() : positiveScore > negativeScore? true: false;
     }
}